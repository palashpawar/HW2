## SDS 365 - Homework 2
## Generates every figure and summary number used on the web site.
## Run from the repository root:  Rscript analysis/generate_site_assets.R

source("R/pca_compress.R")

img_dir  <- "HW2imagefiles"
out_root <- "docs/assets"
extra_k  <- c(15, 20, 30, 50)   # extra ranks shown for the "how large should k be" discussion

## Write a matrix as a grayscale PNG at native resolution (1 pixel per entry).
## `lim` fixes the black/white mapping so every rank-k panel is directly
## comparable to the original.
write_gray_png <- function(M, file, lim) {
  M <- pmin(pmax((M - lim[1]) / diff(lim), 0), 1)
  png(file, width = ncol(M), height = nrow(M))
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  rasterImage(as.raster(M), 0, 0, 1, 1, interpolate = FALSE)
  dev.off()
}

summaries <- list()

for (i in 1:4) {
  name <- sprintf("image%d", i)
  cat("\n====", name, "====\n")
  X <- as.matrix(read.csv(file.path(img_dir, paste0(name, ".csv"))))
  n <- nrow(X); p <- ncol(X)
  out <- file.path(out_root, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  ## Display range: clip the extreme 2% at each end so that a handful of very
  ## bright pixels cannot flatten the rest of the image to a uniform black.
  ## The same mapping is reused for every rank-k panel of this image.
  lim <- as.numeric(quantile(X, c(0.02, 0.98)))
  fit <- pca_decompose(X)

  ## ---- 1. original + rank-k approximations -------------------------------
  write_gray_png(X, file.path(out, "orig.png"), lim)
  for (k in c(1:10, extra_k[extra_k <= p])) {
    write_gray_png(pca_rank_k(fit, k)$approx, file.path(out, sprintf("k%02d.png", k)), lim)
  }

  ## ---- 2. approximation error for every possible k ------------------------
  ks   <- 1:p
  errs <- vapply(ks, function(k) pca_rank_k(fit, k)$error, numeric(1))
  ## cross-check against the closed form ||X - X^(k)||_F^2 = sum_{j>k} lambda_j
  closed <- sqrt(pmax(rev(cumsum(rev(c(fit$values, 0))))[-1], 0))
  cat("max |loop - closed form| =", signif(max(abs(errs - closed)), 3), "\n")
  ## and against the stand-alone pca_compress() function at a few ranks
  spot <- unique(pmin(c(1, 2, 5, 10, p), p))
  cat("pca_compress() agrees at k =", spot, ":",
      isTRUE(all.equal(vapply(spot, function(k) pca_compress(X, k)$error, numeric(1)),
                       errs[spot])), "\n")

  fro_now <- norm(X, type = "F")
  cumv_now <- cumsum(fit$values) / sum(fit$values)
  kmark <- c(which(cumv_now >= 0.95)[1], which(cumv_now >= 0.99)[1], which(cumv_now >= 0.999)[1])

  png(file.path(out, "error.png"), width = 2000, height = 900, res = 170)
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3.0, 1.2))
  plot(ks, errs, type = "l", lwd = 2, col = "#2563eb",
       xlab = "rank k", ylab = expression(paste("|| X - X"^"(k)", " ||"[F])),
       main = "Frobenius error vs. k")
  points(1:10, errs[1:10], pch = 19, cex = 0.6, col = "#dc2626")
  grid(col = "#e5e7eb")
  legend("topright", c("all k", "k = 1..10"), lwd = c(2, NA), pch = c(NA, 19),
         col = c("#2563eb", "#dc2626"), bty = "n", cex = 0.8)

  rel  <- 100 * errs / fro_now
  keep <- which(rel > 1e-3)          # drop the numerically-zero tail
  plot(ks[keep], rel[keep], type = if (length(keep) > 3) "l" else "b", lwd = 2,
       pch = 19, col = "#2563eb", log = "y",
       xlab = "rank k", ylab = "relative error  100 x ||X - X^(k)|| / ||X||   (%)",
       main = "Relative error (log scale)")
  grid(col = "#e5e7eb")
  abline(v = kmark, lty = c(2, 3, 4), col = "#6b7280")
  legend("topright", paste0("k = ", kmark, c(" (95% var.)", " (99% var.)", " (99.9% var.)")),
         lty = c(2, 3, 4), col = "#6b7280", bty = "n", cex = 0.75)
  dev.off()

  ## ---- 3. first eigenvector and first PC score ---------------------------
  u1 <- fit$U[, 1]; z1 <- fit$Z[, 1]
  ## The sign of an eigenvector is arbitrary. Orient u1 so that its loadings are
  ## mostly positive; then a negative score means "this row is darker than the
  ## pattern u1 describes", which is easier to read.
  if (sum(u1) < 0) { u1 <- -u1; z1 <- -z1 }

  png(file.path(out, "pc1.png"), width = 1900, height = 900, res = 180)
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1))
  plot(seq_len(p), u1, type = "h", col = "#2563eb",
       xlab = "column of X (pixel column)", ylab = expression(u[1]),
       main = "First eigenvector")
  abline(h = 0, col = "grey60")
  plot(seq_len(n), z1, type = "l", col = "#b45309",
       xlab = "row of X (pixel row)", ylab = expression(z[1]),
       main = "First PC score")
  abline(h = 0, col = "grey60")
  dev.off()

  ## ---- 4. summary numbers ------------------------------------------------
  fro   <- norm(X, type = "F")
  cumv  <- cumsum(fit$values) / sum(fit$values)
  k_for <- function(prop) which(cumv >= prop)[1]
  summaries[[name]] <- list(
    n = n, p = p, frobenius = fro,
    err = errs, rel = errs / fro, cumvar = cumv,
    k95 = k_for(0.95), k99 = k_for(0.99), k999 = k_for(0.999),
    storage_ratio = (n + p) / (n * p),
    lim = lim, extra_k = extra_k[extra_k <= p]
  )
  cat(sprintf("||X||_F = %.1f | rel. error at k=1: %.4f  k=5: %.4f  k=10: %.4f\n",
              fro, errs[1]/fro, errs[5]/fro, errs[10]/fro))
  cat(sprintf("cum. variance %% at k=1,2,5,10: %s\n",
              paste(sprintf("%.2f", 100*cumv[c(1,2,5,10)]), collapse = ", ")))
  cat(sprintf("k for 95%%/99%%/99.9%% of variance: %d / %d / %d\n",
              k_for(0.95), k_for(0.99), k_for(0.999)))
}

saveRDS(summaries, "analysis/summaries.rds")

## Tidy table of every (image, k) pair, used to build the web pages.
tab <- do.call(rbind, lapply(names(summaries), function(nm) {
  s <- summaries[[nm]]
  data.frame(image = nm, n = s$n, p = s$p, frobenius = s$frobenius,
             k = seq_len(s$p), error = s$err, rel = s$rel, cumvar = s$cumvar,
             k95 = s$k95, k99 = s$k99, k999 = s$k999)
}))
write.csv(tab, "analysis/error_table.csv", row.names = FALSE)

## Flat text dump of the numbers that go into the write-ups.
sink("analysis/summary_table.txt")
for (nm in names(summaries)) {
  s <- summaries[[nm]]
  cat("==", nm, sprintf("(%d x %d)", s$n, s$p), "==\n")
  cat("||X||_F:", round(s$frobenius, 2), "\n")
  cat("k : error : relative error : cumulative variance %\n")
  for (k in c(1:10, s$extra_k, s$p)) {
    cat(sprintf("%4d : %12.3f : %8.5f : %8.4f\n", k, s$err[k], s$rel[k], 100 * s$cumvar[k]))
  }
  cat("k for 95/99/99.9% variance:", s$k95, s$k99, s$k999, "\n")
  cat("storage: rank-k needs k*(n+p) numbers vs n*p =", s$n * s$p,
      "| rank-1 ratio:", round(s$storage_ratio, 4), "\n\n")
}
sink()
cat("\nDone. Assets in", out_root, "\n")
