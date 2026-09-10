# SDS 365 — Homework 2: PCA for Image Compression

Rank-k reconstruction of four images using principal components analysis, plus a static site and an
interactive Shiny app.

**Site:** https://palashpawar.github.io/HW2/
**App design mockup:** https://palashpawar.github.io/HW2/mockup/
**Shiny app:** https://palashpawar.shinyapps.io/pca-image-compressor/

## Layout

| path | what it is |
| --- | --- |
| `R/pca_compress.R` | the compression function: takes an image matrix and a rank `k`, returns the rank-k approximation and its Frobenius error |
| `analysis/generate_site_assets.R` | rebuilds every figure and table used on the site |
| `analysis/error_table.csv`, `summary_table.txt` | error and variance for every `(image, k)` pair |
| `docs/` | the published site (GitHub Pages source: `main` branch, `/docs` folder) |
| `docs/mockup/` | the app design mockup and layout rationale |
| `shiny_app/` | the Shiny implementation (`pca_compress.R` is duplicated here so the deployment is self-contained) |
| `HW2imagefiles/` | the four images from Canvas, as CSV matrices |

## The function

```r
source("R/pca_compress.R")

X   <- as.matrix(read.csv("HW2imagefiles/image3.csv"))
out <- pca_compress(X, k = 30)

out$error       # 38.16  =  ||X - X^(30)||_F
dim(out$approx) # 376 345
```

`V = X'X` is eigendecomposed with `eigen()`, the scores are `Z = XU`, and the rank-k approximation is
`Z[, 1:k] %*% t(U[, 1:k])`. No centering — the identity `X = ZU'` only holds for the uncentered
decomposition.

`pca_decompose()` / `pca_rank_k()` in the same file do the decomposition once so that many values of
`k` can be evaluated without repeating it; they return identical numbers, which
`analysis/generate_site_assets.R` verifies at several ranks and against the closed form
`||X - X^(k)||_F^2 = λ_{k+1} + … + λ_p`.

## Results

| image | size | k = 1 variance | k for 95% | k for 99% | chosen k | relative error there |
| --- | --- | --- | --- | --- | --- | --- |
| image1 | 1000 × 100 | 80.2% | 63 | 91 | 1 | 44.5% (noise is not compressible) |
| image2 | 1000 × 100 | 100.0% | 1 | 1 | 1 | 0.006% (the image is rank 1) |
| image3 | 376 × 345 | 61.8% | 10 | 33 | 30 | 10.6% |
| image4 | 213 × 119 | 40.5% | 11 | 28 | 20 | 13.3% |

## Reproducing

```sh
Rscript analysis/generate_site_assets.R   # regenerate docs/assets
R -e 'shiny::runApp("shiny_app")'         # run the app locally
```

Requires R with `shiny`, `jpeg` and `png` (the analysis script itself only needs base R).

## Deploying the app

```r
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
rsconnect::deployApp("shiny_app", appName = "pca-image-compressor")
```
