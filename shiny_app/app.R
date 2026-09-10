## SDS 365 - Homework 2
## PCA Image Compressor: upload an image, choose a rank k, compare the
## rank-k reconstruction against the original.

library(shiny)

source("pca_compress.R")

MAX_SIDE   <- 900         # uploaded images are subsampled to at most this many pixels per side
MAX_UPLOAD <- 12          # MB
options(shiny.maxRequestSize = MAX_UPLOAD * 1024^2)

## ---------------------------------------------------------------- helpers

## Read an upload into a numeric matrix. Photos are converted to grayscale
## with the usual luminance weights; a CSV is taken as a matrix as-is.
read_image_matrix <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext == "csv") {
    M <- as.matrix(utils::read.csv(path))
    if (!is.numeric(M)) stop("That CSV does not contain a purely numeric matrix.")
    dimnames(M) <- NULL
    return(M)
  }
  arr <- switch(ext,
    "jpg"  = jpeg::readJPEG(path),
    "jpeg" = jpeg::readJPEG(path),
    "png"  = png::readPNG(path),
    stop("Unsupported file type: .", ext))
  if (length(dim(arr)) == 2) return(arr)
  if (dim(arr)[3] >= 3) {
    return(0.299 * arr[, , 1] + 0.587 * arr[, , 2] + 0.114 * arr[, , 3])
  }
  arr[, , 1]
}

## Keep the eigendecomposition affordable: stride-subsample very large images.
shrink <- function(M, max_side = MAX_SIDE) {
  f <- max(1, ceiling(max(dim(M)) / max_side))
  if (f == 1) return(M)
  M[seq(1, nrow(M), by = f), seq(1, ncol(M), by = f), drop = FALSE]
}

## Draw a matrix as a grayscale image, using a fixed black/white mapping so
## that the original and the reconstruction are directly comparable.
draw <- function(M, lim) {
  S <- pmin(pmax((M - lim[1]) / diff(lim), 0), 1)
  op <- par(mar = c(0, 0, 0, 0), bg = "white"); on.exit(par(op))
  plot.new()
  plot.window(c(0, 1), c(0, 1), asp = nrow(M) / ncol(M))
  rasterImage(as.raster(S), 0, 0, 1, 1, interpolate = FALSE)
}

fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)

## ---------------------------------------------------------------- ui

ui <- fluidPage(
  title = "PCA Image Compressor",
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif; }
    .topbar { background:#1d4ed8; color:#fff; padding:14px 20px; border-radius:8px; margin-bottom:18px; }
    .topbar h2 { margin:0; font-size:20px; font-weight:700; }
    .topbar p  { margin:2px 0 0; font-size:12.5px; color:#c7d2fe; }
    .step { font-weight:700; font-size:13px; margin:14px 0 6px; }
    .step span { display:inline-block; width:18px; height:18px; border-radius:9px; background:#1d4ed8;
                 color:#fff; text-align:center; font-size:11px; line-height:18px; margin-right:6px; }
    .statrow { display:flex; justify-content:space-between; font-size:13px; padding:3px 0;
               border-bottom:1px solid #f3f4f6; }
    .statrow b { font-variant-numeric: tabular-nums; }
    .good { color:#15803d; } .bad { color:#b91c1c; }
    .panelhead { font-size:13px; font-weight:700; margin-bottom:4px; }
    .panelhead span { font-weight:400; color:#9ca3af; }
    .hint { font-size:12px; color:#6b7280; }
    .imgbox { border:1px solid #e5e7eb; border-radius:6px; padding:4px; background:#fff; }
    .imgbox.rec { border-color:#1d4ed8; }
  "))),

  div(class = "topbar",
      h2("PCA Image Compressor"),
      p("Rebuild any image from its first k principal components")),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      div(class = "step", span("1"), "Upload an image"),
      fileInput("file", NULL, accept = c(".jpg", ".jpeg", ".png", ".csv"),
                buttonLabel = "Browse", placeholder = "no file selected"),
      div(class = "hint", sprintf(".jpg / .png / .csv matrix, up to %d MB. Colour images are
          converted to grayscale; anything larger than %d pixels a side is subsampled so the
          eigendecomposition stays quick.", MAX_UPLOAD, MAX_SIDE)),
      uiOutput("filenote"),

      div(class = "step", span("2"), "Choose the rank k"),
      sliderInput("k", NULL, min = 1, max = 100, value = 10, step = 1, width = "100%"),
      div(class = "hint", "left = more compression, right = more detail"),
      br(),
      div(class = "hint", "presets (target relative error 25% / 10% / 2%)"),
      div(actionButton("p_tiny", "tiny", class = "btn-sm"),
          actionButton("p_bal", "balanced", class = "btn-sm"),
          actionButton("p_exact", "near-exact", class = "btn-sm")),

      hr(),
      uiOutput("stats"),
      br(),
      downloadButton("dl", "Download compressed image", style = "width:100%")
    ),

    mainPanel(
      width = 8,
      conditionalPanel("output.loaded != 'yes'",
        div(class = "hint", style = "padding:30px 0",
            "Upload an image to begin. No file handy? Any of the",
            tags$code("HW2imagefiles/*.csv"), "matrices from the homework works,
            and so does any photo on your machine.")),
      conditionalPanel("output.loaded == 'yes'",
        fluidRow(
          column(6, div(class = "panelhead", "Original ", textOutput("dim_o", inline = TRUE)),
                    div(class = "imgbox", plotOutput("orig", height = "330px"))),
          column(6, div(class = "panelhead", textOutput("rec_head", inline = TRUE)),
                    div(class = "imgbox rec", plotOutput("rec", height = "330px")))
        ),
        br(),
        div(class = "panelhead", "Error against k ", tags$span("- the red marker is where you are")),
        plotOutput("curve", height = "230px", click = "curve_click"),
        div(class = "hint", "Click the curve to jump to that rank.")
      )
    )
  )
)

## ---------------------------------------------------------------- server

server <- function(input, output, session) {

  ## Read + decompose once per upload; every value of k reuses this.
  img <- reactive({
    req(input$file)
    M <- tryCatch(shrink(read_image_matrix(input$file$datapath, input$file$name)),
                  error = function(e) {
                    showNotification(conditionMessage(e), type = "error", duration = 8)
                    NULL
                  })
    req(M)
    if (nrow(M) < 2 || ncol(M) < 2) {
      showNotification("That file is too small to decompose.", type = "error")
      return(NULL)
    }
    fit  <- pca_decompose(M)
    errs <- sqrt(pmax(rev(cumsum(rev(c(fit$values, 0))))[-1], 0))   # ||X - X^(k)||_F for all k
    cumv <- cumsum(fit$values) / sum(fit$values)
    list(M = M, fit = fit, errs = errs, cumv = cumv,
         fro = norm(M, type = "F"),
         lim = as.numeric(stats::quantile(M, c(0.02, 0.98))))
  })

  ## Smallest k whose relative error is at or below `tol`. Relative error is a
  ## better handle than variance explained here: the images are not centered, so
  ## the first component alone already accounts for most of the "variance" of any
  ## all-positive image while still looking heavily blurred.
  k_for_error <- function(d, tol) {
    hit <- which(d$errs / d$fro <= tol)
    if (length(hit)) hit[1] else length(d$errs)
  }

  observeEvent(img(), {
    d <- img()
    updateSliderInput(session, "k", min = 1, max = ncol(d$M), value = k_for_error(d, 0.10))
  })

  preset <- function(tol) {
    d <- req(img())
    updateSliderInput(session, "k", value = k_for_error(d, tol))
  }
  observeEvent(input$p_tiny,  preset(0.25))
  observeEvent(input$p_bal,   preset(0.10))
  observeEvent(input$p_exact, preset(0.02))

  observeEvent(input$curve_click, {
    d <- req(img())
    updateSliderInput(session, "k", value = max(1, min(ncol(d$M), round(input$curve_click$x))))
  })

  k <- reactive({
    d <- req(img())
    max(1, min(as.integer(input$k), ncol(d$M)))
  })

  approx <- reactive(pca_rank_k(req(img())$fit, k())$approx)

  output$loaded <- reactive(if (is.null(input$file)) "no" else "yes")
  outputOptions(output, "loaded", suspendWhenHidden = FALSE)

  output$filenote <- renderUI({
    d <- req(img())
    div(class = "hint", style = "margin-top:-6px",
        sprintf("%s - using %d x %d = %s numbers", input$file$name,
                nrow(d$M), ncol(d$M), format(length(d$M), big.mark = ",")))
  })

  output$dim_o    <- renderText(sprintf("(%d x %d)", nrow(req(img())$M), ncol(req(img())$M)))
  output$rec_head <- renderText(sprintf("Compressed - k = %d", k()))

  output$orig <- renderPlot(draw(req(img())$M, req(img())$lim))
  output$rec  <- renderPlot(draw(approx(), req(img())$lim))

  output$stats <- renderUI({
    d <- req(img()); kk <- k()
    n <- nrow(d$M); p <- ncol(d$M)
    store <- kk * (n + p); full <- n * p
    ratio <- store / full
    div(
      div(class = "statrow", span("rank k"), tags$b(kk)),
      div(class = "statrow", span("Frobenius error"), tags$b(sprintf("%.3f", d$errs[kk]))),
      div(class = "statrow", span("relative error"), tags$b(fmt_pct(d$errs[kk] / d$fro))),
      div(class = "statrow", span("variance kept"), tags$b(fmt_pct(d$cumv[kk]))),
      div(class = "statrow", span("numbers stored"), tags$b(format(store, big.mark = ","))),
      div(class = "statrow", span("size vs original"),
          tags$b(class = if (ratio <= 1) "good" else "bad", fmt_pct(ratio)))
    )
  })

  output$curve <- renderPlot({
    d <- req(img()); kk <- k(); p <- ncol(d$M)
    rel <- 100 * d$errs / d$fro
    op <- par(mar = c(4, 4.4, 0.6, 1)); on.exit(par(op))
    plot(seq_len(p), rel, type = "l", lwd = 2, col = "#1d4ed8",
         xlab = "rank k", ylab = "relative error (%)")
    grid(col = "#eef0f3")
    abline(v = kk, col = "#9ca3af", lty = 2)
    points(kk, rel[kk], pch = 19, col = "#dc2626", cex = 1.2)
  })

  output$dl <- downloadHandler(
    filename = function() {
      base <- tools::file_path_sans_ext(basename(input$file$name))
      sprintf("%s_pca_k%d.png", base, k())
    },
    content = function(file) {
      d <- req(img())
      S <- pmin(pmax((approx() - d$lim[1]) / diff(d$lim), 0), 1)
      png::writePNG(S, file)
    }
  )
}

shinyApp(ui, server)
