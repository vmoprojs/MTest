####################################################
### File name: pairwiseKStets.r
####################################################


pairwiseKStest <- function(
    X,
    alternative = c("greater","less","two.sided"),
    use = c("asis","pairwise.complete.obs"),
    exact = NULL
){
  alternative <- match.arg(alternative)
  use <- match.arg(use)
  
  if (!is.data.frame(X) && !is.matrix(X)) stop("X must be a data.frame or matrix.")
  X <- as.data.frame(X)
  is_num <- vapply(X, is.numeric, logical(1))
  if (!any(is_num)) stop("X has no numeric columns.")
  if (!all(is_num)) X <- X[, is_num, drop = FALSE]
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  
  cn <- colnames(X); n <- ncol(X)
  P <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(cn, cn))
  
  for (i in seq_len(n)) {
    xi_full <- X[[i]]
    for (j in seq_len(n)) {
      xj_full <- X[[j]]
      if (i == j) { P[i, j] <- 1; next }
      if (use == "pairwise.complete.obs") {
        ok <- stats::complete.cases(xi_full, xj_full)
        xi <- xi_full[ok]; xj <- xj_full[ok]
      } else {
        xi <- xi_full; xj <- xj_full
      }
      if (length(xi) < 1 || length(xj) < 1) { P[i, j] <- NA_real_; next }
      kt <- suppressWarnings(try(
        stats::ks.test(xi, xj, alternative = alternative, exact = exact),
        silent = TRUE
      ))
      P[i, j] <- if (inherits(kt, "try-error")) NA_real_ else as.numeric(kt$p.value)
    }
  }
  
  if (alternative == "less") {
    mes <- "alternative hypothesis: the CDF of x lies below that of y. Rows are x and Columns are y"
    ss <- colSums(P, na.rm = TRUE); names(ss) <- cn; sug <- sort(ss, decreasing = TRUE)
  } else if (alternative == "greater") {
    mes <- "alternative hypothesis: the CDF of x lies above that of y. Rows are x and Columns are y"
    ss <- rowSums(P, na.rm = TRUE); names(ss) <- cn; sug <- sort(ss, decreasing = TRUE)
  } else {
    mes <- "alternative hypothesis: two-sided"
    sug <- "No suggestions"
  }
  
  out <- list(
    KSpwMatrix  = P,
    alternative = mes,
    Suggestion  = sug
  )
  class(out) <- c("pairwiseKStest","list")
  out
}



print.pairwiseKStest <- function(x,
                                 digits = max(3, getOption("digits") - 3),
                                 show_rowmin = TRUE,   # resumen adicional util
                                 ...) {
  op <- options(digits = digits); on.exit(options(op), add = TRUE)
  
  P <- x$KSpwMatrix
  if (is.null(P) || !is.matrix(P)) stop("Objeto invalido: falta 'KSpwMatrix' (matriz de p-values).")
  rn <- rownames(P); cn <- colnames(P); n <- ncol(P)
  
  cat("\n============================================================\n")
  cat("pairwiseKStest summary\n")
  cat("------------------------------------------------------------\n")
  cat("Columns compared     : ", n, "\n", sep = "")
  cat("Alternative (as text): ", x$alternative, "\n", sep = "")
  cat("Note: rows are 'x' and columns are 'y' in ks.test(x, y).\n")
  
  # 1) Matriz de p-values (siempre)
  cat("------------------------------------------------------------\n")
  cat("P-value matrix (rows = x, cols = y):\n")
  Pfmt <- matrix(formatC(P, digits = digits, format = "fg", flag = "#"),
                 nrow = nrow(P), dimnames = dimnames(P))
  print(noquote(Pfmt))
  
  # 2) Suggestion (tu misma logica)
  cat("------------------------------------------------------------\n")
  cat("Suggestion:\n")
  if (is.character(x$Suggestion)) {
    cat(x$Suggestion, "\n")
  } else {
    topk <- x$Suggestion
    lineas <- paste0(names(topk), " = ",
                     formatC(as.numeric(topk), digits = digits, format = "fg", flag = "#"))
    print(noquote(lineas))
  }
  
  # 3) Resumen: minimos por fila (opcional)
  if (isTRUE(show_rowmin)) {
    cat("------------------------------------------------------------\n")
    cat("Row-wise minima (for each x, the y with smallest p):\n")
    # excluir diagonal al buscar minimos
    offdiag <- matrix(TRUE, nrow(P), ncol(P)); diag(offdiag) <- FALSE
    df <- data.frame(x = rn, y = NA_character_, p = NA_real_, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(P))) {
      vals <- P[i, ]
      vals[!offdiag[i, ]] <- NA_real_  # quitar diagonal
      if (all(is.na(vals))) next
      j <- which.min(vals)
      df$y[i] <- cn[j]
      df$p[i] <- vals[j]
    }
    df$p <- formatC(df$p, digits = digits, format = "fg", flag = "#")
    print(df, row.names = FALSE)
  }
  
  cat("============================================================\n")
  invisible(x)
}




