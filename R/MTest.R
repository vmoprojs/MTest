####################################################
### File name: MTest.r
####################################################


MTest <- function(object, nboot = 100, nsam = NULL, trace = FALSE, seed = NULL,
                  valor_vif = 0.9) {
  stopifnot(inherits(object, "lm"))
  
  # Marco del modelo y matriz de diseño (congeladas)
  mf <- stats::model.frame(object)
  y  <- stats::model.response(mf)
  X  <- stats::model.matrix(object, mf)
  n  <- NROW(X)
  if (is.null(nsam)) nsam <- n
  if (!is.null(seed)) set.seed(seed)
  
  # Nombres de predictores (excluye intercepto)
  if (colnames(X)[1] != "(Intercept)") {
    stop("Se espera que la primera columna de X sea el intercepto.")
  }
  term_names <- colnames(X)[-1L]
  p <- length(term_names)
  
  # Funciones para R^2
  r2_fit <- function(Xm, yv) {
    fit <- stats::lm.fit(Xm, yv)
    rss <- sum(fit$residuals^2)
    tss <- sum((yv - mean(yv))^2)
    if (is.finite(rss) && tss > 0) 1 - rss/tss else NA_real_
  }
  r2_aux_j <- function(Xm, j) {
    xj <- Xm[, j, drop = FALSE]
    Z  <- Xm[, setdiff(seq_len(ncol(Xm)), j), drop = FALSE]
    fit <- stats::lm.fit(Z, xj)
    rss <- sum(fit$residuals^2)
    tss <- sum((xj - mean(xj))^2)
    if (is.finite(rss) && tss > 0) 1 - rss/tss else NA_real_
  }
  
  # Métricas del modelo original
  R2_global_tot <- summary(object)$r.squared
  R2j_tot <- vapply(seq_along(term_names), function(k) {
    jcol <- k + 1L
    r2_aux_j(X, jcol)
  }, numeric(1))
  names(R2j_tot) <- term_names
  VIF_tot <- 1 / (1 - R2j_tot)
  VIF_tot[!is.finite(VIF_tot)] <- NA_real_  # evitar Inf cuando R2j_tot ~ 1
  R.tot   <- c(global = R2_global_tot, R2j_tot)
  
  # Contenedores bootstrap
  Bvals   <- matrix(NA_real_, nrow = nboot, ncol = p + 1,
                    dimnames = list(NULL, c("global", term_names)))
  VIFvals <- matrix(NA_real_, nrow = nboot, ncol = p,
                    dimnames = list(NULL, term_names))
  
  if (trace) pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
  
  for (b in seq_len(nboot)) {
    idx <- sample.int(n, nsam, replace = TRUE)
    yb  <- y[idx]
    Xb  <- X[idx, , drop = FALSE]
    
    Bvals[b, "global"] <- r2_fit(Xb, yb)
    for (k in seq_along(term_names)) {
      jcol <- k + 1L
      R2j  <- r2_aux_j(Xb, jcol)
      Bvals[b, term_names[k]] <- R2j
      VIFvals[b, term_names[k]] <- if (is.na(R2j) || R2j >= 1) NA_real_ else 1 / (1 - R2j)
    }
    if (trace) utils::setTxtProgressBar(pb, b)
  }
  if (trace) close(pb)
  
  # --- p-values (con nombres correctos) ---
  stopifnot(colnames(Bvals)[1] == "global")
  pred_names <- colnames(Bvals)[-1L]
  
  # VIF: Pr(R^2_j > valor_vif)
  pval_vif <- colMeans(Bvals[, -1, drop = FALSE] > valor_vif, na.rm = TRUE)
  names(pval_vif) <- pred_names
  
  # Klein: Pr(R^2_global < R^2_j)
  g <- Bvals[, "global"]
  pval_klein <- sapply(pred_names, function(nm) {
    mean(g < Bvals[, nm], na.rm = TRUE)
  })
  names(pval_klein) <- pred_names
  
  # Armar objeto
  res <- list(
    Bvals      = Bvals,
    VIFvals    = VIFvals,
    pval_vif   = pval_vif,
    pval_klein = pval_klein,
    vif.tot    = VIF_tot,
    R.tot      = R.tot,
    nsam       = nsam,
    nboot      = nboot,
    call       = match.call()
  )
  class(res) <- "MTest"
  res
}







print.MTest <- function(x, digits = max(3, getOption("digits") - 3), sort_by = c("klein", "vif"),
                        show_counts = TRUE, ...) {
  stopifnot(is.list(x), !is.null(x$Bvals))
  sort_by <- match.arg(sort_by)
  
  # Restaurar opciones al salir
  op <- options(digits = digits)
  on.exit(options(op), add = TRUE)
  
  B <- x$Bvals
  Rtot <- x$R.tot
  VIFtot <- x$vif.tot
  p_vif <- x$pval_vif
  p_klein <- x$pval_klein
  nboot <- nrow(B)
  nsam <- x$nsam %||% NA_integer_
  
  # Helpers
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  fmt <- function(v, d = digits) formatC(v, digits = d, format = "fg", flag = "#")
  
  # Columnas (primera = global, resto predictores)
  cn <- colnames(B)
  if (is.null(cn) || length(cn) < 2L) {
    stop("Bvals debe tener columnas: 'global' y al menos un predictor.")
  }
  pred_names <- cn[-1L]
  
  # Resúmenes bootstrap (promedios) con NA-safe
  Rest <- colMeans(B, na.rm = TRUE)
  
  # Efectivos por columna (cuántos no-NA)
  eff_counts <- colSums(!is.na(B))
  eff_pred_counts <- eff_counts[-1L]
  # Conteo efectivo global (para R^2 global)
  eff_global <- eff_counts[1L]
  
  # Armar tabla de p-values y orden
  tab_p <- data.frame(
    predictor   = names(p_vif),
    p_vif       = as.numeric(p_vif[names(p_vif)]),
    p_klein     = as.numeric(p_klein[names(p_vif)]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Orden por severidad
  if (sort_by == "klein") {
    tab_p <- tab_p[order(-tab_p$p_klein, -tab_p$p_vif, tab_p$predictor), , drop = FALSE]
  } else {
    tab_p <- tab_p[order(-tab_p$p_vif, -tab_p$p_klein, tab_p$predictor), , drop = FALSE]
  }
  
  # Encabezado
  cat('\n##################################################################')
  cat('\nMTest: a nonparametric bootstrap test for multicollinearity\n')
  cat('------------------------------------------------------------------\n')
  cat('Bootstrap samples requested : ', nboot, '\n', sep = '')
  cat('Effective boots (global R^2): ', eff_global, '\n', sep = '')
  cat('Sample size per boot        : ', nsam, '\n', sep = '')
  cat('------------------------------------------------------------------\n')
  
  # P-values
  cat('Estimated p-values (sorted by ', sort_by, '):\n', sep = '')
  pmat <- cbind(
    predictor = tab_p$predictor,
    `vif p-value` = fmt(tab_p$p_vif),
    `klein p-value` = fmt(tab_p$p_klein)
  )
  print(noquote(pmat), right = TRUE)
  
  cat('Hypothesis:\n')
  cat(' - VIF Ho: R^2_j > threshold from VIF rule\n')
  cat(' - Klein Ho: R^2_global < R^2_j \n')
  
  cat('------------------------------------------------------------------\n')
  
  # R^2 observado (global y auxiliares)
  cat('Observed R^2 (global & per-predictor from VIF):\n')
  # Rtot viene con nombres: "global" y los predictores
  print(noquote(paste0(names(Rtot), " = ", fmt(Rtot))), right = FALSE)
  cat('------------------------------------------------------------------\n')
  
  # R^2 bootstrap (promedios)
  cat('Bootstrap mean R^2 (global & per-predictor):\n')
  # Rest tiene "global" + predictores
  print(noquote(paste0(names(Rest), " = ", fmt(Rest))), right = FALSE)
  cat('------------------------------------------------------------------\n')
  
  # VIF observado
  if (!is.null(VIFtot)) {
    cat('Observed VIF (per predictor):\n')
    vif_line <- paste0(names(VIFtot), " = ", fmt(VIFtot))
    print(noquote(vif_line), right = FALSE)
    cat('------------------------------------------------------------------\n')
  }
  
  # Conteos efectivos por predictor (útil si hubo NA en boots)
  if (isTRUE(show_counts)) {
    cat('Effective bootstrap counts per predictor:\n')
    cnt <- eff_pred_counts[names(p_vif)] %||% eff_pred_counts
    cnt_line <- paste0(names(cnt), " = ", cnt, " / ", nboot)
    print(noquote(cnt_line), right = FALSE)
    if (any(cnt < nboot)) {
      cat('\nNote: Some predictors had fewer effective resamples due to singularities.\n')
    }
    cat('------------------------------------------------------------------\n')
  }
  
  cat('Legend:\n')
  cat(' - vif p-value   = Pr( R^2_j > threshold from VIF rule )\n')
  cat(' - klein p-value = Pr( R^2_global < R^2_j )\n')
  cat('##################################################################\n')
  invisible(x)
}





plot.MTest <- function(x,type = 1,plotly = FALSE,...)
{
  if(!inherits(x,"MTest"))       stop("Enter an object obtained from the function MTest\n")
  ind <- values <- NULL
  boot.sol <- x
  boot.sol <- stack(data.frame(boot.sol$Bvals))
  boot.global <- boot.sol[boot.sol["ind"]=="global",]
  var =
    boot.aux <- boot.sol[boot.sol["ind"]!="global",]

    # ****** ECDF
  g_ecdf_global <- ggplot2::ggplot(boot.global,
                                   ggplot2::aes(values,color = ind)) +
    ggplot2::stat_ecdf(geom = "step")
  g_ecdf_aux <- ggplot2::ggplot(boot.aux,
                                ggplot2::aes(values,color = ind)) +
    ggplot2::stat_ecdf(geom = "step")
  g_ecdf_sol <- ggplot2::ggplot(boot.sol,
                                ggplot2::aes(values,color = ind)) +
    ggplot2::stat_ecdf(geom = "step")

  # ******* Density
  g_dens_sol <- ggplot2::ggplot(boot.sol,
                                ggplot2::aes(x=values,color = ind)) +
    ggplot2::geom_density()

  if(type == 2 & plotly == FALSE)
  {
    print(g_ecdf_sol)
  }
  if(type == 1 & plotly == FALSE)
  {
    print(g_dens_sol)
  }
  if(plotly==TRUE & type == 1)
  {
    g_dens_sol <- plotly::ggplotly(g_dens_sol)
    print(g_dens_sol)
  }
  if(plotly==TRUE & type == 2)
  {
    g_ecdf_sol <- plotly::ggplotly(g_ecdf_sol)
    print(g_ecdf_sol)
  }
}
