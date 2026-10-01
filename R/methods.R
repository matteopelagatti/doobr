#' @export
print.doob_fit <- function(x, ...) {
  s <- x$spec
  cat("Doob scenario-generation fit\n")
  cat(sprintf("  series        : %d modelling observations%s%s\n",
              length(s$w),
              if (s$transform != "none") paste0(", transform = ", s$transform) else "",
              if (s$difference) ", differenced" else ""))
  cat(sprintf("  mean lags     : %s\n", paste(s$mlags, collapse = ", ")))
  cat(sprintf("  scale model   : %s\n",
              if (!s$has_var) "constant (out-of-fold residual variance)"
              else sprintf("lags %s on %s, target = %s",
                           paste(s$vlags, collapse = ", "),
                           if (s$vlag_source == "y") "the series" else "squared residuals",
                           s$var_target)))
  cat(sprintf("  cross-fitting : %d folds of length %d, horizon %d\n",
              x$n_folds, x$fold_length, x$horizon))
  cat(sprintf("  bank          : %d atoms, mean %.2e, variance %.6f\n",
              length(x$bank), mean(x$bank), mean(x$bank^2)))
  cat(sprintf("  scale const   : c_hat = %.4f%s\n", x$c_hat,
              if (is.finite(x$kappa_hat)) sprintf(", kappa_hat = %+.4f", x$kappa_hat) else ""))
  if (x$n_scale_failures > 0L)
    cat(sprintf("  scale failures: %d state(s) excluded from the bank\n",
                x$n_scale_failures))
  if (!is.null(x$md$step1))
    cat(sprintf("  MDH p-values  : step 1 = %.4f, step 2 = %.4f\n",
                x$md$step1, x$md$step2))
  if (!is.null(x$eval))
    cat(sprintf("  mean scores   : ES = %.4f, VS = %.4f (%d paths per fold)\n",
                x$eval$mean_es, x$eval$mean_vs, x$eval$nsim))
  invisible(x)
}

#' @export
summary.doob_fit <- function(object, ...) {
  out <- list(fit = object,
              calibration = if (!is.null(object$eval))
                calibration_tests(object) else NULL)
  structure(out, class = "summary.doob_fit")
}

#' @export
print.summary.doob_fit <- function(x, ...) {
  print(x$fit)
  if (!is.null(x$calibration)) {
    cat("\n")
    print(x$calibration)
  }
  invisible(x)
}

#' @export
print.doob_calibration <- function(x, ...) {
  cat(sprintf("Calibration of the simulated paths (%d folds, %d paths each)\n",
              x$n_folds, x$nsim))
  cat(sprintf("  rank histogram            p = %.4f%s\n", x$rank$pval,
              if (x$rank$reject) "  *" else ""))
  cat(sprintf("  permutation, energy score p = %.4f%s  (Fisher)\n",
              x$perm_es$fisher_pval, if (x$perm_es$reject) "  *" else ""))
  cat(sprintf("  permutation, variogram    p = %.4f%s  (Fisher)\n",
              x$perm_vs$fisher_pval, if (x$perm_vs$reject) "  *" else ""))
  cat(sprintf("  band depth                p = %.4f%s  (Fisher)\n",
              x$band$fisher_pval, if (x$band$reject) "  *" else ""))
  cat(sprintf("  * significant at alpha = %.2f\n", x$alpha))
  invisible(x)
}

#' @export
print.doob_paths <- function(x, ...) {
  cat(sprintf("Simulated sample paths: %d steps x %d paths%s\n",
              x$horizon, x$nsim,
              if (x$transform != "none") sprintf(" (transform = %s)", x$transform) else ""))
  q <- stats::quantile(x$paths[x$horizon, ], c(.05, .5, .95))
  cat(sprintf("  terminal 5%%/50%%/95%%: %.4g / %.4g / %.4g\n", q[1], q[2], q[3]))
  if (x$transform != "none")
    cat(sprintf("  smearing at the final step: %.4f\n", x$smearing[x$horizon]))
  invisible(x)
}

#' @export
print.doob_tree <- function(x, ...) {
  cat(sprintf("Scenario tree: %d stages, branching %s, %d leaves, %d nodes\n",
              x$horizon, paste(x$branching, collapse = " x "),
              x$n_leaves, nrow(x$nodes)))
  invisible(x)
}

#' Plot simulated sample paths
#'
#' @param x A `"doob_paths"` object.
#' @param quantiles Band quantiles to shade.
#' @param main,xlab,ylab Passed to [graphics::plot()].
#' @param ... Further graphical parameters.
#' @return `x`, invisibly.
#' @export
plot.doob_paths <- function(x, quantiles = c(0.05, 0.25, 0.75, 0.95),
                            main = "Simulated sample paths",
                            xlab = "steps ahead", ylab = "value", ...) {
  P  <- x$paths
  st <- seq_len(nrow(P))
  qs <- apply(P, 1L, stats::quantile, probs = sort(quantiles))
  med <- apply(P, 1L, stats::median)
  graphics::plot(st, med, type = "n", ylim = range(qs, med),
                 main = main, xlab = xlab, ylab = ylab, ...)
  nb <- nrow(qs) %/% 2L
  for (i in seq_len(nb))
    graphics::polygon(c(st, rev(st)), c(qs[i, ], rev(qs[nrow(qs) - i + 1L, ])),
                      col = grDevices::adjustcolor("steelblue", 0.18 + 0.16 * (i - 1)),
                      border = NA)
  graphics::lines(st, med, col = "steelblue", lwd = 2)
  invisible(x)
}
