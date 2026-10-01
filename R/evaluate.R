## Per-fold evaluation: simulate over each fold from the models that never saw
## it, using a leave-one-fold-out bank with its own scale constant, and score
## the realised path against the ensemble.
##
## Called from doob_fit(); not exported.
evaluate_folds <- function(object, nsim, xreg_hist, z_raw, w, y_t) {
  spec     <- object$spec
  K        <- object$n_folds
  bank_idx <- object$bank_idx

  fold_bank <- function(fold) {
    keep <- !(bank_idx %in% fold)
    if (sum(keep) < 2L) return(list(bank = object$bank, c = spec$c_hat))
    zk <- z_raw[keep] - mean(z_raw[keep])
    ck <- sqrt(mean(zk^2))
    if (!is.finite(ck) || ck <= 0) return(list(bank = object$bank, c = spec$c_hat))
    list(bank = zk / ck, c = ck)
  }

  window_sims       <- vector("list", K)
  window_sims_model <- vector("list", K)
  window_es <- numeric(K); window_vs <- numeric(K)

  for (k in seq_len(K)) {
    fold  <- object$folds[[k]]
    start <- fold[1L]
    eff_h <- length(fold)
    fb    <- fold_bank(fold)
    raw <- run_recursion(spec, object$models$fold_mean[[k]],
                         object$models$fold_scale[[k]],
                         w[seq_len(start - 1L)],
                         if (spec$vlag_source == "sq_res")
                           spec$sq_series[seq_len(start - 1L)] else NULL,
                         fb$bank, eff_h, nsim,
                         if (is.null(xreg_hist)) NULL
                         else xreg_hist[start:(start + eff_h - 1L), , drop = FALSE],
                         fb$c)
    ## reference level: the last transformed-scale value before the fold
    ref <- if (spec$difference) y_t[start] else NA_real_
    both <- to_original(spec, raw, ref)
    actual_model <- if (spec$difference) y_t[(start + 1L):(start + eff_h)]
                    else w[start:(start + eff_h - 1L)]
    actual <- spec$trf$inv(actual_model)

    window_sims_model[[k]] <- cbind(actual_model, both$model)
    window_sims[[k]]       <- cbind(actual, both$original)
    window_es[k] <- es_score(actual, both$original)
    window_vs[k] <- vs_score(actual, both$original)
  }

  list(window_sims = window_sims, window_sims_model = window_sims_model,
       window_es = window_es, window_vs = window_vs,
       mean_es = mean(window_es), mean_vs = mean(window_vs),
       nsim = nsim)
}
