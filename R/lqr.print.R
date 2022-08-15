lqr.print <- function (model, digits = 3, ...) 
{
  if (!inherits(model,"qr"))
    stop(paste("Class of object", class(model), "not recognized.",  sep = " "))
  cat("\n")
  cat("Fit:",model$dist , "\n")
  cat("Quantile:", model$p, "\n")
  cat("\n")
  cat("Estimates:\n")
  cat("\n")
  print(model$table)
  cat("---\n")
  cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
  cat("\n")
  cat("sigma =", model$sigma, 
      "\n")            
  if (model$dist == "normal" || model$dist == "laplace") {cat("\n")}
  if (model$dist == "t" || model$dist == "slash") {
    cat("nu    =", model$nu, "\n")
    cat("\n")
  }
  if (model$dist == "cont") {
    cat("nu    =", model$nu, "\n")
    cat("gamma =", model$gamma, "\n")
    cat("\n")
  }
  cat('Model selection criteria:\n')
  cat('\n')
  critFin <- c(model$loglik, model$AIC, model$BIC, model$HQ)
  critFin <- round(t(as.matrix(critFin)),digits)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
  print(critFin)
  cat('\n')
  cat("Processing time =",model$time,units(model$time))
  cat('\n')
  
  invisible(model)
}

#lqr.print(model, digits = 3)
