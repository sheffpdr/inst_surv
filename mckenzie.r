#mckenzie method - need to check what assumptions are required for this to be sensible
# also could do with knowing how to incorporate other covariates

trust_props <- aggregate(reg_inst$SurgeryRule2, by = list(reg_inst$Trust), FUN = mean, na.rm=T)
reg_inst$trust_prop <- trust_props[match(reg_inst$Trust, trust_props[,1]),2]
reg_inst <- reg_inst[which(!is.na(reg_inst$SurgeryRule2)),]
reg_inst2 <- reg_inst[which(reg_inst$StageSimple!="IV"),]

library (rootSolve)
IVHR <- function (Time, Status, X, W) {
  # Status is 1 if event occurred , 0 for censored
  # Time is time–to–event Status=l, and time of censoring otherwise
  # X is treatment or exposure of interest
  # W is the instrument
  n <- length (Time)
  ord <- order(-Time)
  Time <- Time [ord]
  Status <- Status[ord]
  X <- X[ord]
  W <- W[ord]
  Est.Equat <- function (beta) {
    HR <- exp(beta *X)
    S.0 <- cumsum(HR)
    S.W1 <- cumsum(W*HR)
    sum(Status * (W - S.W1/S.0))
  }
  out.solution <- multiroot (Est.Equat, start=0)
  beta.hat <- ifelse (out.solution$estim.precis < 0.00001, out.solution$root , NA)
  HR <- exp(beta.hat *X)
  S.0 <-cumsum(HR)
  S.W1 <- cumsum(W*HR)
  S.X1 <- cumsum(X*HR)
  S.W1X1 <- cumsum(W*X*HR)
  Var.E.E <- sum(Status * (W - S.W1/S.0)^2)
  Deriv <- -sum(Status * (S.W1X1/S.0 - (S.W1*S.X1)/S.0^2))
  SE <- sqrt(Var.E.E/Deriv ^ 2)
  list (Est.log.HR=beta.hat , SE=SE, SE.2=SE^2)
}