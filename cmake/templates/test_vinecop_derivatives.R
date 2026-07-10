#!/usr/bin/env Rscript
# Oracle for the analytic full-gradient tests: evaluates RVineGrad and
# RVineHessian on a fixed 5-dimensional D-vine and writes the results
# reordered into vinecopulib's parameter order (trees, then edges, then
# parameters; vinecopulib pair copula (t, e) corresponds to VineCopula
# matrix cell [d - t, e + 1]) with the 90/270-degree sign convention fixed
# (VineCopula negates the parameter for those rotation codes, so gradient
# entries flip sign once per parameter differentiation).
#
# Empirically, RVineGrad's output vector iterates the sub-diagonal cells
# column by column from the left, top to bottom within a column (not the
# order stated in its documentation), with the t-copulas' par2 entries
# appended in the same cell order.
#
# NOTE: VineCopula's RVineGrad/RVineHessian are inconsistent with finite
# differences of RVineLogLik for models containing 90/270-rotated pair
# copulas (wrong entries when the rotated edge feeds deeper trees, and
# garbage output with a rotated family in the last tree), so this model
# uses only 0/180 rotations; 90/270 coverage of the analytic cascade comes
# from the finite-difference tests on the C++ side. The self-checks below
# guard both the order decoding and this assumption.
#
# The model here must stay in sync with the C++ side
# (test_vinecop_class.hpp, VinecopDerivatives tests).
args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[1])

if (!("VineCopula" %in% rownames(installed.packages()))) {
  install.packages("VineCopula", repos = "http://cran.rstudio.com/")
}
suppressMessages(library(VineCopula))

d <- 5
Matrix <- D2RVine(1:d, family = rep(0, d * (d - 1) / 2),
                  par = rep(0, d * (d - 1) / 2))$Matrix
fam <- matrix(0, d, d)
par <- matrix(0, d, d)
par2 <- matrix(0, d, d)
# tree k lives in row d - k + 1; columns 1..d-k
fam[5, 1] <- 1; par[5, 1] <- 0.6 # gaussian
fam[5, 2] <- 3; par[5, 2] <- 2.5 # clayton
fam[5, 3] <- 14; par[5, 3] <- 1.8 # gumbel 180 (survival)
fam[5, 4] <- 6; par[5, 4] <- 2.2 # joe
fam[4, 1] <- 5; par[4, 1] <- 4.0 # frank
fam[4, 2] <- 2; par[4, 2] <- 0.5; par2[4, 2] <- 6 # student
fam[4, 3] <- 13; par[4, 3] <- 1.2 # clayton 180 (survival)
fam[3, 1] <- 4; par[3, 1] <- 1.4 # gumbel
fam[3, 2] <- 1; par[3, 2] <- 0.2 # gaussian
fam[2, 1] <- 3; par[2, 1] <- 2.5 # clayton
RVM <- RVineMatrix(Matrix, family = fam, par = par, par2 = par2)

set.seed(0)
u <- matrix(0.01 + 0.98 * runif(n * d), n, d)

# RVineGrad's cell order (see the note above)
cells <- NULL
for (cc in 1:(d - 1)) {
  for (rr in (cc + 1):d) {
    cells <- rbind(cells, c(rr, cc))
  }
}
grad_index <- function(r, c) which(cells[, 1] == r & cells[, 2] == c)
dd <- nrow(cells)
par2_cells <- which(fam[cells] == 2)

# vinecopulib order: trees, then edges, then parameters of the pair copula,
# with the sign flip for 90/270 rotation codes
reorder_grad <- function(g) {
  out <- NULL
  for (t in 0:(d - 2)) {
    for (e in 0:(d - 2 - t)) {
      r <- d - t
      c <- e + 1
      s <- if (fam[r, c] %in% c(23, 24, 26, 33, 34, 36)) -1 else 1
      out <- c(out, s * g[grad_index(r, c)])
      if (fam[r, c] == 2) {
        out <- c(out, g[dd + which(par2_cells == grad_index(r, c))])
      }
    }
  }
  out
}

# self-check 1: the decoded RVineGrad must match finite differences of
# RVineLogLik cell by cell (guards the order decoding and the exactness of
# RVineGrad for this model)
g <- RVineGrad(u, RVM)$gradient
h <- 1e-5
g_fd <- rep(0, dd + length(par2_cells))
for (i in seq_len(dd)) {
  r <- cells[i, 1]
  c <- cells[i, 2]
  pp <- RVM$par
  pm <- RVM$par
  pp[r, c] <- pp[r, c] + h
  pm[r, c] <- pm[r, c] - h
  g_fd[i] <- (RVineLogLik(u, RVM, par = pp, separate = FALSE)$loglik -
    RVineLogLik(u, RVM, par = pm, separate = FALSE)$loglik) / (2 * h)
}
for (j in seq_along(par2_cells)) {
  i <- par2_cells[j]
  r <- cells[i, 1]
  c <- cells[i, 2]
  pp <- RVM$par2
  pm <- RVM$par2
  pp[r, c] <- pp[r, c] + h
  pm[r, c] <- pm[r, c] - h
  g_fd[dd + j] <- (RVineLogLik(u, RVM, par = RVM$par, par2 = pp,
    separate = FALSE)$loglik -
    RVineLogLik(u, RVM, par = RVM$par, par2 = pm,
      separate = FALSE)$loglik) / (2 * h)
}
if (max(abs(g - g_fd) / pmax(1, abs(g_fd))) > 1e-4) {
  stop("RVineGrad order decoding is inconsistent with finite differences")
}

grad <- reorder_grad(g)
grad1 <- reorder_grad(RVineGrad(u[1, , drop = FALSE], RVM)$gradient)

# Hessian oracle: central finite differences of RVineGrad. (VineCopula's
# RVineHessian output is inconsistent with finite differences of its own
# RVineGrad -- by O(1) relative errors, even for models without rotations --
# so the differentiated gradient is the reliable VineCopula-anchored
# reference for the Hessian of the log-likelihood.)
H_fd <- matrix(0, dd + length(par2_cells), dd + length(par2_cells))
for (i in seq_len(dd)) {
  r <- cells[i, 1]
  c <- cells[i, 2]
  pp <- RVM$par
  pm <- RVM$par
  pp[r, c] <- pp[r, c] + h
  pm[r, c] <- pm[r, c] - h
  RVMp <- RVineMatrix(Matrix, family = fam, par = pp, par2 = par2)
  RVMm <- RVineMatrix(Matrix, family = fam, par = pm, par2 = par2)
  H_fd[, i] <- (RVineGrad(u, RVMp)$gradient -
    RVineGrad(u, RVMm)$gradient) / (2 * h)
}
for (j in seq_along(par2_cells)) {
  i <- par2_cells[j]
  r <- cells[i, 1]
  c <- cells[i, 2]
  pp <- par2
  pm <- par2
  pp[r, c] <- pp[r, c] + h
  pm[r, c] <- pm[r, c] - h
  RVMp <- RVineMatrix(Matrix, family = fam, par = RVM$par, par2 = pp)
  RVMm <- RVineMatrix(Matrix, family = fam, par = RVM$par, par2 = pm)
  H_fd[, dd + j] <- (RVineGrad(u, RVMp)$gradient -
    RVineGrad(u, RVMm)$gradient) / (2 * h)
}
npars <- length(grad)
perm <- reorder_grad(seq_len(dd + length(par2_cells)))
signs <- sign(perm)
perm <- abs(perm)
Hout <- matrix(0, npars, npars)
for (i in seq_len(npars)) {
  for (j in seq_len(npars)) {
    Hout[i, j] <- signs[i] * signs[j] * H_fd[perm[i], perm[j]]
  }
}

write.table(u, file = "temp_vderiv_data", col.names = FALSE, row.names = FALSE)
write.table(RVM$Matrix, file = "temp_vderiv_matrix", col.names = FALSE,
            row.names = FALSE)
write.table(cbind(grad, grad1), file = "temp_vderiv_grad", col.names = FALSE,
            row.names = FALSE)
write.table(Hout, file = "temp_vderiv_hess", col.names = FALSE,
            row.names = FALSE)
