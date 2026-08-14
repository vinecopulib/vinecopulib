#!/usr/bin/env Rscript
# Oracle for the analytic full-gradient / joint-Hessian tests. Evaluates
# VineCopula's RVineGrad and RVineHessian directly (no finite differences) on
# a fixed 5-dimensional D-vine and writes them for the C++ side to compare
# against, element-wise.
#
# Requires VineCopula >= 2.6.2, which fixes RVineGrad/RVineHessian for models
# with rotated families and the RVineHessian sample-size bug (tnagler/
# VineCopula #101, #102). The #101 fix is what lets us use RVineHessian
# *directly* here (it now depends on the data) instead of finite-differencing
# RVineGrad.
#
# Rotations at the vine level: this model deliberately uses only 0/180
# families. vinecopulib and VineCopula do not share the same structure /
# swapped-rotation conventions, so a vine with 90/270-rotated families cannot
# be constructed identically in both libraries (their log-likelihoods, and
# hence gradients/Hessians, then differ -- this is a model-construction
# mismatch, not a derivative-math one). 0/180 families are symmetric under the
# argument swap that distinguishes the conventions, so the two libraries build
# the *same* model and an element-wise comparison is meaningful. Rotated
# families are covered elsewhere: at the bicop level against BiCopDeriv* (all
# rotations; test_bicop_derivatives.R) and at the vine level against
# brute-force finite differences of the whole-vine log-likelihood (mixed
# 0/90/180/270; VinecopDerivatives.hessian_matches_brute_force).
#
# Both the gradient (RVineGrad) and the joint Hessian (RVineHessian) are
# reordered into vinecopulib's parameter order (trees, edges, params; pair
# copula (t, e) <-> VineCopula matrix cell [d - t, e + 1]). The 90/270 sign
# convention (vinecopulib keeps positive natural parameters while VineCopula
# negates the parameter for the rotated codes, so a gradient entry flips sign
# once per parameter differentiation and a Hessian entry (i, j) by s_i * s_j)
# is kept in the reordering below for completeness but is inert here (no
# rotated codes). Ordering (from the VineCopula sources): RVineGrad returns
# the raw C vector reversed (RVineGrad.R: grad1[dd:1], grad2[tt:1]), which is
# the cell order used below; RVineHessian returns the raw C matrix unreversed
# (RVineHessian.R), so it is RVineGrad's order with the theta block and the
# par2 block each reversed -- undone here by the index c(rev(1:dd), dd +
# rev(1:tt)). The C++ side compares element-wise (analytic vs analytic), which
# also confirms this ordering.
#
# The model here must stay in sync with the C++ side (test_vinecop_class.hpp,
# VinecopDerivatives tests).
args <- commandArgs(trailingOnly = TRUE)

# Output goes to the directory the test allocated for this run, so that
# concurrent test processes cannot overwrite each other's files.
outdir <- Sys.getenv("VINECOPULIB_TEST_OUTDIR", unset = ".")
out <- function(f) file.path(outdir, f)

n <- as.numeric(args[1])

need <- "2.6.2"
have <- tryCatch(as.character(packageVersion("VineCopula")),
                 error = function(e) "0")
if (utils::compareVersion(have, need) < 0) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "http://cran.rstudio.com/")
  }
  remotes::install_github("tnagler/VineCopula", upgrade = "never")
}
suppressMessages(library(VineCopula))

d <- 5
Matrix <- D2RVine(1:d, family = rep(0, d * (d - 1) / 2),
                  par = rep(0, d * (d - 1) / 2))$Matrix
fam <- matrix(0, d, d)
par <- matrix(0, d, d)
par2 <- matrix(0, d, d)
# tree k lives in row d - k + 1. Only 0/180 families (see the header note on
# the structure / swapped-rotation convention mismatch); a mix of one- and
# two-parameter families across all trees, including a Student edge feeding a
# deeper tree so the par2 cascade is exercised.
fam[5, 1] <- 1; par[5, 1] <- 0.6 # gaussian
fam[5, 2] <- 3; par[5, 2] <- 2.5 # clayton
fam[5, 3] <- 14; par[5, 3] <- 1.8 # gumbel 180
fam[5, 4] <- 6; par[5, 4] <- 2.2 # joe
fam[4, 1] <- 5; par[4, 1] <- 4.0 # frank
fam[4, 2] <- 2; par[4, 2] <- 0.5; par2[4, 2] <- 6 # student
fam[4, 3] <- 13; par[4, 3] <- 1.2 # clayton 180
fam[3, 1] <- 4; par[3, 1] <- 1.4 # gumbel
fam[3, 2] <- 1; par[3, 2] <- 0.2 # gaussian
fam[2, 1] <- 3; par[2, 1] <- 2.5 # clayton
RVM <- RVineMatrix(Matrix, family = fam, par = par, par2 = par2)

set.seed(0)
u <- matrix(0.01 + 0.98 * runif(n * d), n, d)

# RVineGrad's output order: sub-diagonal cells column by column from the
# left, top to bottom; t-copula par2 entries appended in the same cell order.
cells <- NULL
for (cc in 1:(d - 1)) {
  for (rr in (cc + 1):d) {
    cells <- rbind(cells, c(rr, cc))
  }
}
grad_index <- function(r, c) which(cells[, 1] == r & cells[, 2] == c)
dd <- nrow(cells)
par2_cells <- which(fam[cells] == 2)
rot90 <- c(23, 24, 26, 33, 34, 36)

# maps a length-(dd + #par2) vector in RVineGrad order to vinecopulib order
# (trees, edges, params), with the 90/270 parameter sign flip
reorder_grad <- function(g) {
  out <- NULL
  for (t in 0:(d - 2)) {
    for (e in 0:(d - 2 - t)) {
      r <- d - t
      c <- e + 1
      s <- if (fam[r, c] %in% rot90) -1 else 1
      out <- c(out, s * g[grad_index(r, c)])
      if (fam[r, c] == 2) {
        out <- c(out, g[dd + which(par2_cells == grad_index(r, c))])
      }
    }
  }
  out
}

grad <- reorder_grad(RVineGrad(u, RVM)$gradient)
grad1 <- reorder_grad(RVineGrad(u[1, , drop = FALSE], RVM)$gradient)

# RVineHessian in vinecopulib order: undo the raw-C-order reversal, then apply
# the same cell -> vinecopulib mapping (with signs) as the gradient
npars <- length(grad)
to_grad_order <- c(rev(seq_len(dd)),
                   if (length(par2_cells) > 0)
                     dd + rev(seq_along(par2_cells)) else integer(0))
Hgo <- RVineHessian(u, RVM)$hessian[to_grad_order, to_grad_order]
perm <- reorder_grad(seq_len(dd + length(par2_cells)))
signs <- sign(perm)
perm <- abs(perm)
Hout <- matrix(0, npars, npars)
for (i in seq_len(npars)) {
  for (j in seq_len(npars)) {
    Hout[i, j] <- signs[i] * signs[j] * Hgo[perm[i], perm[j]]
  }
}

write.table(u, file = out("temp_vderiv_data"), col.names = FALSE, row.names = FALSE)
write.table(RVM$Matrix, file = out("temp_vderiv_matrix"), col.names = FALSE,
            row.names = FALSE)
write.table(cbind(grad, grad1), file = out("temp_vderiv_grad"), col.names = FALSE,
            row.names = FALSE)
write.table(Hout, file = out("temp_vderiv_hess"), col.names = FALSE,
            row.names = FALSE)
