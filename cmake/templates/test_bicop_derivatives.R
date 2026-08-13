#!/usr/bin/env Rscript
# Oracle for the analytic-derivative parity tests: evaluates VineCopula's
# BiCopDeriv/BiCopDeriv2/BiCopHfuncDeriv/BiCopHfuncDeriv2 on a fixed grid.
# `family` is the rotation-encoded VineCopula code and `par` is already
# negated for the 90/270 codes (see ParBicopTest). Columns that do not apply
# (second-parameter selectors for one-parameter families) are written as 0
# so that the file stays numeric (tools_eigen::read_matxd cannot parse NA).
args <- commandArgs(trailingOnly = TRUE)

# Output goes to the directory the test allocated for this run, so that
# concurrent test processes cannot overwrite each other's files.
outdir <- Sys.getenv("VINECOPULIB_TEST_OUTDIR", unset = ".")
out <- function(f) file.path(outdir, f)


n <- as.numeric(args[1])
family <- as.numeric(args[2])
par <- as.numeric(args[3])
par2 <- as.numeric(args[4])

# VineCopula >= 2.6.2 fixes the 90/270 second parameter derivatives
# (BiCopDeriv2 / BiCopHfuncDeriv2; tnagler/VineCopula #102).
have <- tryCatch(as.character(packageVersion("VineCopula")),
                 error = function(e) "0")
if (utils::compareVersion(have, "2.6.2") < 0) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "http://cran.rstudio.com/")
  }
  remotes::install_github("tnagler/VineCopula", upgrade = "never")
}

set.seed(0)
# interior points: derivative magnitudes explode near the boundary
u1 <- 0.01 + 0.98 * runif(n)
u2 <- 0.01 + 0.98 * runif(n)

st <- (family == 2) # par2 selectors exist for the t-copula only
zero <- rep(0, n)

# oracles for the first h-function h1(u2|u1): VineCopula's BiCopHfuncDeriv
# differentiates h2(u1|u2), and h1 of a 90-degree code equals h2 of the
# corresponding 270-degree code at swapped arguments (and vice versa)
dec <- family %/% 10
fam_swap <- if (dec == 2) family + 10 else if (dec == 3) family - 10 else family

d1 <- function(deriv, log = FALSE) {
  VineCopula::BiCopDeriv(u1, u2, family, par, par2, deriv = deriv, log = log)
}
d2 <- function(deriv) {
  VineCopula::BiCopDeriv2(u1, u2, family, par, par2, deriv = deriv)
}
h2d1 <- function(deriv) {
  VineCopula::BiCopHfuncDeriv(u1, u2, family, par, par2, deriv = deriv)
}
h2d2 <- function(deriv) {
  VineCopula::BiCopHfuncDeriv2(u1, u2, family, par, par2, deriv = deriv)
}
h1d1 <- function(deriv) {
  VineCopula::BiCopHfuncDeriv(u2, u1, fam_swap, par, par2, deriv = deriv)
}
h1d2 <- function(deriv) {
  VineCopula::BiCopHfuncDeriv2(u2, u1, fam_swap, par, par2, deriv = deriv)
}

results <- cbind(
  u1, u2, # 1:2
  d1("par"), # 3
  if (st) d1("par2") else zero, # 4
  d1("u1"), # 5
  d1("u2"), # 6
  d1("par", log = TRUE), # 7
  if (st) d1("par2", log = TRUE) else zero, # 8
  d2("par"), # 9
  if (st) d2("par2") else zero, # 10
  d2("u1"), # 11
  d2("u2"), # 12
  if (st) d2("par1par2") else zero, # 13
  d2("par1u1"), # 14
  if (st) d2("par2u1") else zero, # 15
  d2("par1u2"), # 16
  if (st) d2("par2u2") else zero, # 17
  h2d1("par"), # 18
  if (st) h2d1("par2") else zero, # 19
  h2d1("u2"), # 20
  h2d2("par"), # 21
  if (st) h2d2("par2") else zero, # 22
  h2d2("u2"), # 23
  if (st) h2d2("par1par2") else zero, # 24
  h2d2("par1u2"), # 25
  if (st) h2d2("par2u2") else zero, # 26
  h1d1("par"), # 27
  if (st) h1d1("par2") else zero, # 28
  h1d1("u2"), # 29 (the swapped call's "u2" is our u1)
  h1d2("par"), # 30
  if (st) h1d2("par2") else zero, # 31
  h1d2("u2"), # 32
  if (st) h1d2("par1par2") else zero, # 33
  h1d2("par1u2"), # 34
  if (st) h1d2("par2u2") else zero # 35
)

write.table(results, file = out("temp_deriv"), col.names = FALSE, row.names = FALSE)
