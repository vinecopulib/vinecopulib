#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Output goes to the directory the test allocated for this run, so that
# concurrent test processes cannot overwrite each other's files.
outdir <- Sys.getenv("VINECOPULIB_TEST_OUTDIR", unset = ".")
out <- function(f) file.path(outdir, f)


n <- as.numeric(args[1])
family <- as.numeric(args[2])
par <- as.numeric(args[3])
par2 <- as.numeric(args[4])

if (!("VineCopula" %in% rownames(installed.packages()))) {
  install.packages("VineCopula", repos = "http://cran.rstudio.com/")
}

set.seed(0)
u1 <- runif(n)
u2 <- runif(n)
taildep <- VineCopula::BiCopPar2TailDep(family, par, par2)
results <- cbind(
  VineCopula::BiCopPar2Tau(family, par, par2), u1, u2,
  VineCopula::BiCopPDF(u1, u2, family, par, par2),
  VineCopula::BiCopCDF(u1, u2, family, par, par2),
  VineCopula::BiCopHfunc1(u1, u2, family, par, par2),
  VineCopula::BiCopHfunc2(u1, u2, family, par, par2),
  VineCopula::BiCopHinv1(u1, u2, family, par, par2),
  VineCopula::BiCopHinv2(u1, u2, family, par, par2),
  VineCopula::BiCopPar2Beta(family, par, par2),
  taildep$lower, taildep$upper
)

write.table(results, file = out("temp"), col.names = FALSE, row.names = FALSE)
