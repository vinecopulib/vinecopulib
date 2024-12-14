#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

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
results <- cbind(
  VineCopula::BiCopPar2Tau(family, par, par2), u1, u2,
  VineCopula::BiCopPDF(u1, u2, family, par, par2),
  VineCopula::BiCopCDF(u1, u2, family, par, par2),
  VineCopula::BiCopHfunc1(u1, u2, family, par, par2),
  VineCopula::BiCopHfunc2(u1, u2, family, par, par2),
  VineCopula::BiCopHinv1(u1, u2, family, par, par2),
  VineCopula::BiCopHinv2(u1, u2, family, par, par2)
)

write.table(results, file = "temp", col.names = FALSE, row.names = FALSE)
