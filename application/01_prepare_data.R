################################################################################
##  application/01_prepare_data.R
##  Builds the analysis data set of Section 6 from curatedOvarianData.
##
##  Input   : six ExpressionSets from the Bioconductor package curatedOvarianData
##            (Ganzfried et al., 2013), or the six *_eset.rda files in
##            application/curatedOvarianData/data/ or in the folder COD_DATA_DIR
##            (see the README, "Real data", for how to obtain them).
##  Steps   : overall survival (days -> years); subjects with positive follow-up;
##            genes measured in all six cohorts; outcome-blind filter keeping the
##            p = 500 genes with the largest average within-cohort IQR rank.
##  Output  : results/analysis_data.rds, results/table_cohorts.csv  (Table 7)
################################################################################
suppressPackageStartupMessages({ library(Biobase); library(survival) })
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {   # folder of this script
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))     # works under Rscript
  if (length(f)) dirname(normalizePath(f[1])) else "." })                          # from any directory
out <- file.path(here, "results"); dir.create(out, showWarnings = FALSE)
P_KEEP <- 500

cohorts <- c(A = "TCGA_eset", B = "GSE32062.GPL6480_eset", C = "GSE9891_eset",
             D = "GSE26712_eset", E = "GSE17260_eset", F = "GSE30161_eset")
## Data source, in this order: (1) COD_DATA_DIR, a folder of *_eset.rda files;
## (2) application/curatedOvarianData/data/ if present (see the README);
## (3) the installed Bioconductor package (recent versions download on first use).
cod_dir <- Sys.getenv("COD_DATA_DIR", file.path(here, "curatedOvarianData", "data"))
load_eset <- function(nm) {
  f <- file.path(cod_dir, paste0(nm, ".rda"))
  if (file.exists(f)) { e <- new.env(); load(f, envir = e); return(e[[nm]]) }
  if (!requireNamespace("curatedOvarianData", quietly = TRUE))
    stop("curatedOvarianData is not installed and ", f, " does not exist; see the README")
  if ("curatedOvarianData" %in% getNamespaceExports("curatedOvarianData"))
    return(curatedOvarianData::curatedOvarianData(nm))           # package >= 1.51
  e <- new.env(); utils::data(list = nm, package = "curatedOvarianData", envir = e); e[[nm]]
}

studies <- lapply(names(cohorts), function(s) {
  es <- load_eset(cohorts[[s]]); pd <- pData(es)
  time <- as.numeric(pd$days_to_death) / 365.25
  status <- as.integer(pd$vital_status == "deceased")
  keep <- !is.na(time) & !is.na(status) & time > 0
  ex <- t(exprs(es)[, keep, drop = FALSE])
  ex <- ex[, colSums(is.na(ex)) == 0, drop = FALSE]
  list(study = s, dataset = sub("_eset$", "", cohorts[[s]]), time = time[keep],
       status = status[keep], expr = ex)
})
common <- Reduce(intersect, lapply(studies, function(z) colnames(z$expr)))
rk <- sapply(studies, function(z) rank(-apply(z$expr[, common], 2, IQR)))
genes <- common[order(rowMeans(rk))][seq_len(P_KEEP)]

X <- do.call(rbind, lapply(studies, function(z) z$expr[, genes]))
dat <- data.frame(study = rep(names(cohorts), sapply(studies, function(z) length(z$time))),
                  time = unlist(lapply(studies, `[[`, "time")),
                  status = unlist(lapply(studies, `[[`, "status")))
saveRDS(list(dat = dat, X = X, genes = genes, n_common = length(common)),
        file.path(out, "analysis_data.rds"))

tab <- do.call(rbind, lapply(studies, function(z) {
  km <- survfit(Surv(z$time, z$status) ~ 1)
  data.frame(Study = z$study, Dataset = z$dataset, N = length(z$time),
             Events = sum(z$status), Censoring = round(100 * mean(z$status == 0)),
             MedianOS = round(unname(summary(km)$table["median"]), 1),
             MedianFU = round(unname(summary(survfit(Surv(z$time, 1 - z$status) ~ 1))$table["median"]), 1))
}))
write.csv(tab, file.path(out, "table_cohorts.csv"), row.names = FALSE)
print(tab); cat("common genes:", length(common), "| kept p =", P_KEEP, "\n")
