################################################################################
##  scripts/audit_consistency.R
##  Code <-> manuscript consistency checks. Run from the repository root:
##      Rscript scripts/audit_consistency.R                      # code side only
##      Rscript scripts/audit_consistency.R path/to/main.tex     # + manuscript
##  The manuscript is not part of this repository; point the script at your
##  local copy (e.g. downloaded from Overleaf, with references.bib next to it).
##  Exit status is non-zero if any check fails.
################################################################################
args <- commandArgs(trailingOnly = TRUE)
ok <- TRUE
check <- function(name, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", name))
  if (!isTRUE(cond)) ok <<- FALSE
}
setup <- new.env(); sys.source("simulation/00_setup.R", envir = setup)
R_PAPER <- 50

## ---- code side -----------------------------------------------------------------
tabs <- c("tab_scenarios", "tab_sim_base", "tab_sim_coverage", "tab_sim_all",
          "tab_oracle_selection", "tab_oracle_inference", "tab_nearcancel", "tab_cohorts",
          "tab_app_methods", "tab_app_genes", "tab_app_frailty", "tab_app_loso")
figs <- c("sim_mcc", "sim_fdr_tpr", "sim_mse_theta", "sim_deviations",
          "oracle_selection", "oracle_qq",
          "app_km", "app_cv_surface", "app_forest", "app_stability")
for (t in tabs) check(paste0("output/tables/", t, ".tex exists"), file.exists(file.path("output/tables", paste0(t, ".tex"))))
for (f in figs) check(paste0("output/figures/", f, ".png exists"), file.exists(file.path("output/figures", paste0(f, ".png"))))
raw <- list.files("simulation/results/raw", "^s[0-9]{2}_r[0-9]{3}\\.csv$", full.names = TRUE)
check(sprintf("%d replicate result files (%d expected)", length(raw), R_PAPER * length(setup$SCENARIOS)),
      length(raw) == R_PAPER * length(setup$SCENARIOS))
if (length(raw)) {
  X <- do.call(rbind, lapply(raw, read.csv))
  check(sprintf("every replicate has all %d methods", length(setup$METHODS)),
        all(tapply(X$method, paste(X$scenario, X$rep), function(m) length(unique(m))) == length(setup$METHODS)))
  reps <- tapply(X$rep, X$scenario, function(r) length(unique(r)))
  check(sprintf("replicates per scenario: min %d, max %d (R = %d required)", min(reps), max(reps), R_PAPER),
        all(reps == R_PAPER))
  check("no hierarchy violations in any fit", sum(X$hier_violation) == 0)
}
orc <- list.files("simulation/results/oracle", "^n[0-9]{4}_r[0-9]{3}\\.csv$", full.names = TRUE)
check(sprintf("%d oracle replicate files (%d expected: 4 sizes x 200)", length(orc), 800),
      length(orc) == 800)
ncf <- list.files("simulation/results/nearcancel", "^r[0-9]{3}\\.csv$", full.names = TRUE)
check(sprintf("%d near-cancellation replicate files (100 expected)", length(ncf)), length(ncf) == 100)
if (length(raw)) {
  st <- list.files("output/tables", "^tab_sim", full.names = TRUE)
  check("simulation tables are newer than the raw results (re-run 02_make_tables_figures.R)",
        length(st) > 0 && max(file.info(raw)$mtime) <= min(file.info(st)$mtime))
}

## ---- manuscript side ---------------------------------------------------------------
if (length(args)) {
  texf <- args[1]; dir <- dirname(texf)
  tex <- paste(readLines(texf, warn = FALSE), collapse = "\n")
  grab <- function(pattern) regmatches(tex, gregexpr(pattern, tex, perl = TRUE))[[1]]
  used <- sub(".*\\{figures/(.*)\\}", "\\1", grab("\\\\includegraphics(\\[[^]]*\\])?\\{figures/[^}]+\\}"))
  for (f in used) check(paste("figure", f, "is produced in output/figures"), file.exists(file.path("output/figures", f)))
  for (f in figs) check(paste("figure", f, "is used in the manuscript"), paste0(f, ".png") %in% used)
  T <- setup$TUNE; B <- setup$BASE
  check("V = 5 folds", grepl("recommend \\$V=5\\$ folds", tex) && T$nfolds == 5)
  rg <- paste(format(T$ratio, trim = TRUE, drop0trailing = TRUE), collapse = ",")
  check(sprintf("ratio grid {%s} stated as in 00_setup.R", rg),
        grepl(paste0("r\\in\\{", rg, "\\}"), tex, fixed = TRUE))
  check("M = 20 lambda values", grepl("\\$M=20\\$", tex) && T$nlambda == 20)
  check("kappa = 0.05", grepl("\\\\kappa=0.05", tex) && T$lambda.min.ratio == 0.05)
  check("d_max = 100", grepl("d_\\{\\\\max\\}=100", tex) && T$dfmax == 100)
  check("base design K = 5, n_k = 200, p = 100", B$p == 100 && B$K == 5 && B$n == 200 &&
          grepl("\\$K=5\\$ studies of \\$n_k=200\\$", tex) && grepl("\\$p=100\\$ covariates", tex))
  check(sprintf("R = %d stated in the manuscript", R_PAPER), grepl(sprintf("\\$R=%d\\$ times", R_PAPER), tex))
  bibf <- file.path(dir, "references.bib")
  if (file.exists(bibf)) {
    bib <- paste(readLines(bibf, warn = FALSE), collapse = "\n")
    keys <- unique(unlist(strsplit(gsub("\\s", "", sub(".*\\{(.*)\\}", "\\1",
              grab("\\\\cite[pt]?\\*?(\\[[^]]*\\])?\\{[^}]+\\}"))), ",")))
    bibkeys <- sub("@\\w+\\{", "", regmatches(bib, gregexpr("@\\w+\\{[^,]+", bib))[[1]])
    miss <- setdiff(keys, bibkeys)
    check(sprintf("all %d citation keys resolve%s", length(keys),
                  if (length(miss)) paste0(" (missing: ", paste(miss, collapse = ", "), ")") else ""),
          length(miss) == 0)
  }
  ntodo <- length(grab("TODO|FIXME|XXX|[Pp]laceholder"))
  check(sprintf("no TODO comments left (%d remaining)", ntodo), ntodo == 0)
  ntbd <- length(grab("\\\\TBD"))
  check(sprintf("no TBD placeholders left (%d remaining)", ntbd), ntbd == 0)
}
cat(if (ok) "\nAll checks passed.\n" else "\nSome checks FAILED (see above).\n")
quit(status = if (ok) 0 else 1)
