################################################################################
##  scripts/audit_consistency.R
##  Mechanical code <-> manuscript consistency checks. Run from the repo root:
##      Rscript scripts/audit_consistency.R
##  Exit status is non-zero if any check fails.
################################################################################
tex <- paste(readLines("manuscript/main.tex", warn = FALSE), collapse = "\n")
bib <- paste(readLines("manuscript/references.bib", warn = FALSE), collapse = "\n")
setup <- new.env(); sys.source("simulation/00_setup.R", envir = setup)
ok <- TRUE
check <- function(name, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", name))
  if (!isTRUE(cond)) ok <<- FALSE
}
grab <- function(pattern) regmatches(tex, gregexpr(pattern, tex, perl = TRUE))[[1]]

## 1. every \input table and \includegraphics figure exists
inp <- sub(".*\\{(.*)\\}", "\\1", grab("\\\\input\\{[^}]+\\}"))
for (f in inp) check(paste("table exists:", f), file.exists(file.path("manuscript", f)))
fig <- sub(".*\\{(.*)\\}", "\\1", grab("\\\\includegraphics(\\[[^]]*\\])?\\{[^}]+\\}"))
for (f in fig) check(paste("figure exists:", f), file.exists(file.path("manuscript", f)))

## 2. every generated output is used by the manuscript (no orphan results)
for (f in list.files("manuscript/tables")) check(paste("table used:", f), any(grepl(f, inp, fixed = TRUE)))
for (f in list.files("manuscript/figures")) check(paste("figure used:", f), any(grepl(f, fig, fixed = TRUE)))

## 3. tuning constants stated in the text equal those in the code
T <- setup$TUNE
check("V = 5 folds", grepl("recommend \\$V=5\\$ folds", tex) && T$nfolds == 5)
check("ratio grid {0.5,1,2}", grepl("r\\\\in\\\\\\{0.5,1,2\\\\\\}", tex) && identical(T$ratio, c(0.5, 1, 2)))
check("M = 20 lambda values", grepl("\\$M=20\\$", tex) && T$nlambda == 20)
check("kappa = 0.05", grepl("\\\\kappa=0.05", tex) && T$lambda.min.ratio == 0.05)
check("d_max = 100", grepl("d_\\{\\\\max\\}=100", tex) && T$dfmax == 100)
B <- setup$BASE
check("base p = 100, K = 5, n_k = 200", B$p == 100 && B$K == 5 && B$n == 200 &&
        grepl("\\$K=5\\$ studies of \\$n_k=200\\$", tex) && grepl("\\$p=100\\$ covariates", tex))
check("13 scenarios", length(setup$SCENARIOS) == 13)
check("replicates R = 25 stated", grepl("\\$R=25\\$ times", tex))
raw <- list.files("simulation/results", "^sim_raw_.*csv$", full.names = TRUE)
if (length(raw)) {
  reps <- sapply(raw, function(f) length(unique(read.csv(f)$rep)))
  check(sprintf("replicates on disk (min %d) equal R = 25", min(reps)), all(reps == 25))
}
## 4. citations resolve
keys <- unique(unlist(strsplit(gsub("\\s", "", sub(".*\\{(.*)\\}", "\\1",
          grab("\\\\cite[pt]?\\*?(\\[[^]]*\\])?\\{[^}]+\\}"))), ",")))
bibkeys <- sub("@\\w+\\{([^,]+),.*", "\\1", grab("@\\w+\\{[^,]+,"), perl = TRUE)
bibkeys <- sub("@\\w+\\{", "", regmatches(bib, gregexpr("@\\w+\\{[^,]+", bib))[[1]])
miss <- setdiff(keys, bibkeys)
check(sprintf("all %d citation keys resolve%s", length(keys),
              if (length(miss)) paste0(" (missing: ", paste(miss, collapse = ", "), ")") else ""),
      length(miss) == 0)
cat(if (ok) "\nAll checks passed.\n" else "\nSome checks FAILED.\n")
quit(status = if (ok) 0 else 1)
