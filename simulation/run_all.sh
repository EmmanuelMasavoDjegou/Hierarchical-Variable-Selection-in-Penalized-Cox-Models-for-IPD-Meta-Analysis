#!/bin/sh
# Runs both simulation experiments of Section 5. Resumable; set NCORES to run
# replicates in parallel, e.g.  NCORES=8 ./run_all.sh
cd "$(dirname "$0")"
Rscript 01_run_simulation.R 50 >> results/run_log.txt 2>&1     # Sections 5.1-5.4
Rscript 03_oracle_study.R 200 >> results/run_log.txt 2>&1      # Section 5.5
