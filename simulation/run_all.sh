#!/bin/sh
# Runs the three simulation experiments of Section 5 and then builds their tables
# and figures (Tables 1-6, A1 and Figures 1-6 in ../output/). Resumable: re-run
# after an interruption and completed replicates are skipped.
#
#   sh simulation/run_all.sh                 # from the repository root, one core
#   NCORES=8 sh simulation/run_all.sh        # replicates in parallel (Linux/macOS)
#
# Single-core run times: about 4 h (01), 15 h (03) and 7 min (05); divide by NCORES.
# Progress is appended to simulation/results/run_log.txt.
set -e
cd "$(dirname "$0")"
mkdir -p results
Rscript 01_run_simulation.R 50       >> results/run_log.txt 2>&1   # Sections 5.1-5.4, R = 50
Rscript 03_oracle_study.R 200        >> results/run_log.txt 2>&1   # Section 5.5, oracle, R = 200
Rscript 05_near_cancellation.R 100   >> results/run_log.txt 2>&1   # Section 5.5, near-cancellation, R = 100
Rscript 02_make_tables_figures.R     >> results/run_log.txt 2>&1   # Tables 1-3, A1; Figures 1-4
Rscript 04_oracle_tables_figures.R   >> results/run_log.txt 2>&1   # Tables 4-6; Figures 5-6
echo "done: tables in output/tables, figures in output/figures (see simulation/results/run_log.txt)"
