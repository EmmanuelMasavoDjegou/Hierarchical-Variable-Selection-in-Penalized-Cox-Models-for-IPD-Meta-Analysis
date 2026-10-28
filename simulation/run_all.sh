#!/bin/sh
# Runs every simulation experiment of Section 5 and then builds their tables and
# figures in ../output/. Resumable: re-run after an interruption and completed
# replicates are skipped.
#
#   sh simulation/run_all.sh                 # from the repository root, one core
#   NCORES=8 sh simulation/run_all.sh        # replicates in parallel (Linux/macOS)
#
# Single-core run times: about 45 h (01), 15 h (03, Gaussian design), 0.5 h (03,
# bounded design), 15 min (05, both designs) and 2.5 h (07); divide by NCORES.
# Progress is appended to simulation/results/run_log.txt.
set -e
cd "$(dirname "$0")"
mkdir -p results
Rscript 01_run_simulation.R 500                 >> results/run_log.txt 2>&1   # Sections 5.1-5.4, R = 500
DESIGN=bounded  Rscript 03_oracle_study.R 200   >> results/run_log.txt 2>&1   # Section 5.5, bounded design
DESIGN=gaussian Rscript 03_oracle_study.R 200   >> results/run_log.txt 2>&1   # Section 5.5, Gaussian design
DESIGN=bounded  Rscript 05_near_cancellation.R 100 >> results/run_log.txt 2>&1 # Section 5.5, near-cancellation
DESIGN=gaussian Rscript 05_near_cancellation.R 100 >> results/run_log.txt 2>&1
Rscript 07_sensitivity.R 20                     >> results/run_log.txt 2>&1   # Section 5.6, sensitivity
Rscript 02_make_tables_figures.R                >> results/run_log.txt 2>&1   # Tables 1-3, A1; Figures 1-4
DESIGN=bounded  Rscript 04_oracle_tables_figures.R >> results/run_log.txt 2>&1 # Tables 4 (top), 5, 6; Figures 5-6
DESIGN=gaussian Rscript 04_oracle_tables_figures.R >> results/run_log.txt 2>&1 # Table 4 (bottom)
Rscript 08_sensitivity_tables.R                 >> results/run_log.txt 2>&1   # Tables 7-8 (Section 5.6)
echo "done: tables in output/tables, figures in output/figures (see simulation/results/run_log.txt)"
