#!/bin/sh
# Runs the simulation of Section 5 (R = 50 replicates per scenario). Resumable.
# Set NCORES to run replicates in parallel, e.g.  NCORES=8 ./run_all.sh
cd "$(dirname "$0")"
Rscript 01_run_simulation.R 50 >> results/run_log.txt 2>&1
