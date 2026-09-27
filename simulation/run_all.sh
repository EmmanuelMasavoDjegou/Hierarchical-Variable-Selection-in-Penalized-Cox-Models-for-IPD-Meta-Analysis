#!/bin/sh
# Runs the simulation of Section 5 (R = 25 replicates per scenario). Resumable.
cd "$(dirname "$0")"
Rscript 01_run_simulation.R 25 >> results/run_log.txt 2>&1
