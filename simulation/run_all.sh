#!/bin/sh
# Runs the simulation of Section 5 (R = 50 replicates per scenario). Resumable:
# replicates already in results/ are skipped, so re-issue after any interruption.
cd "$(dirname "$0")"
Rscript 01_run_simulation.R 50 >> results/run_log.txt 2>&1
