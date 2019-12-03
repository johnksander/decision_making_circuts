#!/usr/bin/env bash

#meant for submitting nets sims, make it a lil easier
#submits a job file & its guest-compute counterpart

S=$1

sbatch --array=1-500 $S
sbatch --array=1-1000 guest_$S
