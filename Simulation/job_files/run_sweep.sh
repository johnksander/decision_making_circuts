#!/usr/bin/env bash


sbatch --array=1-110 D2t_very_slow_pref.sh
sbatch --array=111-550 guest_D2t_very_slow_pref.sh

#do like, 20% paul-compute and 80% guest 
