#!/bin/bash
## SLURM REQUIRED SETTINGS  <--- two hashtags are a comment in Slurm
#SBATCH --partition=dauhajre_nodes
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1

## SLURM reads %x as the job name and %j as the job ID
#SBATCH --job-name=idealshelf
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#Job to run
#mpiexec -np 8 ./roms roms.in
ncjoin -d shelf_avg*
ncjoin -d shelf_grd*
zslice $(seq -1 -1 -500) shelf_grd.nc shelf_avg*




