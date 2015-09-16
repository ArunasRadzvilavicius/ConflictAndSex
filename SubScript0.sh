#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=0:10:0 
#$ -l mem=20M
#$ -l tmpfs=2M
#$ -t 1-200
#$ -N EvoS
#$ -wd /home/ucbprad/scratch/FusionFixation
N=100
sigma=1
fname=$(printf "$R-$M-$mu-%.4f.evofix" $k)
echo $fname
./unisex $N $M $mu $k 0 $sigma $R >> "$fname"
rm EvoS*
