```
module purge
module load slurm seqtk/1.2 blat/35

ls -1 *.fastq | parallel 'seqtk seq -a {} > $(basename -s .fastq {}).fa'
ls -1 *.fa | parallel 'blat {} target.fa $(basename -s .fa {}).psl'
ls -1 *.psl | parallel 'tail -n +6 {} > $(basename -s .psl {})_blat_results.txt'
```

```
module purge
module load slurm R bioconductor

sbatch --wrap="Rscript --vanilla parseReads.R Sansam_d5s_1_3B1_1D2_Linear_991"
```
