```
module purge
module load slurm seqtk/1.2 blat/35

ls -1 *.fastq | parallel 'seqtk seq -a {} > $(basename -s .fastq {}).fa'
ls -1 *.fa | parallel 'blat {} flanking.fa $(basename -s .fa {}).psl'
ls -1 *.psl | parallel 'tail -n +6 {} > $(basename -s .psl {})_blat_results.txt'
```

```
module purge
module load slurm R bioconductor

ls -1 *.fastq | sed 's/.fastq//g' | parallel 'sbatch --wrap="Rscript --vanilla parseReads.R {}"'
```
