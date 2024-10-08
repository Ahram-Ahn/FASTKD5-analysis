#!/bin/bash
#BSUB -J dreem_test
#BSUB -P dreem
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -W 320:00
#BSUB -q bigmem
#BSUB -n 16
#BSUB -R "rusage[mem=15000]"
#BSUB -u axa3528@miami.edu
#BSUB -N axa3528@miami.edu

conda init bash
. ~/.bashrc
conda activate seismic

seismic -v wf --force /scratch/projects/dreem/michele/sequences/reference/human_mt_unprocessed_version2.fasta \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/WT1_R1.fq.gz \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/WT1_R2.fq.gz \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/WT2_R1.fq.gz \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/WT2_R2.fq.gz \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/fastkdko1_R1.fq.gz \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/fastkdko1_R2.fq.gz \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/fastkdko2_R1.fq.gz \
 -x /scratch/projects/dreem/michele/sequences/FASTKD5/fastkdko2_R2.fq.gz \
 --max-fmut-read 0.15 --sep-strands -k 2 --fold -q 0.95 --fold-md 350 -w 100 --svg --pdf --png --graph-giniroll --graph-aucroll
