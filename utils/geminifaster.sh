#!/bin/bash
#SBATCH -p fast
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -o /home/umr8227/gv/dgoudenege/log/gemini.%N.%j.out
#SBATCH -e /home/umr8227/gv/dgoudenege/log/gemini.%N.%j.err
# Load required modules
module load python/3.9
module load hmmer/3.2.1
module load diamond/2.0.9
module load trimmomatic/0.39
module load spades/3.15.2
module load blast/2.9.0
module load interproscan/5.51-85.0
module load trnascan-se/2.0
module load emboss/6.6.0
module load eggnog-mapper/2.1.3
module load muscle/3.8.1551
module load fasttree/2.1.10
module load imagemagick/7.1.0_5
module load entrez-direct/11.0
module load prodigal/2.6.3
module load ppanggolin/1.0.1
module load circos/0.69.6
module load r/3.6.3
module load snippy/4.3.3
module load mafft/7.407
# Launch with script to keep ANSI color output
export isslurm="true"
export geminicpu="4"
export geminimem="8"
export PYTHONPATH=/home/umr8227/gv/dgoudenege/script/python_lib:/home/umr8227/gv/dgoudenege/.local/lib/python3.9/site-packages/
# :/shared/software/miniconda/envs/python-pytorch-tensorflow-3.9-1.11.0-2.6.2/lib/python3.9/site-packages/
# script --flush --quiet --return /dev/null --command "python3.9 /home/umr8227/gv/dgoudenege/script/gemini/gemini.py $@"
script --flush --command "python3.9 /home/umr8227/gv/dgoudenege/script/gemini/gemini.py $@"
# Example
# geminifaster "phage_annotation -i /shared/projects/dynamic/input/phages_FNA -o /shared/projects/dynamic/finalresult/gemini_phage_annotation"