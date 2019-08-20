#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --account=rrg-chauvec
#SBATCH --output /home/chauvec/wg-anoph/DIVERSITY_LABRI/SVN/exp/04-06-2019-savage-vgflow-kallisto/traces/run_kallisto_1.trace
#SBATCH --error  /home/chauvec/wg-anoph/DIVERSITY_LABRI/SVN/exp/04-06-2019-savage-vgflow-kallisto/traces/run_kallisto_1.error

#source activate py37
module load kallisto/0.44.0
cd /home/chauvec/wg-anoph/DIVERSITY_LABRI/SVN/exp/04-06-2019-savage-vgflow-kallisto/

mkdir -p results/run1_1_kallisto_HS20
mkdir -p results/run1_1_kallisto_HS25
kallisto quant -i data/ST150.idx -t 8 -o results/run1_1_kallisto_HS20 -b 100  results/run1_1_HS20_1.fq results/run1_1_HS20_2.fq
kallisto quant -i data/ST150.idx -t 8 -o results/run1_1_kallisto_HS25 -b 100  results/run1_1_HS25_1.fq results/run1_1_HS25_2.fq

mkdir -p results/run1_2_kallisto_HS20
mkdir -p results/run1_2_kallisto_HS25
kallisto quant -i data/ST150.idx -t 8 -o results/run1_2_kallisto_HS20 -b 100  results/run1_2_HS20_1.fq results/run1_2_HS20_2.fq
kallisto quant -i data/ST150.idx -t 8 -o results/run1_2_kallisto_HS25 -b 100  results/run1_2_HS25_1.fq results/run1_2_HS25_2.fq

mkdir -p results/run1_3_kallisto_HS20
mkdir -p results/run1_3_kallisto_HS25
kallisto quant -i data/ST150.idx -t 8 -o results/run1_3_kallisto_HS20 -b 100  results/run1_3_HS20_1.fq results/run1_3_HS20_2.fq
kallisto quant -i data/ST150.idx -t 8 -o results/run1_3_kallisto_HS25 -b 100  results/run1_3_HS25_1.fq results/run1_3_HS25_2.fq

mkdir -p results/run1_4_kallisto_HS20
mkdir -p results/run1_4_kallisto_HS25
kallisto quant -i data/ST150.idx -t 8 -o results/run1_4_kallisto_HS20 -b 100  results/run1_4_HS20_1.fq results/run1_4_HS20_2.fq
kallisto quant -i data/ST150.idx -t 8 -o results/run1_4_kallisto_HS25 -b 100  results/run1_4_HS25_1.fq results/run1_4_HS25_2.fq

mkdir -p results/run1_5_kallisto_HS20
mkdir -p results/run1_5_kallisto_HS25
kallisto quant -i data/ST150.idx -t 8 -o results/run1_5_kallisto_HS20 -b 100  results/run1_5_HS20_1.fq results/run1_5_HS20_2.fq
kallisto quant -i data/ST150.idx -t 8 -o results/run1_5_kallisto_HS25 -b 100  results/run1_5_HS25_1.fq results/run1_5_HS25_2.fq
