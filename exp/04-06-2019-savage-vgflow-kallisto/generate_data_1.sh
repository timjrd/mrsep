#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --account=rrg-chauvec
#SBATCH --output /home/chauvec/wg-anoph/DIVERSITY_LABRI/SVN/exp/04-06-2019-savage-vgflow-kallisto/traces/generate_data_1.trace
#SBATCH --error  /home/chauvec/wg-anoph/DIVERSITY_LABRI/SVN/exp/04-06-2019-savage-vgflow-kallisto/traces/generate_data_1.error

# source activate py37
cd /home/chauvec/wg-anoph/DIVERSITY_LABRI/SVN/exp/04-06-2019-savage-vgflow-kallisto/src
python generate_diversity.py ../data/ST150.fasta 5 3 20 ../results run1_1
python generate_sequence_data.py ../results/run1_1.fasta 100 100 250 HS20 ../results run1_1
python generate_sequence_data.py ../results/run1_1.fasta 100 100 250 HS25 ../results run1_1

python generate_diversity.py ../data/ST150.fasta 5 3 30 ../results run1_2
python generate_sequence_data.py ../results/run1_2.fasta 100 100 250 HS20 ../results run1_2
python generate_sequence_data.py ../results/run1_2.fasta 100 100 250 HS25 ../results run1_2

python generate_diversity.py ../data/ST150.fasta 5 3 40 ../results run1_3
python generate_sequence_data.py ../results/run1_3.fasta 100 100 250 HS20 ../results run1_3
python generate_sequence_data.py ../results/run1_3.fasta 100 100 250 HS25 ../results run1_3

python generate_diversity.py ../data/ST150.fasta 5 3 50 ../results run1_4
python generate_sequence_data.py ../results/run1_4.fasta 100 100 250 HS20 ../results run1_4
python generate_sequence_data.py ../results/run1_4.fasta 100 100 250 HS25 ../results run1_4

python generate_diversity.py ../data/ST150.fasta 5 3 60 ../results run1_5
python generate_sequence_data.py ../results/run1_5.fasta 100 100 250 HS20 ../results run1_5
python generate_sequence_data.py ../results/run1_5.fasta 100 100 250 HS25 ../results run1_5

