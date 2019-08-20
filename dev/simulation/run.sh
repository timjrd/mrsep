#! /usr/bin/env nix-shell
#! nix-shell -i bash
set -e

function unlines {
  for arg in "$@"
  do echo "$arg"
  done
}

if test $# -lt 3 -o $# -gt 10
then
  echo "usage: $0 input.fasta aligned-input.fasta out-dir" \
       "[nb-alleles] [min-copy-nb] [max-copy-nb] [ends-length]" \
       "[take-ref] [seed] [art-flags]"
  exit 1
fi

input=$1
aligned_input=$2
out_dir=$3
nb_alleles=${4:-5}
min_copy_nb=${5:-1}
max_copy_nb=${6:-10}
ends_length=${7:-1000}
take_ref=${8:-"false"}
seed=${9:-$(date +%s)}
art_flags=${10:-"--seqSys HS20 --len 100 --fcov 10"}

art_input=$out_dir/art-input.fasta
art_output=$out_dir/art-output

mkdir -p $out_dir
unlines $input $aligned_input $out_dir \
	$nb_alleles $min_copy_nb $max_copy_nb $ends_length \
	$take_ref $seed "$art_flags" > $out_dir/args.txt
cp $input $out_dir/input.fasta
cp $aligned_input $out_dir/aligned-input.fasta

python src/sample.py $input $nb_alleles $min_copy_nb $max_copy_nb \
       $ends_length $seed > $art_input

$ART_ILLUMINA --in $art_input --out $art_output \
	      --rndSeed $seed $art_flags > /dev/null

python src/input.py $aligned_input $art_output.aln $ends_length $take_ref \
       > $out_dir/output.txt 2> $out_dir/stats.txt

echo $out_dir
