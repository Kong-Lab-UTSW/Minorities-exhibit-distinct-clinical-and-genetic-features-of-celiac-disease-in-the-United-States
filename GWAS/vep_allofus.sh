#!/bin/bash

source ~/.bashrc

conda activate WES

python vep_prep.py $1


conda activate VEP

ANNO=$(echo $1 | sed 's/.csv//' )

vep --cache \
 --dir_cache $workspace/opt/vep \
 --dir_plugins $workspace/opt/vep/Plugins \
 --fasta $workspace/opt/vep/homo_sapiens/108_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
 --fork 8 \
 --assembly=GRCh38 \
 --offline \
 --tab \
 --plugin Downstream \
 --plugin LoF,loftee_path:$workspace/opt/vep/Plugins/loftee \
 --plugin CADD,$workspace/opt/vep/Plugins/whole_genome_SNVs.tsv.gz,$workspace/opt/vep/Plugins/gnomad.genomes.r3.0.indel.tsv.gz \
 --everything \
 --terms SO \
 --pick \
 --transcript_version \
 -i ${ANNO}.tsv \
 -o ${ANNO}.vep.tsv \
 --force_overwrite

conda deactivate 

python vep_merge.py ${ANNO}.tsv ${ANNO}.vep.tsv
