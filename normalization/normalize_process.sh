#!/bin/bash


script_dir=pwd
solr_add=
ab3p_dir=
main_dir=

temp_caption=$main_dir/temp_process_caption
mkdir $temp_caption

caption_dir=$main_dir/caption
fulltext_dir=$main_dir/fulltext
output_dir=$main_dir/normalization

mkdir $output_dir

ori_file=$temp_caption/caption.txt.ori

Ab3P=$temp_caption/Ab3P
mkdir $Ab3P


cd $script_dir
python process_solr.py -f $script_dir/solr/final_gene.tsv.gz -s gene -d True -a $solr_add
python process_solr.py -f $script_dir/solr/protein_no_gene.tsv.gz -s protein -d False -a $solr_add
python write_BioC.py -x $caption_dir/ -o $output_dir/ -a $ori_file -t $temp_caption/documents/ -d $temp_caption/ -f $fulltext_dir/ -p $ab3p_dir/ -s $solr_add



