#!/bin/bash


script_dir=pwd
ner_dir=
genia_dir=
main_dir=

fulltext_dir=$main_dir/fulltext
caption_dir=$main_dir/caption

models_dir=$script_dir/models

temp_fulltext=$main_dir/temp_process_fulltext
mkdir $temp_fulltext
mkdir $temp_fulltext/documents

temp_caption=$main_dir/temp_process_caption
mkdir $temp_caption
mkdir $temp_caption/documents
mkdir $temp_caption/sentences
mkdir $temp_caption/standoff

mkdir $main_dir/recognition/



cd $script_dir
python extract_BioC_text.py -i $fulltext_dir/ -o $temp_fulltext/documents/ 

cd $temp_fulltext/documents
ls *.xml > $documents.lst

mkdir $temp_fulltext/sentences
mkdir $temp_fulltext/standoff

for f in `cat $documents.lst`
do
    in_f=$temp_fulltext/documents/$f
    out_f=$temp_fulltext/sentences/$f
    std_f=$temp_fulltext/standoff/$f
    cd $genia_dir
    ./geniass $in_f $out_f.gss
    ruby sentence2standOff.rb  $in_f $out_f.gss $std_f
    perl geniass-postproc.pl $out_f.gss > $out_f.per
done

cd $script_dir
python combine_sentences.py -i $temp_fulltext/sentences/ -f .per -o $temp_fulltext/fulltext.per

cd $ner_dir/src/tokenizer
./nersuite_tokenizer -multidoc [SEP] < $temp_fulltext/fulltext.per > $temp_fulltext/fulltext.tok

cd $ner_dir/src/gtagger
./nersuite_gtagger -multidoc [SEP] -d $ner_dir/models/gtagger/ < $temp_fulltext/fulltext.tok > $temp_fulltext/fulltext.gtag




cd $script_dir
python extract_BioC_text.py -i $caption_dir/ -o $temp_caption/documents/ 

cd $temp_caption/documents
ls *.xml > $documents.lst

for f in `cat $documents.lst`
do
    in_f=$temp_caption/documents/$f
    out_f=$temp_caption/sentences/$f
    std_f=$temp_caption/standoff/$f
    cd $genia_dir
    ./geniass $in_f $out_f.gss
    ruby sentence2standOff.rb  $in_f $out_f.gss $std_f
    perl geniass-postproc.pl $out_f.gss > $out_f.per
done

cd $script_dir
python combine_sentences.py -i $temp_caption/sentences/ -f .per -o $temp_caption/caption.per

cd $ner_dir/src/tokenizer
./nersuite_tokenizer -multidoc [SEP] < $temp_caption/caption.per > $temp_caption/caption.tok

cd $ner_dir/src/gtagger
./nersuite_gtagger -multidoc [SEP] -d $ner_dir/models/gtagger/ < $temp_caption/caption.tok > $temp_caption/caption.gtag

cd $script_dir
python pre_tokenize.py -g $temp_caption/caption.gtag -t $temp_fulltext/fulltext.tok -o $temp_caption/caption.txt.tok
python fix_offset.py -i $temp_caption/documents/ -f $temp_caption/caption.txt.tok

cd $ner_dir/src/gtagger
./nersuite_gtagger -multidoc [SEP] -d $ner_dir/models/gtagger/ < $temp_caption/caption.txt.off > $temp_caption/caption.txt.gtag


cell_tag=$temp_caption/caption.txt.cell
tiss_tag=$temp_caption/caption.txt.tiss
ggp_tag=$temp_caption/caption.txt.ggp
subc_tag=$temp_caption/caption.txt.subc
chem_tag=$temp_caption/caption.txt.chem
org_tag=$temp_caption/caption.txt.org

cd $ner_dir/nersuite_install/bin
./nersuite tag -C 0.5 -b B-cell:1.5 -multidoc [SEP] -m $models_dir/cell_0.5_1.5.m < $temp_caption/caption.txt.gtag > $cell_tag
./nersuite tag -C 0.5 -b B-tiss:2.5 -multidoc [SEP] -m $models_dir/tiss_0.5_2.5.m < $temp_caption/caption.txt.gtag > $tiss_tag
./nersuite tag -C 0.5 -b B-ggp:1.0 -multidoc [SEP] -m $models_dir/ggp_0.5_1.0.m < $temp_caption/caption.txt.gtag > $ggp_tag
./nersuite tag -C 2 -b B-subc:1.5 -multidoc [SEP] -m $models_dir/subc_2_1.5.m < $temp_caption/caption.txt.gtag > $subc_tag
./nersuite tag -C 0.25 -b B-chem:2.5 -multidoc [SEP] -m $models_dir/chem_0.25_2.5.m < $temp_caption/caption.txt.gtag > $chem_tag
./nersuite tag -C 0.03125 -b B-org:2.5 -multidoc [SEP] -m $models_dir/org_0.03125_2.5.m < $temp_caption/caption.txt.gtag > $org_tag

python extract_CoNLL_recognition.py -x $main_dir/caption/ -o $main_dir/recognition/ -i $tiss_tag,$ggp_tag,$cell_tag,$subc_tag,$chem_tag,$org_tag -a $temp_caption/caption.txt.ori
