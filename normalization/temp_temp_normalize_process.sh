



main_dir=/home/sukaew/Biocreative_database
temp_dir=$main_dir/temp_temp_process_BioIDtest_caption

input_dir=$main_dir/temp_temp_nodict_independent_BioIDtest_caption/xml
output_dir=$main_dir/temp_temp_BioIDtest_with_NER
mkdir $output_dir

ori_file=$temp_dir/BioIDtest.txt.ori

Ab3P=$temp_dir/Ab3P
mkdir $Ab3P

scripts=$main_dir/scripts


cd $scripts

python write_BioC.py -x $input_dir/ -o $output_dir/ -a $ori_file -t $temp_dir/documents/ -d $temp_dir/

cd ../BioID_scorer_1_0_3/scripts
python bioid_score.py --force --test_equivalence_classes ../../test_corpus_20170804/test_equivalence_classes.json ../../score/temp_BioIDtest_with_NER/ ../../test_corpus_20170804/caption_bioc/ ../../temp_BioIDtest_with_NER/



