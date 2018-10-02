#!/usr/bin/python3
# -*- coding: utf-8 -*-

import glob
import argparse


def combine_file(in_folder, file_ext, out_file):
    data = glob.glob(in_folder + '*' + file_ext)
    data.sort()
    all_text = ''
    doc_set = set([])
    for item in data:
        doc_id = item.split('/')[-1]
        if doc_id not in doc_set:
            all_text += '[SEP]' + doc_id + '\n\n'
            doc_set.add(doc_id)
        with open(item, 'rb') as f:
            read_text = f.read()
            if read_text != '':
                all_text += read_text.strip('\n') + '\n'
            all_text += '\n'
    with open(out_file, 'wb') as f:
        f.write(all_text)

                
def argument_parser():
    parser = argparse.ArgumentParser(description="combine file for training and tagging")
    parser.add_argument("-i", "--in_folder", type=str, default= 'sentences/', help="folder where files want to combine locate")
    parser.add_argument("-f", "--file_ext", type=str, default= '.ann', help="file extension we combine")
    parser.add_argument("-o", "--out_file", type=str, default= '/media/suwisa/Biocreative/result_2/devel_original.txt.tag', help="output filename")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    combine_file(args.in_folder, args.file_ext, args.out_file)
