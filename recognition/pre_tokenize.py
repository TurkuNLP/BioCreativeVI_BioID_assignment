#!/usr/bin/env python
# -*- coding: utf-8 -*-        

import argparse

        
def repeat_tokenize(in_gtag, in_tok, out_tok):
    # tax_set = set([])
    # for dictname in ['chebi', 'taxonomy', 'go', 'cell_ontology', 'uberon', 'cellosaurus', 'gene']: 
    #     with open(in_folder + '{}.txt'.format(dictname), 'rb') as f:
    #         for item in f:
    #             tax_set.add(item.split('\t')[0].split())
    new_text = ''
    # token_dict = {}
    tax_set = set([])
    with open(in_tok, 'rb') as f:
        for line in f:
            if line != '\n':
                if "[SEP]" in line:
                    doc_id = line.split("_")[0].split('[SEP]')[1]
                else:
                    lines = line.split('\n')[0].split('\t')
                    tax_set.add(lines[2])
                    # token_dict.setdefault(doc_id, set([])).add(lines[2])
    with open(in_gtag, 'rb') as f:
        for line in f:
            if '[SEP]' in line:                    
                # doc_id = line.split("_")[0].split('[SEP]')[1]
                # tax_set = token_dict[doc_id]
                new_text += line
            elif line[0:5] != '[SEP]' and line != '\n':
                token_list = []
                lines = line.split('\t')
                token = lines[2]
                chunk = lines[5]
                offB = int(lines[0])
                offE = int(lines[1])                    
                if token.isdigit() == False and chunk in ['B-NP\n', 'I-NP\n']:
                    offB = int(lines[0])
                    offE = int(lines[1])
                    i = len(token)
                    while i > 0:
                        if token[0:i] in tax_set:
                            found = token[0:i]
                            token_list.append(found)
                            offE = offB + i
                            new_text += '\t'.join([str(offB), str(offE), found]) + '\n'
                            offB = offE
                            offE = int(lines[1])
                            token = token[len(found):len(token)]
                            i = len(token)
                        i -= 1
                if token != '':
                    new_text += '\t'.join([str(offB), str(offE), token]) + '\n'
            else:                    
                new_text += line
    with open(out_tok, 'wb') as f:
        f.write(new_text)


def argument_parser():
    parser = argparse.ArgumentParser(description="repeat tokenization")
    parser.add_argument("-g", "--in_gtag", type=str, help="in gtag to fix")
    parser.add_argument("-t", "--in_tok", type=str, help="original tokenized full-text file")
    parser.add_argument("-o", "--out_tok", type=str, help="out tokenized file after fixing")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    
    repeat_tokenize(args.in_gtag, args.in_tok, args.out_tok)
    # retokenize(args.dict_folder, args.in_gtag, args.out_tok)


