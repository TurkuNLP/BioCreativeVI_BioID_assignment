#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import os
# import re
import glob
# import htmlentitydefs
import argparse
# from extract_BioC_text import *


def runCommand(cmd_list):
    try:
        subprocess.check_call(cmd_list, shell=True)
    except subprocess.CalledProcessError:
        print 'errors in the called executable', cmd_list
        pass # handle errors in the called executable
    except OSError:
        print 'executable not found', cmd_list
        pass # 

    
def write_temp_file(in_dict, in_file):
    text = ''
    for k, v in in_dict.iteritems():        
        text += 'pmc_' + k + '\n'
        text += v
    with open(in_file, 'wb') as f:
        f.write(text)        

    
def run_ab3p(ab3p_folder, in_file, out_file):
    try:
        cmd_line = 'cd {} && ./identify_abbr {} > {}'.format(ab3p_folder, in_file, out_file)
        runCommand(cmd_line)
    except:
        print('no Ab3P installed/found, normalization performance will be likely low') 
        with open(out_file, 'wb') as f:
            f.write('')
            
        
def old_process_output(out_file, PubMed, fulltext_mapping, in_dict):
    ab3p_mapping = dict()
    with open(out_file, 'rb') as f:
        if PubMed == 'PubMed':
            docs = f.read().split('\n\n')
        else:            
            docs = ['pmc_' + item for item in f.read().split('pmc_')[1:]]    
    for doc in docs:
        lines = doc.split('\n')
        doc_id = lines[0].split('\n')[0]        
        for line in lines[2:]:
            split_line = line.split('|')
            if line[0:2] == "  " and len(split_line) == 3:
                short_form, long_form, conf = split_line
                ab3p_mapping.setdefault(doc_id.split('pmc_')[1], {}).setdefault(short_form[2:], long_form)
    doc_map = {}
    new_mapping = {}
    for k in in_dict.iterkeys():
        pmcid = k.split('_')[0]        
        doc_map.setdefault(k, pmcid)
        try:
            for i, j in ab3p_mapping[k].iteritems():
                new_mapping.setdefault(pmcid, {}).setdefault(i, j)
        except:
            pass
        try:
            for x, y in fulltext_mapping[pmcid].iteritems():
                new_mapping.setdefault(pmcid, {}).setdefault(x, y)        
        except:
            pass            
    new_ab3p_mapping = {}
    for captions, pmcid in doc_map.iteritems():
        try:
            new_ab3p_mapping.setdefault(captions, new_mapping[pmcid])
        except:
            new_ab3p_mapping.setdefault(captions, {})        
    return new_ab3p_mapping


def process_output(out_file):
    ab3p_mapping = dict()
    with open(out_file, 'rb') as f:
        docs = ['pmc_' + item for item in f.read().split('pmc_')[1:]]
    for doc in docs:
        lines = doc.split('\n')
        doc_id = lines[0].split('\n')[0].replace('pmc_', '')
        for line in lines[2:]:
            split_line = line.split('|')
            if line[0:2] == "  " and len(split_line) == 3:
                short_form, long_form, conf = split_line
                ab3p_mapping.setdefault(doc_id, {}).setdefault(short_form[2:], long_form)
    return ab3p_mapping


def combine_caption_fulltext(caption, fulltext):
    caption_docs = {k: k.split('_')[0] for k in caption.keys()}
    for doc_id, abb in caption.iteritems():
        for match_doc in [doc for doc in fulltext.iterkeys() if doc in doc_id]:
            fulltext[match_doc].update(abb)
    for doc_id, abb in fulltext.iteritems():
        for match_doc in [doc for doc in caption.iterkeys() if doc_id in doc]:
            caption[match_doc].update(abb)
    return caption


def ab3p_process(ab3p_folder, data_folder, text_dict, fulltext_dict):
    text_file_in = data_folder + 'Ab3P/text.txt'
    text_file_out = text_file_in.replace('.txt', '.o')
    
    write_temp_file(text_dict, text_file_in)
    run_ab3p(ab3p_folder, text_file_in, text_file_out)
    text_ab3p_dict = process_output(text_file_out)
    
    full_file_in = data_folder + 'Ab3P/fulltext.txt'
    full_file_out = full_file_in.replace('.txt', '.o')
    write_temp_file(fulltext_dict, full_file_in)
    run_ab3p(ab3p_folder, full_file_in, full_file_out)
    fulltext_ab3p_dict = process_output(full_file_out)

    if fulltext_ab3p_dict == {}:        
        return text_ab3p_dict
    else:
        ab3p_mapping = combine_caption_fulltext(text_ab3p_dict, fulltext_ab3p_dict)
        return ab3p_mapping
