#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import re
from process_simstring import match_simstring
import simstring
import gzip
import random


def getNormform(synonym):
    """
    """
    return re.sub("[^a-z0-9]", "", synonym.lower())

    
def cell_combine_matching(predict_list, anno_dict):
    norm_cel_temp = {}
    combine_dict = {}
    for doc_id, annos in anno_dict.iteritems():
        for i, anno in enumerate(annos):
            all_predict = {}
            for j, pred in enumerate(predict_list):
                for concept in pred[doc_id][i]:
                    threshold = concept[4]
                    for cel_id in concept[3]:
                        all_predict.setdefault(threshold, {}).setdefault(j, set([])).add(cel_id)
            if all_predict != {}:
                max_t = max(all_predict.keys())
                dict_id = min(all_predict[max_t].keys())
                cell_ids = list(all_predict[max_t][dict_id])
                new_cell_ids = [item for item in cell_ids if 'CVCL_' in item]
                if new_cell_ids == []:
                    new_cell_ids = [item for item in cell_ids if 'CL:' in item]
                    if new_cell_ids == []:
                        new_cell_ids = cell_ids                    
                # predict_id = cell_ids[0]
                predict_id = random.choice(new_cell_ids)                
                combine_dict.setdefault(doc_id, []).append([anno[0], anno[1], anno[2], 'cell', predict_id])
                norm_cel_temp.setdefault(doc_id, []).append([anno[0], anno[1], anno[2], predict_id, 'cell', max_t])        
            else:
                combine_dict.setdefault(doc_id, []).append([anno[0], anno[1], anno[2], 'cell', 'cell:'+ anno[2]])
                norm_cel_temp.setdefault(doc_id, []).append([anno[0], anno[1], anno[2], '', 'cell', max_t])   
    return norm_cel_temp, combine_dict


def load_cellline_mapping(in_file, col_1, col_2, normform):
    map_dict = dict()
    with gzip.open(in_file, "rb") as f:
        for line in f:
            item = line.strip("\n").split("\t")
            if normform:
                nf = getNormform(item[col_2])
            else:
                nf = item[col_2]
            map_dict.setdefault(nf, set([])).add(item[col_1])            
    return map_dict


def cell_match_simstring(ss_folder, annoDict, ab3p_dict, ms, mapping_dict, db_file, normform, min_threshold=0.01):
    distance_matrix = {0: "exact", 1: "dice", 2: "cosine", 3: "jaccard", 4: "overlap"}
    db = simstring.reader(ss_folder + db_file)
    db.measure = ms
    predict_dict = dict()
    for doc_id, items in annoDict.iteritems():
        predict_dict.setdefault(doc_id, [])
        for j, item in enumerate(items):
            predict_dict[doc_id].append([])
            db.threshold = 1.0
            try:
                cell = ab3p_dict[doc_id][item[2]]
            except:
                cell = item[2]
            if normform:
                mention = getNormform(cell)
            else:
                mention = cell
            match_concept = db.retrieve(mention)            
            while match_concept == () and db.threshold > min_threshold:
                db.threshold = db.threshold - 0.01
                match_concept = db.retrieve(mention)
            if len(match_concept) == 0:
                predict_dict[doc_id][j].append([item[0], item[1], item[2], set([]), 0])
            else:
                for k, concept in enumerate(list(set(match_concept))):
                    predict_dict[doc_id][j].append([item[0], item[1], concept, set([]), db.threshold])
                    for concept_id in mapping_dict[concept]:
                        predict_dict[doc_id][j][k][3].add(concept_id)
    return predict_dict


def normalize_cel(anno_dict, ab3p_dict):    
    ss_folder = os.getcwd() + '/simstring/'
    in_dict = load_cellline_mapping(ss_folder + 'cell.tsv.gz', 0, 1, False)
    db_file = 'cell/cell.db'
    predict_list = []
    predict_list.append(cell_match_simstring(ss_folder, anno_dict, ab3p_dict, simstring.cosine, in_dict, db_file, False, 0.01))
    norm_cel_temp, combine_dict = cell_combine_matching(predict_list, anno_dict)
    return norm_cel_temp, combine_dict


