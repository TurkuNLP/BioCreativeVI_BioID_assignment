#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysolr
import glob
import re
import process_Ab3P as ab3p
import normalize_cellline as norm_cel
import normalize_chemical as norm_che
import normalize_tissue as norm_tis
import normalize_gocc as norm_goc
import normalize_organism as norm_org
import normalize_geneprot as norm_ggp
import canonicalize_entity as canon_ent
import xml.etree.cElementTree as ET


def get_entity_dict(ent_dict):    
    cel_dict = {}
    che_dict = {}
    ggp_dict = {}
    go_dict = {}
    tiss_dict = {}
    org_dict = {}    
    for k, v in ent_dict.iteritems():
        for item in v:            
            if item[3] == 'cell':
                cel_dict.setdefault(k, []).append(item)
            elif item[3] == 'chem':
                che_dict.setdefault(k, []).append(item)
            elif item[3] == 'ggp':
                ggp_dict.setdefault(k, []).append(item)
            elif item[3] == 'subc':
                go_dict.setdefault(k, []).append(item)
            elif item[3] == 'tiss':
                tiss_dict.setdefault(k, []).append(item)
            elif item[3] == 'org':
                org_dict.setdefault(k, []).append(item)            
    return cel_dict, che_dict, ggp_dict, go_dict, tiss_dict, org_dict


def combine_ent_dict(ent_dict, norm_dict_list, ori_dict_list):
    norm_ent_dict = {}
    for k, v in ent_dict.iteritems():
        norm_ent_dict.setdefault(k, [])
        for i, norm_dict in enumerate(norm_dict_list):            
            try:
                anno_dict = norm_dict[k]
                for j, anno in enumerate(anno_dict):
                    text = ori_dict_list[i][k][j][2:3]                    
                    norm_ent_dict[k].append(anno[0:2] + text + anno[3:])
            except:
                pass        
    return norm_ent_dict


def normalize_all(data_folder, fulltext_folder, ab3p_folder, ori_dict, ent_dict, text_dict, offset_dict, fulltext_dict, solr_obj):
    ab3p_dict = ab3p.ab3p_process(ab3p_folder, data_folder, text_dict, fulltext_dict)
    cel_dict, che_dict, ggp_dict, go_dict, tis_dict, org_dict = get_entity_dict(ent_dict)
    cel_ori, che_ori, ggp_ori, go_ori, tis_ori, org_ori = get_entity_dict(ori_dict)

    canon_ggp = canon_ent.canonicalize_ggp_org(ggp_dict)        

    norm_org_temp, norm_org_dict, all_rank_count = norm_org.normalize_org(org_dict, ab3p_dict)
    norm_cel_temp, norm_cel_dict = norm_cel.normalize_cel(cel_dict, ab3p_dict)
    norm_che_dict = norm_che.normalize_che(che_dict, ab3p_dict)
    norm_go_dict = norm_goc.normalize_goc(go_dict, ab3p_dict)
    norm_tis_temp, norm_tis_dict = norm_tis.normalize_tis(tis_dict, ab3p_dict)    

    norm_ggp_dict = norm_ggp.normalize_ggp(offset_dict, norm_cel_temp, norm_org_temp, norm_tis_temp, ggp_dict, canon_ggp, ab3p_dict, text_dict, all_rank_count, solr_obj)
    
    norm_dict_list = [norm_cel_dict, norm_tis_dict, norm_go_dict, norm_che_dict, norm_org_dict, norm_ggp_dict]
    ori_dict_list = [cel_ori, tis_ori, go_ori, che_ori, org_ori, ggp_ori]
    
    norm_ent_dict = combine_ent_dict(ori_dict, norm_dict_list, ori_dict_list)
    
    return norm_ent_dict
