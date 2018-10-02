#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import copy
import re
import gzip
import simstring
from collections import Counter
import random
import numpy as np
from difflib import SequenceMatcher
import sys
import socket


def getNormform(synonym):
    """
    get norm form with no space
    """
    return re.sub("[^a-z0-9]", "", synonym.lower())

        
def load_mapping(in_file, col_1, col_2):
    map_dict = dict()
    with gzip.open(in_file, "rb") as f:
        for line in f:
            item = line.strip("\n").split("\t")
            map_dict.setdefault(item[col_1], set()).add(item[col_2])
    return map_dict


def load_multiple_mapping(in_file, col_1, col_2):
    map_dict = dict()
    with gzip.open(in_file, "rb") as f:
        for line in f:
            item = line.strip("\n").split("\t", 1)
            map_dict.setdefault(item[col_1], set()).add(item[col_2])
    return map_dict


def replace_abbreviation(abb_dict, anno_dict):
    element_dict = {'Ac': 'Actinium', 'Ag': 'Silver', 'Al': 'Aluminum', 'Am': 'Americium', 'Ar': 'Argon', 'As': 'Arsenic', 'At': 'Astatine', 'Au': 'Gold', 'B': 'Boron', 'Ba': 'Barium', 'Be': 'Beryllium', 'Bh': 'Bohrium', 'Bi': 'Bismuth', 'Bk': 'Berkelium', 'Br': 'Bromine', 'C': 'Carbon', 'Ca': 'Calcium', 'Cd': 'Cadmium', 'Ce': 'Cerium', 'Cf': 'Californium', 'Cl': 'Chlorine', 'Cm': 'Curium', 'Co': 'Cobalt', 'Cr': 'Chromium', 'Cs': 'Cesium', 'Cu': 'Copper', 'Db': 'Dubnium', 'Ds': 'Darmstadtium', 'Dy': 'Dysprosium', 'Er': 'Erbium', 'Es': 'Einsteinium', 'Eu': 'Europium', 'F': 'Fluorine', 'Fe': 'Iron', 'Fm': 'Fermium', 'Fr': 'Francium', 'Ga': 'Gallium', 'Gd': 'Gadolinium', 'Ge': 'Germanium', 'H': 'Hydrogen', 'He': 'Helium', 'Hf': 'Hafnium', 'Hg': 'Mercury', 'Ho': 'Holmium', 'Hs': 'Hassium', 'I': 'Iodine', 'In': 'Indium', 'Ir': 'Iridium', 'K': 'Potassium', 'Kr': 'Krypton', 'La': 'Lanthanum', 'Li': 'Lithium', 'Lr': 'Lawrencium', 'Lu': 'Lutetium', 'Md': 'Mendelevium', 'Mg': 'Magnesium', 'Mn': 'Manganese', 'Mo': 'Molybdenum', 'Mt': 'Meitnerium', 'N': 'Nitrogen', 'Na': 'Sodium', 'Nb': 'Niobium', 'Nd': 'Neodymium', 'Ne': 'Neon', 'Ni': 'Nickel', 'No': 'Nobelium', 'Np': 'Neptunium', 'O': 'Oxygen', 'Os': 'Osmium', 'P': 'Phosphorus', 'Pa': 'Protactinium', 'Pb': 'Lead', 'Pd': 'Palladium', 'Pm': 'Promethium', 'Po': 'Polonium', 'Pr': 'Praseodymium', 'Pt': 'Platinum', 'Pu': 'Plutonium', 'Ra': 'Radium', 'Rb': 'Rubidium', 'Re': 'Rhenium', 'Rf': 'Rutherfordium', 'Rg': 'Roentgenium', 'Rh': 'Rhodium', 'Rn': 'Radon', 'Ru': 'Ruthenium', 'S': 'Sulfur', 'Sb': 'Antimony', 'Sc': 'Scandium', 'Se': 'Selenium', 'Sg': 'Seaborgium', 'Si': 'Silicon', 'Sm': 'Samarium', 'Sn': 'Tin', 'Sr': 'Strontium', 'Ta': 'Tantalum', 'Tb': 'Terbium', 'Tc': 'Technetium', 'Te': 'Tellurium', 'Th': 'Thorium', 'Ti': 'Titanium', 'Tl': 'Thallium', 'Tm': 'Thulium', 'U': 'Uranium', 'Uub': 'Ununbium', 'Uuh': 'Ununhexium', 'Uuo': 'Ununoctium', 'Uup': 'Ununpentium', 'Uuq': 'Ununquadium', 'Uus': 'Ununseptium', 'Uut': 'Ununtrium', 'Uuu': 'Ununium', 'V': 'Vanadium', 'W': 'Tungsten', 'Xe': 'Xenon', 'Y': 'Yttrium', 'Yb': 'Ytterbium', 'Zn': 'Zinc', 'Zr': 'Zirconium'}
    new_dict = copy.deepcopy(anno_dict)
    for k, v in anno_dict.iteritems():
        doc_id = k.split('_')[0]
        found_dict = dict()
        new_dict.setdefault(k, v[0:3])
        for i, anno in enumerate(v):
            new_dict[k][i].append("") 
            new_dict[k][i].append([])
            try:
                new_dict[k][i][4] = found_dict[anno[2]]
            except:
                new_text = anno[2]
                try:
                    new_text = abb_dict[doc_id][anno[2]]
                except:
                    pass
                try:
                    new_text = element_dict[anno[2]]
                except:
                    pass
                split_text = re.split(r',\s|/|\band\b|\bor\b', new_text)    
                new_split_text = [new_e.strip() for new_e in split_text]
                new_dict[k][i][4] = new_split_text
                found_dict.setdefault(anno[2], new_split_text)
    return new_dict


def socket_check_value(sock, dict_name, entrezgene_id):
    try:
        message = '"{}", "{}"'.format(dict_name, entrezgene_id)
        sock.sendall(message)
        value = sock.recv(2000)
        fam_id = set([item for item in value.split("'") if item.isdigit() == True or 'CHEBI:' in item]) # .replace("set(['", '').replace("'])", '').split("', '"))
    except:
        fam_id = set([])
    return fam_id


def chem_match_simstring(ss_folder, annoDict, db_file, umls_chemical_symbol):

    db = simstring.reader(ss_folder + db_file)
    db.measure = simstring.cosine
    pred_dict = copy.deepcopy(annoDict)
    found_dict = dict()
    for doc_id, items in annoDict.iteritems():
        for j, item in enumerate(items):
            pred_dict[doc_id][j].append([])
            pred_dict[doc_id][j].append([])
            try:
                pred_dict[doc_id][j][5:8] = found_dict[item[2]]
            except:
                for entity in item[4]:
                    db.threshold = 1.0
                    mention = getNormform(entity)
                    if len(mention) < 5:
                        threshold = 1.0
                    elif len(mention) in range(5, 11):
                        threshold = 0.8
                    else:
                        threshold = 0.6
                    match_concept = db.retrieve(mention)
                    while match_concept == () and db.threshold > threshold:
                        db.threshold = db.threshold - 0.01
                        match_concept = db.retrieve(mention)
                    if len(match_concept) != 0:
                        for concept in match_concept:
                            if concept not in pred_dict[doc_id][j][5]:
                                pred_dict[doc_id][j][5].append(concept)
                                pred_dict[doc_id][j][7].append(db.threshold)
                                # concepts = socket_check_value(sock, 'umls_chemical_symbol', concept)
                                concepts = umls_chemical_symbol[concept]
                                for concept_id in concepts:
                                    pred_dict[doc_id][j][6].append(concept_id)
                found_dict.setdefault(item[2], pred_dict[doc_id][j][5:8])
    return pred_dict


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def ranking_score(item_list):
    nf_list, id_list, score_list = item_list
    new_dict = dict()
    for i, score in enumerate(score_list):
        new_dict.setdefault(score, dict()).setdefault(nf_list[i], id_list[i])
    return new_dict


def combine_prediction(anno_dict, pred_dict, umls_chemical_symbol):#, ctd_symbol_dict): #, umls_chemical_detail):#, tty_class):

    # docs = anno_dict.keys()
    # docs.sort()
    found_dict = dict()
    new_dict = copy.deepcopy(anno_dict)
    for doc, v in pred_dict.iteritems():
        for i, item in enumerate(v):
            try:
                new_dict[doc][i][3:5] = found_dict[item[2]]
            except:
                mesh_set = set([])
                nf_set = set([])
                new_dict[doc][i][4] = ''
                for name in item[4]:
                    try:
                        # exact_match = list(ctd_symbol_dict[getNormform(name)])
                        concepts = umls_chemical_symbol[concept]
                        # concept = socket_check_value(sock, 'umls_chemical_symbol', getNormform(name))
                        exact_match = list(concept)
                        exact_match.sort(reverse=True)
                        if len(set(exact_match)) == 1:# step 1: take exact match from normform (if match term = 1)
                            new_dict[doc][i][4] += exact_match[0] + '|'
                        elif len(set(exact_match)) > 1:# step 2: take RANDOM exact match from normform (if match term > 1)
                            new_dict[doc][i][4] += random.choice(list(exact_match)) + '|'
                    except:
                        score_dict = ranking_score(item[5:])
                        sorted_score = sorted(score_dict.keys(), reverse=True)
                        if sorted_score != []:
                            max_pred = len(item[4])
                            l = 0
                            cum_score = []
                            nf_set = set([])
                            while l < len(sorted_score) and l < max_pred:
                                max_score = sorted_score[l]
                                nf_set = nf_set | set([nf for nf, cui in score_dict[max_score].iteritems()])    
                                mesh_set = mesh_set | set([cui for nf, cui in score_dict[max_score].iteritems()])
                                cum_score.append(max_score)
                                l += 1
                            if len(mesh_set) == max_pred:
                                for scui in mesh_set:
                                    new_dict[doc][i][4] += scui + '|'
                            else:
                                if mesh_set != set([]):
                                    new_dict[doc][i][4] += random.choice(list(mesh_set)) + '|'
                        else:
                            if new_dict[doc][i][4] == '':
                                try:
                                    sim_dict = dict()
                                    for x, norm_id in enumerate(new_dict[doc][0:i]):
                                        if norm_id[3] != '':
                                            sim_score = similar(norm_id[2], new_dict[doc][i][2])
                                            sim_dict.setdefault(sim_score, dict()).setdefault(x, norm_id[3].strip('|'))
                                    match_cui = sim_dict[max(sim_dict.keys())]
                                    new_dict[doc][i][4] += max(match_cui.keys()) + '|'
                                except:
                                    pass
        for k, v in new_dict.iteritems():
            for anno in v:
                anno[3] = '|'.join(set(anno[3].split('|'))).strip('|')
        found_dict.setdefault(anno[2], new_dict[doc][i][3:5])
    return new_dict


def normalize_che(chem_test, ab3p_test):
    data_folder = os.getcwd()
    ss_folder = data_folder + "/simstring/"

    import time
    start_time = time.time()
    umls_file = "norm_chemical.tsv.gz"
    umls_chemical_symbol = load_mapping(ss_folder + umls_file, 1, 0)
        
    split_test = replace_abbreviation(ab3p_test, chem_test)
                
    cosine_umls_Test = chem_match_simstring(ss_folder, split_test, 'pubchem/pubchem.db', umls_chemical_symbol)
    combine_chemical_Test = combine_prediction(split_test, cosine_umls_Test, umls_chemical_symbol) 

    final_dict = {}
    for k, v in combine_chemical_Test.iteritems():
        for anno in v:
            if '|' in anno[4]:
                anno_id = anno[4].strip('|')
                if 'CHEBI' in anno_id:
                    final_dict.setdefault(k, []).append([anno[0], anno[1], anno[2], anno[3], anno_id])
                else:
                    final_dict.setdefault(k, []).append([anno[0], anno[1], anno[2], anno[3], 'PubChem:'+anno_id])    
            else:
                final_dict.setdefault(k, []).append([anno[0], anno[1], anno[2], anno[3], 'molecule:'+anno[2]])
    return final_dict
