#!/usr/bin/env python
# -*- coding: utf-8 -*

import os
import re
import gzip
import random
import simstring
import subprocess
import cPickle as pickle
from difflib import SequenceMatcher
    

def load_org_mapping(in_file, col_1, col_2, dtype):
    map_dict = dict()
    with gzip.open(in_file, "rb") as f:
        for line in f:
            item = line.strip("\n").split("\t")
            if dtype == 'set':
                map_dict.setdefault(item[col_1], set()).add(item[col_2])
            elif dtype == 'lst':
                map_dict.setdefault(item[col_1], []).append(item[col_2])
            elif dtype == 'int':
                map_dict.setdefault(item[col_1], int(item[col_2]))
            elif dtype == 'str':
                map_dict.setdefault(item[col_1], item[col_2])
    return map_dict


def create_abbreviation(item):
    sci_name = set([])
    items = item.split()
    if len(items) >= 2:
        for splitter in set([".", ". ", " "]):
            if item[0].isupper() and item[1].islower():
                sci_name.add(items[0][0] + splitter + items[1]) 
                sci_name.add(items[0][0:2] + splitter + items[1])
    return sci_name


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

    
def take_popular(match_species, all_rank_count):#, ssp_rank_count, no_rank_count):
    pmc_count = dict()
    for item in match_species:
        try:
            pmc_count.setdefault(int(all_rank_count[item]), []).append(item)
        except:
            pass
    new_pmc = sorted(pmc_count.keys(), reverse=True)    
    try:
        return random.choice(pmc_count[new_pmc[0]])
    except:
        return ""

        
def check_long_form(doc, alice_mapping):
    """
    usage: get the mapping b/w long_form: short form for each document 
    input: take the document in
    output: return the matching short_form: long_form
    """
    try:
        long_name = alice_mapping[doc]
    except:
        long_name = dict()
    return long_name

    
def check_normform(long_name, short_name):
    try:
        return long_name[short_name]
    except:
        return short_name
    

def load_normalize_org_mapping(data_folder):
    """
    load all the mapping files used in normalization of organisms
    """
    map_folder = data_folder + "map_files/"     
    ncbi_symbol_dict = load_org_mapping(map_folder + 'new_all_org.tsv.gz', 0, 1, 'lst')
    id_rank_dict = load_org_mapping(map_folder + "map_organism2rank.tsv.gz", 0, 1, 'str')
    tax_tree = pickle.load(gzip.open(data_folder + "pickle/tax_tree.pkl.gz", "rb"))
    with gzip.open(map_folder + "spp_only.txt.gz", "rb") as f:
        all_under_spp = set([line.strip("\n") for line in f])
    lower_rank_map = load_org_mapping(map_folder + "above_spp.tsv.gz", 0, 1, 'set')
    all_rank_count = load_org_mapping(map_folder + "taxid_pmc_count.tsv.gz", 0, 1, 'int')
    with gzip.open(map_folder + "model_org_wiki_taxid.txt.gz", "rb") as f:
        model_org = set([line.strip("\n") for line in f if line.strip("\n") != ""])
    return ncbi_symbol_dict, id_rank_dict, tax_tree, all_under_spp, lower_rank_map, all_rank_count, model_org


def check_org_simstring(ss_folder, db, min_threshold, mention):
    db.threshold = 1.0
    match_concept = db.retrieve(mention)
    while match_concept == () and db.threshold > min_threshold:
        db.threshold = db.threshold - 0.01
        match_concept = db.retrieve(mention)
    return match_concept, db.threshold


def take_model_org(maxtax, id_rank_dict, tax_tree, model_org):
    match_pid = set([])
    for query_taxid in maxtax:
        query_lineage = tax_tree[query_taxid]
        for model_id in model_org:
            model_lineage = tax_tree[model_id]
            common_lineage = set(model_lineage) & set(query_lineage)
            for item in common_lineage:
                if id_rank_dict[item] == 'genus' or id_rank_dict[item] == 'subgenus':
                    match_pid.add(model_id)
    return list(match_pid)


def prev_exact(new_pred, maxtax, i):
    match_pid = ""
    source = ""
    j = i - 1
    while j >= 0:
        if new_pred[j][3] != '' and new_pred[j][3] in maxtax:
            match_pid = new_pred[j][3]
            source = 'prev mention, same species'
            return match_pid, source
        j = j - 1
    return match_pid, source


def prev_abbreviation(new_pred, maxtax, i):
    match_pid = ""
    source = ""
    j = i - 1
    while j >= 0:
        if new_pred[j][3] != "": #previous mention
            abbrev_name = create_abbreviation(new_pred[j][2]) # create abbreviation names for previous mention
            if new_pred[i][2] in abbrev_name:# and new_pred[j][3] in maxtax:
                match_pid = new_pred[j][3]        
                source = "abbreviation match"
                return match_pid, source
                j = -1
        j = j - 1
    return match_pid, source    


def prev_acronym(new_pred, maxtax, i):
    match_pid = ""
    source = ""
    j = i - 1
    while j >= 0:
        if new_pred[j][3] != "":
            acronym = create_acronym(new_pred[j][2])            
            if new_pred[i][2].upper() in acronym:
                match_pid = new_pred[j][3]        
                source = "acronym match"
                return match_pid, source
            else:
                for acro in acronym:
                    if acro in new_pred[i][2] or new_pred[i][2] in acro:
                        match_pid = new_pred[j][3]        
                        source = "acronym match"
                        return match_pid, source                        
        j = j - 1
    return match_pid, source    


def prev_genus(new_pred, maxtax, i, tax_tree, id_rank_dict):
    match_pid = ""
    source = ""
    j = i - 1
    upper_rank = {}
    for item in maxtax:
        match_tax = [tax_tree[item][x+1] for x, t_rank in enumerate(tax_tree[item]) if id_rank_dict[t_rank] == 'species']
        for match_t in match_tax:
            upper_rank.setdefault(match_t, set([])).add(item)
        # upper_rank.setdefault(match_tax[0], set([])).add(item)
    while j >= 0:
        if new_pred[j][3] != '':
            prev_tax = [tax_tree[new_pred[j][3]][k+1] for k, t_rank in enumerate(tax_tree[new_pred[j][3]]) if id_rank_dict[t_rank] == 'species']
            try:
                shared_genus = list(upper_rank[prev_tax[0]])
                if len(shared_genus) == 1:
                    return shared_genus[0], 'shared genus'
            except:
                pass
        j = j - 1
    return match_pid, source    


def take_prev_mention(new_pred, i, maxtax, mention, tax_tree, id_rank_dict):
    source = ""
    match_pid = ""    
    match_pid, source = prev_exact(new_pred, maxtax, i)
    if match_pid == "":
        match_pid, source = prev_abbreviation(new_pred, maxtax, i)
    if match_pid == "":
        match_pid, source = prev_acronym(new_pred, maxtax, i)
    if match_pid == "":
        match_pid, source = prev_genus(new_pred, maxtax, i, tax_tree, id_rank_dict)
    return match_pid, source


def take_lower_rank(maxtax, lower_rank_map, id_rank_dict):
    spp_id = []
    for tax_id in maxtax:
        if tax_id in lower_rank_map.keys():
            try:
                taxset = lower_rank_map[tax_id] - set([tax_id])
                spp_id += list(taxset)                
            except:
                pass
    return spp_id


def take_species(maxtax, all_under_spp, tax_tree):
    match_pid = dict()
    common = set([])
    # step 1: take species ranking
    for pid in maxtax:
        if pid in all_under_spp:
            match_pid.setdefault(pid, tax_tree[pid])
    all_pid = match_pid.keys()
    for pid in all_pid:
        for k, v in match_pid.iteritems():
            if pid != k and pid in v:
                common.add(pid)
    return set(match_pid.keys()) - common


def take_rank(ss_folder, id_rank_dict, tax_tree, maxtax, new_pred, i, curr_pred, all_rank_count, model_org, all_under_spp, mention, lower_rank_map):
    new_id = ""
    source = ""
    match_species = take_species(maxtax, all_under_spp, tax_tree)
    if len(match_species) == 1:
        new_id = list(match_species)[0]
        source = "match species"
    else:
        if len(match_species) == 0:
            match_species = take_lower_rank(maxtax, lower_rank_map, id_rank_dict)
        new_id, source = take_prev_mention(new_pred, i, match_species, mention, tax_tree, id_rank_dict)
        if new_id == "":
            match_model_org = take_model_org(match_species, id_rank_dict, tax_tree, model_org)
            if len(match_model_org) == 1:
                new_id = list(match_model_org)[0]
                source = "match species --> model organism"
            else:
                new_id = take_popular(match_species, all_rank_count)
                if new_id != '':
                    source = "match species --> previous mention --> popular studies"
    return new_id, source


def short_name_similar(curr_pred, mention):
    new_id = ''
    source = ''
    score = 0
    sim_dict = dict()
    for item in curr_pred:
        sim_score = similar(item[2], mention)
        sim_dict.setdefault(sim_score, [])
        for taxid in item[3]:
            if taxid != "":
                sim_dict[sim_score].append(taxid)
    score = sim_dict.keys()
    score.sort()
    if score != []:
        new_id = random.choice(sim_dict[score[-1]])
        source = 'short name similar'
        score = score[-1]
    return new_id, source, score

        
def take_most_similar(new_pred, i, mention):
    new_id = ''
    source = ''
    score = 0
    sim_dict = dict()
    for item in new_pred:
        if item[3] != '':
            sim_score = similar(item[2], mention)
            sim_dict.setdefault(sim_score, item[3])
    score = sim_dict.keys()
    score.sort()
    if score != []:
        new_id = sim_dict[score[-1]]
        source = 'most similar match'
        score = score[-1]
    return new_id, source, score


def create_acronym(item):
    acronym = set([])
    split_name = re.split("(virus|viruses)", item)[0]
    if len(split_name) >1:
        try:
            name_1 = "".join([name[0].upper() for name in re.split("[^A-Za-z0-9]", split_name) if name != ""]) + 'V'
            acronym.add(name_1)
            name_3 = name_1.strip('V') + ' VIRUS'
            acronym.add(name_3)
        except:
            pass
    return acronym


def prev_common_str(new_pred, i, mention, extra_name):
    match_pid = ""
    score_dict = dict()
    curr_name = ' '.join([name.strip().lower() for name in re.split(extra_name, mention) if name != ''])
    curr_name = set(re.split("[^A-Za-z0-9]", curr_name)) - set([''])
    for j, item in enumerate(new_pred):
        if item[3] != '':
            prev_name = ' '.join([name.strip().lower() for name in re.split(extra_name, new_pred[j][2]) if name != ''])        
            prev_name = set(re.split("[^A-Za-z0-9]", prev_name)) - set([''])
            common =  curr_name & prev_name
            divider = min([len(curr_name), len(prev_name)])
            if divider == 0:
                score = 0
            else:
                score = len(common) / divider
            score_dict.setdefault(score, set([])).add(item[3])
    if score_dict != {}:
        maxscore = score_dict.keys()
        maxscore.sort()
        if maxscore[-1] > 0.5:
            match_pid = random.choice(list(score_dict[maxscore[-1]]))
            return match_pid, 'common token', score    
    return match_pid, '', 0    


def prev_substring(new_pred, i, mention, extra_name):
    match_pid = ""
    score_dict = dict()
    score = 0
    for j, item in enumerate(new_pred):
        if item[3] != '':
            if item[2] in mention:
                score = len(item[2]) / len(mention)*1.0
            elif mention in item[2]:
                score = len(mention) / len(item[2])*1.0
            score_dict.setdefault(score, set([])).add(item[3])
    if score_dict != {}:
        maxscore = score_dict.keys()
        maxscore.sort()
        if maxscore[-1] > 0.5:
            match_pid = random.choice(list(score_dict[maxscore[-1]]))
    return match_pid, '', 0    


def match_org_simstring(ss_folder, annoDict, mapping_dict, db, min_threshold, alice_mapping, all_under_spp):
    extra_name = re.compile('spec\.\s|nov\.\s|n\.\s|f\.\s|gen\.\s|comb\.\s|nom\.\s|approb\.\s|rev\.\s|nov\.\s|gen\.\s|spec\.\s|ord\.\s|nov\.\s|pv\.\s|var\.\s|variety\s|spp\.\s|aff\.\s|nr\.\s|ex\.\s|s\.l\.\s|subsp\.\s|strain\s|Strain\s|n\.\s|sp[.]\s')
    predict_dict = {}
    found_dict = {}
    for doc, items in annoDict.iteritems():
        long_name = check_long_form(doc, alice_mapping)
        predict_dict.setdefault(doc, [])
        for j, item in enumerate(items):
            predict_dict[doc].append([])
            try:
                match_concept, f_threshold = found_dict[item[2]]
            except:
                new_name = check_normform(long_name, item[2])
                mention = ' '.join([name.strip() for name in re.split(extra_name, new_name) if name != '']).lower()
                if mention == ' ':
                    predict_dict[doc][j].append([item[0], item[1], item[2], set([]), 0])
                else:
                    match_concept, f_threshold = check_org_simstring(ss_folder, db, min_threshold, mention)
            if match_concept == ():
                predict_dict[doc][j].append([item[0], item[1], item[2], set([]), 0])
            else:
                common = set([comm for comm in match_concept if set(mapping_dict[comm]) & all_under_spp != set([])])
                if common != set([]):
                    match_concept = list(common)
                else:
                    match_concept = list(match_concept)
                for k, concept in enumerate(match_concept):
                    predict_dict[doc][j].append([item[0], item[1], concept, set([]), f_threshold])
                    for concept_id in mapping_dict[concept]:
                        predict_dict[doc][j][k][3].add(concept_id)
            found_dict.setdefault(item[2], [match_concept, f_threshold])                
    return predict_dict


def combine_prediction(ss_folder, anno, pred_list, id_rank_dict, tax_tree, ncbi_symbol_dict, all_rank_count, model_org, all_under_spp, lower_rank_map, alice_mapping):
    """
    1. match the tokens using all three dictionaries, NCBI, UMLS, abbreviation
    2. combine alll three matchings to take only maximum matchs using db.threshold
    3. among the maximum, take the match that is "species" in taxonomic ranking
    4. if none of the max match is "species" ranking --> look into the taxonomic tree  
    """
    extra_name = re.compile('spec\.\s|nov\.\s|n\.\s|f\.\s|gen\.\s|comb\.\s|nom\.\s|approb\.\s|rev\.\s|nov\.\s|gen\.\s|spec\.\s|ord\.\s|nov\.\s|pv\.\s|var\.\s|variety\s|spp\.\s|aff\.\s|nr\.\s|ex\.\s|s\.l\.\s|subsp\.\s|strain\s|Strain\s|n\.\s|sp[.]\s')
    docs = anno.keys()
    new_predict = {}
    for doc in docs:
        mentions = []
        known_predict = {}
        long_name = check_long_form(doc, alice_mapping)
        new_predict.setdefault(doc, [])
        for i, item in enumerate(anno[doc]):
            mention = ' '.join([name.strip() for name in re.split(extra_name, item[2]) if name != '']).lower()
            mentions.append(mention)
            if mention in known_predict.keys():
                new_predict[doc].append(item[0:3] + known_predict[mention])
            else:
                new_predict[doc].append([item[0], item[1], item[2], "", "", ""])
                n_name = check_normform(long_name, item[2])
                try: # candidate taxid is from the exact string matching
                    exact_string = list(ncbi_symbol_dict[mention])
                    new_predict[doc][i][3], new_predict[doc][i][4] = take_rank(ss_folder, id_rank_dict, tax_tree, exact_string, new_predict[doc], i, pred_list[doc], all_rank_count, model_org, all_under_spp, mention, lower_rank_map)
                    new_predict[doc][i][5] = 1.0
                except: # candidate taxid from fuzzy matching
                    if pred_list[doc][i][0][3] != set([]):
                        maxtax = set([])
                        curr_score = pred_list[doc][i][0][4]
                        for x in pred_list[doc][i]:
                            for y in x[3]:
                                maxtax.add(y)
                        new_predict[doc][i][3:5] = take_rank(ss_folder, id_rank_dict, tax_tree, maxtax, new_predict[doc], i, pred_list[doc], all_rank_count, model_org, all_under_spp, mention, lower_rank_map)
                        new_predict[doc][i][5] = curr_score 
                if new_predict[doc][i][3] != "":
                    known_predict.setdefault(mention, new_predict[doc][i][3:])
        for i, item in enumerate(anno[doc]):            
            new_name = mentions[i]
            if new_predict[doc][i][3] == '':
                new_predict[doc][i][3:] = prev_common_str(new_predict[doc], i, new_name, extra_name)
            if new_predict[doc][i][3] == "":
                new_predict[doc][i][3:] = prev_substring(new_predict[doc], i, new_name, extra_name)  
            if new_predict[doc][i][3] == '':
                new_predict[doc][i][3:] = take_most_similar(new_predict[doc], i, new_name)
            known_predict.setdefault(new_name, new_predict[doc][i][3:])
    return new_predict


def normalize_org(org_dict, abbr_dict):
    print 'load organism mapping files ...'
    data_folder = os.getcwd() + '/'
    ncbi_symbol_dict, id_rank_dict, tax_tree, all_under_spp, lower_rank_map, all_rank_count, model_org = load_normalize_org_mapping(data_folder)
    min_threshold = 0.6
    ss_folder = data_folder + "simstring/"
    ncbi_db = 'new_all_org/new_all_org.db'
    db = simstring.reader(ss_folder + ncbi_db)
    db.measure = simstring.cosine
    pred_tests = match_org_simstring(ss_folder, org_dict, ncbi_symbol_dict, db, min_threshold, abbr_dict, all_under_spp)
    combined_tests = combine_prediction(ss_folder, org_dict, pred_tests, id_rank_dict, tax_tree, ncbi_symbol_dict, all_rank_count, model_org, all_under_spp, lower_rank_map, abbr_dict)
    final_dict = {}
    temp_dict = {}
    for k, v in combined_tests.iteritems():
        for anno in v:
            temp_dict.setdefault(k, []).append([anno[0], anno[1], anno[2], anno[3], anno[5]])
            if anno[3] == '':
                final_dict.setdefault(k, []).append([anno[0], anno[1], anno[2], 'org', 'organism:' + anno[2]])
            else:
                final_dict.setdefault(k, []).append([anno[0], anno[1], anno[2], 'org', 'NCBI taxon:' + anno[3]])
    return temp_dict, final_dict, all_rank_count # combined_tests
