#!/usr/bin/env python
# -*- coding: utf-8 -*-        

import os
from collections import Counter
import random
import copy
import roman
import re
import gzip
import simstring
import subprocess
import cPickle as pickle
from difflib import SequenceMatcher


def extract_entity_abstract(data_folder, text_dict, norm_ent_dict, threshold):
    """
    extract extra mentioned of organisms that NER missed 
    """
    for doc_id in norm_ent_dict.iterkeys():
        org_name = {item[2]: item[3:] for item in norm_ent_dict[doc_id] if item[-1] >= threshold}
        org_offB = [item[0] for item in norm_ent_dict[doc_id] if item[-1] >= threshold]
        full_text = text_dict[doc_id]
        
        for org_str in org_name: #loop through all organism names
            norm_org = remove_brackets(org_str)
            try:
                offB_index = [m.start() for m in re.finditer(norm_org, full_text)]
            except:
                offB_index = []
            for offB in offB_index:
                if str(offB) not in org_offB: # check if that organism is already recognized
                    offE = offB + len(norm_org)
                    new_list = [str(offB), str(offE), norm_org] + org_name[org_str]
                    norm_ent_dict[doc_id].append(new_list)
        norm_ent_dict[doc_id].sort(key=lambda x: int(x[0]))
    
    return norm_ent_dict


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


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def runCommand(cmd_list):
    try:
        subprocess.check_call(cmd_list, shell=True)
    except subprocess.CalledProcessError:
        # print cmd_list
        pass # handle errors in the called executable
    except OSError:
        # print cmd_list
        pass # executable not found

    
def getNormform(synonym):
    return re.sub("[^a-z0-9]", "", synonym.lower())


def getPureform(synonym):
    return re.sub("[\[\]\(\)\{\}]", "", synonym)


def getNormform_space(synonym):
    """
    get norm form with space
    """
    a = re.sub("[^a-z0-9]", " ", synonym.lower())
    return " ".join(a.split())


def turn_gold_pred(gold_dict):
    new_dict = dict()
    for doc, annos in gold_dict.iteritems():
        new_dict.setdefault(doc, [])
        for anno in annos:
            new_dict[doc].append([anno[0], anno[1], anno[2], '', [], [], [anno[2]], 0.0])
    return new_dict


def modify_match_dis_normform(data_folder, annoDict, db_file, len_text, min_threshold=0.01):
    print data_folder + "simstring/" + db_file
    db = simstring.reader(data_folder + "simstring/" + db_file)
    db.measure = simstring.cosine
    predict_dict = dict()
    for doc_id, items in annoDict.iteritems():
        j = 0
        for item in items:
            db.threshold = 1.0
            mention = getNormform_space(item[2])
            match_concept = db.retrieve(mention)
            while match_concept == () and db.threshold > min_threshold:
                db.threshold = db.threshold - 0.01                
                match_concept += db.retrieve(mention)
            if len(match_concept) == 1 and len(mention) > len_text:
                predict_dict.setdefault(doc_id, [])
                for concept in set(match_concept):
                    predict_dict[doc_id].append([item[0], item[1], concept, '9606', db.threshold])
                j += 1
    return predict_dict


def modify_match_tis_normform(data_folder, annoDict, db_file, len_text, min_threshold=0.01):
    predict_dict = dict()
    for doc_id, items in annoDict.iteritems():
        j = 0
        for item in items:
            mention = getNormform_space(item[2])
            predict_dict.setdefault(doc_id, []).append([item[0], item[1], mention, '9606', 0.0])
            j += 1
    return predict_dict


def match_cel_org(cel_dict, acc_org):
    new_dict = copy.deepcopy(cel_dict)
    for k, v in new_dict.iteritems():
        for anno in v:
            if 'CVCL' in anno[3]:# and anno[4] != '':
                anno[4] = acc_org[anno[3]]
    
    return new_dict


def remove_brackets(synonym):
    a = re.sub("[\{\}\[\]\(\)]", " ", synonym)
    return " ".join(a.split())


def combine_disease_organism_cellline(doc_id, norm_dis, norm_cel, norm_org):
    spp = []
    if doc_id in norm_dis.keys():
        spp += norm_dis[doc_id]
    if len(spp) < 3:
        spp = []
    if doc_id in norm_org.keys():
        spp += norm_org[doc_id]        
    if doc_id in norm_cel.keys():
        for anno in norm_cel[doc_id]:
            if anno[3] == '10029': #expression vector for human genes
                spp.append([anno[0], anno[1], anno[2], '9606', anno[4]])
            elif anno[3].isdigit == True:
                spp.append(anno)
    sort_spp = sorted(spp, key=lambda x: int(x[0]))
    return sort_spp


def divide_parentheses(s):
    """
    form varieties of entity
    by replacing the text outside the parentheses
    with the text inside the parentheses
    """
    new_s = ''
    text = [m.group() for m in re.compile('\(.*?\)').finditer(s)]        
    spans = [m.start() for m in re.compile('\(.*?\)').finditer(s)]
    if text != [] and spans != []:
        offE = 0
        for i, item in enumerate(spans):
            offE = item + len(text[i])
            new_s += s[item+1: offE-1] + ' '
        new_s += s[offE+1: len(s)]
        for term in text:
            s = s.replace(term, '')
        return [s, new_s]
    else:
        return [s]


def roman_to_int_repl(match):
    try:
        return str(roman.fromRoman(match.group(0)))
    except:
        return str(match.group(0))

    
def roman2arabic(s):
    regex = re.compile(r'\b(?=[MDCLXVI]+\b)M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})\b')
    return regex.sub(roman_to_int_repl, s.upper()).lower()


def greek_end(mention):
    greek_dict = {'alpha': 'a', 'beta': 'b'}
    greek_set = set(['omega', 'upsilon', 'omicron', 'lambda', 'kappa', 'iota', 'theta', 'zeta', 'epsilon', 'delta', 'gamma', 'beta', 'alpha'])
    pat_2 = re.compile('(omega|upsilon|omicron|lambda|kappa|iota|theta|zeta|epsilon|delta|gamma|beta|alpha)', re.X)
    t_alpha = re.split(pat_2, mention)
    alpha = ' '.join([item.strip() for item in t_alpha if item != '']).strip()
    no_greek = [item.strip() for item in t_alpha if item.strip() not in greek_set and item != '']
    with_greek = [item.strip() for item in t_alpha if item.strip() in greek_set and item != '']
    beta = ' '.join(no_greek + with_greek)
    try:
        gamma = ' '.join(no_greek + [greek_dict[with_greek[0]]])
    except:
        gamma = ''
    return alpha, beta, gamma


def modify_match_ggp_normform(data_folder, pred_dict):
    """
    generate the combinations of names
    1. normform
    2. parenthesis form
    3. roman --> arabic numbers
    4. move first greek letter to the last
    """
    in_paren = ''
    found_dict = dict()
    for doc_id, annos in pred_dict.iteritems():
        for j, item in enumerate(annos):
            try:
                pred_dict[doc_id][j][6] = found_dict[item[2]]
            except:
                in_paren = divide_parentheses(item[2])
                pred_dict[doc_id][j][6] += [text for text in set(in_paren) if text not in pred_dict[doc_id][j][6] and text != '']
                combine = []
                for old_ggp in pred_dict[doc_id][j][6]:
                    arabic = roman2arabic(old_ggp)
                    alpha, beta, gamma = greek_end(old_ggp)
                    combine += [arabic, alpha, beta, gamma]
                pred_dict[doc_id][j][6] += [ggp for ggp in set(combine) if ggp not in pred_dict[doc_id][j][6] and ggp != '']
                found_dict.setdefault(item[2], pred_dict[doc_id][j][6])
            pred_dict[doc_id][j][7] = 1.0
    return pred_dict    


def add_abbreviation(abbr_dict, in_dict):
    docs = set(abbr_dict.keys()) & set(in_dict.keys())
    for k in docs:
        for i, anno in enumerate(in_dict[k]):
            long_form = ""
            temp_dict = {b.lower(): a for a, b in abbr_dict[k].iteritems()}            
            try:
                long_form = abbr_dict[k][anno[2]]
            except:
                try:
                    long_form = temp_dict[anno[2].lower()]
                except:
                    try:
                        pattern = '|'.join(abbr_dict[k].keys())
                        match = re.match(r'({})'.format(pattern), anno[2]) 
                        long_form = anno[2].replace(match.group(1), abbr_dict[k][match.group(1)])
                    except:
                        pass
            if long_form != '' and long_form not in in_dict[k][i][6]:
                in_dict[k][i][6].append(long_form)
    return in_dict


def add_canonical(canon_dict, in_dict):
    new_dict = copy.deepcopy(in_dict)
    for k, v in canon_dict.iteritems():
        for i, anno in enumerate(v):
            new_dict[k][i][6] += [item for item in anno[4]]
            if ' genes' in anno[2] and 'gen' in new_dict[k][i][6]:
                new_dict[k][i][6].remove('gen')
    return new_dict


def extract_org_ggp(new_dict, pred_org_id_dict, sci_name, threshold):
    """
    extract mentions of organisms inside the mention of ggp
    only take matching organism that exceed the threshold
    create all mentions of organisms
    """
    one_dict = {'9606': 'h', '10090': 'm', '10116': 'r'}
    three_dict = {'9606': 'Hum', '10090': 'm', '10116': 'r'}
    found_dict = dict()
    for doc, v in new_dict.iteritems():
        # full_dict = {getPureform(item[2]): item[3] for item in pred_org_id_dict[doc] if item[3] != '' and item[-1] >= threshold]
        if doc in pred_org_id_dict.keys():# and full_dict != {}:
            short_dict = {}
            full_dict = {}
            for item in pred_org_id_dict[doc]:
                
                if item[-1] >= threshold:
                    full_dict.setdefault(getPureform(item[2]), item[3])
                    if item[3] in sci_name.keys():
                        sci_split = sci_name[item[3]].split()
                        if len(sci_split) >= 2 and sci_split[1] != 'sp.':
                            full_dict.setdefault(sci_split[0] + ' ' + sci_split[1], item[3])
                            full_dict.setdefault(sci_split[0][0] + '. ' + sci_split[1], item[3])
                            short_dict.setdefault(sci_split[0][0] + sci_split[1][0:2], item[3])
                            short_dict.setdefault(sci_split[0][0] + sci_split[1][0], item[3])
                            if item[3] in one_dict.keys():
                                short_dict.setdefault(one_dict[item[3]], item[3])
                                short_dict.setdefault(three_dict[item[3]], item[3])
            full_names = "|".join([(item) for item in full_dict.keys() if item[0] != '['])
            short_names = "|".join([(item) for item in short_dict.keys() if item[0] != '['])
            for i, annos in enumerate(v):
                try:
                    new_dict[doc][i][4], new_dict[doc][i][6] = found_dict[annos[2]]
                except:
                    full_forms = annos[6]
                    for anno in full_forms:
                        try:
                            a = re.split(r'\b({})\b'.format(full_names), anno, re.IGNORECASE)
                        except:
                            a = [anno]
                        new_anno = [item.strip() for item in a if item not in full_dict.keys()]
                        spp_text = [item for item in a if item in full_dict.keys()]
                        if spp_text != [] and full_names != '':
                            new_ent = ' '.join(new_anno)
                            spp = spp_text[0]
                            if full_dict[spp] not in new_dict[doc][i][4]: 
                                new_dict[doc][i][4].append(full_dict[spp])
                            if new_ent not in new_dict[doc][i][6]: 
                                new_dict[doc][i][6].append(new_ent.strip())
                        new_anno = re.split(r'^({})(?=[A-Z])'.format(short_names), anno)
                        if len(new_anno) > 1 and short_names != '':
                            new_ent = ''.join(new_anno[2:])
                            spp = new_anno[1]
                            if short_dict[spp] not in new_dict[doc][i][4]: 
                                new_dict[doc][i][4].append(short_dict[spp])
                            if new_ent not in new_dict[doc][i][6]:
                                new_dict[doc][i][6].append(new_ent.strip())
                    found_dict.setdefault(annos[2], [new_dict[doc][i][4], new_dict[doc][i][6]])
    
    return new_dict


def solr_search(solr_obj, symbol, taxid, gene_set, cond, field, order):
    entrezgene_id = ""
    # params = {'fl': "{}".format(field), 'sort': "{}".format(order), 'wt': "csv", 'rows': "20"}
    params = {'fl': "{}".format(field), 'sort': "{}".format(order), 'wt': "csv", 'rows': "20"}
    # print symbol, cond, taxid
    results = solr_obj.search("symbol:{}{} AND ncbitax_id:({})".format(symbol, cond, taxid), **params)
    match = []
    if results:
        if field == 'entrezgene_id':
            match = [result[field] for result in results if field in result.keys()]
        else:
            match = [result[field][0].decode('utf-8').encode('ascii') for result in results if field in result.keys()]
    try:
        common = gene_set & set(match)
        if common != set([]):
            match = list(common)
        entrezgene_id = str(random.choice(match))
    except:
        entrezgene_id = ""
    return entrezgene_id


def map_normform_entrez(pred_ggp_train, cond, field, order, org_under_spp, solr_obj):    
    origin = copy.deepcopy(pred_ggp_train)
    gene_dict = {}
    for doc, v in pred_ggp_train.iteritems():
        for i, anno in enumerate(v):
            gene_set = set([])
            j = 0
            res = ''            
            anno[5].sort(key=lambda x: len(x))
            combine = anno[5][::-1]
            while res == '' and j < len(combine) and anno[3] == '':
                for spp in anno[4]:
                    try:
                        res = gene_dict[spp][combine[j]]
                    except:
                        res = solr_search(solr_obj, combine[j], spp, gene_set, cond, field, order)
                        gene_dict.setdefault(spp, dict()).setdefault(combine[j], res)
                        # res = check_DB_mapping(oriform, combine[j], spp, gene_set)
                    if res == '':
                        try:
                            new_spp = ' '.join(org_under_spp[spp])
                            res = solr_search(solr_obj, combine[j], new_spp, gene_set, cond, field, order)
                            gene_dict.setdefault(spp, dict()).setdefault(combine[j], res)
                        except:
                            pass
                    else:
                        gene_set.add(res)
                        origin[doc][i][3] = res
                        origin[doc][i][4] = spp
                        j = len(combine)
                        break
                j += 1
    return origin


def map_organism_ncbi(pred_ggp_train, pred_org_id_train, pred_dis_train, pred_cel_id_train, offset_dict):
    origin = copy.deepcopy(pred_ggp_train)
    for doc, v in pred_ggp_train.iteritems():
        combined_list = combine_disease_organism_cellline(doc, pred_dis_train, pred_cel_id_train, pred_org_id_train)
        for j, anno in enumerate(v):
            check_spp = []
            if anno[4] == []:                
                sent = [offs for offs in offset_dict[doc] if int(offs[0]) <= int(anno[0]) and int(offs[1]) >= int(anno[1])][0]
                spp_list = []
                spp_offI = []
                spp_offB = []
                spp_offA = []
                spp_offP = []
                for item in combined_list:
                    spp_offA.append(item[3])
                    if int(item[0]) >= int(sent[0]) and int(item[1]) <= int(sent[1]):
                        if int(item[0]) >= int(anno[0]) and int(item[1]) >= int(anno[1]): # same sentence, inside entity
                            spp_offI.append(item[3])
                        elif int(item[1]) <= int(anno[0]): # same sentence, before entity
                            spp_offB.append(item[3])
                    elif int(item[1]) < int(sent[0]):
                        spp_offP.append(item[3])
                for item in spp_offI[::-1] + spp_offB[::-1]:#, (spp_offE, 0),]:
                    spp_list.append(item)
                check_spp += [spp_id for spp_id in spp_list if spp_id not in check_spp and spp_id != '']
                check_spp += [str(org_id[0]) for org_id in Counter(spp_offP).most_common() if org_id[0] not in check_spp and org_id[0] != '']
                check_spp += [spp_id for spp_id in spp_offP[::-1] if spp_id not in check_spp and spp_id != '']
                check_spp += [str(org_id[0]) for org_id in Counter(spp_offA).most_common() if org_id[0] not in check_spp and org_id[0] != '']
                check_spp += [spp_id for spp_id in spp_offA[::-1] if spp_id not in check_spp and spp_id != '']                        
                origin[doc][j][4] = check_spp
                if origin[doc][j][4] == [] or origin[doc][j][4] == ['']:
                    origin[doc][j][4] = ['9606']
            all_nf = list(set([getNormform(item) for item in set(origin[doc][j][6])]))
            origin[doc][j][5] = all_nf
            if '' in origin[doc][j][4]:
                origin[doc][j][4].remove('')
            if '' in origin[doc][j][5]:
                origin[doc][j][5].remove('')
    return origin


def normalize_ggp(offset_dict, norm_cel_dict, norm_org_dict, norm_dis_dict, ggp_dict, canon_dict, abbr_dict, text_dict, all_rank_count, solr_obj):
    # acc_org = {}
    data_folder = os.getcwd() + '/'
    org_under_spp_prot = pickle.load(gzip.open(data_folder + 'pickle/org_under_spp_protein.pkl.gz', 'rb'))
    org_under_spp_gene = pickle.load(gzip.open(data_folder + 'pickle/org_under_spp_gene.pkl.gz', 'rb'))
    sci_name_dict = pickle.load(gzip.open(data_folder + "pickle/sci_name_gene.pkl.gz", "rb"))
    with gzip.open(data_folder + "map_files/new_acc_ncbitax_id.txt.gz", "rb") as f:
        acc_org = {line.strip('\n').split('\t')[0]: line.strip('\n').split('\t')[1] for line in f}
    norm_cel_org = match_cel_org(norm_cel_dict, acc_org)
    ext_norm_cel = extract_entity_abstract(data_folder, text_dict, norm_cel_org, 1.0)
    
    # normalize tissue
    umls_db_dis = 'dis_UMLS/dis_UMLS.db'
    dis_len_text = 3
    norm_dis_org = modify_match_tis_normform(data_folder, norm_dis_dict, umls_db_dis, dis_len_text, 1.0)
    ext_norm_dis = extract_entity_abstract(data_folder, text_dict, norm_dis_org, 1.0)
    
    # normalize organism
    ext_org_dict = copy.deepcopy(norm_org_dict)
    ext_norm_org = extract_entity_abstract(data_folder, text_dict, ext_org_dict, 1.0)

    # print sci_name_dict
    ori_ggp_dict = turn_gold_pred(ggp_dict)
    abb_ggp_dict = add_abbreviation(abbr_dict, ori_ggp_dict)
    can_ggp_dict = add_canonical(canon_dict, abb_ggp_dict)
    org_ggp_dict = extract_org_ggp(can_ggp_dict, norm_org_dict, sci_name_dict, 1.0)
    pred_ggp_dict = modify_match_ggp_normform(data_folder, org_ggp_dict)
    norm_ggp_dict = map_organism_ncbi(pred_ggp_dict, ext_norm_org, ext_norm_dis, ext_norm_cel, offset_dict)
    
    norm_ggp_dict_1 = map_normform_entrez(norm_ggp_dict, "", 'entrezgene_id', 'type ASC,entrezgene_id DESC', org_under_spp_gene, solr_obj)
    norm_ggp_dict_2 = map_normform_entrez(norm_ggp_dict_1, "", 'PID', 'type ASC,entrezgene_id ASC', org_under_spp_prot, solr_obj)
    
    norm_ggp_dict_3 = {}
    for k, v in norm_ggp_dict_2.iteritems():
        norm_ggp_dict_3.setdefault(k, [])
        for item in v:            
            if item[3] == '':
                norm_ggp_dict_3[k].append(item[0:3] + ['ggp', 'protein:' + item[2]])
            else:
                if item[3].isdigit() == True:
                    norm_ggp_dict_3[k].append(item[0:3] + ['ggp', 'NCBI gene:' + item[3]])    
                else:
                    norm_ggp_dict_3[k].append(item[0:3] + ['ggp', 'Uniprot:' + item[3]])
    return norm_ggp_dict_3
