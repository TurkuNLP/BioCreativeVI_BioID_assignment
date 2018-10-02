#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import argparse
import simstring
import subprocess
import cPickle as pickle



def runCommand(cmd_list):
    try:
        subprocess.check_call(cmd_list, shell=True)
    except subprocess.CalledProcessError:
        print cmd_list
        pass # handle errors in the called executable
    except OSError:
        print cmd_list
        pass # executable not found

        
def create_ss_db(ss_folder, in_txt):
    runCommand("simstring -b -d {} < {}".format(ss_folder + in_txt.replace(".txt", ".db"), ss_folder + in_txt))

    
def getNormform(synonym):
    """
    """
    return re.sub("[^a-z0-9]", "", synonym.lower()[:100])


def getNormform_space(synonym):
    """
    """
    return re.sub("[^a-z0-9]", " ", synonym.lower())

    
def match_simstring(ss_folder, annoDict, ms, mapping_dict, db_file, normform, min_threshold=0.01):
    distance_matrix = {0: "exact", 1: "dice", 2: "cosine", 3: "jaccard", 4: "overlap"}
    db = simstring.reader(ss_folder + db_file)
    db.measure = ms
    predict_dict = dict()
    for doc_id, items in annoDict.iteritems():
        predict_dict.setdefault(doc_id, [])
        for j, item in enumerate(items):
            predict_dict[doc_id].append([])
            db.threshold = 1.0
            if normform:
                mention = getNormform(item[2])
            else:
                mention = item[2]
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


def argument_parser():
    parser = argparse.ArgumentParser(description="extract title from pubmed central documents")
    parser.add_argument("-a", "--anno_dict", type=str, help="annotation file in Python dictionary format")    
    parser.add_argument("-o", "--out_folder", type=str, help="out folder for prediction")    
    parser.add_argument("-s", "--ss_folder", type=str, help="simstring folder")    
    parser.add_argument("-b", "--db_file", type=str, help="simstring database file name")
    parser.add_argument("-d", "--target_mapping", type=str, help="dictionary of what to be mapped with")
    parser.add_argument("-m", "--ms", default="cosine", type=str, help="similarity measure (exact, dice, cosine, jaccard, overlap)")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    with open(args.anno_dict, "rb") as f:
        annoDict = pickle.load(f)
    with open(args.target_mapping, "rb") as f:
        mapping_dict = pickle.load(f)
    print mapping_dict.keys()[9]
    predict_dict = match_simstring(args.ss_folder, annoDict, ms, mapping_dict, args.db_file)
    with open(out_folder + "predict_{}_{}.pkl".format(ms, db.split('db')[0]), "wb") as f:
        pickle.dump(predict_dict, f)
    
