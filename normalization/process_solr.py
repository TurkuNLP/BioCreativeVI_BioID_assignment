#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import pysolr
import gzip
import argparse
import re
import requests
import cPickle as pickle
import os

stype_order = {'default_symbol': 1, 'Allergen': 2, 'CD_antigen': 3, 'locustag': 4, 'preferred_names': 5, 'names': 6, 'synonym': 7, 'Synonym': 7, 'description': 8, 'Description': 8}


def getNormform(synonym):
    """
    """
    a = re.sub("[^a-z0-9]", "", synonym[0:255].lower())
    return " ".join(a.split())


def get_gene(filename):
    quotes = []
    pred_dict = {}
    done = set([])
    with gzip.open(filename, "rb") as s:
        f = s.readlines()
    for line in f:            
        entrezgene_id, symbol, stype, ncbitax_id = line.strip('\n').split('\t')
        symbol = getNormform(symbol)
        text = '\t'.join([entrezgene_id, symbol, ncbitax_id])
        pred_dict.setdefault('\t'.join([entrezgene_id, symbol, ncbitax_id]), []).append(stype_order[stype])
    f = open(filename.replace('.tsv.gz', '_dict.pkl'), 'wb')
    pickle.dump(quotes, f)
    f.close()
    for j, v in pred_dict.iteritems():
        k = j.split('\t')
        v.sort()
        if args.symbol_type == 'gene' and k[1] != '':
            quotes.append({"entrezgene_id" : int(k[0]), "symbol": k[1], 'ncbitax_id': int(k[2]), 'type': v[0]})
        elif args.symbol_type == 'protein' and k[1] != '':
            quotes.append({"PID" : k[0], "symbol": k[1], 'ncbitax_id': int(k[2]), 'type': v[0]})
    f = open(filename.replace('.tsv.gz', '_list.pkl'), 'wb')
    pickle.dump(quotes, f)
    f.close()
    return quotes


def index_solr(solr_address, quotes, del_record):
            
    solr = pysolr.Solr(solr_address)    
    print("Example quote:", quotes[0])
    print("Indexing quotes...")
    if del_record:
        print del_record
        solr.delete(q='*:*')
    i = 0
    while i < len(quotes):
        stime = time.time()
        solr.add(quotes[i: i + 1000000], commit =True)
        i += 1000000
    solr.add(quotes[i: len(quotes)], commit =True)
    print("%d quotes indexed to Solr." % len(quotes))


def clean_up(quotes, filename):
    print 'before deleteing', len(quotes)
    del_list = set([])
    for i, item in enumerate(quotes):
        for k, v in item.iteritems():
            if v == "":
                del_list.add(i)
    del_idx = list(del_list)
    del_idx.sort(reverse=True)
    for idx in del_idx:
        quotes.pop(idx)
    print 'after deleting', len(quotes)
    f = open(filename.replace('.tsv.gz', '_list.pkl'), 'wb')
    pickle.dump(quotes, f)
    f.close()
    return quotes


def argument_parser():
    parser = argparse.ArgumentParser(description="processing the data for the solr")
    parser.add_argument("-f", "--filename", type=str, default='solr/final_gene.tsv.gz', help="input file to be put into solr")
    parser.add_argument("-s", "--symbol_type", type=str, default='gene', help="gene or protein")
    parser.add_argument("-d", "--del_record", default=True, type=str, help="delete records in solr or not")
    parser.add_argument("-a", "--solr_address", type=str, help="solr core address where gene protein mapping ")    
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    quotes = get_gene(args.filename)
    index_solr(args.solr_address, quotes, args.del_record)
