#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import codecs
import re
import glob
# import xml.etree.cElementTree as ET
import xml.etree.ElementTree as ET
import argparse
import htmlentitydefs
import normalize_entity as norm_ent
import process_Ab3P as ab3p 
import cPickle as pickle
import time



            
def convert_html_ascii(old_s, unicode_ascii):
    temp ={u'\x93': '"', u'\x94': '"', '&gt;': '>', '&lt;': '<'}
    s = old_s
    matches = re.findall("&#\d+;", s)
    if len(matches) > 0:
        hits = set(matches)
        for hit in hits:
            name = hit[2:-1]
            try:
                entnum = int(name)
                a = unicode_ascii[unichr(entnum)]
                s = s.replace(hit, a)
            except KeyError:
                entnum = int(name)
                a = temp[unichr(entnum)]
                s = s.replace(hit, a)
            except ValueError:
                pass        
    matches = re.findall("&#[xX][0-9a-fA-F]+;", s)
    if len(matches) > 0:
        hits = set(matches)
        for hit in hits:
            hex = hit[3:-1]
            try:
                entnum = int(hex, 16)
                a = unicode_ascii[unichr(entnum)]
                s = s.replace(hit, a)
            except KeyError:
                entnum = int(hex, 16)
                a = temp[unichr(entnum)]
                s = s.replace(hit, a)
            except ValueError:
                pass        
    matches = re.findall("&\w+;", s)
    hits = set(matches)
    amp = "&amp;"
    if amp in hits:
        hits.remove(amp)
    for hit in hits:
        name = hit[1:-1]
        if htmlentitydefs.name2codepoint.has_key(name):
            a = unicode_ascii[unichr(htmlentitydefs.name2codepoint[name])]
            s = s.replace(hit, a)
    a = '&'
    s = s.replace(amp, a)
    return s


def read_mapping(f, fn="mapping data"):
    """
    Reads in mapping from Unicode to ASCII from the given input stream
    and returns a dictionary keyed by Unicode characters with the
    corresponding ASCII characters as values. The expected mapping
    format defines a single mapping per line, each with the format
    CODE\tASC where CODE is the Unicode code point as a hex number and
    ASC is the replacement ASCII string ("\t" is the literal tab
    character). Any lines beginning with "#" are skipped as comments.
    """
    # read in the replacement data
    linere = re.compile(r'^([0-9A-Za-z]{4,})\t(.*)$')
    mapping = {}
    for i, l in enumerate(f):
        # ignore lines starting with "#" as comments
        if len(l) != 0 and l[0] == "#":
            continue
        m = linere.match(l)
        assert m, "Format error in %s line %s: '%s'" % (fn, i+1, l.replace("\n","").encode("utf-8"))
        c, r = m.groups()
        c = unichr(int(c, 16))
        assert c not in mapping or mapping[c] == r, "ERROR: conflicting mappings for %.4X: '%s' and '%s'" % (ord(c), mapping[c], r)
        # exception: literal '\n' maps to newline
        if r == '\\n':
            r = '\n'
        mapping[c] = r
    return mapping


def extract_annotation(in_xml, unicode_ascii):
    entity = {'Uniprot': 'ggp', 'protein': 'ggp', 'NCBI gene': 'ggp',
              u'protein': 'ggp', u'gene': 'ggp', 'Rfam': 'ggp',
              'CL': 'cell', 'cellline': 'cell', 'cell': 'cell',
              'molecule': 'chem', 'CHEBI': 'chem', 'PubChem': 'chem',
              'Uberon': 'tiss', 'tissue': 'tiss',
              'GO': 'subc', 'subcellular': 'subc',
              'NCBI taxon': 'org', 'organism': 'org'}    
    anno_dict = {}
    ori_dict = {}
    tree = ET.parse(open(in_xml))
    docs = tree.getroot() # get root element
    for doc in docs.findall('.//document'):        
        doc_id = doc.findall('id')[0].text.replace(' ', '_').strip('.')
        for passage in doc.findall('passage'):
            for anno in passage.findall('annotation'):
                string = ''
                for infon in anno.findall('infon'):
                    if infon.attrib['key'] == 'type':
                        anno_type = infon.text.split('|')
                        anno_type_split = anno_type[0].split(':', 1)
                        if len(anno_type_split) != 2:
                            db_name, anno_id = 'cellline', anno_type
                        else:
                            db_name, anno_id = anno_type_split
                for texts in anno.findall('text'):
                    string += texts.text + " "
                for location in anno.findall('location'):
                    offB = int(location.attrib['offset'])
                    offE = offB + int(location.attrib['length'])
                string = string.strip()
                if db_name in entity.keys():
                    if anno_id != string:
                        anno_dict.setdefault(doc_id, []).append([offB, offE, convert_html_ascii(string, unicode_ascii), entity[db_name], anno_id])
                        ori_dict.setdefault(doc_id, []).append([offB, offE, string, entity[db_name], anno_id])
                    else:
                        anno_dict.setdefault(doc_id, []).append([offB, offE, convert_html_ascii(string, unicode_ascii), entity[db_name], anno_id])
                        ori_dict.setdefault(doc_id, []).append([offB, offE, string, entity[db_name], anno_id])
    return ori_dict, anno_dict



def all_entity_extract(xml_folder):
    cwd = os.getcwd()
    with open(os.path.join(cwd, 'entities.dat'), 'rb') as f:
        unicode_ascii = read_mapping(f)
    ascii_dict = {}
    xml_dict = {}
    for in_xml in glob.glob(xml_folder + '*'):
        doc_id = in_xml.split('/')[-1]        
        ori_dict, anno_dict = extract_annotation(in_xml, unicode_ascii)
        
        ascii_dict.update(anno_dict)
        xml_dict.update(ori_dict)
    return xml_dict, ascii_dict



xml_dict_curr, ascii_dict_curr = all_entity_extract('/home/sukaew/Biocreative_test/test_input/output/')

xml_dict_subm, ascii_dict_subm = all_entity_extract('/home/sukaew/Biocreative_database/BioIDtest_with_NN/')


for k, v in xml_dict_subm.iteritems():
    for i, item in enumerate(v):
        anno = xml_dict_curr[k][i]
        if item != anno:
            print item
            print anno
            print ''

            
