#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysolr
import os
import codecs
cwd = os.getcwd()
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


def rewrite_entity(new_ent):
    ent_line = new_ent.split("\n")
    offset = []
    tokens = []
    entity = []
    for line in ent_line:
        if line != '':
            col = line.split('\t')
            offset.append([col[0], col[1]])
            tokens.append(col[2])
            entity.append(col[3].split('-')[1])
    off_str = []
    off_str = []
    txt_str = ''
    for i, item in enumerate(offset):
        if i == 0:
            off_str.append([item[0], item[1]])
            txt_str += tokens[i]
        else:
            dist = int(item[0]) - int(off_str[-1][1])
            if dist == 1:
                txt_str += ' ' + tokens[i]
                off_str[-1][1] = item[1]
            elif dist == 0:
                txt_str += tokens[i]
                off_str[-1][1] = item[1]
            else:
                off_str.append([item[0], item[1]])
                txt_str += ' ' + tokens[i]
    off_text = ''
    for item in off_str:
        off_text += " ".join(item) + ";"
    off_B, off_E = off_text.split(";", 1)[0].split(' ')
    return off_B, off_E, txt_str, entity[0]


def load_prediction(pred_file, off_file):
    ent_list = []
    new_ent = ''
    result = ''
    for idx, line in enumerate(pred_file):
        if line != '\n':
            lines = line.split('\t')
            if 'B-' in lines[-1]:
                b_type = lines[-1].split('B-')[-1]
                if new_ent != '': # encounter the non-B and non-I, insert result to ent_list
                    result = rewrite_entity(new_ent)
                    ent_list.append(result)
                new_ent = off_file[idx] + '\t' + lines[-1] + '\n'
            elif 'I-' in lines[-1]:
                i_type = lines[-1].split('I-')[-1]                
                try: 
                    if i_type != b_type: # type of entity in I-tag is not the same as B-tag, insert the B-tag to the ent_list
                        result = rewrite_entity(new_ent)
                        ent_list.append(result)                        
                        new_ent += off_file[idx] + '\t' + lines[-1].replace('I-', 'B-') + '\n'# then start new entity with B-tag
                    else: # I-tag continue with the previous B-tag
                        new_ent += off_file[idx] + '\t' + lines[-1] + '\n'
                except UnboundLocalError: # this is for the case where the B-tag is missing, we change i-tag to be b-tag
                    b_type = lines[-1].split('B-')[-1]
                    if new_ent != '':
                        result = rewrite_entity(new_ent)
                        ent_list.append(result)
                    new_ent = off_file[idx].split('\n')[0] + '\t' + lines[-1]                    
    if new_ent != '':
        result = rewrite_entity(new_ent)
        # result = [off_B, off_E, txt_str, entity]
        ent_list.append(result)
    ent_list = list(ent_list)
    return ent_list


def write_result(xml_folder, out_folder, norm_ent_dict, unicode_ascii, unanno_folder):
    doc_set = set([item.split('_')[0] for item in norm_ent_dict.keys() if item != ''])
    for in_xml in glob.glob(xml_folder + '*.xml'):
        if in_xml.split('/')[-1].split('.xml')[0] in doc_set:
            write_annotation(in_xml, norm_ent_dict, out_folder, unicode_ascii, unanno_folder)

            
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


def get_doctext(text_folder, unicode_ascii):
    text_dict = dict()
    for in_xml in glob.glob(text_folder + '*.xml'):
        with open(in_xml, 'rb') as f:
            doc_id = in_xml.split('/')[-1].split('.xml')[0]
            doc_text = f.read()            
            new_text = convert_html_ascii(doc_text, unicode_ascii)
            text_dict.setdefault(doc_id, new_text)
    return text_dict


def get_fulltext(text_folder, unicode_ascii):
    # offset_dict = {}
    text_dict = {}
    for in_xml in glob.glob(text_folder + '*.xml'):
        tree = ET.parse(open(in_xml))
        docs = tree.getroot() # get root element
        for doc in docs.findall('.//document'):
            doc_id = doc.findall('id')[0].text.replace(" ", "_").strip('.')
            text_list = []
            for passage in doc.findall('passage'):
                for text in passage.findall('text'):
                    doc_text = ET.tostring(text)                
                    ext_text = re.findall('<text>(.*?)</text>', doc_text, re.DOTALL)
                    if len(ext_text) == 1:
                        text_list.append(ext_text[0])
            text_dict.setdefault(doc_id, convert_html_ascii(' '.join(text_list), unicode_ascii))
    return text_dict
                    

def get_entity(in_pred, off_file):
    pred_dict = {}
    off_dict = {}
    ent_dict = {}
    with open(in_pred, 'rb') as f:
        pred_lines = f.read().split('[SEP]')
    with open(off_file, 'rb') as f:
        off_lines = f.read().split('[SEP]')
    for pred, off in zip(pred_lines, off_lines):
        pred_split = pred.split('\n')
        off_split = off.split('\n')
        doc_id = pred_split[0].split('.xml.per')[0]        
        pred_dict.setdefault(doc_id, pred_split[1:])
        off_dict.setdefault(doc_id, off_split[1:])
    for doc_id, content in pred_dict.iteritems():
        pred_list = content
        off_list = off_dict[doc_id]
        ent_list = load_prediction(pred_list, off_list)
        ent_dict.setdefault(doc_id, ent_list)        
    return ent_dict


def get_offset(in_off):
    off_dict = dict()
    with open(in_off, 'rb') as f:
        docs = f.read().split('[SEP]')
        for doc in docs[1:]:
            sentences = doc.split('\n\n') 
            doc_id = sentences[0].strip('.xml.per')
            for sentence in sentences[1:]:
                tokens = sentence.split('\n')
                if tokens != [''] and tokens != ['', '']:
                    offB = tokens[0].split('\t')[0]
                    offE = tokens[-1].split('\t')[1]                    
                    off_dict.setdefault(doc_id, []).append([int(offB), int(offE)])
    return off_dict


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
                        anno_dict.setdefault(doc_id, []).append([offB, offE, convert_html_ascii(string, unicode_ascii), entity[db_name]])
                        ori_dict.setdefault(doc_id, []).append([offB, offE, string, entity[db_name]])
                    else:
                        anno_dict.setdefault(doc_id, []).append([offB, offE, convert_html_ascii(string, unicode_ascii), entity[db_name]])
                        ori_dict.setdefault(doc_id, []).append([offB, offE, string, entity[db_name]])
    return ori_dict, anno_dict


def write_annotation(in_xml, norm_ent_dict, out_folder, unicode_ascii, unanno_folder):
    ent_dict = {'cell': 'cell:', 'chem': 'molecule:', 'tiss': 'tissue:', 'ggp': 'protein:', 'subc': 'subcellular:', 'org': 'organism:'}
    root = ET.parse(open(unanno_folder + in_xml.split('/')[-1]))
    docs = root.getroot() # get root element
    count = 0
    for doc in docs.findall('.//document'):
        old_doc_id = doc.findall('id')[0].text.replace(' ', '_')
        doc_id = convert_html_ascii(old_doc_id, unicode_ascii)
        if doc_id in norm_ent_dict.keys():
            ent_list = norm_ent_dict[doc_id]
            ent_list.sort()
            for passage in doc.findall('passage'):
                for idx, ent in enumerate(ent_list):
                    anno = ET.Element('annotation')
                    anno.set('id', str(idx +1))                    
                    infon = ET.Element('infon')
                    infon.set('key', 'type')
                    ent_type = ent_dict[ent[3]]
                    try:
                        infon.text = ent[4]
                    except:
                        infon.text = ent_type + ent[2]
                        count += 1
                    anno.append(infon)
                    passage.append(anno)                    
                    location = ET.SubElement(anno, 'location')
                    location.set('offset', str(ent[0]))
                    location.set('length', str(int(ent[1]) - int(ent[0])))
                    texts = ET.SubElement(anno, 'text')
                    texts.text = ent[2]                    
    root = ET.ElementTree(docs)
    # print out_folder + in_xml.split('/')[-1]
    root.write(out_folder + in_xml.split('/')[-1])
    return count


def combine_ent_dict(ent_dict, norm_dict_list):
    norm_ent_dict = {}
    for k, v in ent_dict.iteritems():
        norm_ent_dict.setdefault(k, [])
        for norm_dict in norm_dict_list:            
            try:
                anno_dict = norm_dict[k]
                for anno in anno_dict:
                    if len(anno) == 5:
                        norm_ent_dict[k].append(anno)
                    elif len(anno) == 4:
                        norm_ent_dict[k].append(anno + anno[3])
            except:
                pass
    return norm_ent_dict


def argument_parser():
    """
    assuming in the released code, 
    fulltext_folder is the same as input folder unless it is specified separately 
    """
    parser = argparse.ArgumentParser(description="inject entity prediction into BioC XML")
    parser.add_argument("-x", "--xml_folder", type=str, help="folder where input TEES xml resides")
    parser.add_argument("-o", "--out_folder", type=str, help="folder where output xml will reside")
    parser.add_argument("-a", "--off_file", type=str, help='folder where .txt.ori are')
    parser.add_argument("-t", "--text_folder", type=str, help='folder where .txt.ori are')
    parser.add_argument("-d", "--data_folder", type=str, help='folder where .txt.ori are')    
    parser.add_argument("-f", "--fulltext_folder", type=str, help='full-text folder')    
    parser.add_argument("-p", "--ab3p_folder", type=str, help='folder where Ab3P installed')
    parser.add_argument("-s", "--solr_address", type=str, help='folder where Ab3P installed')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    
    solr_obj = pysolr.Solr(args.solr_address)    
    with open(os.path.join(cwd, 'entities.dat'), 'rb') as f:
        unicode_ascii = read_mapping(f)

    ascii_dict = {}
    xml_dict = {}
    for in_xml in glob.glob(args.xml_folder + '*'):
        
        doc_id = in_xml.split('/')[-1]        
        ori_dict, anno_dict = extract_annotation(in_xml, unicode_ascii)
        
        ascii_dict.update(anno_dict)
        xml_dict.update(ori_dict)

    offset_dict = get_offset(args.off_file)
    
    text_dict = get_doctext(args.text_folder, unicode_ascii)
    
    fulltext_dict = get_fulltext(args.fulltext_folder, unicode_ascii)        
    
    norm_ent_dict = norm_ent.normalize_all(args.data_folder, args.text_folder, args.ab3p_folder, xml_dict, ascii_dict, text_dict, offset_dict, fulltext_dict, solr_obj)
    
    write_result(args.xml_folder, args.out_folder, norm_ent_dict, unicode_ascii, args.xml_folder)
