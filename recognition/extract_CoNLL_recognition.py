#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import glob
import xml.etree.cElementTree as ET
import argparse
import htmlentitydefs


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
                    b_type = lines[-1].split('-')[-1]
                    if new_ent != '':
                        result = rewrite_entity(new_ent)
                        ent_list.append(result)
                    new_ent = off_file[idx].split('\n')[0] + '\t' + lines[-1].replace('I-', 'B-') + '\n'
                    
    if new_ent != '':
        result = rewrite_entity(new_ent) # result = [off_B, off_E, txt_str, entity]
        ent_list.append(result)
    ent_list = list(ent_list)
    return ent_list

def write_annotation(in_xml, norm_ent_dict, out_folder, unicode_ascii):
    ent_dict = {'cell': 'cell:', 'chem': 'molecule:', 'tiss': 'tissue:', 'ggp': 'protein:', 'subc': 'subcellular:', 'org': 'organism:'}
    root = ET.parse(open(in_xml))
    docs = root.getroot() # get root element
    count = 0
    for doc in docs.findall('.//document'):
        old_doc_id = doc.findall('id')[0].text.replace(' ', '_')
        doc_id = convert_html_ascii(old_doc_id, unicode_ascii)

        if '3988959' in doc_id:
            doc_id = doc_id.encode("utf8")
        if doc_id in norm_ent_dict.keys():
            ent_list = norm_ent_dict[doc_id]
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
    # root = ET.ElementTree(docs)
    root.write(out_folder + in_xml.split('/')[-1])
    return count

def write_result(xml_folder, out_folder, norm_ent_dict, unicode_ascii):
    doc_set = set([item.split('_')[0] for item in ent_dict.keys() if item != ''])    
    for in_xml in glob.glob(xml_folder + '*.xml'):
        if in_xml.split('/')[-1].split('.xml')[0] in doc_set:
            write_annotation(in_xml, norm_ent_dict, out_folder, unicode_ascii)

            
def convert_html_ascii(old_s, unicode_ascii):
    temp ={u'\x93': '"', u'\x94': '"', '&gt;': '>', '&lt;': '<', u'\xa0':'_'}
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
    if s == ' ':
        return '-'
    return s

def get_doctext(in_pred, unicode_ascii):
    text_dict = dict()
    for in_xml in glob.glob(in_pred + '*.xml'):
        with open(in_xml, 'rb') as f:            
            old_doc_id = in_xml.split('/')[-1].split('.xml')[0]
            doc_id = convert_html_ascii(old_doc_id, unicode_ascii)
            doc_text = f.read()
            new_text = convert_html_ascii(doc_text, unicode_ascii)
            text_dict.setdefault(doc_id, new_text)
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
        old_doc_id = pred_split[0].split('.xml.per')[0]
        doc_id = convert_html_ascii(old_doc_id, unicode_ascii)
        pred_dict.setdefault(doc_id, pred_split[1:])
        off_dict.setdefault(doc_id, off_split[1:])
    for doc_id, content in pred_dict.iteritems():
        pred_list = content
        off_list = off_dict[doc_id]
        ent_list = load_prediction(pred_list, off_list)
        ent_dict.setdefault(doc_id, ent_list)        
         
        # print pred.split('\n')[1]
        # pred_dict.setdefault(doc_id, 
    # for filesection in all_files[1:]:
    #     doc_id = filesection.split('\n\n')[0].split('_')[0]
    #     ent_list = load_prediction(filesection, off_folder, file_ext)
    #     ent_dict.setdefault(doc_id, ent_list)
    return ent_dict
            
def argument_parser():
    """
    assuming in the released code, 
    fulltext_folder is the same as input folder unless it is specified separately 
    """
    parser = argparse.ArgumentParser(description="inject entity prediction into BioC XML")
    parser.add_argument("-x", "--xml_folder", type=str, default='/home/suwisa/Biocreative_database/BioIDdevel/caption_bioc_unannotated/', help="folder where TEES xml resides")
    parser.add_argument("-o", "--out_folder", type=str, default='/home/suwisa/Biocreative_database/output/', help="folder where output xml will reside")
    parser.add_argument("-i", "--in_preds", type=str, default='/home/suwisa/Biocreative_database/temp_process_BioIDdevel_caption/BioIDdevel.txt.tag', help="prediction file, if multiple files, separate by comma with no space")    
    parser.add_argument("-a", "--off_file", type=str, default='/home/suwisa/Biocreative_database/temp_process_BioIDdevel_caption/BioIDdevel.txt.ori', help='folder where .txt.ori are')
    args = parser.parse_args()
    return args
    
if __name__ == "__main__":
    args = argument_parser()
    doc_set = set([])
    with open('entities.dat', 'rb') as f:
        unicode_ascii = read_mapping(f)
    ent_dict = {}
    in_preds = args.in_preds.split(',')
    for in_pred in in_preds:
        dict_1 = get_entity(in_pred, args.off_file)
        for k, v in dict_1.iteritems():
            ent_dict.setdefault(k, [])
            for item in v:
                ent_dict[k].append(item)
    
    write_result(args.xml_folder, args.out_folder, ent_dict, unicode_ascii)
    
