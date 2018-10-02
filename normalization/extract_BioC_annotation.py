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


def write_no_anno(filename, gold):
    new_doc = []
    for lines in all_lines:
        if lines == '':
            new_doc.append(lines)
        else:
            if gold == 'begin':
                new_doc.append('O\t' + lines)
            elif gold == 'last':
                new_doc.append(lines + '\tO')
            else:
                new_doc.append(lines)
    doc_text = "\n".join(new_doc)
    return doc_text


def adjust_offsets(annos):
    new_anno = []
    for offset, text, ent, eid in annos:
        new_off = []
        i = 0
        while i < len(offset):
            new_off += range(offset[i], offset[i+1])
            i += 2
        new_anno.append([new_off, text, ent, eid])
    return new_anno


def add_annotation(all_lines, filename, term_dict, doc_id, entity, unicode_ascii):
    b_tag = 'B-{}'
    i_tag = 'I-{}'    
    new_doc = []
    gold_doc = []
    ner_doc = []
    checking = []
    # print [item for item in term_dict.keys() if '3569655' in item]
    try:
        annos = term_dict[doc_id.strip('.xml.per')]        
        annos = sorted(annos, key=lambda x:x[0])
        new_anno = adjust_offsets(annos)
    except:
        annos = []
    if annos == []:
        doc_text = write_no_anno(all_lines, 'last')
        gold_text = write_no_anno(all_lines, 'begin')
        ner_text = write_no_anno(all_lines, 'no_label')
    else:
        match_ent = []
        done = True
        takeB = False
        for lidx, lines in enumerate(all_lines):
            if lines != '':
                # tokens = lines.split('\t')
                # token = convert_html_ascii(tokens[2], unicode_ascii)
                # lines = '\t'.join(tokens[0:2] + [token] + tokens[3:])
                new_line = lines + '\tO'
                gold_line = 'O\t' + lines
                line = lines.split('\t')
                off_range = set(range(int(line[0]), int(line[1])))
                if done == True:
                    match_ent = [anno for anno in new_anno if off_range  & set(anno[0]) != set([])]
                    match_ent.sort(key=lambda x:len(x[0]))
                    string = ''
                if match_ent != []:
                    ent_type = entity[match_ent[-1][2]].lower()
                    if takeB == False:
                        new_line = lines + '\t' + b_tag.format(ent_type)
                        gold_line = b_tag.format(ent_type) + '\t' + lines
                        string += line[0]
                        takeB = True
                    else:
                        new_line = lines + '\t' +  i_tag.format(ent_type)
                        gold_line = i_tag.format(ent_type) + '\t' + lines
                    if int(line[1]) > match_ent[-1][0][-1]:
                        string += '\t' + line[1]
                        done == True
                        takeB = False
                        checking.append(string)
                ner_doc.append(lines)
                new_doc.append(new_line)
                gold_doc.append(gold_line)
            else:
                ner_doc.append(lines)
                new_doc.append(lines)
                gold_doc.append(lines)
        doc_text = "\n".join(new_doc)
        gold_text = "\n".join(gold_doc)
        ner_text = "\n".join(ner_doc)
    return doc_text, gold_text, ner_text

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
    if s == ' ':
        return '-'
    return s


def extract_annotation(in_xml, entity, unicode_ascii):
    """
        <infon key="type">gene:HSP70</infon>
        <infon key="sourcedata_figure_annot_id">4</infon>
        <infon key="sourcedata_article_annot_id">275</infon>
        <location offset="134" length="5"/>
        <text>HSP70</text>
    """
    doc_dict = {}
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
                # print '\t'.join([doc_id, offB, offE, string, db_name, anno_id, fig_anno_id, doc_anno_id])
                string = string.strip()
                if db_name in entity.keys():
                    if anno_id != string:
                        doc_dict.setdefault(doc_id, []).append([[offB, offE], string, db_name, anno_id])
                    else:
                        doc_dict.setdefault(doc_id, []).append([[offB, offE], string, db_name, '0'])
    return doc_dict
             
                
def argument_parser():
    parser = argparse.ArgumentParser(description="extract entities from BioC XML")
    parser.add_argument("-i", "--in_folder", type=str, default= '/media/suwisa/Biocreative/BioIDtraining_2/caption_bioc/', help="folder where annotated xml resides")
    parser.add_argument("-o", "--out_folder", type=str, default= '/media/suwisa/Biocreative/result_2/pre_train/', help="folder where we want to add annotation")
    parser.add_argument("-f", "--filename", type=str, default= '/media/suwisa/Biocreative/final/devel.txt.gtag', help="gtag file that already fixed offset")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    entity = {'Uniprot': 'ggp', 'protein': 'ggp', 'NCBI gene': 'ggp',
              u'protein': 'ggp', u'gene': 'ggp', # 'Rfam': 'ggp', 'Corum': 'ggp',
              'CL': 'cell', 'cellline': 'cell', 'cell': 'cell',
              'molecule': 'chem', 'CHEBI': 'chem', 'PubChem': 'chem',
              'Uberon': 'tiss', 'tissue': 'tiss',
              'GO': 'subc', 'subcellular': 'subc',
              'NCBI taxon': 'org', 'organism': 'org'}
    with open('entities.dat', 'rb') as f:
        unicode_ascii = read_mapping(f)
    args = argument_parser()
    doc_dict = {}
    text_dict = {}
    for in_xml in glob.glob(args.in_folder + '*'):
        doc_id = in_xml.split('/')[-1]        
        anno_dict = extract_annotation(in_xml, entity, unicode_ascii)
        doc_dict.update(anno_dict)
    new_doc = ''
    new_gold = ''
    new_ner = ''
    
    with open(args.filename, 'rb') as f:
        all_docs = f.read().split('[SEP]')
    for doc in all_docs[1:]:
        lines = doc.split('\n')
        doc_id = lines[0].split('.xml.per')[0]
        
        if '3569655' in doc_id:
            
            print doc_id
        all_lines = lines[2:]
        
        doc_text, gold_text, ner_text = add_annotation(all_lines, args.filename, doc_dict, doc_id, entity, unicode_ascii)
        
        new_doc += '[SEP]' + doc_id + '\n\n' + doc_text
        new_gold += '[SEP]' + doc_id + '\n\n' + gold_text
        new_ner += '[SEP]' + doc_id + '\n\n' + ner_text        
        
    with open(args.filename.replace('.gtag', '.ann'), 'wb') as a:
        a.write(new_doc.replace('\n\n\n', '\n\n'))
    with open(args.filename.replace('.gtag', '.mod'), 'wb') as a:
        a.write(new_gold.replace('\n\n\n', '\n\n'))
    with open(args.filename.replace('.gtag', '.ner'), 'wb') as a:
        a.write(new_ner.replace('\n\n\n', '\n\n'))
        
    with open(args.filename.replace('.gtag', '.cell'), 'wb') as a:
        a.write(new_gold.replace('\n\n\n', '\n\n').replace('B-tiss', 'O').replace('I-tiss', 'O').replace('B-subc', 'O').replace('I-subc', 'O').replace('B-ggp', 'O').replace('I-ggp', 'O').replace('B-org', 'O').replace('I-org', 'O').replace('B-chem', 'O').replace('I-chem', 'O'))

    with open(args.filename.replace('.gtag', '.subc'), 'wb') as a:
        a.write(new_gold.replace('\n\n\n', '\n\n').replace('B-cell', 'O').replace('I-cell', 'O').replace('B-tiss', 'O').replace('I-tiss', 'O').replace('B-ggp', 'O').replace('I-ggp', 'O').replace('B-org', 'O').replace('I-org', 'O').replace('B-chem', 'O').replace('I-chem', 'O'))

    with open(args.filename.replace('.gtag', '.ggp'), 'wb') as a:
        a.write(new_gold.replace('\n\n\n', '\n\n').replace('B-cell', 'O').replace('I-cell', 'O').replace('B-subc', 'O').replace('I-subc', 'O').replace('B-tiss', 'O').replace('I-tiss', 'O').replace('B-org', 'O').replace('I-org', 'O').replace('B-chem', 'O').replace('I-chem', 'O'))

    with open(args.filename.replace('.gtag', '.org'), 'wb') as a:
        a.write(new_gold.replace('\n\n\n', '\n\n').replace('B-cell', 'O').replace('I-cell', 'O').replace('B-subc', 'O').replace('I-subc', 'O').replace('B-ggp', 'O').replace('I-ggp', 'O').replace('B-tiss', 'O').replace('I-tiss', 'O').replace('B-chem', 'O').replace('I-chem', 'O'))

    with open(args.filename.replace('.gtag', '.chem'), 'wb') as a:
        a.write(new_gold.replace('\n\n\n', '\n\n').replace('B-cell', 'O').replace('I-cell', 'O').replace('B-subc', 'O').replace('I-subc', 'O').replace('B-ggp', 'O').replace('I-ggp', 'O').replace('B-org', 'O').replace('I-org', 'O').replace('B-tiss', 'O').replace('I-tiss', 'O'))

    with open(args.filename.replace('.gtag', '.tiss'), 'wb') as a:
        a.write(new_gold.replace('\n\n\n', '\n\n').replace('B-cell', 'O').replace('I-cell', 'O').replace('B-subc', 'O').replace('I-subc', 'O').replace('B-ggp', 'O').replace('I-ggp', 'O').replace('B-org', 'O').replace('I-org', 'O').replace('B-chem', 'O').replace('I-chem', 'O'))




