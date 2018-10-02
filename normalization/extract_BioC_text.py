#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import htmlentitydefs
import re
# import xml.etree.cElementTree as ET
import xml.etree.ElementTree as ET # use this one for the ET.tostring
import glob


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

def mapchar(c, mapping):
    if c in mapping:
        return mapping[c]
    else:
        # make a note of anything unmapped
        global missing_mappings, options
        missing_mappings.add("%.4X" % ord(c))
        # remove missing by default, output codepoint as hex as an option
        if not options.hex:
            return ''
        else:
            return "<%.4X>" % ord(c)
        
def replace_mapped(root, mapping):
    # TODO: inefficient, improve
    s = ""
    o = ""
    u_count = 0
    for i, c in enumerate(root):
        o += '\t'.join([str(i+(u_count*4)), str(len(s))])  + '\n'
        if ord(c) >= 128:
            new_char = mapchar(c, mapping)
            s += new_char
            u_count += 1
        else:
            s += c
    return s, o

def convert_html_entities(old_s, unicode_ascii):
    temp ={u'\x93': '"', u'\x94': '"', '&gt;': '>', '&lt;': '<', u'\xa0': '_'}
    map_dict = {}
    s = old_s
    map_idx = ''
    matches = re.findall("&#\d+;", s)
    if len(matches) > 0:
        hits = set(matches)
        for hit in hits:
            name = hit[2:-1]
            try:
                entnum = int(name)
                a = unicode_ascii[unichr(entnum)]
                map_dict.setdefault(hit, [a, len(hit)])
                s = s.replace(hit, a)
            except KeyError:
                a = temp[unichr(entnum)]
                map_dict.setdefault(hit, [a, len(hit)])
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
                map_dict.setdefault(hit, [a, len(hit)])
                s = s.replace(hit, a)
            except KeyError:
                a = temp[unichr(entnum)]
                map_dict.setdefault(hit, [a, len(hit)])
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
            map_dict.setdefault(hit, [a, len(hit)])
            s = s.replace(hit, a)
    a = '&'
    map_dict.setdefault(amp, [a, len(amp)])
    s = s.replace(amp, a)
    all_list = {}
    for k, v in map_dict.iteritems():
        html_list = {html_s.start(): [html_s.end(), len(v[0]), v[1]] for html_s in re.finditer(re.compile(k), old_s)} #html_s.end(), html_s.expand(), html_s.group(), html_s.groupdict(), html_s.groups(), html_s.span(), html_s.start()
        all_list.update(html_list)
    for k, v in all_list.iteritems():
        map_idx += '\t'.join([str(k), str(v[0]), str(v[1]), str(v[2])]) + '\n'
    return s, map_idx

def temp_get_doctext(in_xml, unicode_ascii):
    # offset_dict = {}
    text_dict = {}
    tree = ET.parse(open(in_xml))
    docs = tree.getroot() # get root element
    # doc_xml = in_xml.split('/')[-1].split('.xml')[0]
    for doc in docs.findall('.//document'):
        doc_id = doc.findall('id')[0].text.replace(" ", "_").strip('.')
        for passage in doc.findall('passage'):
            # off_text = ET.tostring(passage.find('offset'))            
            # offset = int(re.findall('<offset>(.*?)</offset>', off_text, re.DOTALL)[0])
            try:
                doc_text = ET.tostring(passage.find('text'))                
                ext_text = re.findall('<text>(.*?)</text>', doc_text, re.DOTALL)
                if len(ext_text) == 1:
                    text_dict.setdefault(doc_id, []).append(ext_text[0])
            except:
                pass
    return text_dict


def argument_parser():
    parser = argparse.ArgumentParser(description="extract entities from BioC XML")
    parser.add_argument("-i", "--in_folder", type=str, default= '../input/', help="full path of xml file resides")
    parser.add_argument("-o", "--out_folder", type=str, default= '../temp_process/documents/', help="folder where text files resides")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    text_dict = dict()
    with open('entities.dat', 'rb') as f:
        unicode_ascii = read_mapping(f)
    for in_xml in glob.glob(args.in_folder + '*.xml'):
        doc_id = in_xml.split('/')[-1].strip('.xml')        
        text_dict = temp_get_doctext(in_xml, unicode_ascii)
        for doc_id, all_text in text_dict.iteritems():
            text = '\n'.join(all_text)
            doc_text, map_idx = convert_html_entities(text, unicode_ascii)
            doc_id_text, map_id_idx = convert_html_entities(doc_id, unicode_ascii)
            with open(args.out_folder + doc_id_text + '.xml', 'w') as f:
                f.write(text + '\n')
            # with open(args.out_folder + doc_id_text + '.txt', 'w') as f:
                # f.write(doc_text + '\n')
            with open(args.out_folder + doc_id_text + '.xml.idx', 'w') as f:
                f.write(map_idx + '\n')
