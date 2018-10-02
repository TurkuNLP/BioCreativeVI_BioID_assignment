#!/usr/bin/python3
# -*- coding: utf-8 -*-


import re
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
    if s == ' ' or s == '\n':
        return '-'
    return s


def fix_gtag_offset(in_folder, filename, all_lines, unicode_ascii):
    byte_dict = {8: 3, 7: 3, 6: 2, 5: 1, 4: 1}
    new_doc = []
    old_doc = []
    mapping = {}
    with open(in_folder + filename.replace('.per', '.idx'), 'r') as f:
        for l1 in f:
            l2 = l1.split('\n')[0].split('\t')
            if l2 != ['']:
                mapping.setdefault(int(l2[0]), [int(l2[1]), int(l2[2]), int(l2[3])])
    if mapping != {}:
        all_keys = mapping.keys()
        j = 1
        new_doc.append(all_lines[0])
        old_doc.append(all_lines[0])
        prev_line = new_doc[-1].split('\t', 2)
        prevE = int(prev_line[1])
        x = 0
        while j < len(all_lines):
            if all_lines[j-1] != '':
                x = j -1
            else:
                x = j -2
            old_line = all_lines[x].split('\t', 2) 
            line = all_lines[j]
            if line != '':
                oldE = int(old_line[1])
                lines = line.split('\t')                
                k = j
                offB = int(lines[0]) 
                offE = int(lines[1])
                gap = offB - oldE
                new_text = lines[2]
                new_offB = int(prevE) + gap
                new_offE = new_offB + len(new_text)
                if offB in all_keys:
                    old_offE, len_text, len_byte = mapping[offB]                    
                    while int(offE) < int(old_offE):
                        next_line = all_lines[k+1].split('\t')
                        new_text += next_line[2]
                        offE = next_line[1]
                        k += 1
                    new_offE = new_offB + byte_dict[len_byte]
                    j = k
                new_doc.append('\t'.join([str(new_offB), str(new_offE), convert_html_ascii(new_text, unicode_ascii)]))
                old_doc.append('\t'.join([str(new_offB), str(new_offE), new_text]))
                prevE = new_offE
            else:
                new_doc.append(line)
                old_doc.append(line)
            j += 1
    else:
        new_doc = all_lines
        old_doc = all_lines
    all_texts = '\n'.join(new_doc)
    old_texts = '\n'.join(old_doc)
    return all_texts, old_texts


def argument_parser():
    parser = argparse.ArgumentParser(description="fix offsets")
    parser.add_argument("-i", "--in_folder", type=str, help="main input folder")
    parser.add_argument("-f", "--filename", type=str, help="tokenized file that is already fixed")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    with open('entities.dat', 'rb') as f:
        unicode_ascii = read_mapping(f)
    new_text = ''
    old_text = ''
    with open(args.filename, 'rb') as f:
        all_docs = f.read().split('[SEP]')
    for doc in all_docs[1:]:
        lines = doc.split('\n')
        doc_id = lines[0]
        all_lines = lines[2:]        
        new_str, old_str = fix_gtag_offset(args.in_folder, doc_id, all_lines, unicode_ascii)
        new_text += '[SEP]' + doc_id + '\n\n' + new_str
        old_text += '[SEP]' + doc_id + '\n\n' + old_str
    with open(args.filename.replace('.tok', '.off'), 'wb') as f:
        f.write(new_text)
    with open(args.filename.replace('.tok', '.ori'), 'wb') as f:
        f.write(old_text)

    
