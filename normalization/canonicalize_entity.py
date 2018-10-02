#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import re
import copy
import os
import sys
import strip_contexts_v3 as sc3


def canonicalize_ggp_org(anno_dict):
    cwd = os.getcwd()
    # canon_path = cwd + '/canonicalization/'
    # sys.path.append(canon_path)
    # os.chdir(canon_path)
    symbDBs, orgDB = sc3.init()

    canon_dict = copy.deepcopy(anno_dict)
    for k, v in anno_dict.iteritems():
        for i, item in enumerate(v):
            ggp_list = set([])
            org_list = set([])
            (pref,symb,suff,(B,E)), matchType = sc3.stripSymbolLIB(item[2])
            if matchType=="contains":
                B=0
                E=len(item[2]) #contains means we didn't get stripped deep enough
            elif matchType=="OVERLAP" and symb=="xxxorgxxx":
                B=0
                E=len(item[2])
            canon = item[2][B:E].strip()
            ggp_list = set([item[2][b:bX].strip() for (b,bX,e),dbName in sc3.simstringMatches(item[2],symbDBs).iteritems()])
            org_list = set([item[2][b:bX].strip() for (b,bX,e) in sc3.simstringMatches(item[2],(("-org-",orgDB),))])
            ggp = list(ggp_list - set([canon]))
            ggp.sort(key=lambda x: len(x), reverse=True)
            org = list(org_list)
            ggp.insert(0, canon)
            canon_dict[k][i].append('')
            canon_dict[k][i].append(ggp)
    os.chdir(cwd)
    print(canon_dict['5048354_Figure_7-C-D'])
    return canon_dict

