#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import simstring
from produce_contexts import simstringMatches, best_matches, bestIntervals, normform
import re

tokBoundaryRe=re.compile(r"^[^a-zA-Z0-9]+$")

def longest_prefix(symbol,db):
    matches=list(simstringMatches(symbol,(("-xxx-",db),)))
    matches.sort()
    if not matches:
        return "", symbol, [], []
    after=[None for x in range(len(matches))] #list of matches that go after this one as token boundary
    for idx,(b,eX,e) in enumerate(matches):
        for idx2 in range(idx+1,len(matches)):
            b2,eX2,e2=matches[idx2]
            if b2>e:
                if tokBoundaryRe.match(symbol[e:b2]):
                    after[idx]=b2
                break
    value=[(None,None) for x in range(len(symbol))]
    for revIdx in range(len(matches)-1,-1,-1):
        b,eX,e=matches[revIdx]
        ownValue=e-b
        nextOne=after[revIdx]
        backtrace=[]
        if nextOne:
            assert value[nextOne]!=(None,None)
            ownValue+=value[nextOne][0]
            backtrace=value[nextOne][1]
        if value[b]==(None,None) or ownValue>value[b][0]:
            value[b]=ownValue,[revIdx]+backtrace
    #assemble the result
    firstB=matches[0][0]
    if firstB==0 or tokBoundaryRe.match(symbol[:firstB]):
        lastIntervalIdx=value[firstB][1][-1]
        b,eX,e=matches[lastIntervalIdx]
        return symbol[:e],symbol[e:],matches,value[firstB][1]
    else:
        return "", symbol, [], []
    
reverseDict={"default_symbol":10,"official_symbol":10,"official_name":10,"synonym":7,"homologene_id":6,"ensemblcluster_id":5,"ensembl_id":4,"locustag":3,"cd_antigen":2,"description":2,"locusname":1,"alergen":0}
def matchSort(((b1,eX1,e1),type1),((b2,eX2,e2),type2)):
    if type1!=type2:
        return cmp(reverseDict[type1],reverseDict[type2]) #higher priority last, will sort with reverse
    l1=e1-b1
    l2=e2-b2
    if l1!=l2:
        return cmp(l1,l2) #longer hit last (will sort with reverse)
    s1=e1-eX1
    s2=e2-eX2
    if s1!=s2:
        return -1*cmp(s1,s2) #hit without plural last (will sort with reverse)
    return 0

def hitsByLen(symbol,dbs):
    matches=simstringMatches(symbol,dbs)
    matches=sorted(matches.items(),cmp=matchSort,reverse=True)
    return matches
    
def subsetMatches(matches,(no_b,no_e)): 
    #Filter from matches intervals that overlap with [no_b,no_e) and rewrite matches above no_e to 0-based
    rP=[]
    rS=[]
    for (b,eX,e) in matches:
        if e<=no_b:
            rP.append((b,eX,e))
        elif b>=no_e:
            rS.append((b-no_e,eX-no_e,e-no_e))
    return rP,rS

def formatORGContext(s,matches):
    # matches is a sorted list of (b,eX,e) of organism matches
    mapBack=[] #list of indices into the original string s from the replaced string
    contextpieces=[]
    lastIdx=0
    for (b,eX,e) in matches:
        if b>lastIdx:
            contextpieces.append(s[lastIdx:b])
            mapBack.extend(range(lastIdx,b))
        contextpieces.append("xxxorgxxx")
        mapBack.extend([b for x in range(len("xxxorgxxx"))])
        lastIdx=e
    else:
        contextpieces.append(s[lastIdx:len(s)])
        mapBack.extend(range(lastIdx,len(s)))
    res="".join(contextpieces)
    assert len(mapBack)==len(res)
    return res, mapBack

atLeastOneCharRe=re.compile(r"^.*[a-zA-Z]")
def readAffixFile(fName,reverse):
    res={} #normform:count
    f=open(fName,"r")
    for line in f:
        line=line.strip()
        if not line:
            continue
        line=line.replace("-org-","xxxorgxxx")
        count,s=line.split(" ",1)
        if not atLeastOneCharRe.match(s):
            continue
        s_n=normform(s)
        if len(s_n)<3:
            continue
        if reverse:
            s_n=s_n[::-1]
        res[s_n]=res.get(s_n,0)+int(count)
    return res

def init():
    global orgDB, symbDBs, prefixDB, suffixDB, prefixCounts, suffixCounts
    orgDB=simstring.reader("data/all_species_db")
    orgDB.measure=simstring.overlap
    orgDB.threshold=1
    symbDBs=[]    
    for dbName in ("default_symbol","official_symbol","official_name","synonym","homologene_id","ensemblcluster_id","ensembl_id","locustag","cd_antigen","description","locusname","alergen"):
        dbFile="data/all_symbols_%s_db"%dbName
        try:
            simDB=simstring.reader(dbFile)
        except IOError:
            continue
        simDB.measure=simstring.overlap
        simDB.threshold=1
        symbDBs.append((dbName,simDB))
    prefixDB=simstring.reader("data/all_prefixes_db")
    prefixDB.measure=simstring.overlap
    prefixDB.threshold=1
    suffixDB=simstring.reader("data/all_suffixes_rev_db")
    suffixDB.measure=simstring.overlap
    suffixDB.threshold=1   
    prefixCounts=readAffixFile("src_data/prefixes2.txt",False)
    suffixCounts=readAffixFile("src_data/suffixes2.txt",True)
    print >> sys.stderr, "INIT DONE"
    return symbDBs, orgDB


def stripSymbol(symbol,symbDBs,orgDB,prefixDB,suffixDB,prefixCounts,suffixCounts):
    symbolMatches=hitsByLen(symbol,symbDBs)
    orgMatches=list(simstringMatches(symbol,(("-org-",orgDB),)))
    if not symbolMatches: #No known symbol here, find the longest prefix/suffix and strip
        orgs=bestIntervals(orgMatches)
        repl,mapBack=formatORGContext(symbol,orgs)
        pref,pref_remain,prefix_matches,prefix_matches_subset=longest_prefix(repl,prefixDB)
        suff,suff_remain,suffix_matches,suffix_matches_subset=longest_prefix(repl[::-1],suffixDB)
        if len(normform(pref))+len(normform(suff))==0:
            return ("",symbol,"",(0,len(symbol))), "NOAFFIX"
        elif len(normform(pref))+len(normform(suff))<len(normform(repl)): #Something remains, yay!
            B,E=len(pref),len(repl)-len(suff)
            B,E=mapBack[B],mapBack[E-1]+1
            while B<E and symbol[B].isspace():
                B+=1
            while B<E-1 and symbol[E-1].isspace():
                E-=1
            return (pref, symbol[B:E], suff[::-1],(B,E)), "PRSUFF"
        else:
            strips=[]
            #We have an overlap, one way or another, revert to the original count-based algorithm
            for idx in prefix_matches_subset:
                b,eX,e=prefix_matches[idx]
                try:
                    count=prefixCounts[normform(repl[b:eX])]
                except:
                    count=0
                strips.append((b,e,count))
            for idx in suffix_matches_subset[::-1]:
                b,eX,e=suffix_matches[idx]
                try:
                    count=suffixCounts[normform(repl[::-1][b:eX])]
                except:
                    count=0
                strips.append((len(repl)-e,len(repl)-b,count))
            B=0
            E=len(repl) 
            while True:
                if not strips:
                    break
                if strips[0][2]>strips[-1][2]:
                    newB=strips[0][1]
                    newE=E
                    strips.pop(0)
                else:
                    newE=strips[-1][0]-1
                    newB=B
                    strips.pop(-1)
                if newB>=newE:
                    break
                B=newB
                E=newE
            # print [B, E, mapBack, repl, symbol]
            return ("",repl,"",(mapBack[B],mapBack[E-1]+1)),"OVERLAP"
    #We have a symbol match
    stripped_versions=[]
    for (b,eX,e),symbType in symbolMatches: #In order, longest to shortest single GGP match in the string
        pref=symbol[:b]
        suff=symbol[e:]
        if normform(pref)=="" and normform(suff)=="": #Perfect match
            return ("",symbol[b:eX],"",(b,eX)), "match"
        orgP,orgS=subsetMatches(orgMatches,(b,e)) #so these are all organism matches compatible with this symbol
        orgPMatches_resolved=bestIntervals(orgP) #widest-spanning subset of orgs in prefix
        orgSMatches_resolved=bestIntervals(orgS) #widest-spanning subset of orgs in suffix
        fullPrefix,mapBackPref=formatORGContext(symbol[:b],orgPMatches_resolved)
        fullSuffix,mapBackSuff=formatORGContext(symbol[e:],orgSMatches_resolved)
        pref,pref_remain,delme,delme2=longest_prefix(fullPrefix,prefixDB)
        suff,suff_remain,delme,delme2=longest_prefix(fullSuffix[::-1],suffixDB)
        suff,suff_remain=suff[::-1],suff_remain[::-1]
        if normform(pref_remain)=="" and normform(suff_remain)=="":
            return ("",symbol[b:eX],"", (b,eX)), "match" #stripped down to a perfect match
        else:
            stripped_versions.append((pref_remain,symbol[b:e],suff_remain, (b-len(pref),e+len(suff))))
    else:
        if stripped_versions:
            stripped_versions.sort(key=lambda (pref,sym,suff,be): len(pref)+len(sym)+len(suff))
            pref,sym,suff,be=stripped_versions[0]
            return ("",(pref+sym+suff),"",be), "contains"
        else:
            assert False

def stripSymbolLIB(symbol):
    global symbDBs,orgDB,prefixDB,suffixDB,prefixCounts,suffixCounts
    return stripSymbol(symbol,symbDBs,orgDB,prefixDB,suffixDB,prefixCounts,suffixCounts)
