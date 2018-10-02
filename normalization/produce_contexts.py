import sys
import simstring
import gzip

import MySQLdb as db
import os.path
import re

MYDIR=os.path.dirname(os.path.abspath(__file__))
#this is where the various blacklists, etc are
DATADIR=os.path.join(MYDIR,"data")

### WARNING WARNING WARNING
###
### Because simstring does not support Unicode, everything below assumes ASCII to work properly
###
### PARAMETERS

MINLEN=3 #Minimum length of a match to really count it for this purpose

#beginRe=re.compile(r"(?<=^|[^a-zA-Z0-9])[a-zA-Z0-9]")
#endRe=re.compile(r"[a-zA-Z0-9](?=$|[^a-zA-Z0-9]|e?s[ -]|e?s$|e?s$)")

wordRe=re.compile(r"[a-zA-Z0-9]+")

replRe=re.compile(r"[^a-z0-9]+")
def normform(entityName):
    lower=entityName.lower()
    return replRe.sub("",lower)[:100]

def isalnum(c):
    return c.isalpha() or c.isdigit()

openParRe=re.compile(r"^.*([\(\{\[])[^\]\)\}]+$")
pluralRe=re.compile(r"^s|es|\(s\)")

def simstringMatches(s,dbs):
    """s is a string to look in, dbs is a list of (dbname,db) in descending order of priority. All valid matches (on token boundary) will be found and returned as {(beg,end):dbname}"""
    hits={} #{(b,eX,e):dict} #eX is the non-suffix end of the entity
    normformChars=[]
    normformOffsets=[] #for each normform index, stores index into s
    for idx in xrange(len(s)):
        if isalnum(s[idx]):
            normformChars.append(s[idx].lower())
            normformOffsets.append(idx)
    normform="".join(normformChars)
    assert len(normform)==len(normformOffsets)
    if not normform or len(normform)<MINLEN:
        return hits
    for dbName,db in dbs:
        response=db.retrieve(normform)
        for r in response:
            if len(r)<MINLEN:
                continue
            normBeg=0
            while True: #Go through all hits of r
                normBeg=normform.find(r,normBeg)
                if normBeg==-1: #no more hits
                    break
                end=normformOffsets[normBeg+len(r)-1]+1 #end-1 because end-1 is the last character of the match and it will be in the normformOffsets, then +1 to get back into the Python-style right-open intervals
                beg=normformOffsets[normBeg]
                normBeg+=len(r) #this is where we search from next
                #Now we need to see if that hit is valid
                #1) beg should be beginning of string, or preceded by non-alphanum
                if beg>0 and isalnum(s[beg-1]): #not a valid hit
                    continue
                #Now the end of the string is worse, we need to check for plurals and other stuff
                if end==len(s): #we end at the end of the string -> done
                    if (beg,end,end) not in hits:
                        hits[(beg,end,end)]=dbName
                        continue
                #Take care of expanding parentheses
                if end<len(s) and s[end] in (")","}","]"): #check for un-closed parentheses
                    match=openParRe.match(s[beg:end])
                    if match and (match.group(1),s[end]) in (("(",")"),("{","}"),("[","]")):
                        end+=1
                endX=end
                match=pluralRe.match(s[end:end+3])
                if match: #we may have a plural here
                    end+=len(match.group())
                #Now end points to the end of the match
                if end==len(s) or not isalnum(s[end]): #we have a token-boundary match
                    if (beg,endX,end) not in hits:
                        hits[beg,endX,end]=dbName
    return hits
        
def bestIntervals(ints):
    """ ints is a list of [a,bX,b) intervals where a,b are
    integers. If b>bX, then there is plural ending at [b:bX]. Returns
    a list with the longest non-overlapping subset. Longest in the
    sense of the sum of interval lengths. When unsure, single long
    interval is preferred over several short. When equivalent
    solutions are present, one is chosen arbitrarily"""
    #Algorithm:
    #http://pages.cs.wisc.edu/~shuchi/courses/787-F09/scribe-notes/lec3.pdf
    ints.sort() #sort by beginning interval
    JS=[] #A list where for each i we store the first j that comes after i ends, if any. JS is at least one shorter than ints, since the last interval(s) have nothing after them
    for i in xrange(len(ints)):
        ends=ints[i][2]
        for j in xrange(i+1,len(ints)):
            if ints[j][0]>=ends:
                JS.append(j)
                break
        else:
            JS.append(None)
    value=[-1 for x in range(len(ints))]
    #Fill in the values for all but the last interval(s)
    for i in xrange(len(value)-1,-1,-1):
        wi=ints[i][2]-ints[i][0]-0.5
        if JS[i]==None:
            value_j=0
        else:
            value_j,tmp=value[JS[i]] #tmp is the stored trace-back info we don't need at this point
        if i==len(value)-1:
            value_next=0
        else:
            value_next,tmp=value[i+1]
        if wi+value_j>=value_next:
            value[i]=(wi+value_j,i) #interval i is part of the optimal solution
        elif wi+value_j<value_next:
            value[i]=(value_next,None) #interval i is NOT part of the optional solution None is the backtrace
    #Now just reassemble the solution:
    selected=[]
    i=0
    while i<len(ints):
        w,next=value[i]
        if next==None: #i was not selected, continue to the next
            i+=1
        elif next==i: #i was selected, continue to JS[i] or quit if None
            selected.append(ints[i])
            if JS[i]==None:
                break
            else:
                i=JS[i]
    return selected

def best_matches(s,dbs):
    """s is a string to look in, dbs is a dictionary {dbname,db}
       returns the best non-overlapping matches from s"""
    hits=simstringMatches(s,dbs) #{(b,e):dictname}
    intervals=list(hits.keys())
    best=bestIntervals(intervals)
    result=[] #listof (beg,end,dict)
    for b,eX,e in best:
        result.append((b,eX,e,hits[(b,eX,e)]))
    return result
