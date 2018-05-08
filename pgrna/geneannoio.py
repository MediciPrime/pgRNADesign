# exon annotation

from __future__ import print_function
import re;

from pgrna.CoverageOverlap import *;
from pgrna.ClassDef import *;

import collections;
import bisect;
import sys;
import copy;

def getgeneannofromgtf(filename,targetgenename,initexoncvg=None):
  '''
  Read gene annotation from gtf file
  '''
  # reading exon coverage
  targetgenecvg={}; # exons of target genes
  targettranscriptcvg={}; # transcripts of target gene
  if initexoncvg is None:
    exoncvg={}; # exons of all gene
  else:
    exoncvg=initexoncvg.copy();
    for esk in exoncvg.keys():
      targetgenecvg[esk]=[]
      targettranscriptcvg[esk]=[]
  
  geneidre=re.compile('gene_id "(\S+)";');
  nexon=0;
  for line in open(filename):
    field=line.strip().split('\t');
    if len(field)<9:
      continue;
    #if field[1]!='protein_coding': # and field[1]!='lincRNA':
    #  continue;
    #if field[2]!='exon' and field[2]!='transcript': # add exon and transcript TSS into consideration
    #  continue;
    gname=geneidre.findall(field[8]);
    if len(gname)==0:
      #print(gname[0]);
      continue;
    # mark the record
    chrname='chr'+field[0];
    if chrname not in exoncvg:
      exoncvg[chrname]=[];
      targetgenecvg[chrname]=[];
      targettranscriptcvg[chrname]=[];
    exonstart=int(field[3]);
    exonend=int(field[4]);
    exonori=field[6];
    if gname[0]=='ENSG00000256982':
      print(gname[0]+'  '+chrname+':'+str(exonstart)+'-'+str(exonend));
      #sys.exit(0)
    # first set up the targettranscriptcvg
    if field[2]=='transcript' and gname[0] in targetgenename:
        targettranscriptcvg[chrname]+=[(exonstart,exonend)];
        # pass;
    # change the transcript to TSS
    if field[2]=='transcript':
      if exonori=='+':
        exonend=exonstart;
        exonstart-=5000;
      else:
        exonstart=exonend;
        exonend+=5000;
    else:
      # extend exon region to make sure it does not go too close
      if field[1]=='protein_coding':
        exonstart-=200
        exonend+=200
    if gname[0] in targetgenename and field[2]=='exon':
        targetgenecvg[chrname]+=[(exonstart,exonend)];
        #pass;
    #if field[2]!='transcript':
    if field[1] == 'protein_coding' and (field[2] == 'exon' or field[2] == 'transcript'):
      exoncvg[chrname]+=[(exonstart,exonend)];
    nexon+=1;
  
  exoncvg2={};
  for (ec,ev) in exoncvg.items():
    exoncvg2[ec]=CvgObj();
    exoncvg2[ec].add_record(ev);
  
  exoncvg=None;
  exoncvg=exoncvg2;
  
  exoncvg2={};
  for (ec,ev) in targetgenecvg.items():
    exoncvg2[ec]=CvgObj();
    exoncvg2[ec].add_record(ev);
  targetgenecvg=None;
  targetgenecvg=exoncvg2;
  
  exoncvg2={};
  for (ec,ev) in targettranscriptcvg.items():
    exoncvg2[ec]=CvgObj();
    exoncvg2[ec].add_record(ev);
  targettranscriptcvg=None;
  targettranscriptcvg=exoncvg2;
  exoncvg2={};
  
  return (exoncvg,targetgenecvg,targettranscriptcvg);


