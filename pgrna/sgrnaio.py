# reading sgrna list
# here, enhancer = gene

from __future__ import print_function

from pgrna.ClassDef import *;

import collections;
import bisect;
import sys;

class FIELD_INDEX:
  # required field
  chr_i=-1;
  start_i=-1;
  end_i=-1;
  symbol_i=-1;
  strand_i=-1;
  short_seq=-1;
  # optional field
  hg_0=-1;
  hg_1=-1;
  hg_01=-1;
  hg_10=-1;
  score_i=-1;




def getfieldindex(field):
  '''
  Get the index of the corresponding columns
  '''
  fic=FIELD_INDEX();
  for fid in range(len(field)):
    tf=field[fid]
    if tf == 'chr':
      fic.chr_i= fid;
    if tf == 'start':
      fic.start_i= fid;
    if tf == 'end':
      fic.end_i= fid;
    if tf == 'gene_symbol':
      fic.symbol_i= fid;
    if tf == 'strand':
      fic.strand_i= fid;
    if tf == 'short_seq':
      fic.short_seq = fid;
    if tf == 'score':
      fic.score_i=fid;
    if tf == 'hg19_19_0':
      fic.hg_0=fid;
    if tf == 'hg19_19_1':
      fic.hg_1=fid;
    if tf == 'hg19_19_01':
      fic.hg_01=fid;
    if tf == 'hg19_19_10':
      fic.hg_10=fid;
  if fic.chr_i==-1 or fic.start_i==-1 or fic.end_i==-1 or fic.strand_i==-1 or fic.symbol_i==-1 or fic.short_seq==-1:
    print('Error: cannot find sequence field of sgRNAs. They should have the column header "short_seq".');
    sys.exit(-1)
  return fic;

def readsgrnainfo(filename,debug=False,enhancer=False,effcutoff=0.0,hg19_19_0_cutoff=0,hg19_19_10_cutoff=0,chromosome='all'):
  '''
  reading possible sgRNA list
  filter unnecessary sgRNAs
  '''
  allsg={};
  allsgrec=0;
  targetgenename=[];
  n=0;
  allsgids={};
  fic=None;
  for line in open(filename):
    field=line.strip().split();
    n+=1;
    if n==1:
      fic=getfieldindex(field);
      continue;
    if n%100000 == 1:
      print('.',end='');
    #if debug and n>100000:
    #  break;
    gid=field[fic.symbol_i];
    chrv=field[fic.chr_i];
    if debug and chromosome !='all':
      if chrv!=chromosome:
        continue
    try:
      vstart=int(field[fic.start_i]);
      vend=int(field[fic.end_i]);
    except ValueError:
      continue
    vseq=field[fic.short_seq];
    ori=field[fic.strand_i];
    # remove if duplicated
    thissgid=chrv+':'+str(vstart)+':'+ori;
    if thissgid in allsgids:
      # continue;
      pass;
    allsgids[thissgid]=0;
    #
    if fic.score_i==-1:
      effs=100;
    else:
      effs=float(field[fic.score_i]);
    if effs<effcutoff:
      continue;
    if fic.hg_0==-1:
      hg19_19_0=0
    else:
      hg19_19_0=int(field[fic.hg_0]);
    if fic.hg_1==-1:
      hg19_19_10=0
    else:
      hg19_19_10=int(field[fic.hg_1]);
    if hg19_19_0 >hg19_19_0_cutoff or hg19_19_10>hg19_19_10_cutoff:
      if gid == 'breast_prostate3':
        print('skipped for off target ...');
      continue;
    # count g/G
    if vseq.upper().count('G')>7:
      if gid == 'breast_prostate3':
        print('skipped for >7 Gs ...');
      continue;
    # remove ESP3I (BSMBI) cutting sides
    # if vseq.upper().count('CGTCTC')>0 or vseq.upper().count('GAGACG')>0:
    if vseq.upper().count('GGTCTC')>0 or vseq.upper().count('GAGACC')>0:
      continue;
    tif=SGRNA_CLASS();
    tif.chr=chrv; tif.start=vstart; tif.end=vend; tif.seq=vseq; tif.ori=ori;tif.effscore=effs;
    # 
    allsgrec+=1;
    gsubid=gid.split('_');
    uniqgsid={};
    for gs in gsubid:
      gsname=gs.split('_');
      # print('gsname:'+gsname[1]);
      if enhancer:
        #identifier=gsname[0]+'_'+gsname[1];  # enhancer ID: using DNase peak as ID
        identifier=gsname[0];  # gene ID
      else:
        identifier=gsname[0];  # gene ID
      if identifier not in targetgenename:
        # print('gene name:'+identifier);
        targetgenename+=[identifier]; # set gene name into white list
      if identifier not in uniqgsid:
        uniqgsid[identifier]=[];
      uniqgsid[identifier]+=[gs];
    # for each different gene id, try to add this record separately
    for (gk, gv) in uniqgsid.items():
      if gk not in allsg:
        allsg[gk]=[];
      allsg[gk]+=[tif];
  print();
  print('sgRNA: '+str(n)+', remaining sgRNA:'+str(allsgrec));
  print('target genes:'+str(len(targetgenename)));
  print('target genes:'+str(targetgenename));
  targetgenename=set(targetgenename);
  print('all sgRNA ids: '+' # '.join([ k+':'+str(len(v))  for (k,v) in allsg.items()]));
  return (allsg,targetgenename);


