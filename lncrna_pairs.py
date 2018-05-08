#!/usr/bin/env python3

from __future__ import print_function

import sys;
import random;
import re;
import argparse;
import collections;
import bisect;
import pgrna;
from pgrna import *;
from pgrna.CoverageOverlap import *;
from pgrna.geneinfoio import *;
from pgrna.sgrnaio import *;
from pgrna.geneannoio import *;
from pgrna.sgrna_selection import *;
from pgrna.getnegcontrol_nontarget import *;
from pgrna.ClassDef import *;



# argparse [BR]:  All relevant arguments called by the user
parser=argparse.ArgumentParser('lncrna_pairs: get sgRNA pairs (pgRNAs) for knocking out lncRNA promoters and bodies.');
parser.add_argument('-o','--output-file',help='Prefix of the output file, default results/sample1.',default='results/sample1');
parser.add_argument('-n','--max-pairs',type=int,help='Maximum number of sgRNA pairs per gene, default -1 (no limits)',default=-1);
parser.add_argument('--sgrna',required=True,help='sgRNA sequence file.',default='');
parser.add_argument('--debug',action='store_true',help='sgRNA sequence file.',default='');
parser.add_argument('--chromosome',help='chromosome to search. default is all',default='all');
#parser.add_argument('--design-version',choices=['v1','v2','v3','v8'],default='v8');
parser.add_argument('--gene-annotation',required=True, help='The annotation of genes and lncRNAs in gtf format');
parser.add_argument('--gene-info',required=True, help='The list of target lncRNAs');
parser.add_argument('--blackout-region',help='The list of blackout region (in bed file). If specified, none of the pgRNAs will span over the blackout region.');
parser.add_argument('stype',choices=['promoter', 'gene','aavs1','nontarget'],
                    help='Types of pgRNAs: promoter (lncRNA promoters), gene (lncRNA promoter+gene boday), aavs1 (aavs1 control), nontarget (nontarget control).'); # design pattern: promoters (200bp-5kbp), gene (5k-100k, including promoters) 
args=parser.parse_args();

# cvgTest()

args.output_file+='.'+args.stype;
if args.chromosome != 'all':
  args.output_file+='.'+args.chromosome

exonannofile=args.gene_annotation;



if args.stype=='nontarget':
  generateNonTarget(args.sgrna,args.output_file+'.design.txt',400);
  sys.exit(-1);
    


# reading possible sgRNA list
print('Reading sgRNAs...')
if args.stype=='promoter' or args.stype=='gene':
  (allsg, targetgenename)=readsgrnainfo(args.sgrna,debug=args.debug,chromosome=args.chromosome);
elif args.stype=='aavs1':
  (allsg, targetgenename)=readsgrnainfo(args.sgrna,debug=args.debug,effcutoff=-10000);

# processing blackout region
blackout={}
if args.blackout_region  is not None:
  for lines in open(args.blackout_region):
    field=lines.strip().split()
    fchr=field[0]
    fstart=int(field[1])
    fend=int(field[2])
    if fchr not in blackout:
      blackout[fchr]=[]
    blackout[fchr]+=[(fstart,fend)]

# reading exon coverage
print('Checking exon coverage...')
(exoncvg,targetgenecvg,targettranscriptcvg)=getgeneannofromgtf(exonannofile,targetgenename,initexoncvg=blackout);

# 

##################
# reading gene information: 

print('Reading gene information...')
if args.stype=='promoter' or args.stype=='gene':
  # read refFlat gene promoters
  allgene=readgeneinfo(args.gene_info,targetgenename);
elif args.stype=='aavs1':
  pseudogene=TRANSCRIPT_CLASS();
  pseudogene.symbol='AAVS1';
  allgene={'AAVS1':pseudogene};
  



###########

print('Generating pairs...')
# output files
# now, for each gene, iterate the sgRNAs
osgfile=open(args.output_file+'.design.txt','w');
ogenefile=open(args.output_file+'.gene.txt','w');
osgbedfile=open(args.output_file+'.design.bed','w');

osggoodsgbed=open(args.output_file+'.goodsgrna.bed','w');

nindex=0
for (gid,ngene) in allgene.items():
  if ngene.symbol not in allsg:
    print('Warning: gene '+ngene.str()+' not in sgRNA list.',file=sys.stderr);
  goodsgrna=[];
  for sg in allsg[ngene.symbol]:
    if isgoodsgrna(sg,args.stype,ngene,exoncvg)==1:
      goodsgrna+=[sg];
      print('\t'.join(sg.bedfield()),file=osggoodsgbed);
  nindex+=1
  print(str(nindex)+'/'+str(len(allgene))+': '+ gid+' '+ngene.symbol+', good sgrna:'+str(len(goodsgrna))+' from '+str(len(allsg[ngene.symbol]))+' candidates..');
  #
  #
  print(ngene.str());

  # if args.stype=='promoter' or args.stype=='promoter-pos' or args.stype=='enhancer':
   # enumerate pairs
  (selpairs,ncandidate)=enumeratepairs(goodsgrna,ngene,exoncvg,targetgenecvg,args);
    
  # output candidate pairs
  print(str(len(selpairs))+' chosen.');
  if args.stype=='promoter' or args.stype=='promoter-pos':
    gid=ngene.symbol+'_p'+str(ngene.tss());
  else:
    gid=ngene.symbol;
  for selpair in selpairs:
    vl1=selpair.printfield(gid);
    vl1+=[args.stype];
    print('\t'.join(vl1),file=osgfile);
    vl1bed=selpair.bedfield(gid);
    print('\t'.join(vl1bed),file=osgbedfile);
  # output gene information
  gfield=[gid,str(ncandidate),str(len(selpairs))];
  gfield+=[ngene.strtab()];
  print('\t'.join(gfield),file=ogenefile);
    
  

osgbedfile.close();
ogenefile.close();
osgfile.close();
osggoodsgbed.close();



