#!/usr/bin/env python3
"""
Identify unique barcodes for each sgRNA
"""

from __future__ import print_function

import sys;
import os;
import re;


def switchsgs(field):
  sg1f=[t for t in field[4:10]];
  sg2f=[t for t in field[10:16]];
  field[4:10]=sg2f;
  field[10:16]=sg1f;
  return field;

barcode={};

nremove=0;

sgpool=[];
sgpool2=[];
sglibrary={};

# save the sgpool
nline=0;
for line in open(sys.argv[1]):
  field=line.strip().split('\t');
  sgseq1=field[7].upper();
  sgseq2=field[13].upper();
  field[7]=sgseq1;
  field[13]=sgseq2;
  
  sgpool+=[field]
  if sgseq1 not in sglibrary:
    sglibrary[sgseq1]=0;
  sglibrary[sgseq1]+=1;
  if sgseq2 not in sglibrary:
    sglibrary[sgseq2]=0;
  sglibrary[sgseq2]+=1;
  nline+=1;

print('Records: '+str(nline));



while True:
  # keep those with unique sgrna (in current selection) and without conflicting with barcodes
  nline=0;
  for field in sgpool:
    sgseq1=field[7];
    sgseq2=field[13];
    #print(str(sgseq1));
    # if the 1st sg is unique: save directly
    if sglibrary[sgseq1]==1  and sgseq1 not in barcode:
      barcode[sgseq1]=field;
    elif sglibrary[sgseq2]==1  and sgseq2 not in barcode:
      # switch sgrna 1 and sgrna2
      field=switchsgs(field);
      barcode[sgseq2]=field;
    else:
      # non-unique, save for later
      sgpool2+=[field];
    nline+=1;
  
  # the remaining one, search for the unique barcodes
  sglibrary={};
  for field in sgpool2:
    sgseq1=field[7];
    sgseq2=field[13];
    if sgseq1 not in sglibrary:
      sglibrary[sgseq1]=0;
    sglibrary[sgseq1]+=1;
    if sgseq2 not in sglibrary:
      sglibrary[sgseq2]=0;
    sglibrary[sgseq2]+=1;
  print('Initial pool: '+str(len(sgpool))+', Remaining pool:'+str(len(sgpool2)));
  if len(sgpool)==len(sgpool2):
    break;
  sgpool=sgpool2;
  sgpool2=[];


print('Trying to add one of two sgRNAs to library...');

while True:
  sglibrary={};
  # create a library
  for field in sgpool:
    sgseq1=field[7];
    sgseq2=field[13];
    if sgseq1 not in sglibrary:
      sglibrary[sgseq1]=0;
    sglibrary[sgseq1]+=1;
  
  nline=0;
  sgpool2=[];
  for field in sgpool:
    sgseq1=field[7];
    sgseq2=field[13];
    #print(str(sgseq1));
    # if the 1st sg is unique: save directly
    if sglibrary[sgseq1]==1  and sgseq1 not in barcode:
      barcode[sgseq1]=field;
    elif sgseq2 not in sglibrary  and sgseq2 not in barcode:
      # switch sgrna 1 and sgrna2
      field=switchsgs(field);
      barcode[sgseq2]=field;
    else:
      # non-unique, save for later
      sgpool2+=[field];
    nline+=1;
  print('Initial pool: '+str(len(sgpool))+', Remaining pool:'+str(len(sgpool2)));
  if len(sgpool)==len(sgpool2):
    break;
  sgpool=sgpool2;
  sgpool2=[];
 
print('Trying to add two of two sgRNAs to library...');

while True:
  sglibrary={};
  # create a library
  for field in sgpool:
    sgseq1=field[7];
    sgseq2=field[13];
    if sgseq2 not in sglibrary:
      sglibrary[sgseq2]=0;
    sglibrary[sgseq2]+=1;
  
  nline=0;
  sgpool2=[];
  for field in sgpool:
    sgseq1=field[7];
    sgseq2=field[13];
    #print(str(sgseq1));
    # if the 1st sg is unique: save directly
    if sglibrary[sgseq2]==1  and sgseq2 not in barcode:
      # switch sgrna 1 and sgrna2
      field=switchsgs(field);
      barcode[sgseq2]=field;
    elif sgseq1 not in sglibrary  and sgseq1 not in barcode:
      barcode[sgseq1]=field;
    else:
      # non-unique, save for later
      sgpool2+=[field];
    nline+=1;
  print('Initial pool: '+str(len(sgpool))+', Remaining pool:'+str(len(sgpool2)));
  if len(sgpool)==len(sgpool2):
    break;
  sgpool=sgpool2;
  sgpool2=[];
 

print('Final round:');

nline=0;
sgpool2=[];
for field in sgpool:
  sgseq1=field[7];
  sgseq2=field[13];
  #print(str(sgseq1));
  # if the 1st sg is unique: save directly
  if  sgseq1 not in barcode:
    barcode[sgseq1]=field;
  elif  sgseq2 not in barcode:
    # switch sgrna 1 and sgrna2
    field=switchsgs(field);
    barcode[sgseq2]=field;
  else:
    # non-unique, save for later
    sgpool2+=[field];
  nline+=1;

print('Initial pool: '+str(len(sgpool))+', Remaining pool:'+str(len(sgpool2)));

# write the improved design
keeppoolfd=open(sys.argv[1]+'.unique.csv','w');

for (k,field) in barcode.items():
  print(','.join(field),file=keeppoolfd);
  

keeppoolfd.close();


