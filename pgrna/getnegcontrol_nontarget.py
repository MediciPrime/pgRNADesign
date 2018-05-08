#!/usr/bin/env python

from __future__ import print_function

import sys;
import random;


def generateNonTarget(ngctrlfile,outrecordfile,npair):
  '''
  Generating non-targeting control pairs
  '''
  #ngctrlfile='GeCKOv2NonTarget.csv'
  #outrecordfile='negative_control.nontarget.txt';
  
  allsgrecord=[];
  nl=0;
  for line in open(ngctrlfile):
    nl+=1;
    if nl==1:
      continue;
    field=line.strip().split('\t');
    allsgrecord+=[field];
  
  # npair=100;
  nrec=len(allsgrecord);
  print('Reading '+str(nrec)+' records.');
  
  orfhd=open(outrecordfile,'w');
  seldict={};
  for i in range(npair):
    while True:
      sind=sorted(random.sample(range(0,nrec),2));
      print('Sampling '+str(sind[0])+' and '+str(sind[1]));
      sgx=allsgrecord[sind[0]];
      sgy=allsgrecord[sind[1]];
      sgxseq=sgx[-3].upper();
      sgyseq=sgy[-3].upper();
      if str(sind) not in seldict:
        seldict[str(sind)]=0;
        break;
      # if sgxseq.count('CGTCTC')==0 and sgxseq.count('GAGACG')==0 and sgyseq.count('CGTCTC')==0 and sgyseq.count('GAGACG')==0 :
      #  break;
    printid=sgx[3]+'_'+sgy[3];
    # write the record file
    print('\t'.join([printid, '0','0','','NA','0','0',sgx[-3],'.',sgx[-1],'NA','0','0',sgy[-3],'.',sgy[-1],'nontarget']),file=orfhd);
  orfhd.close();
