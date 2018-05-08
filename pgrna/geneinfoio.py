# gene info definition
# note that, enhancer = gene in this script

from __future__ import print_function
import sys;

from pgrna.ClassDef import *;

def readgeneinfo(file,targetgenename):
  '''
  read refFlat information
  '''
  allgene={};
  allgenerec=0;
  for line in open(file):
    allgenerec+=1;
    if allgenerec<=1:
      #continue;
      pass;
    field=line.strip().split('\t');
    gid=field[0]+'_'+field[1];
    gsymbol=field[0];
    if gsymbol not in targetgenename:
      continue;
    vstart=int(field[4]);
    vend=int(field[5]);
    gchr='chr'+field[2];
    ori=field[3];
    tss=0;
    cdsstart=0;
    if ori =='+':
      tss=vstart;
      cdsstart=int(field[6]);# CDS start
    elif ori == '-':
      tss=vend;
      cdsstart=int(field[7]);# CDS start
    else:
      if allgenerec>0:
        print('Error: unrecognized string for orientation. Skip this line.',file=sys.stderr);
        sys.exit(-1);
    if allgenerec>0:
      tssid=gchr+':'+str(tss);
      ngene=TRANSCRIPT_CLASS();
      ngene.chr=gchr;ngene.symbol=gsymbol; ngene.start=vstart; ngene.end=vend; ngene.ori=ori;ngene.cdsstart=cdsstart;
      if tssid in allgene:
        if allgene[tssid].len() < ngene.len(): # choose the longest one
          allgene[tssid]=ngene;
      else:
        allgene[tssid]=ngene;
    # end if
  return allgene;


def readenhancerinfo(file,targetgenename):
  '''
  read potential enhancer information
  '''
  allgene={};
  allgenerec=0;
  for line in open(file):
    allgenerec+=1;
    if allgenerec<=1:
      pass;
      # continue;
    field=line.strip().split('\t');
    gchr=field[0];
    enhancername_long=field[3];
    enhf=enhancername_long.split('_');
    if len(enhf)<2:
      print('Error: cannot recognize name field.');
      sys.exit(-1);
    gid=enhf[0]+'_'+enhf[1];
    gsymbol=gid;
    if gsymbol not in targetgenename:
      continue;
    try:
      vstart=int(field[1]);
      vend=int(field[2]);
      ori=field[5];
      cdsstart=int(enhf[1]);
    except ValueError:
      continue;
    if ori =='+':
      pass;
    elif ori == '-':
      pass;
    else:
      if allgenerec>0:
        print('Error: unrecognized string for orientation. Skip this line.',file=sys.stderr);
        sys.exit(-1);
    if allgenerec>0:
      ngene=TRANSCRIPT_CLASS();
      ngene.chr=gchr;ngene.symbol=gsymbol; ngene.start=vstart; ngene.end=vend; ngene.ori=ori;ngene.cdsstart=cdsstart;
      if gid in allgene:
        print('Warning: duplicated enhancers: '+ngene.str());
      else:
        allgene[gid]=ngene;
    # end if
  return allgene;



def readenhanceradhocinfo(file):
  '''
  read a list of ad-hoc enhancer information
  '''
  allgene={};
  allgenerec=0;
  enhancerids=[];
  for line in open(file):
    allgenerec+=1;
    if allgenerec<=1:
      # continue;
      pass;
    field=line.strip().split('\t');
    gchr=field[0];
    enhancername_long=field[3];
    enhf=enhancername_long;
    gid=enhf;
    gsymbol=gid;
    vstart=int(field[1]);
    vend=int(field[2]);
    ori=field[5];
    cdsstart=0;
    if ori =='+':
      pass;
    elif ori == '-':
      pass;
    else:
      if allgenerec>0:
        print('Error: unrecognized string for orientation. Skip this line.',
              file=sys.stderr);
        sys.exit(-1);
    if allgenerec>0:
      ngene=TRANSCRIPT_CLASS();
      ngene.chr=gchr;ngene.symbol=gsymbol; ngene.start=vstart;
      ngene.end=vend; ngene.ori=ori; ngene.cdsstart=cdsstart;
      ngene.containingids=[];
      print('gene: '+str(gsymbol));
      for xid in gsymbol.split(';'):
        xidfield=xid.split('_');
        if xid.startswith('chr'):
          eid=xidfield[0]+'_'+xidfield[1];
        else:
          eid=xid;
        print('-->'+eid);
        ngene.containingids+=[eid];
        enhancerids+=[eid];
      if gid in allgene:
        print('Warning: duplicated enhancers: '+ngene.str());
      else:
        allgene[gid]=ngene;
    # end if
  return (allgene,enhancerids);


def cleansgbyadhocenhancer(enhancerids,allgene,allsg):
  '''
  Only keep sgRNAs within ad-hoc enhancers
  '''
  namemap={};
  nrec=0;
  for (seid,eid) in allgene.iteritems():
    for k in eid.containingids:
      namemap[k]=seid;
      print(k+'->'+seid);
      nrec+=1;
  print('SE ID mapped to enhancer ID:'+str(len(namemap))+' ->'+str(nrec));
  #
  newsg={};
  nsg=0;
  print('genes appear in sgRNA list: '+' # '.join(allsg.keys()));
  for (sgid,sglist) in allsg.iteritems():
    if sgid in namemap:
      newname=namemap[sgid];
      if newname not in newsg:
        newsg[newname]=[];
      for x in sglist:
        if x not in newsg[newname]:
          newsg[newname]+=[x];
          nsg+=1;
  print('New sg list: '+str(len(newsg))+', sgRNAs:'+str(nsg));
  print('--'+' # '.join(newsg.keys()));
  return newsg;
  
  

