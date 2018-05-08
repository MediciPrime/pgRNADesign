#!/usr/bin/env python

from __future__ import print_function

import sys;
import random;
import re;
import argparse;
import collections;
import bisect;

class CvgObj:
  '''
    A coverage object
  '''
  cvgcore=None;
  cvgindex=[];
  def __init__(self):
    self.cvgcore=collections.OrderedDict();
  def add_record(self,cvglist,debug=False):
    '''
    Adding records from a list of ranges
    All ranges must be added using this function once.
    '''
    self.cvgcore=collections.OrderedDict();
    for z in cvglist:
      if len(z)<2:
        raise Exception('Error: the input must be a list of ranges');
      self.cvgcore[z[0]]=0;
      self.cvgcore[z[1]]=0;
    self.cvgindex=sorted(self.cvgcore.keys());
    for z in cvglist:
      ix=bisect.bisect_left(self.cvgindex,z[0]);
      iy=bisect.bisect_left(self.cvgindex,z[1]);
      for d in range(ix,iy):
        if d<0 or d >=len(self.cvgindex):
          continue;
        ddv=self.cvgindex[d];
        self.cvgcore[ddv]+=1;
    if debug==True:
      for z in range(len(self.cvgindex)):
        zk=self.cvgindex[z];
        zv=self.cvgcore[zk];
        print(str(zk)+':'+str(zv)+',',end='');
        print();
  def check_range(self,x,y,debug=False):
    '''
    Check if range [x,y] is in the selected region
    '''
    nsi=bisect.bisect_left(self.cvgindex,x)-1;
    nsj=bisect.bisect_right(self.cvgindex,y);
    if debug:
      print('index range:'+str(nsi)+'-'+str(nsj));
    for di in range(nsi,nsj):
      if di<0 or di>=len(self.cvgindex):
        continue;
      ekey=self.cvgindex[di];
      if debug:
        print('index:'+str(di)+',key:'+str(ekey)+',value:'+str(self.cvgcore[ekey]));
      if self.cvgcore[ekey]>0:
        return True;
    return False;
    

class BinCvgObj:
  '''
    A binary coverage obj
  '''
  cvgindex=[];
  start=0;
  end=0;
  def __init__(self, nx,ny):
    self.start=nx;
    self.end=ny;
    slen=ny-nx;
    self.cvgindex=[0]*slen;
  def add_record(self,px,py):
    '''
    Adding records from range [px,py]
    Multiple calling is allowed
    '''
    if px<self.start or py>self.end:
      print('Error: added range ('+str(px)+','+str(py)+') outside boundary ('+str(self.start)+','+str(self.end)+').');
      return -1;
    for k in range(px-self.start,py-self.start):
      self.cvgindex[k]+=1;
    return 0;
  def check_range(self,x,y):
    '''
    Check if range [x,y] is in the selected region
    '''
    max_x=max(self.start,x);
    min_y=min(self.end,y);
    if max_x >= min_y:
      return False;
    for t in range(max_x,min_y):
      tid=t-self.start;
      if self.cvgindex[tid]>0:
        return True;
    return False;
  def get_avg_cvg(self,x,y):
    '''
    Get the average coverage 
    '''
    max_x=max(self.start,x);
    min_y=min(self.end,y);
    if max_x >= min_y:
      return 0.0;
    sumcvg=0;
    for t in range(max_x,min_y):
      tid=t-self.start;
      sumcvg+= self.cvgindex[tid];
    return sumcvg/(y-x);




if __name__ == '__main__':
  '''
  test case
  '''
  print('d0:');
  d0=CvgObj();
  d0.add_record([(3,6),(5,9),(15,23)],debug=True);
  print(d0.check_range(1,10)==True);
  print(d0.check_range(1,3)==True);
  print(d0.check_range(9,13)==True);
  print(d0.check_range(8,12)==True);
  print(d0.check_range(-1,4)==True);
  print(d0.check_range(21,29)==True);
  print(d0.check_range(13,18)==True);
  print(d0.check_range(-1,1)==False);
  print(d0.check_range(25,300)==False);
  print(d0.check_range(10,12)==False);
  
  print('d1:');
  d1=CvgObj();
  d1.add_record([(3,6),(6,9),(15,23)],debug=True);
  print(d0.check_range(1,5)==True);
  print(d0.check_range(2,3)==True);
  print(d0.check_range(4,5)==True);

  print('d2:');
  d2=CvgObj();
  print(d2.check_range(1,3)==False);

def cvgTest():
  print('Test...');


