#
from __future__ import print_function

import sys;
import random;
import re;
import argparse;
import collections;
import bisect;
import math;
from pgrna.CoverageOverlap import *;
from pgrna.ClassDef import *;
from pgrna.geneinfoio import *;
from pgrna.sgrnaio import *;
from pgrna.geneannoio import *;


def isgoodsgrna(sg,sgtype,ngene,exoncvg):
  '''
  Whether a sgrna is a good one
  '''
  tss=ngene.tss();  ori=ngene.ori;
  if sgtype=='promoter' or sgtype=='promoter-pos':
    if sg.ori!=ori:
      pass;
      # return 0;
    if ori=='+':
      if sg.start<tss-5500: 
        return 0;
    else:
      if sg.end>tss+5500:
        return 0;
    if sg.chr not in exoncvg:
      print('Warning: '+sg.chr+' not present in exon coverage.');
    else:
      checkstart=sg.start;
      checkend=sg.end;
      if checkstart<checkend and  exoncvg[sg.chr].check_range(checkstart,checkend):
        return 0;
    return 1;
  elif sgtype=='gene' or sgtype=='enhancer-pos':
    # enhancers/genes
    #if sg.start>coding+2500 or sg.end<coding-2500:
    #  return 0;
    if ori=='+':
      if sg.start<tss-5000: 
        return 0;
    else:
      if sg.end>tss+5000:
        return 0;
    if sg.chr not in exoncvg:
      print('Warning: '+sg.chr+' not present in exon coverage.');
    else:
      checkstart=sg.start;
      checkend=sg.end;
      if checkstart<checkend and  exoncvg[sg.chr].check_range(checkstart,checkend):
        return 0;
  elif sgtype=='enhancer-adhoc':
    if sg.end<ngene.start or sg.start>ngene.end:
      return 0;
  elif sgtype=='aavs1':
    return 1;
  else:
    print('Error: type not recognized.');
    sys.exit(-1);
  return 1;



def weighted_choice(choices):
    weights = choices;
    total = 0
    cum_weights = []
    for w in weights:
        total += w
        cum_weights.append(total)
    x = random.random() * total
    i = bisect.bisect(cum_weights, x)
    return i;


def samplebyfrequency(pairlist,n,reweight=2.0,coverageweight=False):
  '''
  Random sample by considering the frequency of each sgRNA.
  if the sgRNA score is chosen, then it's more likely that it is picked up again.
  '''
  pairweight=[256.0]*len(pairlist);
  chosenlist=[];
  # range coverage
  rangex=min([pl.sg1.start for pl in pairlist]);
  rangey=max([pl.sg2.end for pl in pairlist]);
  bcg=BinCvgObj(rangex,rangey);
  for i in range(n):
    # print('Weight:'+str(pairweight));
    chindex=weighted_choice(pairweight);
    chosenpair=pairlist[chindex];
    chosenlist+=[chosenpair];
    pairweight[chindex]=0.0;
    # re-weight
    for j in range(len(pairlist)):
      if pairlist[j].sg1== chosenpair.sg1 or pairlist[j].sg2==chosenpair.sg1:
        pairweight[j]*=reweight;
      if pairlist[j].sg1==chosenpair.sg2 or pairlist[j].sg2==chosenpair.sg2:
        pairweight[j]*=reweight;
    if coverageweight:
      bcg.add_record(chosenpair.sg1.start,chosenpair.sg2.end);
      # reweight 
      for j in range(len(pairlist)):
        tgwt=bcg.get_avg_cvg(pairlist[j].sg1.start,pairlist[j].sg2.end);
        pairweight[j]*=math.exp(-0.1*tgwt);
    for j in range(len(pairlist)):
      if pairweight[j]>0 and pairweight[j]<1e-5:
        pairweight[j]=1e-5;
  # end for  
  return chosenlist;
  

def pickuppairs(pairlist,binsize,iteminbin,args):
  '''
  Choose a set of sgRNA pairs 
  '''
  #  first, sort for score
  pairlist_filter=[x for x in pairlist if x.score>0]
  print('Filtered pairs:'+str(len(pairlist_filter)))
  sgrnapairs_sorted=sorted(pairlist_filter, key=lambda x: x.score, reverse=True);
  bindict={};
  binscore={};
  for sgpair in sgrnapairs_sorted:
    bin=int(sgpair.dist/binsize);
    if bin<0:
      bin=0;
    if bin not in bindict:
      bindict[bin]=[];
      binscore[bin]=100000000;
    if len(bindict[bin])>iteminbin*1.2 and sgpair.score<binscore[bin]:
      continue;
    bindict[bin]+=[sgpair];
    binscore[bin]=sgpair.score;
  goodsgpair=[];
  
  if args.stype=='enhancer-adhoc':
    rew=0.3;
  else:
    rew=0.3;
  for (bin,pl) in bindict.items():
    print('-bin:'+str(bin)+', score cutoff:'+str(binscore[bin])+', size:'+str(len(pl)));
    ## random choice
    # cchosen=random.sample(pl,min(iteminbin,len(pl)));
    ## random weighted choice
    if args.stype=='enhancer-adhoc':
      cchosen=samplebyfrequency(pl,min(iteminbin,len(pl)),reweight=rew,coverageweight=True); 
    else:
      cchosen=samplebyfrequency(pl,min(iteminbin,len(pl)),reweight=rew); 
    goodsgpair+=cchosen;
  return goodsgpair;



def getpairscores(sgp1,sgp2,ngene,dist,exoncvg,args):
  '''
  Get the sgRNA pair score and the associated pair information for promoter targeting pairs
  '''
  pairscore=0.0;
  paircode='';
  # adjust the TSS a little bit for positive control genes
  gtss=ngene.tss();
  # transcripts: enlarge the distance between tss and cds start a little bit
  gcdss=ngene.cdsstart
  if abs(gtss-ngene.cdsstart)<700:
    if ngene.ori=='+':
      #gtss-=100;
      gcdss+=700
    else:
      #gtss+=100;
      gcdss-=700
  # score for low efficiency
  if sgp1.effscore<0 or sgp2.effscore<0:
    paircode+='LOWEF;';
    pairscore-=1000;
  # score for orientation 
  if sgp1.ori ==ngene.ori and sgp2.ori==ngene.ori:
    pairscore+=10;
    paircode+='+ORI;';
  # make sure distance is within the range, and before the cds
  if  dist>=1000 and dist <=3000:
    pairscore+=100;
    paircode+='+GOODDIST;';
  if  dist<1000: 
    pairscore-=(1000-dist)/20;
    paircode+='NDIST;';
  if  dist>2000: 
    pairscore-=(dist-2000)/1000*100;
    paircode+='NDIST;';
  if dist<500 or dist >=3000:
    pairscore=float("-inf"); paircode+="NNNDIST";
    return (pairscore,paircode);
  if ngene.ori=='+':
    if sgp2.start>gtss-1 and sgp2.end<gcdss-10:
      pairscore+=100;
      paircode+='+5UTR;';
    if sgp2.end>=gcdss-10:
      #pairscore-=100;
      #paircode+='CDS;';
      pass
  else:
    if sgp1.end<gtss+1 and sgp1.start>gcdss+10:
      pairscore+=100;
      paircode+='+5UTR;';
    if sgp1.start<=gcdss-10:
      #pairscore-=100;
      #paircode+='CDS;';
      pass
  # end if 
  # check sgrna spanning across the TSS
  if sgp1.end<gtss and gtss < sgp2.start:
    pairscore+=100;
    paircode+='+TSS;';
  else: # if not, make it as close as possible
    mindist=min(abs(sgp1.end-gtss),abs(sgp2.start-gtss));
    if mindist<50:
      mindist=50;
    if mindist>1:
      pairscore=float("-inf"); paircode+="TOOFAR";
      return (pairscore,paircode);
    else:
      pairscore+=20.0-mindist/100;
      paircode+='+TSSNEARBY;';
  # do not overlap with any exons or promoters of refseq.
  if sgp1.chr not in exoncvg:
    print('Warning: '+sgp1.chr+' not present in exon coverage.');
  else:
    checkstart=sgp1.start;
    checkend=sgp2.end
    if checkstart<checkend and  exoncvg[sgp1.chr].check_range(checkstart,checkend):
      pairscore-=1000;
      paircode+='EXON;';
    if ngene.ori=='+':
      checkstart=sgp1.start;
      checkend=ngene.tss()-10;
      #checkend=sgp2.end
    else:
      checkstart=ngene.tss()+10;
      #checkstart=sgp1.start
      checkend=sgp2.end;
    if checkstart<checkend and  exoncvg[sgp1.chr].check_range(checkstart,checkend):
      pairscore-=1000;
      paircode+='EXON;';

  # end if
  return (pairscore,paircode);
  


def getpairscores_enhancer(sgp1,sgp2,ngene,dist,exoncvg,targetgenecvg,args):
  '''
  Get the sgRNA pair score and the associated pair information for promoter+exon targeting pairs
  '''
  pairscore=0.0;
  paircode='';
  # adjust the TSS a little bit for positive control genes
  gtss=ngene.tss();
  # transcripts: enlarge the distance between tss and cds start a little bit
  gcdss=ngene.cdsstart
  if abs(gtss-ngene.cdsstart)<3000:
    if ngene.ori=='+':
      #gtss-=100;
      gcdss+=3000
    else:
      #gtss+=100;
      gcdss-=3000
  # score for low efficiency
  if sgp1.effscore<0 or sgp2.effscore<0:
    paircode+='LOWEF;';
    pairscore-=1000;
  # score for orientation 
  if sgp1.ori ==ngene.ori and sgp2.ori==ngene.ori:
    pairscore+=10;
    paircode+='+ORI;';
  # make sure distance is within the range, and before the cds
  if  dist>=3000 and dist <=5000:
    pairscore+=100;
    paircode+='+GOODDIST;';
  if  dist<3000: 
    pairscore-=(3000-dist);
    paircode+='NDIST;';
  if  dist>4000: 
    pairscore-=(dist-4000)/1000*100;
    paircode+='NDIST;';
  if dist<10 or dist >=5000:
    pairscore=float("-inf"); paircode+="NNNDIST";
    return (pairscore,paircode);
  if ngene.ori=='+':
    if sgp2.start>gtss-1 and sgp2.end<gcdss-10:
      pairscore+=100;
      paircode+='+5UTR;';
    if sgp2.end>=gcdss-10:
      #pairscore-=100;
      #paircode+='CDS;';
      pass
  else:
    if sgp1.end<gtss+1 and sgp1.start>gcdss+10:
      pairscore+=100;
      paircode+='+5UTR;';
    if sgp1.start<=gcdss-10:
      #pairscore-=100;
      #paircode+='CDS;';
      pass
  # end if 
  # check sgrna spanning across the TSS
  if sgp1.end<gtss and gtss < sgp2.start:
    pairscore+=100;
    paircode+='+TSS;';
  else: # if not, make it as close as possible
    mindist=min(abs(sgp1.end-gtss),abs(sgp2.start-gtss));
    if mindist<50:
      mindist=50;
    pairscore+=20.0-mindist/100;
    paircode+='+TSSNEARBY;';
  # do not overlap with any exons or promoters of refseq.
  if sgp1.chr not in exoncvg:
    print('Warning: '+sgp1.chr+' not present in exon coverage.');
  else:
    checkstart=sgp1.start;
    checkend=sgp2.end
    if checkstart<checkend and  exoncvg[sgp1.chr].check_range(checkstart,checkend):
      pairscore-=1000;
      paircode+='EXON;';
    if ngene.ori=='+':
      checkstart=sgp1.start;
      checkend=ngene.tss()-10;
      #checkend=sgp2.end
    else:
      checkstart=ngene.tss()+10;
      #checkstart=sgp1.start
      checkend=sgp2.end;
    if checkstart<checkend and  exoncvg[sgp1.chr].check_range(checkstart,checkend):
      pairscore-=1000;
      paircode+='EXON;';
  if sgp1.chr not in targetgenecvg:
    print('Warning: '+sgp1.chr+' not present in target gene coverage.');
  else:
    checkstart=sgp1.start;
    checkend=sgp2.end
    if checkstart<checkend and  targetgenecvg[sgp1.chr].check_range(checkstart,checkend):
      pairscore+=300;
      paircode+='+TARGETEXON;';
    else:
      pairscore-=1000
      paircode+='NOEXONCVG;';


  # end if
  return (pairscore,paircode);
  





def getpairscores_enhancer_0(sgp1,sgp2,ngene,dist,exoncvg,args):
  '''
  Get the sgRNA pair score and the associated pair information
  '''
  pairscore=0.0;
  paircode='';
  # score for low efficiency
  if sgp1.effscore<0 or sgp2.effscore<0:
    paircode+='LOWEF;';
    pairscore-=1000;
  # make sure distance is within the range, and before the cds
  if  dist>=150 and dist <=300:
    pairscore+=100;
    paircode+='+GOODDIST;';
  if  dist<150: 
    pairscore-=(150-dist);
    paircode+='NDIST;';
  if  dist>300: 
    pairscore-=100-(dist-300)/700*100;
    paircode+='NDIST;';
  if dist<10 or dist >=1000:
    pairscore=float("-inf"); paircode+="NNNDIST";
    return (pairscore,paircode);
  if sgp2.start>ngene.cdsstart and sgp1.end<ngene.cdsstart:
    paircode+='+SUMMIT;';
    pairscore+=100;
  else:
    pairscore=float("-inf"); paircode+="NOVERLAP"; # do not cut the corresponding enhancer
    return (pairscore,paircode);
  if sgp2.start>ngene.end and sgp1.end<ngene.start:
    paircode+='+INRANGE;';
    pairscore+=100;
  # end if 
  # check sgrna spanning across the TSS
  # do not overlap with any exons or promoters of refseq.
  if sgp1.chr not in exoncvg:
    print('Warning: '+sgp1.chr+' not present in exon coverage.');
  else:
    checkstart=sgp1.start;
    checkend=sgp2.end;
    if checkstart<checkend and  exoncvg[sgp1.chr].check_range(checkstart,checkend):
      # shall we delete those overlapping with gene exons?
      #pairscore-=1000;paircode+='EXON;';
      pairscore=float("-inf"); paircode+="EXON;"; # do not cut the corresponding enhancer
    return (pairscore,paircode);
  # end if
  return (pairscore,paircode);
  


def getpairscores_enhancer_adhoc(sgp1,sgp2,ngene,dist,exoncvg,args):
  '''
  Get the sgRNA pair score and the associated pair information
  '''
  pairscore=0.0;
  paircode='';
  # score for low efficiency
  if sgp1.effscore<0 or sgp2.effscore<0:
    paircode+='LOWEF;';
    pairscore-=100;
  # make sure distance is within the range, and before the cds
  if  dist>=150 and dist <=600:
    pairscore+=100;
    paircode+='+GOODDIST;';
  if  dist<150: 
    pairscore-=(150-dist);
    paircode+='NDIST;';
  if  dist>600: 
    pairscore-=100-(dist-600)/1400*100;
    paircode+='NDIST;';
  if dist<10 or dist >=2000:
    pairscore=float("-inf"); paircode+="NNNDIST";
    return (pairscore,paircode);
  if sgp2.start>ngene.end and sgp1.end<ngene.start:
    pairscore+=100;
  # end if 
  # check sgrna spanning across the TSS
  # do not overlap with any exons or promoters of refseq.
  if sgp1.chr not in exoncvg:
    print('Warning: '+sgp1.chr+' not present in exon coverage.');
  else:
    checkstart=sgp1.start;
    checkend=sgp2.end;
    if checkstart<checkend and  exoncvg[sgp1.chr].check_range(checkstart,checkend):
      pairscore-=1000;
      paircode+='EXON;';
  # end if
  return (pairscore,paircode);
  
def enumeratepairs_pairwise(goodsgrna,ngene,exoncvg,targetgenecvg,args,blockedpairs=[]):
  '''
  Enumerate all vs. all pairs
  '''
  sgrnapairs=[];
  blockedpairs_list=[(x.sg1,x.sg2) for x in blockedpairs];
  # now, loop in sgRNAs
  for i in range(len(goodsgrna)):
    sgi=goodsgrna[i];
    for j in range(i+1,len(goodsgrna)):
      sgj=goodsgrna[j];
      # calculate the pair information
      (dist,sgp1,sgp2)=sgDistance(sgi,sgj);
      if args.stype=='promoter' or args.stype=='promoter-pos':
        (pairscore,paircode)=getpairscores(sgp1,sgp2,ngene,dist,exoncvg,args);
      elif args.stype=='gene' or args.stype=='enhancer-pos':
        (pairscore,paircode)=getpairscores_enhancer(sgp1,sgp2,ngene,dist,exoncvg,targetgenecvg,args);
      elif args.stype=='enhancer-adhoc' or args.stype=='aavs1':
        (pairscore,paircode)=getpairscores_enhancer_adhoc(sgp1,sgp2,ngene,dist,exoncvg,args);
        # print('enhancer-adhoc disabled.');
        # sys.exit(-1);
      if pairscore==float("-inf"):
        continue;
      # save pair
      thispair=SGRNA_PAIR();
      thispair.sg1=sgp1;
      thispair.sg2=sgp2;
      thispair.score=pairscore;
      thispair.dist=dist;
      thispair.code=paircode;
      if (sgp1,sgp2) not in blockedpairs_list:
        sgrnapairs+=[thispair];
      else:
        print('Duplicated list, skipped...');
    # end for j
  # end for i
  return sgrnapairs;

def enumeratepairs_nbyn_distfunc(dist):
  '''
  score function for distance
  '''
  score=0;
  if dist<150 and dist > 50:
    score=100;
  elif dist<=50:
    score=100-(50-dist)/50.0*100;
  elif dist>=150:
    score=100-(dist-150)/850.0*100;
  else:
    score=-100;
  if score<-100:
    score=-100;
  return score;

def enumeratepairs_nbyn(goodsgrna,ngene,exoncvg,args,n):
  '''
  choose an n by n sgRNA pairs
  '''
  if args.stype=='promoter' or args.stype=='promoter-pos':
    ksplit=ngene.tss();
  else:
    ksplit=ngene.cdsstart;
  good_left=[];
  good_right=[];
  for gs in goodsgrna:
    if gs.end < ksplit:
      gdist=abs(gs.end-ksplit);
      gscore=enumeratepairs_nbyn_distfunc(gdist);
      good_left+=[(gs,gscore)];
    elif gs.start > ksplit:
      gdist=abs(gs.start-ksplit);
      gscore=enumeratepairs_nbyn_distfunc(gdist);
      good_right+=[(gs,gscore)];
  # 
  print(str(n)+' by '+str(n)+': choosing from '+str(len(good_left))+' by '+str(len(good_right))+' pairs.');
  
  if  len(good_left)==0 or len(good_right)==0:
    return [];
  if  len(good_left)>=n and len(good_right)>=n:
    # if we have enough pairs
    good_left_sort=sorted(good_left,key=lambda t: t[1], reverse=True);
    good_right_sort=sorted(good_right,key=lambda t: t[1], reverse=True);
    print('Left: '+';'.join([','.join(t[0].printfield())+',score '+str(t[1]) for t in good_left_sort]));
    print('Right: '+';'.join([','.join(t[0].printfield())+',score '+str(t[1]) for t in good_right_sort]));
    good_left_sel=good_left_sort[:n];
    good_right_sel=good_right_sort[:n];
  elif len(good_left)<n:
    good_left_sel=good_left;
    nrt=min(int(math.ceil(n*n/len(good_left_sel))),len(good_right));
    good_right_sort=sorted(good_right,key=lambda t: t[1], reverse=True);
    good_right_sel=good_right_sort[:nrt];
  elif len(good_right)<n:
    good_right_sel=good_right;
    nrt=min(int(math.ceil(n*n/len(good_right_sel))),len(good_left));
    good_left_sort=sorted(good_left,key=lambda t: t[1], reverse=True);
    good_left_sel=good_left_sort[:nrt];
  
  # enumerate pairs
  sgrnapairs=[];
  # now, loop in sgRNAs
  for i in range(len(good_left_sel)):
    sgi=good_left_sel[i][0];
    for j in range(len(good_right_sel)):
      sgj=good_right_sel[j][0];
      # calculate the pair information
      (dist,sgp1,sgp2)=sgDistance(sgi,sgj);
      if args.stype=='promoter' or args.stype=='promoter-pos':
        (pairscore,paircode)=getpairscores(sgp1,sgp2,ngene,dist,exoncvg,args);
      elif args.stype=='gene' or args.stype=='enhancer-pos':
        (pairscore,paircode)=getpairscores_enhancer(sgp1,sgp2,ngene,dist,exoncvg,args);
      else:
        pairscore=0;
        paircode='ERROR;';
      ## if args.stype=='enhancer-adhoc' or args.stype=='aavs1':
      ##  (pairscore,paircode)=getpairscores_enhancer_adhoc(sgp1,sgp2,ngene,dist,exoncvg,args);
      if pairscore==float("-inf"):
        print('Skipped. Code: '+paircode);
        continue;
      # save pair
      paircode+="5X5;";
      thispair=SGRNA_PAIR();
      thispair.sg1=sgp1;
      thispair.sg2=sgp2;
      thispair.score=pairscore;
      thispair.dist=dist;
      thispair.code=paircode;
      sgrnapairs+=[thispair];
    # end for j
  # end for i
  return sgrnapairs;




def enumeratepairs(goodsgrna,ngene,exoncvg,targetgenecvg,args):
  '''
  Pick up a few sgRNA pairs by enumerating all possible sgRNA pairs
  '''
  # for enhancers, emply an n by n
  if args.stype=='gene' or args.stype=='enhancer-pos' or args.stype=='promoter' or args.stype=='promoter-pos':
    # sgrnapairs0=enumeratepairs_nbyn(goodsgrna,ngene,exoncvg,args,5);
    sgrnapairs0=[];
  else:
    sgrnapairs0=[];
  
  
  if args.stype=='promoter' or args.stype=='promoter-pos':
    binsize=10000;
    nsel=10;
  elif args.stype=='gene' or args.stype=='enhancer-pos':
    binsize=10000;
    nsel=10;
  elif args.stype=='enhancer-adhoc':
    binsize=10000;
    nsel=450;
  elif args.stype=='aavs1':
    binsize=10000;
    nsel=400;
  else:
    binsize=10000;
    nsel=10;
  
  # first, choose n by n
  if len(sgrnapairs0)>=nsel:
    return (sgrnapairs0,nsel);
  # 
  # the remaining is chosen from all by all enumerate
  # print('Chosen '+str(len(sgrnapairs0))+' from n by n. The rest are chosen by random:');
  sgrnapairs=enumeratepairs_pairwise(goodsgrna,ngene,exoncvg,targetgenecvg,args,blockedpairs=sgrnapairs0);
  
  print('Candidate sgRNA pairs:'+str(len(sgrnapairs)));
  # random sample
  selpairs=pickuppairs(sgrnapairs,binsize,nsel-len(sgrnapairs0),args);
  
  finalpairs=sgrnapairs0+selpairs;
  
  return (finalpairs,len(sgrnapairs));


def tilingpairs(goodsgrna,ngene,exoncvg,args):
  '''
  Make sgRNA pairs tiling
  '''
  # sort sgRNA
  sortsg=sorted(goodsgrna, key=lambda x: x.start);
  


