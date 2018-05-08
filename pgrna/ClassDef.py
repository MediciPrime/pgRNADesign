#

from __future__ import print_function
# class definition

class SGRNA_CLASS:
  '''
  Defining a sgRNA class
  '''
  chr='';
  start=0;
  end=0;
  seq='';
  ori='+';
  effscore=0.0;
  score=0.0;
  otherinfo=[]; # optional field
  def printfield(self):
    v=[self.chr, str(self.start), str(self.end), self.seq, self.ori,str(self.effscore)];
    return v;
  def bedfield(self):
    id=self.seq;
    v=[self.chr, str(self.start), str(self.end), self.seq,str(self.effscore), self.ori];
    return v;
    


class SGRNA_PAIR:
  '''
  Defining a pair of sgRNAs
  '''
  sg1=None;
  sg2=None;
  score=0.0;
  dist=0.0;
  code=''; # define the type code of this pair
  def printfield(self,geneid):
    id=geneid+'_dist'+str(self.dist);
    strfield=[id,str(self.dist),str(self.score),self.code];
    strfield+=self.sg1.printfield();
    strfield+=self.sg2.printfield();
    return strfield;
  def bedfield(self,geneid):
    id=geneid+'_dist'+str(self.dist)+'_'+self.code;
    dist1=self.sg1.end-self.sg1.start;
    dist2=self.sg2.end-self.sg2.start;
    block2start=self.sg2.start-self.sg1.start;
    bedelement=[self.sg1.chr,str(self.sg1.start),str(self.sg2.end),id, \
      str(self.score),self.sg1.ori,str(self.sg1.start),str(self.sg2.end), \
      '0','2',','.join([str(dist1),str(dist2)]),','.join(['0',str(block2start)])]
    return bedelement;
    

class TRANSCRIPT_CLASS:
  '''
  Defining a transcript, or an enhancer
  '''
  symbol='';
  chr=''
  start=0;
  end=0;
  ori='+';
  cdsstart=0;  # used in a transcript for cds, or in an enhancer its center
  containingids=[]; # used in adhoc or super enhancers; the containing enhancer ID
  ######
  # member functions
  def str(self):
    vs=self.symbol+'('+self.chr+':'+str(self.start)+'-'+str(self.end)+';'+self.ori+' CDS='+str(self.cdsstart)+';TSS='+str(self.tss())+')';
    return vs;
  def strtab(self):
    vs='\t'.join([self.symbol,self.chr,str(self.start),str(self.end),self.ori,str(self.cdsstart),str(self.tss())]);
    return vs;
  def tss(self):
    if self.ori=='+':
      return self.start;
    else:
      return self.end;
  def len(self):
    return self.end-self.start;


def sgDistance(sgi,sgj):
  if sgi.end>sgj.end:
    dist=sgi.start-sgj.end;
    sgp1=sgj; # the one with smaller coordinates
    sgp2=sgi; # the one with larger coordinates
  else:
    dist=sgj.start-sgi.end;
    sgp1=sgi;
    sgp2=sgj;
  return (dist,sgp1,sgp2);

