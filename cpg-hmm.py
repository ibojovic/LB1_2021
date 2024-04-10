#!/usr/bin/env python3
import sys
import numpy as np
import pickle

def get_seq(seqfile):
	seq=''
	for line in open(seqfile):
		if line[0]!='>':seq=seq+line.rstrip()
	return seq


def get_ranges(bedfile):
	l=[]
	for line in open(bedfile):
		v=line.split()
		l.append([v[1],int(v[2]),int(v[3])])
	return l


def match_seq(bed,seq):
	# 0-start half open [min,max)
	s=''
	n=len(bed)
	p0=0
	for i in range(n):
		s=s+(bed[i][1]-p0)*'N'
		s=s+(bed[i][2]-bed[i][1])*'Y'
		p0=bed[i][2]
	s=s+(len(seq)-bed[i][2])*'N'
	return s


def count(seq,cpg,nuc='ACGTN',state='NY'):
  n=len(seq)
  k=len(nuc)
  m=len(state)
  t=np.zeros((m,m))
  e=np.zeros((m,k))
  n=len(seq)
  ps0=state.find(cpg[0])
  pn0=nuc.find(seq[0])
  e[ps0][pn0]=e[ps0][pn0]+1
  for i in range(1,n):
    ps=state.find(cpg[i])
    pso=state.find(cpg[i-1])
    pn=nuc.find(seq[i])
    pno=nuc.find(seq[i-1])
    if ps>-1 and pn>-1:
      e[ps][pn]=e[ps][pn]+1
      if pso>-1 and pno>-1:
         t[pso][ps]=t[pso][ps]+1
  for i in range(m):
    e[i,:]=e[i,:]/np.sum(e[i,:])
    t[i,:]=t[i,:]/np.sum(t[i,:])
  return e,t


def print_results(e,t,nuc='ACGTN',state='NY'):
  for i in range(len(state)):
    for j in range(len(state)):
      print (state[i],'->',state[j],'%.6e\t' %t[i][j],end='')
    print ('')
  for i in range(len(state)):
    print ('State:',state[i],end='\t')
    for j in range(len(nuc)):
        print (nuc[j]+':','%.3f\t' %e[i][j],end='')
    print ('')


if __name__ == '__main__':
  if len(sys.argv)<3:
    print ('python cpg-hmm.py bedfile seqfile')
    sys.exit()
  hmm={}
  state='NY'
  nuc='ACGTN'
  bedfile=sys.argv[1]
  seqfile=sys.argv[2]
  seq=get_seq(seqfile)
  bed=get_ranges(bedfile)
  cgp=match_seq(bed,seq)
  seq=seq.upper()
  e,t=count(seq,cgp)
  print_results(e,t)
  hmm['E']=e
  hmm['T']=t
  pickle.dump(hmm,open('../data/cpg-hmm.pik','wb'))
