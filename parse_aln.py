#!/usr/bin/python
import sys
import numpy as np


def get_aln(alnfile):
  d_aln={}
  f=open(alnfile,'r')
  for line in f:
    if line.find('sp|')!=0: continue
    l=line.rstrip().split()
    pid=l[0].split('|')[-1]
    seq=l[1]
    d_aln[pid]=d_aln.get(pid,'')+seq
  return d_aln


def get_profile(d_aln):
  profile=[]
  n=len(list(d_aln.values())[0])
  sids=d_aln.keys()
  for i in range(n):
    aas=[d_aln[j][i] for j in sids]
    vaas=get_iprofile(aas)
    tot=float(vaas[:20].sum())
    vaas[:20]=vaas[:20]/tot
    profile.append(vaas)
  return profile


def get_iprofile(aas,aa_list='ACDEFGHIKLMNPQRSTVWY-'):
	v=np.zeros(len(aa_list))
	for aa in aas:
		pos=aa_list.find(aa)
		if pos>-1: v[pos]=v[pos]+1
	return v


def print_profile(profile,aa_list='ACDEFGHIKLMNPQRSTVWY-'):
	n=len(profile)
	for i in range(n):
		pi=profile[i][:20]
		s=0.0
		for j in range(20):
			if pi[j]>0: s=s-pi[j]*np.log(pi[j])
		pm=pi.argmax()
		print (i+1,aa_list[pm],pi[pm],s,profile[i][20])


if __name__ == '__main__':
  alnfile=sys.argv[1]
  d_aln=get_aln(alnfile)
  profile=get_profile(d_aln)
  print_profile(profile)
