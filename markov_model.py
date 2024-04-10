#!/usr/bin/python
import sys, random
import numpy as np


def get_alphabet(seq):
  list_events=[i for i in seq]
  s=set(list_events)
  alphabet=list(s)
  alphabet.sort()
  return ''.join(alphabet)


def get_tmatrix(seq,alphabet):
	n=len(alphabet)
	tm=np.zeros((n,n))
	l=len(seq)
	for i in range(l-1):
		p1=alphabet.find(seq[i])
		p2=alphabet.find(seq[i+1])
		if p1==-1 or p2==-1: continue
    tm[p1][p2]+=1
	for i in range(n):
		tm[i,:]=tm[i,:]/np.sum(tm[i,:])
	return tm


def get_probility(seq,tm,alphabet):
	p=1.0
	l=len(seq)
	for i in range(l-1):
		p1=alphabet.find(seq[i])
		p2=alphabet.find(seq[i+1])
		if p1==-1 or p2==-1: continue
    # Warning if tm[p1][p2]=0 -> p=0
    p=p*tm[p1][p2]
	return p


def get_logprobility(seq,tm,alphabet):
	logp=0.0
	l=len(seq)
	for i in range(l-1):
		p1=alphabet.find(seq[i])
		p2=alphabet.find(seq[i+1])
		if p1==-1 or p2==-1: continue
    # Warning if tm[p1][p2]=0 -> logp=inf
    logp=logp-np.log10(tm[p1][p2])
	return logp


def get_shuffle(seq):
	l=[i for i in seq]
	random.shuffle(l)
	return ''.join(l)


def print_matrix(m,alphabet):
  print ('Transition Matrix:')
  n=len(alphabet)
  for i in range(n):
    for j in range(n):
      print (alphabet[i],'->',alphabet[j],'%.3f\t' %m[i][j], end = '')
    print ('')


if __name__ == '__main__':
  if len(sys.argv)<2:
    print ('python markov_model.py string_of_events')
    sys.exit()
  seq=sys.argv[1]
  alphabet=get_alphabet(seq)
  tm=get_tmatrix(seq,alphabet)
  p=get_logprobility(seq,tm,alphabet)
  print_matrix(tm,alphabet) 
  print ("\n-Log(Probability):",'%.3f' %p)
  for i in range(5):
    rseq=get_shuffle(seq)
    p=get_logprobility(rseq,tm,alphabet)
    print ("\n-Log(Probability):",'%s %.3f' %(rseq,p))
 
