#!/usr/bin/env python
import sys
import numpy as np
import pickle


def calculate_forward(seq,t,e,pb,pe,nuc='ACGTN',state='NY'):
	n=len(seq)
	ns=len(state)
	m=np.zeros((ns,n))
	fp=0.0
	for i in range(ns):
		m[i][0]=pb[i]*e[i][nuc.find(seq[0])]
	for i in range(1,n):
		for j in range(ns):
			for k in range(ns):
				m[j][i]=m[j][i]+m[k][i-1]*t[k][j]*e[j][nuc.find(seq[i])]
	for i in range(ns):
		fp=fp+m[i][-1]*pe[i]
	return m,fp


def calculate_logforward(seq,t,e,pb,pe,nuc='ACGTN',state='NY'):
  n=len(seq)
  ns=len(state)
  m=np.zeros((ns,n))
  fp=0.0
  for i in range(ns):
    m[i][0]=-np.log10(pb[i]*e[i][nuc.find(seq[0])])
  for i in range(1,n):
    for j in range(ns):
      tmp=0.0
      for k in range(ns):
        tmp=tmp+10**-(m[k][i-1]-m[0][i-1])*t[k][j]*e[j][nuc.find(seq[i])]
      m[j][i]=m[0][i-1]-np.log10(tmp)
  tmp=0.0
  for i in range(ns):
    tmp=tmp+10**-(m[i][-1]-m[0][-1])*pe[i]
  fp=m[0][-1]-np.log10(tmp)
  return m,fp



if __name__ == '__main__':
  if len(sys.argv)<3:
    print ('python forward-hmm.py seqfile picklefile')
    sys.exit()
  state='NY'
  nuc='ACGT'
  seq=sys.argv[1]
  picklefile=sys.argv[2]
  hmm=pickle.load(open(picklefile,'rb'))
  # Is assumed that pbegin is 0.5
  pb=np.array([0.5,0.5])
  # Is assumed that pend is 0.5 - this not correct
  pe=np.array([0.5,0.5])
  t=hmm['T']
  e=hmm['E']
  m,fp=calculate_logforward(seq,t,e,pb,pe)
  print ('-Log10 Probability:','%.3f' %fp)
