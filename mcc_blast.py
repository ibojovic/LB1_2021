#!/usr/bin/python
import sys, math

def get_blast(filename):
  flist=[]
  d={}
  f=open(filename)
  for line in f:
    v=line.rstrip().split()
    d[v[0]]=d.get(v[0],[])
    d[v[0]].append([float(v[1]),int(v[2])])
  for v in d.values():
    v.sort()
    flist.append(v[0])
  return flist


def get_cm(data,th):
  # CM=[[TP,FP],[FN,TN]]
  # 0 = Negatives 1=Positives
  cm=[[0.0,0.0],[0.0,0.0]]
  for i in data:
    if i[0]<th and i[1]==1: 
      cm[0][0]=cm[0][0]+1
    if i[0]>=th and i[1]==1:
      cm[1][0]=cm[1][0]+1
    if i[0]<th and i[1]==0:
      cm[0][1]=cm[0][1]+1
    if i[0]>th and i[1]==0:
      cm[1][1]=cm[1][1]+1
  return cm
		

def get_acc(cm):
  return float(cm[0][0]+cm[1][1])/(sum(cm[0])+sum(cm[1]))


def mcc(m):
  d=(m[0][0]+m[1][0])*(m[0][0]+m[0][1])*(m[1][1]+m[1][0])*(m[1][1]+m[0][1])
  return (m[0][0]*m[1][1]-m[0][1]*m[1][0])/math.sqrt(d)



if __name__ == "__main__":
  filename=sys.argv[1]
  #th=float(sys.argv[2])
  data=get_blast(filename)   
  for i in range(20):
    th=10**-i
    cm=get_cm(data,th)
    print 'TH:',th,'ACC:',get_acc(cm),'MCC:',mcc(cm),cm
