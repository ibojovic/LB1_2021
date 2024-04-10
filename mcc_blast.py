import math 

# Training confusion matrix
mtrain=[[340,290],[1,9710]]

# Training confusion matrix
mtest=[[340,200],[1,10147]]


def accuracy(m):
  return float(m[0][0]+m[1][1])/(sum(m[0])+sum(m[1]))

def mcc(m):
  d=(m[0][0]+m[1][0])*(m[0][0]+m[0][1])*(m[1][1]+m[1][0])*(m[1][1]+m[0][1])
  return (m[0][0]*m[1][1]-m[0][1]*m[1][0])/math.sqrt(d)

print ('ACC Training= %5.3f' %accuracy(mtrain))
print ('MCC Testing= %5.3f\n' %mcc(mtrain))

print ('ACC Training= %5.3f' %accuracy(mtest))
print ('MCC Testing= %5.3f\n' %mcc(mtest))
