#Python script for superimposing two 3D structures

#!/usr/bin/env python
import sys
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np


def get_list(naa):
	vpos=[]
	l1=naa.split(',')
	for i in l1:
		l2=i.split('-')
		if len(l2)==1: 
			vpos=vpos+l2
		elif len(l2)==2:
			try:
				vpos=vpos+[str(j) for j in range(int(l2[0]),int(l2[1])+1)]
			except:
				sys.stderr.write('ERROR: Incorrect input format '+naa)
				sys.exit(1)
		else:
			sys.stderr.write('ERROR: Incorrect input format '+naa)
			sys.exit(1)
	return vpos


def get_ca_atoms(pdbfile,chain,rlist,atom='CA'):
		l_coord=[]
		fpdb=open(pdbfile)
		for line in fpdb:
			if line[:4] != 'ATOM': continue
			if line[21] != chain: continue
			if line[22:26].strip() not in rlist: continue
			if line[12:16].strip() != atom: continue
      if line[16]!=" " and line[16]!="A": continue
			x=float(line[30:38])
			y=float(line[38:46])
			z=float(line[46:54])
			l_coord.append([x,y,z])
		return l_coord


def get_rmsd(coord1,coord2):
		if len(coord1)!=len(coord2):
				sys.stderr.write('ERROR: The sets of coordinates have different size.')
				sys.exit(1)
		svd=SVDSuperimposer()
		svd.set(np.array(coord1),np.array(coord2))
		svd.run()
		rmsd=svd.get_rms()
		rot,tran=svd.get_rotran()
		return rot,tran,rmsd


if __name__ == '__main__': 
	if len(sys.argv)==7:
		pdbfile1=sys.argv[1]
		pdbfile2=sys.argv[2]
		chain1=sys.argv[3]
		chain2=sys.argv[4]
		list1=get_list(sys.argv[5])
		list2=get_list(sys.argv[6])
		l_coord1=get_ca_atoms(pdbfile1,chain1,list1)
		l_coord2=get_ca_atoms(pdbfile2,chain2,list2)
		rot,tran,rmsd=get_rmsd(l_coord1,l_coord2)
		print (f'R= [[ {rot[0][0]:7.3f}, {rot[0][1]:7.3f}, {rot[0][2]:7.3f} ],\n'+\
		       f'    [ {rot[1][0]:7.3f}, {rot[1][1]:7.3f}, {rot[1][2]:7.3f} ],\n'+\
					 f'    [ {rot[2][0]:7.3f}, {rot[2][1]:7.3f}, {rot[2][2]:7.3f} ]]')
		print (f'T= [ {tran[0]:5.3f}, {tran[1]:5.3f}, {tran[2]:5.3f} ]')
		print (f'RMSD= {rmsd:5.3f}')
	else:
		print ("python super_pdb.py pdb1 pdb2 chain1 chain2 residues1 residues2")
