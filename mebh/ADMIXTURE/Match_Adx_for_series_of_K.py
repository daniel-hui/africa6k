#!/usr/bin/python

"""
Match_Adx_for_series_of_K.py is a python program that tries to match the admixture components across admixture outputs for different K values
The input is a file that lists the paths for a series of ADMIXTURE output .Q files. One file path per line.
The code assumes that if you input N files, then there are N consequetive K values included. So if the lowest K value is K1, then assumes the files have K values: K1, K1+1, K1+2,..., K1+(N-1).
The code inspects the input files and sorts them by K value, so they can be input in any order.
The code starts with K1 and K2=K1+1 and tries to align each component in K2 to a component in K1. The component in K2 with the worst matching score is designated as the "new" ancestry introduced when going from K1 to K2. 
This process is iterated, matching the components in the next-highest K to the components in the current K. 

Denote the smallest K value by Ks and the largest K value by Ke=Ks+(N-1).
The output is a file with Ke rows.
Each row represents one "ancestry thread" or "ancestry type". Graphically, these ancestry threads should be given the same color in the admixture barcharts. Each column in a row has the form "K:c", where c is the colum in the .Q file with K components that belongs to this ancestry type. c is 1-indexed, i.e. the first column is c=1, the final column is c=K.
The rows are in order of how many K values contain that ancestry type. So the first Ks rows contain Ke columns. The Ks+1 row contains Ke-1 columns. The K2+2 row contains Ke-2 columns. The Ke row contains 1 column. 
"""
import os
import re
import sys
import argparse
import numpy as np
import scipy.stats as stats

####################################################################
# Functions
def alignAncestries(Qs,Kvals,r):
	# r is the distance measure exponent: d = |u-v|^r
	# Returns comp = dictionary, key = ancestry type index c (0-indexed), value = list of tuples (K,k) where k is the column (1-indexed) that correspond to ancestry type c for admixture run with K components. 
	
	numK = len(Kvals)
	comp = {}
	
	if numK > 1:
		for i1 in range(numK-1):
			K1 = Kvals[i1]
			# initialize comp
			if i1 == 0:
				for k in range(1,K1+1):
					comp[k] = [(K1,k)]
			Q1 = Qs[K1]
			K2 = Kvals[i1+1]
			# align the ancestries for the next K value (K2) to the present K value (K1)
			if( K2 == (K1+1)):
				Q2 = Qs[K2]

				# a21 = dictionary. key = tuple of (K2,k2), value = matching ancestry tuple (K1,k1). If k2 is the new ancestry type introduced at K2, then value = "None".
				a21 = alignTwoQ(K1,Q1,K2,Q2,r) 

				# update the comp dictionary
				cvec = sorted(comp.keys())
				cMax = cvec[-1]
				for t2 in a21.keys():
					t1 = a21[t2]
					if t1 != "None":
						for c in comp.keys():
							if t1 in comp[c]:
								comp[c].append(t2)
					if t1 =="None":
						cNew = cMax + 1
						comp[cNew] = [t2]
	else:
		K1 = Kvals[0]
		for k in range(1,K1+1):
			comp[k] = [(K1,k)]
			
	return comp

def alignTwoQ(K1,Q1,K2,Q2,r):
	# compute R2 distance between each ancestry in Q1, k1, and each ancestry in Q2, k2
	# dist measure is dot product of proportions across all individuals to the power of r, sum_i |Q_ik1 - Q_ik2|^r
	N,K1 = Q1.shape
	N,K2 = Q2.shape

	d = np.zeros((K1,K2))
	for k1 in range(1,K1+1):
		c1 = k1 -1
		for k2 in range(1,K2+1):
			c2 = k2 -1
			d[c1,c2] = np.sum(np.abs(Q1[:,c1]-Q2[:,c2])**r)
	
	# create two lists: one of distances and one of (k1,k2) tuples
	dists = []
	x = []
	for k1 in range(1,K1+1):
		c1 = k1 -1
		for k2 in range(1,K2+1):
			c2 = k2 - 1
			dists.append(d[c1,c2])
			x.append((k1,k2))
	order = np.argsort(dists)
	x = np.array(x,dtype=object)
	xorder = x[order]
	
	# go through the sorted (k1,k2) list and match each k2 to k1, among those that remain
	k2_matched = []
	k1_matched = []
	a21 = {}
	for i,(k1,k2) in enumerate(xorder):
		if (not k1 in k1_matched) and (not k2 in k2_matched):
			k1_matched.append(k1)
			k2_matched.append(k2)
			a21[(K2,k2)] = (K1,k1)
	k2all = range(1,K2+1)
	k2Unmatched = list(set(k2all) - set(k2_matched))
	k2new = k2Unmatched[0]
	a21[(K2,k2new)] = "None"
	return a21

# Do main
if __name__=="__main__":
#	print('Do the main stuff')
	parser = argparse.ArgumentParser(description='Process input for Align_Adx_to_Ref.py')
	parser.add_argument('--in',type=str,required=True,help="path to file that lists the input ADMIXTURE .Q files to align. One file path per line. For N .Q files, assumes there are N consequetive K values represented by these files.")
	parser.add_argument('--r',type=float,default=3,help="exponent to use in matching. sum_i|Q[i,k1]-Q[i,k2]|^r" )
	parser.add_argument('--out',type=str,required=True,help="output file that lists which columns in the .Q files belong to each of the \"ancestry threads\".")

	args = vars(parser.parse_args())

	# Read in the input .Q files and check that there are consequetive K values
	Qs = {}
	with open(args['in'],'r') as f:
		for line in f:
			line = line.strip()
			Q = np.loadtxt(line,dtype=float)
			N,K = Q.shape
			Qs[K] = Q
			print(line,N,K)
			
	Kvals = list(sorted(Qs.keys()))
	nK = len(Kvals)
	Ks = Kvals[0]
	Ke = Kvals[-1]
	if Ke != (Ks + (nK-1)):
		print('Error: not-consequetive K values!')
		exit(0)
	
	# Run the aligner
	comp = alignAncestries(Qs,Kvals,args['r'])
	print(comp)
	tuple_to_comp = {}
	cvec = sorted(comp.keys())
	cmax = cvec[-1]
	for i in range(1,cmax+1):
		for t in comp[i]:
			tuple_to_comp[t] = i
	
	# Write output
	with open(args['out'],'w') as f:
		for c in cvec:
			rowvec = ["%d:%d" % (tup[0],tup[1]) for tup in comp[c]]
			row = "\t".join(rowvec)
			f.write('%s\n' % row)