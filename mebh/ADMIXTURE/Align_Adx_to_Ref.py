#!/usr/bin/python

"""
Align_Adx_to_Ref.py is a python program that aligns the columns in Admixture P and Q files to match the those in a reference Q file.
Assumes that both have the same K value.
The purpose is to match ADMIXTURE results between a reference dataset and a dataset with the same individuals but different sites (or at least not entirely the same set of sites). 

Input files: 
	1. Reference admixture Q file
	2. Base path to input admixture P and Q files
	3. Base path output for aligned P and Q files
"""
import os
import re
import sys
import argparse
import numpy as np
import scipy.stats as stats

####################################################################
# Do main
if __name__=="__main__":
#	print('Do the main stuff')
	parser = argparse.ArgumentParser(description='Process input for Align_Adx_to_Ref.py')
	parser.add_argument('--ref',type=str,help="base path to reference admixture P and Q files")
	parser.add_argument('--refplink',type=str,help="base path to reference plink files")
	parser.add_argument('--in',type=str,help="base path to input admixture P and Q files" )
	parser.add_argument('--inplink',type=str,help="base path to input plink files")
	parser.add_argument('--r',type=float,help="exponent to use in matching. sum_i|Q[i,k1]-Q[i,k2]|^r" )
	parser.add_argument('--m',type=str,help="mode to run:score,full")   # "score" just prints the matching score to stdout. "full" also writes the aligned .Q and .P files
	parser.add_argument('--out',type=str,help="base path to output, aligned, admixture P and Q files")
	
	args = vars(parser.parse_args())
#	print('args',args)
	
	Qref = np.loadtxt(args['ref'] + '.Q',dtype=float)
	Qin = np.loadtxt(args['in'] + '.Q',dtype=float)
	N,K = Qref.shape
	
	
	r = float(args['r'])
	S = np.zeros((K,K),dtype=float)
	
	for k1 in range(K):
		for k2 in range(K):
			S[k1,k2] = (np.sum(np.abs(Qin[:,k1]-Qref[:,k2])**r)/N)**(1/r)

	f = np.zeros(K,dtype=int)
	scores = np.zeros(K,dtype=float)
	valsort = np.sort(S,axis=None)
	k1fit = []
	k2fit = []
	i=0
	while(len(k2fit)<K):
		v = valsort[i]
		idxs = np.where(S==v)

		k1 = idxs[0]  # input k val
		k2 = idxs[1]  # ref k val
		if k1 not in k1fit and k2 not in k2fit:
			f[k2] = k1
			scores[k2] = v
			k1fit.append(k1)
			k2fit.append(k2)
		i+=1
	f = list(f)
	if args['m']=='score':
		print(np.mean(scores))
#		print('Ave_Score\t%f\n' % (np.mean(scores)))
		del Qref
		del Qin
	elif args['m']=='full':
		print('map f:',f)
		print('scores:',scores)
		print('ave score:',np.mean(scores))
		Qout = Qin[:,f]
		print('Qin[0,:]',Qin[0,:])
		print('Qref[0,:]',Qref[0,:])
		print('Qout[0,:]',Qout[0,:])
		np.savetxt(args['out'] + '.Q',Qout,fmt='%f',delimiter=' ')
		del Qref
		del Qin
		del Qout
		Pin = np.loadtxt(args['in'] + '.P',dtype=float)
		Pout = Pin[:,f]
		np.savetxt(args['out'] + '.P',Pout,fmt='%f',delimiter=' ')
	