#!/bin/bash
##############################################################
# This scripts runs ADMIXTURE in the default "unsupervised" mode
# Makes file structure:
# <outdir>/<K>/Rep_<r>/<bedname>.<K>.Q
##############################################################
# Input
bedfile=$1			# .bed file for the input subjects & sites (ideally unrelated and unlinked)
K1=$2				# Will run ADMIXTURE for K1 <= K <= K2
K2=$3				#
reps=$4				# Will run replicates indixed by those listed in reps, as each replicate run is given name indexed by the replicate number. e.g. to run replicates indexed by numbers 1, 3,4 and 5, the input would be a string "1 3 4 5".
outdir=$5			# Path to directory that will hold all output files

# Hardwired paths and constants. Can/should be turned into input
admixturetool=/local3/mhansen/Tools/admixture_linux-1.3.0/admixture
#clumpp=/local3/mhansen/Tools/CLUMPP_Linux64.1.1.2/CLUMPP
numcv=5
numnodes=4
##############################################################
name=${bedfile##*/}
name=${name%.*}
##############################################################

if [ ! -d $outdir ]
then
	mkdir $outdir
fi

for K in `seq $K1 $K2`
do
	outdirK=$outdir/$K
	if [ ! -d $outdirK ]
	then
		mkdir $outdirK
	fi

	# Run R replicates
	for r in $reps
	do
		repdir=$outdirK/Rep_$r
		if [ ! -d $repdir ]
		then
			mkdir $repdir
		fi
		cd $repdir
		if [ ! -e $name.$K.Q ]
		then
			$admixturetool -s time -j$numnodes --cv=$numcv $bedfile $K | tee log.$K
		fi
		cd -
	done
done
