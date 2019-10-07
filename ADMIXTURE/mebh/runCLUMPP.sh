#!/bin/bash
##############################################################
# This scripts runs CLUMPP on set of replicate ADMIXTURE .Q files
##############################################################
#Input
outdir=$1			# Path to the directory that holds all the ADMIXTURE replicates. Expects folder structure: <outdir>/<K>/Rep_<r>/<name>.<K>.Q
K1=$2				# CLUMPP will be run for each K in the range K1 <= K <= K2
K2=$3				#
R=$4				# Number of ADMIXTURE replicates for each K that were run
famfile=$5			# Path to the plink style .fam file for the samples used in the ADMIXTURE .bed file
paramfile=$6		# Path to the base CLUMPP parameters file. For each K, this file gets modified and copied into the appropriate CLUMPP directory for that K value
#adxname=$7				# Name prefix of the ADMIXTURE .Q files. i.e. a .Q file has the name <adxname>.<K>.Q

# Hardwired paths. Can/should be turned into input
clumpptool=/local3/mhansen/Tools/CLUMPP_Linux64.1.1.2/CLUMPP
#famfile=/local3/mhansen/Projects/CTA/CTA_Main/Data/Genotypes/CTAanc.fam
##############################################################
name=${famfile##*/}
name=${name%.*}
##############################################################
# Run CLUMPP
for K in `seq $K1 $K2`
do
	outdirK=$outdir/$K
	clumppdir=$outdirK/CLUMPP
	if [ ! -d $clumppdir ]
	then
		mkdir $clumppdir
	fi
	paramfileK=$clumppdir/CTA.$K.paramfile
	echo -e "DATATYPE 0" > $paramfileK
	echo -e "INDFILE" CTA.$K.infile >> $paramfileK
	echo -e "OUTFILE" CTA.$K.outfile >> $paramfileK
	echo -e "MISCFILE" CTA.$K.misc >> $paramfileK
	echo -e "K" $K >> $paramfileK
	cat $paramfile >> $paramfileK
	clumppout=$clumppdir/CTA.$K.outfile
	if [ ! -e $clumppout ]
	then
		clumppfile=$clumppdir/CTA.$K.infile
		if [ ! -e $clumppfile ]
		then
		#	rm $clumppfile
		#fi
			for r in `seq 1 $R`
			do
				qfile=$outdirK/Rep_$r/$name.$K.Q
				awk '{out=NR " " NR " (x) 1 : ";for(i=1;i<=NF;i++)out=out " " $i; print out}' $qfile >> $clumppfile
				if [ $r -le $R ]
				then
					echo -e "" >> $clumppfile
				fi
			done
		fi
		cd $clumppdir
		$clumpptool CTA.$K.paramfile
		cd -
	fi
	newQ=$clumppdir/$name.$K.Q
	awk '{out=$6;for(i=7;i<=NF;i++){out = out " " $i}; print out}' $clumppdir/CTA.$K.outfile > $newQ
	fullQ=$clumppdir/$name.ids.Q
	paste <(cut -d' ' -f1,2 $famfile) $newQ > $fullQ
done
