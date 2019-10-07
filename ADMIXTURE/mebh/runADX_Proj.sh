#!/bin/bash
##############################################################
# This scripts runs ADMIXTURE in "projection" mode
##############################################################
# Input
bedfile_proj=$1				# .bed file for just the subjects to be projected (must be same sites as in the input reference bed file)
bedfile_ref=$2				# .bed file for the reference subjects & sites (unrelated and unlinked)
K1=$3						# Will run ADMIXTURE for K1 <= K <= K2
K2=$4						#
reps=$5						# Will run replicates indixed by those listed in reps, as each replicate run is given name indexed by the replicate number. e.g. to run replicates indexed by numbers 1, 3,4 and 5, the input would be a string "1 3 4 5". 
refdir=$6					# Path to directory that contains the reference .bed file. Assumed that it does *not* end with "/"
outdir=$7					# Path to directory that will hold all output files
name_comb=$8				# Name to give all final output reference + projected ADMIXTURE .Q and.P files

# Hardwired paths. Can/should be turned into input
admixturetool=/local3/mhansen/Tools/admixture_linux-1.3.0/admixture
#clumpp=/local3/mhansen/Tools/CLUMPP_Linux64.1.1.2/CLUMPP
##############################################################
name_proj=${bedfile_proj##*/}
name_proj=${name_proj%.*}
name_ref=${bedfile_ref##*/}
name_ref=${name_ref%.*}
bimfile_proj=${bedfile_proj%.*}.bim
famfile_proj=${bedfile_proj%.*}.fam
famfile_ref=${bedfile_ref%.*}.fam
numcv=5
numnodes=4
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
		cp $refdir/$K/Rep_$r/$name_ref.$K.P $repdir/$name_proj.$K.P.in
		cp $bedfile_proj $repdir/$name_proj.bed
		cp $bimfile_proj $repdir/$name_proj.bim
		cp $famfile_proj $repdir/$name_proj.fam

		cd $repdir
		if [ ! -e $name_proj.$K.Q ]
		then
			$admixturetool -s time -j$numnodes --cv=$numcv -P $name_proj.bed $K | tee log.$K
		fi
		echo "Done with the projection"
		cd -

		# Now we combine the projected data with the reference data
#		echo "Combinding raw Q file"
		cat $repdir/$name_proj.$K.Q $refdir/$K/Rep_$r/$name_ref.$K.Q > $repdir/$name_comb.$K.Q.raw
#		echo "Combinding raw fam file"
		cat $famfile_proj $famfile_ref > $repdir/$name_comb.fam.raw
		cd $repdir
#		echo "Pasting fam and Q file together, then sorting"
		paste $name_comb.fam.raw $name_comb.$K.Q.raw > tmp
		sort -k1,1 -k2,2 tmp > tmpsort
		cut -f1-6 tmpsort > $name_comb.fam
		cut -f7- tmpsort > $name_comb.$K.Q

		# Clean up extra files
		rm $name_proj.bed
		rm $name_proj.bim
		rm $name_proj.fam
		rm $name_proj.$K.P.in
		rm tmp
		rm tmpsort
		rm $name_comb.$K.Q.raw
		rm $name_comb.fam.raw
		cd -
	done
done
