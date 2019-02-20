#RUN_fineSTRUCTURE.sh

# Software:
# - fineSTRUCTURE v2+ and all acompaning scripts
# - PLINK v1.07+
# - Perl v5.0+
# - R v3.0+

# Required data: (all copied/linked in the same directory) 
# - Binary PLINK dataset (replace FAMID with the POP tag)
# - Genetic recombination map (genetic_map_chr${chr}_combined_b37.txt)

# Usage:
# bash RUN_fineSTRUCTURE.sh [PLINK dataset suffix]

#<START>

# 1. Preparations
echo ""
echo "STARTING> fineSTRUCTURE pipeline for dataset: " $1

# Input files lists:
allphs=""
allrec=""

# 2. Split a binary PLINK dataset into chromosomes and run SHAPEIT to phase them
for chr in $(seq 1 22); do
	echo ""
	if [ -f "$1.chr$chr.bed" ]
	then
		echo "FOUND> PLINK binary data set for chromosome " $chr
		echo "SKIP> Step 1 PLINK"
	else
		echo "STEP1> Extracting chromosome " $chr
		plink --bfile $1 --chr $chr --allow-no-sex --make-bed \
		--out $1.chr$chr ;
	fi
	echo ""
	if [ -f "$1.chr$chr.haps" ]
	then
		echo "FOUND> Phased haplotypes for chromosome " $chr
		echo "SKIP> Step 2 SHAPEIT"
	else
		echo "STEP2> Phasing chromosome " $chr
		shapeit --input-bed $1.chr$chr.bed $1.chr$chr.bim $1.chr$chr.fam \
			 -M genetic_map_chr${chr}_combined_b37.txt \
			 -O $1.chr$chr
	fi

# 3. Transform IMPUTE2 files to fineSTRUCTURE input
	echo ""
	if [ -f "$1.chr$chr.phase" ]
	then
		echo "FOUND> Phased input files for chromosome " $chr
		echo "SKIP> Step 3 IMPUTE2fSTR"
	else
		echo "STEP3> Transforming SHAPEIT output for chromosome " $chr
		impute2chromopainter.pl $1.chr$chr.haps $1.chr$chr
	fi

# 4. Make custom recombination maps from the BIM files to use in fineStructure
	echo ""
	if [ -f "$1.chr$chr.recomb" ]
	then
		echo "FOUND> Recombination map for chromosome " $chr
		echo "SKIP> Step 4 RecombMap"
	else
		echo "STEP4> Building recombination map for chromosome " $chr
		awk 'BEGIN{print"start.pos","recom.rate.perbp"}
			NR==FNR{a[$1]; if(NR>1) print $1,$2; next}
			(($4 in a)==0){print $4,"NA"
			}' genetic_map_chr${chr}_combined_b37.txt \
			$1.chr$chr.bim | sort -n -k 1 > $1.chr$chr.recomb
		echo "Imputing recombination map for chromosome " $chr
		Rscript IMPUTE_recomb_rates.R $1.chr$chr.recomb
		awk 'NR==FNR{a[$4]; next}(FNR==1 || $1 in a){print $0
		}' $1.chr$chr.bim $1.chr$chr.recomb > tmp
		mv tmp $1.chr$chr.recomb
	fi

# 5. Register input files
	echo ""
	echo "STEP5> Adding chromosome $chr to fsSTR list"
	allphs="$allphs $1.chr$chr.phase"
	allrec="$allrec $1.chr$chr.recomb"
done 

# 6. Create sample list

#Extra: if IDs need to be renamed with FAMID+number
#awk 'BEGIN{p="";n=0}(p != $1){p=$1;n=1}{if(n<10)$2=p"0"n;else $2=p""n;print $0;n+=1}' $1.fam > tmp ; mv tmp > $1.fam

awk '{print $2,$1,"1"}' $1.fam > $1.ids

echo ""
echo "RUNNING> fineSTRUCTURE for dataset: " $1
# 7. Start the fineSTRUCTURE analyses
fs $1.cp -idfile $1.ids \
	-phasefiles $allphs \
	-recombfiles $allrec \
	-go

echo ""	
echo "FINISHED> fineSTRUCTURE pipeline for dataset: " $1

#Extra: Move all log files
echo ""
echo "Moving log files to ./logs/"
mkdir logs
mv *.mm ./logs/
mv *.log ./logs/
mv *.nosex ./logs/

#<END>
	