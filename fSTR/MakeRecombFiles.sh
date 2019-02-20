#MakeRecombFiles.sh
#Makes custom recombination maps from BIM files to use in fineStructure

for chr in $(seq 1 22); do
	echo "Processing chromosome" $chr
	awk 'BEGIN{print"start.pos","recom.rate.perbp"}NR==FNR{a[$1]=$2;next}{if($4 in a) print $4,a[$4];else print $4,"0"}' genetic_map_chr${chr}_combined_b37.txt $1.chr${chr}.bim > $1.chr${chr}.rec
done 


awk 'BEGIN{print"start.pos","recom.rate.perbp"}NR==FNR{a[$1]=$2; if(NR>1)print $1,$2;next}{if($4 in a) print $4,a[$4];else print $4,"NA"}' genetic_map_chr22_combined_b37.txt Bot+PACcghp.s.chr22.bim | sort -n -k 1 > BotPAC.chr22.recmap+.recomb