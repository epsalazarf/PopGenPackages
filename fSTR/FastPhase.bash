#FastPhase.sh
#Splits a binary PLINK dataset into chromosomes and runs SHAPEIT to phase them.

for chr in $(seq 1 22); do
     echo "Processing chromosome" $chr
     plink --bfile $1 \
           --chr $chr \
           --make-bed \
           --out $1.chr$chr ;
     shapeit --input-bed $1.chr$chr.bed $1.chr$chr.bim $1.chr$chr.fam \
             -M /data/reference_panels/genetic_map_b37/genetic_map_chr${chr}_combined_b37.txt \
             -O $1.chr${chr}.phs
done 
