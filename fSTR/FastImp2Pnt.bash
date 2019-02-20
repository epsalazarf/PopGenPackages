#FastImp2Pnt.sh
#Transforms IMPUTE2 files to fineSTRUCTURE input.

for chr in $(seq 1 22); do
    echo "Processing chromosome" $chr
    perl impute2chromopainter.pl $1${chr}.phs.haps $1${chr}
done 
