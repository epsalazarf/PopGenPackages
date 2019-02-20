# fineSTRUCTURE Tutorial

*By Pavel Salazar-Fernandez (epsalazarf@gmail.com)*

*Human Population and Evolutionary Genomics Lab | LANGEBIO*

## About
*Documentation: [PaintMyChromosomes.com](http://www.paintmychromosomes.com/)*

This document explains how to manipulate PLINK data for performing fineSTRUCTURE analysis and plotting.

## fineSTRUCTURE
fineSTRUCTURE is a fast and powerful algorithm for identifying population structure using dense sequencing data.  By using the output of ChromoPainter as a (nearly) sufficient summary statistic, it is able to perform model-based Bayesian clustering on large datasets, including full resequencing data, and can handle up to 1000s of individuals. Full assignment uncertainty is given.

For installation instructions, introductory tutorials and detailed manuals, check the software web page:

<http://www.maths.bris.ac.uk/~madjl/finestructure/finestructure.html>

Citation:
>Lawson, Hellenthal, Myers, and Falush (2012), "Inference of population structure using dense haplotype data", PLoS Genetics, 8 (e1002453)

### 1. Preparations

**Requires:**
- [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) (version 1.07+)

1. Get the PLINK binary files (.bed, .bim, .fam) to be analyzed. The dataset must contain only the samples of interest. Preferably, replace the samples' `FAMID`s with their population tags.
2. Link/copy the PLINK data files to a new directory.
3. Link/copy the human genetic map files to the new directory, called `genetic_map_chr[*]_combined_b37.txt`, one for each chromosome (Total: 22). These will be used for posterior data processing. If you don't have them, download the "HapMap phase II b37" file from the [SHAPEIT website](http://www.shapeit.fr/pages/m02_formats/geneticmap.html).

*Note:*
The authors recommend to re-label your samples to *[Population tag]+[number]*, such as `ABC01` for the first sample from the 'ABC' population. If your IDs need to be renamed , and the FAMID contains the POP tag, run to **replace** the FAM file:
`awk 'BEGIN{p="";n=0}(p != $1){p=$1;n=1}{if(n<10)$2=p"0"n;else $2=p""n;print $0;n+=1}' DATA.fam > tmp ; mv tmp > DATA.fam`

### 2. Phasing the data

**Input:**
- PLINK binary dataset: BED, BIM and FAM.
- Human genetic maps (`genetic_map_chr##_combined_b37.txt`).


**Requires:**
- [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) (version 1.07+)
- [SHAPEIT](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) (version 2+)

1. Ready the input files: put the PLINK binary dataset and the genetic maps in the same directory.
2. Split the PLINK dataset into chromosomes, one dataset for each one:
`plink --bfile DATA --chr [##] --make-bed --out DATA.chr##`
3. Start the SHAPEIT program with the following command:
`shapeit --input-bed DATA.chr##.bed DATA.chr##.bim DATA.chr##.fam -M genetic_map_chr##_combined_b37.txt -O DATA.chr##`
4. Repeat for each chromosome.
5. The following script automates the previous steps for the 22 chromosomes. Keep in mind that the phasing process can last hours:

```bash
# FastPhase.sh
# Splits a binary PLINK dataset into chromosomes and runs SHAPEIT to phase them.

for chr in $(seq 1 22); do
     plink --bfile $1 \
           --chr $chr \
           --make-bed \
           --out $1.chr$chr ;
     shapeit --input-bed $1.chr$chr.bed $1.chr$chr.bim $1.chr$chr.fam \
             -M /genetic_map_chr${chr}_combined_b37.txt \
             -O $1.chr${chr}
done 
```

**Command**:
`bash FastPhase.sh DATA`

**Output:**
- `DATA.chr##.haps` (for each chromosome).
- Other SHAPEIT files that will not be used for this pipeline.

### 3. Transform the phased data
**Input:**
- `DATA.chr##.haps` (for each chromosome).

**Requires:**
- [PERL 5.0+](http://www.cpan.org/)
- [impute2chromopainter.pl](http://www.maths.bris.ac.uk/~madjl/finestructure/manualse11.html#x14-2100011)

1. Check that `.haps` have been created correctly for all chromosomes.
2. Run the `impute2chromopainter.pl` script:
`perl impute2chromopainter.pl DATA.chr##.haps DATA.chr##`
3. Automated script:

```bash
#FastImp2fST.sh
#Transforms IMPUTE2 files to fineSTRUCTURE input.

for chr in $(seq 1 22); do
    perl impute2chromopainter.pl $1${chr}.haps $1${chr}
done 
```
**Command**:
`bash FastImp2fST.sh DATA`

**Output:**
- `DATA.chr##.phase` (for each chromosome).

### 4. Create the other input files
#### Sample list
The sample list contains no headers, a sample per row and three columns: sample ID (Label), population (POP) and a boolean value indicating if the sample should be analyzed (1 = yes, 0 =no). It is recommended that your samples are labeled as *[Population tag]+[number]*, such as `ABC01` for the first sample from the 'ABC' population (see example below):

|*Label*|*Population*|*Include?*|
|---|---|---|
|ABC01|ABC|1|
|ABC02|ABC|1|
|DEF01|DEF|1|
|...|...|...|

If your FAM file has been modified to contain the population tags (POP) in the FAMID column, you can easily create the sample list file running:
`awk '{print $2,$1,"1"}' DATA.fam > DATA.ids`

#### Recombination maps
1. Check that the BIM and genetic map files are in the same directory.
2. Generate the recombination file with the following scripts:
`awk 'BEGIN{print"start.pos","recom.rate.perbp"}
  NR==FNR{a[$1]; if(NR>1) print $1,$2; next}
  (($4 in a)==0){print $4,"NA"
  }' genetic_map_chr${chr}_combined_b37.txt \
  DATA.chr##.bim | sort -n -k 1 > DATA.chr##.recomb`
3. Impute the missing data with the following R script:

```r
# IMPUTE_recomb_rates.R
# Linearly imputes recombination rates for SNPs with missing data.

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
x <- read.table(file, header=T)

library(zoo)
x[,2] <- na.approx(x[,2],x[,1], na.rm=F)
x[is.na(x[,2]),2] <- 0

write.table(x,file,row.names = F, quote = F, eol = "\n")
```

**Command**:
`Rscript IMPUTE_recomb_rates.R DATA.chr##.recomb`

4. Filter the recombination map file to enlist only those SNPs present in the dataset:
```bash
awk 'NR==FNR{a[$4]; next}(FNR==1 || $1 in a){print $0
}' DATA.chr##.bim DATA.chr##.recomb > tmp
mv tmp DATA.chr##.recomb
```

5. Automated:
*(Requires Step 2 script to be saved as a file named `IMPUTE_recomb_rates.R`)*

```bash
# MakeRecombFiles.sh
# Extracts and imputes recombination rates for all SNPs in a PLINK dataset.

for chr in $(seq 1 22); do
  awk 'BEGIN{print"start.pos","recom.rate.perbp"}
    NR==FNR{a[$1]; if(NR>1) print $1,$2; next}
    (($4 in a)==0){print $4,"NA"
    }' genetic_map_chr${chr}_combined_b37.txt \
    DATA.chr##.bim | sort -n -k 1 > DATA.chr##.recomb
  echo "Imputing recombination map for chromosome " ##
  Rscript IMPUTE_recomb_rates.R DATA.chr##.recomb
  awk 'NR==FNR{a[$4]; next}(FNR==1 || $1 in a){print $0
  }' DATA.chr##.bim DATA.chr##.recomb > tmp
  mv tmp DATA.chr##.recomb
done
```

**Command**:
`bash MakeRecombFiles.sh DATA`

**Output:**
- `DATA.ids`
- `DATA.chr##.recomb` (for each chromosome).

### 5. Run fineSTRUCTURE
>Note: It is highly recommended that the user completes the examples given within the package to get a better understanding of the program. Check the tutorial here: [fineSTRUCTURE:  How to use this software](http://www.maths.bris.ac.uk/~madjl/finestructure/manualse2.html#x5-40002)

**Input:**
- `DATA.ids`
- `DATA.chr##.phase` (for each chromosome).
- `DATA.chr##.recomb` (for each chromosome).

**Requires:**
- [fineSTRUCTURE v2](<http://www.maths.bris.ac.uk/~madjl/finestructure/finestructure.html>)

1. Download and install fineSTRUCTURE. Detailed instructions and troubleshooting is available at its [webpage](http://www.maths.bris.ac.uk/~madjl/finestructure/finestructure.html).
```bash
tar -xzvf fs-2.0.4.tar.gz # (check correct version)
cd fs-2.0.4
./configure
make
sudo make install # optional; installs the fs executable
```
2. Run fineSTRUCTURE. Note that all phase and recombination files for EACH chromosome must be explicitly declared. This analysis may take hours depending of the data.

`fs DATA.cp -idfile DATA.ids -phasefiles DATA.chr#1.phase DATA.chr#2.phase -recombfiles DATA.chr#1.rec DATA.chr#2.rec -go`

**Output:**
- DATA_linked.chunkcounts.out
- DATA_linked_mcmc.xml
- DATA_linked_tree.xml

### 6. Fully automated script
Provided with a PLINK binary dataset and the genetic maps, and having installed the fineSTRUCTURE package and all other software dependencies, this script will automate all the previously described steps.

NOTE: Be sure that the `IMPUTE_recomb_rates.R` script is also available to the script.

**Input:**
- PLINK binary dataset: BED, BIM and FAM.
- Human genetic maps (`genetic_map_chr##_combined_b37.txt`).- 

**Command:**
`bash RUN_fs_pipeline.sh DATA`

**Output:**
- DATA_linked.chunkcounts.out
- DATA_linked_mcmc.xml
- DATA_linked_tree.xml

## Visualizing fineSTRUCTURE Output
### GUI Application
The GUI application is the easiest and most effective tool for visualizing the fineSTRUCTURE output. It can generate different types of plots and has many easy-to-use customization tools to modify them as needed.

**Download and installation instructions:**
<http://people.maths.bris.ac.uk/~madjl/finestructure/finestructure.html>

**Input:**
- DATA_linked.chunkcounts.out
- DATA_linked_mcmc.xml
- DATA_linked_tree.xml

**Output:**
- Co-ancestry heat map (PNG exportable)
- Pairwise coincidence matrix (PNG exportable)
- PCA plots (PNG exportable)
- PCA values (CSV)

### R-based Plotting Library
R functions for fineSTRUCTURE results are provided by the software creators, with an example script. It requires intermediate R-programming skills to be used effectively, as a customized R script has not been developed yet.

**finestructureR:**
<http://people.maths.bris.ac.uk/~madjl/finestructure/finestructureR.html>

### R Shiny: fineStructure PCA Plot
While the GUI provides a PCA plotting function, the interface does not allow for a more detailed exploration of the results.

`fSTR/app.r`, a customized R Shiny app, adds interactivity to the PCA plot and allows for customization and zooming.

#### Saving a fineSTRUCTURE PCA as a CSV matrix

1. Open the Finestructure GUI.
2. Go to the menu File > Manage files. In the 'Raw data file', click 'Change data file' and select the file with the `.chunkscount.out` suffix. 
3. Click the 'Read data file' to oad the data file. If the file is loaded correctly, a co-ancestry matrix is plotted. 
4. Start the PCA Plot module (Plot > Principal Components Analysis). A new window will show a preview of the PCA data.
5. Click the 'Export CSV Data' button to the left and save the file.
6. Run the `app.R` (more info below). Select the CSV file in the 'Open File' window to load the data. 
7. A plot should appear. You can explore and customize the plot with the tool bar at the left panel.
8. Save the current plot visualization to a PNG file with a right-click on top of the plot area.

#### Running a R Shiny App
To open a Shiny app, you have these alternatives:
* Open [R Studio](http://www.rstudio.com/), load the `app.R` file and then click the button `Run`. A web browser window must open.
* Copy the entire code of `app.R` as text and paste it in an R console. A web browser window containing the app must open.

If the app is unable to open or breaks during operation, close the app window, use the `STOP` button to end all current processes, and restart it.

Once the script loads, a "Open File" window should appear. Select the CSV file and a PCA should plot in a browser window. Use the options in the left panel to customize the plot.

**Features:**
- Plot types: Can select between points or tags for the plot.
- Select Population: Shows only a particular population.
- Emphasize Population: Highlights points or tags for a chosen population.
- Show/Hide Legend: Displays all values from the selected category.
- Interactive Zoom: select an area and double click to zoom in, double click again to zoom out.
- Area Info: brush an area to see the selected points info and coordinates.

**Output:**
- Plots can be saved by right-clicking the image. Default height: 900px.

___
