# RNASEQ Salmon to  Visualization
Contains scripts used to generate Heatmaps for TPM expression for transcritps/loci of interest from filtered paired-end RNA reads 
For this repository, put all filtered fastq files within a single directory containing additionally your species reference trancriptome (cdna_file.fa )

This repository contains two scripts: 
1. rna_salmon.sh : for use within a HPC, generates expression values from filtered fastq files ; can be generated with Step 1 of my [RNA-seq pipeline](https://github.com/gih0004/RNA_Seq_featurecounts).
2. Salmon_Visualization: R script, this is used to generate expression heatmaps for transcript/loci of interest
   
## Rna_salmon.sh
To quantifiy expression, make index using transcriptome(cDNA) fasta file
```ruby
salmon index -t {reference_cdna} -i {index_name}.index
```
Now that an index has been made, we can use this index to quantify the expression of reads in TPM. This step generates files that have other values within (not just TPM), but we will focus on using TPM values for visualization. If you desire to use another value it is also possible with minor modifications to the salmon visualization R script. For more info on other values : https://www.youtube.com/watch?v=TTUrtCY2k-w
```ruby
for fq1 in ./*.bam #<>This should be changed into whichever string you have last in sample names common between all samples 
do

    base="${fq1%.bam}"

    salmon quant -i salmon_ATL.v3.index -l A -1 ${base}.filtered_1.fq -2 ${base}.filtered_2.fq --gcBias -o ./salmon/${base}
done
```
Salmon’s main output is its quantification file. This file is a plain-text, tab-separated file with a single header line (which names all of the columns). This file is named quant.sf and appears at the top-level of Salmon’s output directory.


## Salmon_visualization_Rstudio

Now that we have the desired expresion files <filename>quant.sf within individual directories ./salmon/<filename>, the indiviual quantification files for each of your samples need to be downloaded and all put within one local directory for correct use in Rstudio. 

To use the rest of the Salmon_visualization script, if your organism does not have a reference TxDb within tximport, it is neccesary to create a dataframe similar to: 
|  Transcript_id   | Gene_id          |
|------------------|------------------|
| Transcript_A_1   |      Gene_A      |
| Transcript_B_1   |      Gene_B      |

Alternatively, you can make a TxDb from your GTF file following the instructions found in this [docuumentation](https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/building-a-txdb-objecthttps://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/building-a-txdb-object)

Once that is completed, we can then use this file for our visualization 


We first set the directory containing all files and load the required packages to visualize. 

dir <- ("C:/Users/")  # Replace with your directory path

```ruby
# Load required packages
library("RColorBrewer")
library("tximport")
library("GenomicFeatures")
library("pheatmap")
```

Reading in data for visualization
```ruby
mySalmons <- list.files(pattern=glob2rx("*.sf"), # Create a vector containing the names of your files
                        all.files=T, 
                        full.names=F)


#Gives names to the samples from mysalmons, make sure this in in the same order as files in your directory 
names(mySalmons)<- c("Sample1","Sample2")
#Quality Check 
names(mySalmons)

#Read in the samples metadata
samps <- read.csv("RNAmetadata.csv")  #change with file for metadata
samps

#This can be used to change the metadata variable names
#samps$<Factor> <- gsub("Metadata.Variable", "Arbitrary string name", samps$<Factor>) 

```

This step is where we load in the dataframe created earlier (or imported in from TxDb in Tximport) and  makes an expression dataframe (cts_filteresd) 


Options used for reading salmon outputs:
TxOut is for downstream analysis 
countsfromambundance specifies what value to use character: either "no" (default), "scaledTPM", or "lengthScaledTPM".

```ruby
txg <- read.table("gtf.csv", header = FALSE, sep = ",") #data frame that has two columns : transcript and gene id 
colnames(txg) <- c("transcript_id", "gene_id")
txg # df that one columns has all transcript id and the second one, the related gene_id


#Read into R dataframe the Salmon Output
txi <- tximport(mySalmons, type="salmon", tx2gene=txg , txOut=TRUE, countsFromAbundance="lengthScaledTPM") #TxOut is for downstream analysis and countsfromambundance specifies what value to use 
txi #shpuld be a list of a couple of 4 elements for each of the observations (samples)


#Using the Txi object to create a dataframe 
cts <- txi$counts #making abundance value TPM df from salmon output as a variable 
cts


#keeps only those that in one sample have more than 0 
cts_filtered <- cts[rowSums(cts != 0, na.rm = TRUE) > 0, ]


cts_filtered #expression dataframe with trancript id with more than 0 as a value and in this case without unanchored transcripts ; i took out the unanchored because theres no gen_id for them 

```



the first is the trasncript_id or the feature Salmon used to quanitfy TPM and the second column should be gene_id. This will be used to find transcprits within the samples and associate them with an arbitrary name to later on visualize using a heatmap. 

The output of this should be a heatmap similar to:
![image](https://github.com/gih0004/RNA_salmon_to_visualization/assets/114354096/ba58e5ae-a293-41d4-b7b0-eee26538ac62)
