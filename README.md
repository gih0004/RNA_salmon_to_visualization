# RNA_salmon_to_visualization
Contains scripts used to generate Heatmaps for TPM expression for transcritps/loci of interest from filtered paired-end RNA reads 
For this repository, put all filtered fastq files within a single directory containing additionally your species reference trancriptome (cdna_file.fa )

This repository contains two scripts: 
1. rna_salmon: for use within a HPC, generates expression values from filtered fastq files ; can be generated with my RNA_seq_pipeline
2. Salmon_Visualization: R script, this is used to generate expression heatmaps for loci of interested 
   
### Rna_salmon.sh
To quantifiy expression, make index using transcriptome fasta (used high confidence model)
```ruby
salmon index -t {species reference cdna file} -i {desired_index_name}.index
```
Now that an index has been made, we can use this index to quantify the expression of reads in TPM. This step generates files that have other values within (not just TPM), but we will focus on using TPM values for visualization. If you desire to use another value it is also possible with minor modifications to the salmon visualization R script. For more info on other values : https://www.youtube.com/watch?v=TTUrtCY2k-w
```ruby
for fq1 in ./*.bam #<>This should be changed into whichever string you have last in sample names common between all samples 
do

    base="${fq1%.bam}"

    salmon quant -i salmon_ATL.v3.index -l A -1 ${base}.filtered_1.fq -2 ${base}.filtered_2.fq --gcBias -o ./salmon/${base}
done
```

Now we have the desired expresion files <filename>quant.sf within individual directories ./salmon/<filename>. The indiviual quantification files need to be downloaded and all put within one local directory for correct use in Rstudio. 

### Salmon_visualization_Rstudio

We first set the directory containing all files and load the required packages to visualize. 

dir <- ("C:/Users/")  # Replace with your directory path

```ruby
# Load required packages
library("RColorBrewer")
library("tximport")
library("GenomicFeatures")
library("pheatmap")
```

To use the rest of the Salmon_visualization script it is neccesary to have the correct metadata file for your samples and a csv file (excel) file where you have two columns: the first is the trasncript_id and the second column should be gene_id. This will be used to find transcprits within the samples and associate them with an arbitrary name to later on visualize using a heatmap. 

The output of this should be a heatmap similar to:
![image](https://github.com/gih0004/RNA_salmon_to_visualization/assets/114354096/ba58e5ae-a293-41d4-b7b0-eee26538ac62)
