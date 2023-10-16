###### SALMON TPM VALUES TO HEATMAP GENE IDS OF INTEREST ####### 

dir <- ("C:/Users/")  # Replace with your directory path


# Load required packages
library("RColorBrewer")
library("tximport")
library("GenomicFeatures")
library("pheatmap")

# Create a mapping table using tximport
mySalmons <- list.files(pattern=glob2rx("*.sf"), # Create a vector containing the names of your files
                        all.files=T, 
                        full.names=F)



#Gives names to the samples from mysalmons, make sure this in in the same order as files in your directory 
names(mySalmons)<- c("POT_AA", # Name the elements of m fread, header = TRUE, sep = "\t")yData.
                        "POT_AB","POT_AC","POT_AD","POT_AE","POT_AF","POT_AG","POT_AH","POT_AI",
                        "POT_AJ","POT_AK","POT_AL","POT_AM","POT_AN","POT_AO","POT_AP","POT_AQ",
                        "POT_AR","POT_AS","POT_AT","POT_AU","POT_AV","POT_AW","POT_AX","POT_AY",
                        "POT_AZ","POT_BA","POT_BB","POT_BC","POT_BD","POT_BE","POT_BF","POT_BG",
                        "POT_BH","POT_BI","POT_BJ","POT_BK","POT_BL","POT_BM","POT_BN","POT_BO",
                        "POT_BP","POT_BQ","POT_BR","POT_BS","POT_BT","POT_BU","POT_BV","POT_BW")
#Quality Check 
names(mySalmons)


#Read in the samples metadata
samps <- read.csv("RNAexpressionmetadata.csv")  #change with file for metadata
samps  #Quakity check, should have a dataframe with observation # = to samples and variable # = to metadata columns 


#Now, to use salmon sf outputs, we need to generate a dataframe that links EACH trasncript ids to gene ids. 
#If the gtf file has trasncript id and gene id in *ALL* the rows for features, you can use the txg. I had to create one in ARC using as mine had gaps 
txg <- read.table("gtf.csv", header = FALSE, sep = ",") #data frame that has two columns : transcript and gene id 
colnames(txg) <- c("transcript_id", "gene_id")
txg # df that one columns has all transcript id and the second one, the related gene_id 


#Load Salmon Output
txi <- tximport(mySalmons, type="salmon", tx2gene=txg , txOut=TRUE, countsFromAbundance="lengthScaledTPM") #TxOut is for downstream analysis and countsfromambundance specifies what value to use 
txi #shpuld be a list of a couple of 4 elements for each of the observations (samples)


#Using the Txi object to create a dataframe 
cts <- txi$counts #making abundance value TPM df from salmon output as a variable 
cts
cts_filtered <- cts[rowSums(cts != 0, na.rm = TRUE) > 0, ] #keeps only those that in one sample have more than 0 

cts_filtered <- cts[rowSums(cts != 0, na.rm = TRUE) > 0 & !grepl("^Soltu.Atl_v3.S", rownames(cts)), ] #filters out the counts from transcripts unanchored (the ones that are ATL.v3_s<number<)
cts_filtered #expression dataframe with trancript id with more than 0 as a value and in this case without unanchored transcripts ; i took out the unanchored because theres no gen_id for them 

#Quality Check
cts_filtered["Soltu.Atl_v3.02_4G005370.2",]  # change with a transcript id of interest to view expression value per sample 


names(cts_filtered)<- c("POT_AA", # Name the elements of m fread, header = TRUE, sep = "\t")yData.
                        "POT_AB","POT_AC","POT_AD","POT_AE","POT_AF","POT_AG","POT_AH","POT_AI",
                        "POT_AJ","POT_AK","POT_AL","POT_AM","POT_AN","POT_AO","POT_AP","POT_AQ",
                        "POT_AR","POT_AS","POT_AT","POT_AU","POT_AV","POT_AW","POT_AX","POT_AY",
                        "POT_AZ","POT_BA","POT_BB","POT_BC","POT_BD","POT_BE","POT_BF","POT_BG",
                        "POT_BH","POT_BI","POT_BJ","POT_BK","POT_BL","POT_BM","POT_BN","POT_BO",
                        "POT_BP","POT_BQ","POT_BR","POT_BS","POT_BT","POT_BU","POT_BV","POT_BW")
#Salmon reads : 
range(colSums(cts)/1e6) # column sums are equal to the number of mapped paired end reads per experiment
#3#[1] 12.28212 25.35258 ; hence, the experiment has between 12-25 million paired end reeads that were mapped to the transcriptome using salmon

#generate z scores 
cts_zscore<-t(apply(cts_filtered,1,scale))
colnames(cts_zscore) <- colnames(cts)
cts_zscore  

#Now you have a dataframe with expression values for all genes and can generate a heatmap
#You also have Z score values for all genes


###### Visualization ######


#Create dataframe with annotation info based on metadata
annot_info <- samps[,c('Treatment','Tissue','Timepoint')] #this Vector should be the different columns of interest in metadata
 
rownames(annot_info) <- names(mySalmons)  # gives names to rows instead of sample numbers 


#Load Transcripts you want to graph as a csv gile (can create in excel)
#The csv file should be two columns : 1) transcript_id 2) gene_name  
mapid_Df <- read.csv("C:/Users/guill/Desktop/Research/ATL_DiffExp_Data/featureCounts_Rstudio/Tuber_transcript_id.csv")

#Generate a vector with transcript_id of the genes of interest to look for within cts_filtered df (which contains transcript ids and expression level of each)
tuber_loci <- c(mapid_Df$transcript_id)

#Vector that contains the gene names for the transcript ids to be graphed 
tuber_loci_gene<-c(mapid_Df$gene_name)

# #subsets cts_filtered based on tuber loci 
subset_tuberdf <- cts_filtered[rownames(cts_filtered) %in% tuber_loci,]

# Extract the transcript IDs from mapid_Df and gene names
transcript_ids <- mapid_Df$transcript_id
gene_names <- mapid_Df$gene_name

# Create a vector to store the new row names
new_row_names <- character(length = nrow(subset_tuberdf))

# Loop through transcript IDs and gene names
for (i in 1:length(transcript_ids)) {
  # Find the matching row index in subset_tuberdf
  matching_row_index <- which(rownames(subset_tuberdf) == transcript_ids[i])
  
  # If a match is found, update the new row name
  if (length(matching_row_index) > 0) {
    new_row_names[matching_row_index] <- gene_names[i]
  }
}

# Update row names in subset_tuberdf
rownames(subset_tuberdf) <- new_row_names


# Reorder subset_tuberdf based on the new_row_names
subset_tuberdf <- subset_tuberdf[order(new_row_names), ]


# Create a sequential color palette with 5 colors from the "Blues" palette
my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(9) #this is to map the expression values
my_colors <- colorRampPalette(brewer.pal(9, "BuPu"))(9) # this is eventually used for ann colors

ann_colors = list(
  Treatment = c("ElevatedT"="#FF4500", "AmbientT"= "cyan"),
  Tissue = c("Leaf"="#008000" , "Small.Tuber"="gold" , "Medium.Tuber"="darkgoldenrod2" , "Big.Tuber"= "darkgoldenrod4"), 
  Timepoint = c("30.days"= my_colors[1], "60.days" = my_colors[3] , "90.days" = my_colors[6],"120.days"= my_colors[9])
)


pheatmap(subset_tuberdf,color= my_palette, cluster_rows=FALSE,show_rownames=TRUE, show_colnames = FALSE,
         cluster_cols=FALSE,annotation_col = annot_info,
         cutree_rows = 4 , 
         annotation_colors = ann_colors,
         angle_col ="315",
         gaps_row= c(4,8,11,12,13,16,23,24))
