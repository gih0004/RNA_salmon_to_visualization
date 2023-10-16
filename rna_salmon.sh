module load salmon

# make index using transcriptome fasta (used high confidence model)
salmon index -t ATL_v3.hc_gene_models.repr.cdna.fa.gz -i salmon_ATL.v3.index
mkdir salmon_quant
for fq1 in ./*.bam #<>This should be changed into whichever string you have last in sample names common between all samples 
do

    Seq="${fq1%.bam}"

    salmon quant -i salmon_ATL.v3.index -l A -1 ${Seq}.filtered_1.fq -2 ${Seq}.filtered_2.fq --gcBias -o ./salmon/${Seq}


done
