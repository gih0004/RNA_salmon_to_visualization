module load salmon

# make index using transcriptome fasta (used high confidence model)

mkdir salmon_quant

salmon index -t {species reference cdna file} -i {desired_index_name}.index

for fq1 in ./*.bam #<>This should be changed into whichever string you have last in sample names common between all samples 
do
    base="${fq1%.bam}"
#Generate transcript expression quantification files for each sample 
    salmon quant -i salmon_ATL.v3.index -l A -1 ${base}.filtered_1.fq -2 ${base}.filtered_2.fq --gcBias -o ./salmon/${base}


done
