# RNA-seq-pipeline
This is a RNA seq pipeline developed in snakemake work engine. All the file names here are general and no data is added because of privacy of our collaborators data.Add your data folder with file name as mentioned in metadata.tsv. The organims is Mycobacterium Tuberculosis Strain H37Rv and it is advised to change the names and folders whereever required.
# Additional requirements:
Perform initial quality control check for the reads.
Install snakemake, fastp , bowtie2 , samtools , featureCounts and R
# Command
snakemake --cores 8
