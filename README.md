# GEN_711_Final_Project

## GEN 711 Final Project

## Mystery Bacterial Genome

## Member 

Andrew Ricketts

## Background
Four mystery DNA sequences were provided that represent two different bacterial individuals, two files per bacterium. The data was provided in FASTQ files and will be processed to construct the genome of the bacteria to identify the genus and species. 
## Methods
The data was provided via the mystery bacterial genome tutorial. The files are 250bp paired-end Illumina HiSeq 2500 reads. The analysis was conducted on the RON server. 
### FastQC
FastQC runs a quality control of provided fastq files on both forward and reverse reads separately. The ouput of the program is an html file for each fastq file which can be accessed through a web browser. It provides multiple statistics and blots to show average quality and gc content of reads. 
### Trimmomatic
Trimmomatic trims off low quality bases, typically near the ends of reads, and adapter sequnces from the reads. Fastq files of both forward and reverse reads are provided as the input. The ouput is 4 fastq files: two files for the paired forward and reverse reads; and two unpaired files for both forward and reverse. 
### SPAdes
SPAdes is a genome assembler best used for bacteria. The four trimmed sequence files are inputted into the program which then aligns the reads and attempts to construct a genome. The output is multiple versions of a contigs fasta file which has the assembeled contigs. Multiple contigs are often produced when they cannot link with other contigs or there is contamination. 
### QUAST
QUAST assesses quality of the genome asembly. It looks at genome fragmentation, whether one or a few contigs is prepresentative of the genome or if there are manby small contigs. It intakes the contigs fasta file and provides quality reports in multiple file types. 
### BUSCO
BUSCO assesses the completeness of the assembled genome in the contigs fasta file. The program compares a large set of orthologous genes which 90% of all orgranisms in a certain group contain, the group in this analysis was bacteria. A complete genome would have most of these genes. The output is multiple files which report the completeness. 
### PROKKA
PROKKA is a genome annotation tool for prokaryotes which analyses the genome for where gene sequences are likely. It reports the locations and names of known genes, unknown predicted genes, and different RNAs. The contigs fasta file is inputted and multiple files are produced which each show different data like nucleotide sequence, amino acid sequence, etc.
### Extract Sequences
The extract_sequences program is utilized to find a certain gene or RNA in PROKKA ffn files. It helps extract certain sequences to compare with known organism sequences. The name of the gene or RNA is identified to find the correct sequence. The output can be put into a fasta file. 
### BLAST
#### BLAST Database
makeblastdb is a program which creates a local database of seqeunce data to quickly compare sequences with other org anisms. 
#### Local BLAST
One option for the blastn program is one against the local database. A gene or RNA seqeunce in a fasta file can be provided into the program which will compare the provided sequence with that of other organisms in the database. The output type must be provided, here a tsv file was used. Likely organism IDs are provided in the output.
#### Entire Assembly BLAST
The authors of the tutorial provided the blob_blast.sh program which BLAST's the entire genome against the local database and prepares a file which can later be used in the blob tools program. The input is the full genome fasta file and outputs a megablast.out file. 
### BWA and Samtools
BWA is a read mapping program, meaning reads are aligned to a given genome assembly and analyses the amount of reads (depth) at each position. The input is the reference genome assembly, here it was the assembled genomes, and the trimmed forward and reverse reads. The output is a sam file which is then inputted into samtools to be converted into a bam file, a binary version of a sam file. 
### Blobtools
Blobtools is a group of mutliple programs which can help visualize different parts of the genome like coverage and GC content. THe programs utilize the bam file and local blast file to construct the graphs. 
### Filter Contigs
filter_contigs_by_list.py was a program provided by the authors which filters and gets rid of contigs which don't fit set criteria for the genome. Grep is used to create a list of contigs to keep from one of the blobtools files in a new txt file. The genome fasta file and the txt file are provided into the program to create a new fasta file with only the contigs kept. 
### UniVec BLAST 
UniVec is a database of DNA sequences of various vectors. blastn takes the filtered fasta file and compares it to the UniVec database to check for contamination. The ouptut is a .6 file. 
## Findings

## References
Andrews, S. (2010).  Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics - FastQC a quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170 
Camacho, C. (2024, June 25). Blast+ Release notes. BLAST® Help [Internet]. https://www.ncbi.nlm.nih.gov/books/NBK131777/ 
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2). https://doi.org/10.1093/gigascience/giab008 
Laetsch, D. R., & Blaxter, M. L. (2017). Blobtools: Interrogation of genome assemblies. F1000Research, 6, 1287. https://doi.org/10.12688/f1000research.12232.1 
Li, H., & Durbin, R. (2010). Fast and accurate long-read alignment with Burrows–Wheeler transform. Bioinformatics, 26(5), 589–595. https://doi.org/10.1093/bioinformatics/btp698 
Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., & Zdobnov, E. M. (2021). Busco update: Novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular Biology and Evolution, 38(10), 4647–4654. https://doi.org/10.1093/molbev/msab199 
Mikheenko, A., Prjibelski, A., Saveliev, V., Antipov, D., & Gurevich, A. (2018). Versatile Genome Assembly evaluation with Quast-LG. Bioinformatics, 34(13), i142–i150. https://doi.org/10.1093/bioinformatics/bty266 
Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using Spades de Novo Assembler. Current Protocols in Bioinformatics, 70(1). https://doi.org/10.1002/cpbi.102 
Seemann, T. (2014). Prokka: Rapid Prokaryotic Genome Annotation. Bioinformatics, 30(14), 2068–2069. https://doi.org/10.1093/bioinformatics/btu153 
