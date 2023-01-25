Manipulating VCF files
#Qn 1. Describe the format of the file and the data stored
#The VCF (variant calling format) file is a standard ileformat for storing variation data.
#It is used by large scale variant mapping projects and is also the standard output of variant calling software such as GATK 
#and the standard input for variant analysis tools such as the VEP or for variation archives like EVA.
#It is a text file containing meta-information lines (included after the ## string), a header line (included after #), and then data lines 
#Each data line contains information about a position in the genome. 
#The format can also contain genotype information on samples for each position.
#A single ‘fileformat’ field is always required, must be the first line in the file
# It details the VCF format version number. For example, for VCF version 4.2, as for this specific file

#Qn 2. What does the header section of the file contain
#The header line names the 8 fixed, mandatory columns which are;
#CHROM - chromosome
#POS - reference position
#ID - identifier
#REF - reference bases
#ALT - alternate bases 
#QUAL - quality
#FILTER - filter status
#INFO - additional information


#Qn 3. How many samples are in the file
#this command lists the sample names in the file

bcftools query -l sample.vcf

#this command counts the number of samples in the file
#there are 6 samples in the vcf file
bcftools query -l sample.vcf  | wc -l

#Qn 4. How many variants are in the file
  #using grep
grep -v '^#' sample.vcf | wc -l # there are 398246 variants
  
  #using bcftools
bcftools query -f '%ALT\n' sample.vcf | wc -l
bcftools view -H sample.vcf | wc -l

#Qn 5. How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? Save the output to a tab-delimited file
#this command saves the output in a text file called cpqr.txt
bcftools query -f '%CHROM\t%POS\t%INFO/QD\t%INFO/MQ\n' sample.vcf > cpqr.txt

#6. Extract data that belongs to chromosomes 2,4 and MT
# this command extracts the data and saves it to a text file called extract.txt
awk '$1=="2" || $1=="4" || $1=="MT"' sample.vcf > extract.txt 

#Qn 7. Print out variants that do not belong to chr20:1-30000000
grep -v '^##' sample.vcf > minus#.txt
awk '$1 != "20" || ($1 == "chr20" && ($2 < 1 || $2 > 30000000)) {print $1, $2, $4, $5}' minus#.txt > chr20.txt


#Qn 8. Extract variants that belong to SRR13107019
#this command extracts the variants and saves to a file called SRR7019.tt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -s SRR13107019 sample.vcf > SRR7019.txt

#Qn 9. Filter out variants with a QualByDepth above 7
#this command saves variants with a QualByDepth above 7 in a filve called variants_7.vcf
awk -F '\t' '{if ($6>=7) print $0}' sample.vcf > variants_7.vcf


#Qn 10. How many contigs are referred to in the file. Check the header section
#cheching number of contigs
grep -c "^##contig" sample.vcf #there are 2211 contigs
#or
bcftools view -h sample.vcf | grep "#contig" | wc -l

#11. Comment on the eighth and ninth columns of the file
#The eight column is the additional information column (INFO). 
#INFO fields are encoded as a semicolon-separated series of short keys with optional values in the format: <key>=<data>[,data] 
#The exact format of each INFO sub-field should be specified in the meta-information 

#The ninth column, called FORMAT, is an additional column in a vcf file containing genotype information
#The FORMAT field specifies the data types and order (colon-separated alphanumeric String)
#The first sub-field must always be the genotype (GT) if it is present
#As with the INFO field, there are several common, reserved keywords that are standards across the community
#For example GT for genotype, DP for read depth, etc


#Qn 12. Extract data on the read depth of called variants for sample SRR13107018
#data on read depth for SRR13107018 saved in alleles.vcf
bcftools query -f '%DP\n' -s SRR13107018 sample.vcf > SRR7018.vcf


#Qn 13. Extract data on the allele frequency of alternate alleles. Combine this data with the
#chromosome and position of the alternate allele
#data on allele frequency
bcftools query -f '%CHROM\t%POS\t%AF\n' sample.vcf >alleles.vcf


Manipulating SAM files
#Qn 1. Describe the format of the file and the data stored
#SAM stands for Sequence Alignment/Map format, and it is a TAB-delimited text format 
#It consists a header section, which is optional, and an alignment section. 
#If present, the header must be prior to the alignments. Header lines start with ‘@’, while alignment lines do not. 
#Each alignment line has 11 mandatory fields for essential alignment information such as mapping position
#It can also have a variable number of optional fields for flexible or aligner specific information.
#SAM files are a type of text file format that contains the alignment information of various sequences that are mapped against reference sequences. 
#These files can also contain unmapped sequences. 
#Since SAM files are a text file format, they are more readable by humans and will be used as the examples for this section.


#Qn 2. What does the header section of the file contain
#Each header line begins with the character ‘@’ followed by one of the two-letter header standardized record type codes for example @HD, @SQ, etc. 
#In the header, each line is TAB-delimited and, apart from @CO lines, each data field follows a format ‘TAG:VALUE’ 
#where TAG is a two-character string that defines the format and content of VALUE. 
#Thus header lines match /^@(HD|SQ|RG|PG)(\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ or /^@CO\t.*/. 
#Within each (non-@CO) header line, no field tag may appear more than once and the order in which the fields appear is not significant.


#3. How many samples are in the file
#using samtools
samtools view -H sample.sam | grep -c '@RG'

#Qn 4. How many alignments are in the file
#counting number of alignments in the file
samtools view -c -F 4 sample.sam

#Qn 5. Get summary statistics for the alignments in the file
#this command gets the statistics and saves them to a text file callled stats.txt
samtools flagstat sample.sam > stats.txt


#Qn 6. Count the number of fields in the file
#number of fields in the file. this command gives 18 fields
grep -v "^@" sample.sam | awk -F'\t' '{print NF; exit}' 
#Qn 7. Print all lines in the file that have @SQ and sequence name tag beginning with NT_
#this command saves the lines in a text file names sq.txt
grep '@SQ.*NT_' sample.sam > sq.txt

#Qn 8. Print all lines in the file that have @RG and LB tag beginning with Solexa
#this command saves the lines in a text file called rg.txt
grep '@RG.*LB:Solexa' sample.sam > rg.txt

#Qn 9. Extract primarily aligned sequences and save them in another file
#this command saves the sequences in a sam file called primary_alignment.sam
awk '$1 !~ /^@/ && $2 == "99" || $2 == "83"' sample.sam > primarily_sequences.sam

#this command saves the sequences in a sam file called primary_alignment.bam
samtools view -f 2 -b sample.sam > primary_alignment.bam

#Qn 10. Extract alignments that map to chromosomes 1 and 3. Save the output in BAM
#format
# this command saves the extracted data in a bam file called filtered.bam
grep -E "^[^@]*\t(1|3)\t" sample.sam | samtools view -bS - > filtered.bam

#Qn 11. How would you obtain unmapped reads from the file
#this command saves the reads in a sam file called unmapped_reads.bam
samtools view -f 4 sample.sam > unmapped_reads.sam
#this command converts the unmapped reads sam file to bam format
samtools view -bS unmapped_reads.sam > unmapped_reads.bam

#Qn 12. How many reads are aligned to chromosome 4
grep -c "^[^@]\s*4\s" sample.sam #using this command, it shows that there are no reads aligned to chromosome 4

#13. Comment of the second and sixth column of the file
#The secong column represents the FLAG column. This column contains data which shows the properties of the alignment
#For example; paired or single end, mapped or unmapped, primary or secondary alignment, etc

#The sixth column is the CIGAR String column
#This column has data in form of letters to show how the samples aligned to the reference and at what length
#For example, M means alignment matches, X means sequence miss match, N means skipped region, etc 

#14. Extract all optional fields of the file and save them in “optional_fields.txt”
awk '{for(i=11;i<=NF;i++)print $i}' sample.sam > optional_fields.txt