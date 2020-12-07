hg38GTF = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Annotation/Genes
hg38FASTA = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta
ANNOTATIONS = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations

WHERE =/rds/project/yhbl2/rds-yhbl2-genehunter/SM/Naira2
hg38BOWTIE2INDEX = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index
GATKRESOURCES = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/hg38gatk

ANNOVAR = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/tools/annovar

hg19Fasta = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/build37.2/Sequence/WholeGenomeFasta
hg19Bowtie2Index = /rds/project/yhbl2/rds-yhbl2-genehunter/SM/annotations/Homo_sapiens/NCBI/build37.2/Sequence/Bowtie2Index



S1 = Sample-02_S1
S2 = Sample-22_S3 
S3 = Sample-37_S5  
S4 = Sample-52_S7  
S5 = Sample-80_S9
S6 = Sample-17_S2  
S7 = Sample-30_S4  
S8 = Sample-42_S6  
S9 = Sample-68_S8  
S10 = Sample-90_S10


SAMPLES = ${S1} ${S2} ${S3} ${S4} ${S5} ${S6} ${S7} ${S8} ${S9} ${S10}
R_ID = @NDX550116_59_H7N52BGXG

galore/Sample-90_S10_L001_R1_001.trimmed.fastq.gz: 
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L001_R1_001.fastq.gz ${WHERE}/$(i)_L001_R2_001.fastq.gz; ) 
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L002_R1_001.fastq.gz ${WHERE}/$(i)_L002_R2_001.fastq.gz; ) 
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L003_R1_001.fastq.gz ${WHERE}/$(i)_L003_R2_001.fastq.gz; )
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L004_R1_001.fastq.gz ${WHERE}/$(i)_L004_R2_001.fastq.gz; )




Sample-90_S10_L004.sam: galore/Sample-90_S10_L001_R1_001.trimmed.fastq.gz 
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)_L001_R1_001_val_1.fq.gz -2 galore/$(i)_L001_R2_001_val_2.fq.gz -S $(i)_L001.bowtie2.sam;)
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)_L002_R1_001_val_1.fq.gz -2 galore/$(i)_L002_R2_001_val_2.fq.gz -S $(i)_L002.bowtie2.sam;)
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)_L003_R1_001_val_1.fq.gz -2 galore/$(i)_L003_R2_001_val_2.fq.gz -S $(i)_L003.bowtie2.sam;)
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)_L004_R1_001_val_1.fq.gz -2 galore/$(i)_L004_R2_001_val_2.fq.gz -S $(i)_L004.bowtie2.sam;)


Sample-90_S10.sam: Sample-90_S10_L004.sam
	 $(foreach i, $(SAMPLES), samtools merge -f $(i).bowtie2.sam $(i)_L*.bowtie2.sam ;) 

Sample-90_S10.RG.bowtie2.sam:  Sample-90_S10.sam
	$(foreach i, $(SAMPLES), picard AddOrReplaceReadGroups I=$(i).bowtie2.sam O=$(i).RG.bowtie2.sam SO=coordinate RGID=${R_ID} RGSM=$(i) RGPL=Illumina RGLB=$(i) RGPU=${R_ID}.$(i) ;)

Sample-90_S10.dedupped.bam: Sample-90_S10.RG.bowtie2.sam
	 $(foreach i, $(SAMPLES),  picard MarkDuplicates I=$(i).RG.bowtie2.sam O=$(i).dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$(i).output.metrics;)  

Sample-90_S10.recal.table: Sample-90_S10.dedupped.bam
	 $(foreach i, $(SAMPLES), gatk BaseRecalibrator  -I $(i).dedupped.bam -R ${hg38FASTA}/genome.fa --known-sites ${GATKRESOURCES}/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites ${GATKRESOURCES}/Homo_sapiens_assembly38.known_indels.vcf --known-sites ${GATKRESOURCES}/Mills_and_1000G_gold_standard.indels.hg38.vcf  -O $(i).recal_data.table;) 


Sample-90_S10.recalibrated.bam: Sample-90_S10.recal.table 
	 $(foreach i, $(SAMPLES), gatk ApplyBQSR -I $(i).dedupped.bam -R ${hg38FASTA}/genome.fa --bqsr-recal-file $(i).recal_data.table -O $(i).recalibrated.bam;)


Sample-90_S10.vcf:  Sample-90_S10.recalibrated.bam
	$(foreach i, $(SAMPLES), gatk --java-options "-Xmx4g" HaplotypeCaller  -R  ${hg38FASTA}/genome.fa -I $(i).recalibrated.bam -O $(i).vcf -bamout $(i).bamout.bam;)

Sample-90_S10.GATK.hg38_multianno.txt: Sample-90_S10.vcf  
	$(foreach i, $(SAMPLES), ${ANNOVAR}/table_annovar.pl $(i).vcf.gz  ${ANNOVAR}/humandb/ -buildver hg38 -out $(i).GATK.hg38 -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad_exome,avsnp147,dbnsfp33a,clinvar_20170130,revel,regsnpintron,dbscsnv11 -operation g,g,f,f,f,f,f,f,f,f,f -nastring . -vcfinput;)

human.hg38.excl.tsv:
	wget https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv

Sample-90_S10.delly.vcf: Sample-90_S10.vcf human.hg38.excl.tsv 
	$(foreach i, $(SAMPLES), delly call -x human.hg38.excl.tsv  -o $(i).delly.bcf -g ${hg38FASTA}/genome.fa $(i).recalibrated.bam;)
	$(foreach i, $(SAMPLES), bcftools view $(i).delly.bcf > $(i).delly.vcf;)

Sample-90_S10.delly.annotated.tsv2: Sample-90_S10.delly.vcf 
	$(foreach i, $(SAMPLES), $$ANNOTSV/bin/AnnotSV -SVinputFile ./$(i).delly.vcf -outputFile ./$(i).delly.annotated.tsv -genomeBuild GRCh38;)

clean: 
	rm *L00*.sam 
	rm *avinput  
