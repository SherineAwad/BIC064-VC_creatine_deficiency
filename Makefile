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


galore/Sample-90_S10_L001_R1_001.trimmed.fastq.gz: 
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L001_R1_001.fastq.gz ${WHERE}/$(i)_L001_R2_001.fastq.gz; ) 
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L002_R1_001.fastq.gz ${WHERE}/$(i)_L002_R2_001.fastq.gz; ) 
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L003_R1_001.fastq.gz ${WHERE}/$(i)_L003_R2_001.fastq.gz; )
	$(foreach i,$(SAMPLES), trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir ${WHERE}/fastqc" -o ${WHERE}/galore --paired ${WHERE}/$(i)_L004_R1_001.fastq.gz ${WHERE}/$(i)_L004_R2_001.fastq.gz; )




Sample-90_S10_L004.sam: galore/Sample-90_S10_L001_R1_001.trimmed.fastq.gz 
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)i_L001_R1_001_val_1.fq.gz -2 galore/$(i)_L001_R2_001_val_2.fq.gz -S $(i)_L001.bowtie2.sam;)
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)i_L002_R1_001_val_1.fq.gz -2 galore/$(i)_L002_R2_001_val_2.fq.gz -S $(i)_L002.bowtie2.sam;)
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)i_L003_R1_001_val_1.fq.gz -2 galore/$(i)_L003_R2_001_val_2.fq.gz -S $(i)_L003.bowtie2.sam;)
	$(foreach i, $(SAMPLES),  bowtie2 -x ${hg38BOWTIE2INDEX}/genome -1 galore/$(i)i_L004_R1_001_val_1.fq.gz -2 galore/$(i)_L004_R2_001_val_2.fq.gz -S $(i)_L004.bowtie2.sam;)


Sample-90_S10.sam: Sample-90_S10_L004.sam
	 $(foreach i, $(PATIENTS), samtools merge -f $(i).bowtie2.sam $(i)_L*.bowtie2.sam ;) 

Sample-90_S10.RG.bowtie2.sam:  Sample-90_S10.sam
	picard AddOrReplaceReadGroups I=${S1}.bowtie2.sam O=${S1}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG RGSM=Sample-02 RGPL=Illumina RGLB=Sample-02_S1 RGPU=@NDX550116_59_H7N52BGXG.Sample-02
	
	picard AddOrReplaceReadGroups I=${S2}.bowtie2.sam O=${S2}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG RGSM=Sample-22 RGPL=Illumina RGLB=Sample-22_S3 RGPU=@NDX550116_59_H7N52BGXG.Sample-22
	
	picard AddOrReplaceReadGroups I=${S3}.bowtie2.sam O=${S3}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-37  RGPL=Illumina  RGLB=Sample-37_S5 RGPU=@NDX550116_59_H7N52BGXG.Sample-37
	
	picard AddOrReplaceReadGroups I=${S4}.bowtie2.sam O=${S4}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-52  RGPL=Illumina  RGLB=Sample-52_S7 RGPU=@NDX550116_59_H7N52BGXG.Sample-52
	
	picard AddOrReplaceReadGroups I=${S5}.bowtie2.sam O=${S5}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-80  RGPL=Illumina  RGLB=Sample-80_S9 RGPU=@NDX550116_59_H7N52BGXG.Sample-80
	
	picard AddOrReplaceReadGroups I=${S6}.bowtie2.sam O=${S6}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-17  RGPL=Illumina  RGLB=Sample-17_S2 RGPU=@NDX550116_59_H7N52BGXG.Sample-17
	
	picard AddOrReplaceReadGroups I=${S7}.bowtie2.sam O=${S7}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-30  RGPL=Illumina  RGLB=Sample-30_S4 RGPU=@NDX550116_59_H7N52BGXG.Sample-30
	
	picard AddOrReplaceReadGroups I=${S8}.bowtie2.sam O=${S8}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-42  RGPL=Illumina  RGLB=Sample-42_S6 RGPU=@NDX550116_59_H7N52BGXG.Sample-42
	
	picard AddOrReplaceReadGroups I=${S9}.bowtie2.sam O=${S9}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-68   RGPL=Illumina  RGLB=Sample-68_S8 RGPU=@NDX550116_59_H7N52BGXG.Sample-68
	
	picard AddOrReplaceReadGroups I=${S10}.bowtie2.sam O=${S10}.RG.bowtie2.sam SO=coordinate \
	RGID=@NDX550116_59_H7N52BGXG  RGSM=Sample-90   RGPL=Illumina  RGLB=Sample-90_S10 RGPU=@NDX550116_59_H7N52BGXG.Sample-90
	

Sample-90_S10.dedupped.bam: Sample-90_S10.RG.bowtie2.sam
	 $(foreach i, $(PATIENTS),  picard MarkDuplicates I=$(i).RG.bowtie2.sam O=$(i).dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$(i).output.metrics;)  

Sample-90_S10.recal.table: Sample-90_S10.dedupped.bam
	 $(foreach i, $(PATIENTS), gatk BaseRecalibrator  -I $(i).dedupped.bam -R ${hg38FASTA}/genome.fa --known-sites ${GATKRESOURCES}/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites ${GATKRESOURCES}/Homo_sapiens_assembly38.known_indels.vcf --known-sites ${GATKRESOURCES}/Mills_and_1000G_gold_standard.indels.hg38.vcf  -O $(i).recal_data.table;) 


Sample-90_S10.recalibrated.bam: Sample-90_S10.recal.table 
	 $(foreach i, $(PATIENTS), gatk ApplyBQSR -I $(i).dedupped.bam -R ${hg38FASTA}/genome.fa --bqsr-recal-file $(i).recal_data.table -O $(i).recalibrated.bam;)


Sample-90_S10.vcf:  Sample-90_S10.recalibrated.bam
	$(foreach i, $(PATIENTS), gatk --java-options "-Xmx4g" HaplotypeCaller  -R  ${hg38FASTA}/genome.fa -I $(i).recalibrated.bam -O $(i).vcf -bamout $(i).bamout.bam;)

NNM-80_S4.GATK.hg38_multianno.txt: Sample-90_S10.vcf  
	$(foreach i, $(PATIENTS), ${ANNOVAR}/table_annovar.pl $(i).vcf.gz  ${ANNOVAR}/humandb/ -buildver hg38 -out $(i).GATK.hg38 -remove -protocol refGene,ensGene,cytoBand,exac03,gnomad_exome,avsnp147,dbnsfp33a,clinvar_20170130,revel,regsnpintron,dbscsnv11 -operation g,g,f,f,f,f,f,f,f,f,f -nastring . -vcfinput;)



clean: 
	rm *L00*.sam 
	rm *avinput  
