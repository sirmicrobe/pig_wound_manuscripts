#Code and scripts used to analyze the Porcine full-thickness burn wounds
#Sequence data can be found in BioProjects PRJNA491911 for PA14 and PRJNA633671 for PAO1 

##############################
## Pig pre-processing ##
##############################

#Trim and quality filter all reads from RSCV samples
for i in {1..9}; do trimmomatic PE -phred33 /home/dmux/170123/ErinGloag/012317_0$i/*R1_001.fastq.gz /home/dmux/170123/ErinGloag/012317_0$i/*R2_001.fastq.gz "$i"_forward_paired.fq.gz "$i"_forward_unpaired.fq.gz "$i"_reverse_paired.fq.gz "$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE--pe.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done

for i in {10..96}; do trimmomatic PE -phred33 /home/dmux/170123/ErinGloag/012317_$i/*R1_001.fastq.gz /home/dmux/170123/ErinGloag/012317_$i/*R2_001.fastq.gz "$i"_forward_paired.fq.gz "$i"_forward_unpaired.fq.gz "$i"_reverse_paired.fq.gz "$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE--pe.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done

#Trim reads from non-RSCV samples
for i in {01..84}; do java -jar trimmomatic-0.36.jar PE -threads 4 -phred33 /home/data/dmux/171122/ErinGloag/112217_$i/*R1_001.fastq.gz /home/data/dmux/171122/ErinGloag/112217_$i/*R2_001.fastq.gz "$i"_forward_paired.fq.gz "$i"_forward_unpaired.fq.gz "$i"_reverse_paired.fq.gz "$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/cwm47/build/Trimmomatic-0.36/adapters/NexteraPE--pe.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done

#last batch of non-RSCV samples
module load trimmomatic/trimmomatic-0.36
for i in {77..81}; do trimmomatic PE -threads 4 -phred33 /home/data/dmux/171207/ErinGloag/120717_$i/*R1_001.fastq.gz /home/data/dmux/171207/ErinGloag/120717_$i/*R2_001.fastq.gz /home/cwm47/pig_pa/trim3/"$i"_forward_paired.fq.gz /home/cwm47/pig_pa/trim3/"$i"_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim3/"$i"_reverse_paired.fq.gz /home/cwm47/pig_pa/trim3/"$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/cwm47/build/Trimmomatic-0.36/adapters/NexteraPE--pe.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done


##############################
## Variant calling ##
##############################

#Variant calling with breseq
#running on RSCV isolates with PA14 background
for i in 1 2 3 6 8 9 10 12 14 16 17 20 23 29 36 37 38 40 42 43 45 46 89; do breseq -r ~/ref_genomes/Paeruginosa/PA14_GCF_000014625.1_ASM1462v1_genomic.gbff /home/cwm47/pig_pa/trim/"$i"_forward_paired.fq.gz /home/cwm47/pig_pa/trim/"$i"_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim/"$i"_reverse_paired.fq.gz -o ~/pig_pa/new_ref_breseq/pa14_set/"$i"_breseq_out -j 8;done
#---> bowtie2  :: version 2.2.6 [/usr/bin/bowtie2]
#---> R        :: version 3.2.3

#running on all RSCV isolates with A107/08 pig with PAO1 background
for i in 5 7 11 15 18 21 22 24 25 26 30 31 32 33 34 35 39 41 44 90; do breseq -r ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.gbff /home/cwm47/pig_pa/trim/"$i"_forward_paired.fq.gz /home/cwm47/pig_pa/trim/"$i"_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim/"$i"_reverse_paired.fq.gz -o ~/pig_pa/new_ref_breseq/pao1_set/"$i"_breseq_a107_out -j 8;done
#PAO1 reads vs PAO1 reference (to subtract out false positives)
breseq -r ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.gbff /home/cwm47/pig_pa/trim/90_forward_paired.fq /home/cwm47/pig_pa/trim/90_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim/90_reverse_paired.fq.gz -o ~/pig_pa/new_ref_breseq/pao1_set/90_breseq_a107_out -j 8

#non-RSCVs from PA14
module load breseq/breseq-0.31.0
for i in {32..84}; do breseq -r ~/ref_genomes/Paeruginosa/PA14_GCF_000014625.1_ASM1462v1_genomic.gbff /home/cwm47/pig_pa/trim2/"$i"_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/"$i"_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim2/"$i"_reverse_paired.fq.gz -o /home/cwm47/pig_pa/wt_set/"$i"_breseq_wt_out -j 8;done

module load breseq/breseq-0.31.0
for i in 01 04 10 11 13 14 15 16 28 29 30; do breseq -r ~/ref_genomes/Paeruginosa/PA14_GCF_000014625.1_ASM1462v1_genomic.gbff /home/cwm47/pig_pa/trim2/"$i"_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/"$i"_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim2/"$i"_reverse_paired.fq.gz -o /home/cwm47/pig_pa/wt_set/"$i"_breseq_pa14_a107_out -j 4;done

#non-RSCVs from PAO1
for i in 02 03 05 06 07 08 09 12 13 17 18 19 20 21 22 23 24 25 26 27 31; do breseq -r ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.gbff /home/cwm47/pig_pa/trim2/"$i"_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/"$i"_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim2/"$i"_reverse_paired.fq.gz -o /home/cwm47/pig_pa/wt_set/"$i"_breseq_pa01_a107_out -j 4;done

#last batch fo non-RSCVs from PAO1
for i in {77..81}; do breseq -r ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.gbff /home/cwm47/pig_pa/trim3/"$i"_forward_paired.fq.gz /home/cwm47/pig_pa/trim3/"$i"_forward_unpaired.fq.gz /home/cwm47/pig_pa/trim3/"$i"_reverse_paired.fq.gz -o /home/cwm47/pig_pa/wt_set/wt_pig_a107/wt_152/"$i"_breseq_wt_pa01_out -j 8;done



##############################
## Post-processing ##
##############################

#Did redundant post processing in both R and with gdtools. R script available on github
#subtract out reference sequence
for i in 02 03 05 06 07 08 09 12 13 17 18 19 21 22 23 24 25 26 27; do gdtools SUBTRACT -o "$i"_subtractref_output.gd "$i"_breseq_pa01_a107_out/output/output.gd ~/pig_pa/new_ref_breseq/pao1_set/90_breseq_a107_out/output/output.gd; done
#annotate each gdtools file
for i in 02 03 05 06 07 08 09 12 13 17 18 19 21 22 23 24 25 26 27; do gdtools ANNOTATE -o "$i"_pa01_subtractref.txt -f TSV -r ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.gbff "$i"_subtractref_output.gd; done
#change outputs to sample name - run locally
for i in 02 03 05 06 07 08 09 12 13 17 18 19 21 22 23 24 25 26 27; do sed -i '' "s/output/sample_$i/g" "$i"_pa01_subtractref.txt; done

#non-RSCV mutations
#subtract out reference sequence
for i in {32..84}; do gdtools SUBTRACT -o "$i"_subtractref_output.gd "$i"_breseq_wt_out/output/output.gd ~/pig_pa/new_ref_breseq/pa14_set/89_breseq_out/output/output.gd; done
#annotate each gdtools file
for i in {32..84}; do gdtools ANNOTATE -o "$i"_wt_pa14_subtractref.txt -f TSV -r ~/ref_genomes/Paeruginosa/PA14_GCF_000014625.1_ASM1462v1_genomic.gbff "$i"_subtractref_output.gd; done
# cut off previous analysis and ran the last few high mapping samples
for i in 70 82 83; do gdtools ANNOTATE -o "$i"_wt_pa14_subtractref.txt -f TSV -r ~/ref_genomes/Paeruginosa/PA14_GCF_000014625.1_ASM1462v1_genomic.gbff "$i"_subtractref_output.gd; done
#change outputs to sample name - run locally
for i in 32 33 34 35 36 37 38 39 40 42 43 44 47 48 49 50 51 62 63 64 67 68 70 82 83; do sed -i '' "s/output/sample_$i/g" "$i"_wt_pa14_subtractref.txt; done
#subtract ref
for i in {77..81}; do gdtools SUBTRACT -o "$i"_152_subtractref_pa01_output.gd "$i"_breseq_wt_pa01_out/output/output.gd ~/pig_pa/new_ref_breseq/pao1_set/90_breseq_a107_out/output/output.gd; done
#annotate each gdtools file
for i in {77..81}; do gdtools ANNOTATE -o "$i"_152_wt_pa01_subtractref.txt -f TSV -r ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.gbff "$i"_152_subtractref_pa01_output.gd; done
#change outputs to sample name - run locally
for i in {152..156}; do sed -i '' "s/output/sample_$i/g" "$i"_wt_pa01_subtractref.txt; done

##############################
## Nanopore MinION Sequencing ##
##############################

#Basecalling
#ONT Albacore Sequencing Pipeline Software (version 2.2.5):
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -r -i /home/data/Nanopore/180801/20180801_1916_ChrisPigPA/fast5/pass -t 24 -s /home/cwm47/pig_pa/nanopore/08022018 --barcoding -o fastq -q 0

#ONT Albacore Sequencing Pipeline Software (version 2.2.5):
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -r -i 20180825_1823_ChrisPA/fast5/pass/ -t 24 -s /home/cwm47/pig_pa/nanopore/08252018 --barcoding -o fastq -q 0 
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -r -i /home/data/Nanopore/180824/20180824_1649_ChrisPA/fast5/pass -t 24 -s /home/cwm47/pig_pa/nanopore/08242018 --barcoding -o fastq -q 0

#nanopore B23 and 130: ONT Albacore Sequencing Pipeline Software (version 2.2.5d) #-f FLO-MIN106 -k SQK-LSK109
multi_to_single_fast5 --input_path /home/data/Nanopore/190605/20190604_2205_MN18102_FAK93217_df0744a2/ --save_path /home/data/Nanopore/190605/20190604_2205_MN18102_FAK93217_df0744a2//single_reads --recursive -t 16

read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -r -i /home/data/Nanopore/190605/20190604_2205_MN18102_FAK93217_df0744a2/single_reads/ -t 16 -s /home/dmux/nanopore/190604 --barcoding -o fastq -q 0 # 2 days, 14:05:30


#Concatenating reads from multiple runs
# 2 PA14 
cat ~/pig_pa/nanopore/08022018/workspace/pass/barcode01/*.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode01/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode01/*.fastq > pa14_2_final.fastq
# 5 PAO1
cat pao1_5_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode02/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode02/*.fastq > pao1_5_final.fastq
# 14 PA14
cat pa14_14_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode03/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode03/*.fastq > pa14_14_final.fastq
#15 PAO1
cat pao1_15_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode04/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode04/*.fastq > pao1_15_final.fastq
#16 PA14
cat pa14_16_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode05/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode05/*.fastq > pa14_16_final.fastq
#31 PAO1
cat pao1_31_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode06/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode06/*.fastq > pao1_31_final.fastq
#39 PAO1
cat pao1_39_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode07/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode07/*.fastq > pao1_39_final.fastq
#46 PAO1
cat pao1_46_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode08/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode08/*.fastq > pao1_46_final.fastq
#100 PA14
cat pa14_100_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode09/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode09/*.fastq > pa14_100_final.fastq
#126 PAO1
cat pao1_126_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode10/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode10/*.fastq > pao1_126_final.fastq
#146 PAO1
cat pao1_146_catreads.fastq /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode11/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode11/*.fastq > pao1_146_final.fastq
#S54485
cat ~/pig_pa/nanopore/08022018/workspace/pass/barcode12/*.fastq  /home/cwm47/pig_pa/nanopore/08242018/workspace/pass/barcode12/*.fastq /home/cwm47/pig_pa/nanopore/08252018/workspace/pass/barcode12/*.fastq > PAS54485_12_final.fastq

##############################
## Hybrid Genome Assemblies ##
##############################

#Used Unicycler v0.4.4 to assemble genomes
  Program         Version   Status  
  spades.py       3.13.0    good    
  racon           -         good    
  makeblastdb     2.7.1+    good    
  tblastn         2.7.1+    good    
  bowtie2-build   2.2.6     good    
  bowtie2         2.2.6     good    
  samtools        1.4       good    
  java            11.0.1    good    
  pilon           1.22      good 
#Hybrid assembly sample 2:
unicycler -1 /home/cwm47/pig_pa/trim/2_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim/2_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim/2_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pa14_2_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_2_PA14_assembly -t 16 --mode conservative --keep 3
#assembly stats script
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_2_PA14_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_2_PA14_assembly/stat_02_assembly.txt
#Hybrid assembly sample5:
unicycler -1 /home/cwm47/pig_pa/trim2/02_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/02_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/02_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pao1_5_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_05_PAO1_assembly -t 16 --keep 2 --mode conservative
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_05_PAO1_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_05_PAO1_assembly/stat_05_assembly.txt
#Hybrid assembly sample 14:
unicycler -1 /home/cwm47/pig_pa/trim/14_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim/14_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim/14_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pa14_14_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_14_PA14_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_14_PA14_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_14_PA14_assembly/stat_14_assembly.txt
#samples 15 
unicycler -1 /home/cwm47/pig_pa/trim2/05_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/05_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/05_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pao1_15_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_15_PAO1_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_15_PAO1_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_15_PAO1_assembly/stat_15_assembly.txt
#sample 16
unicycler -1 /home/cwm47/pig_pa/trim/16_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim/16_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim/16_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pa14_16_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_16_PA14_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_16_PA14_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_16_PA14_assembly/stat_16_assembly.txt
#samples 31 
unicycler -1 /home/cwm47/pig_pa/trim2/18_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/18_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/18_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pao1_31_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_31_PAO1_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_31_PAO1_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_31_PAO1_assembly/stat_31_assembly.txt
# 39
unicycler -1 /home/cwm47/pig_pa/trim2/24_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/24_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/24_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pao1_39_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_39_PAO1_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_39_PAO1_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_39_PAO1_assembly/stat_39_assembly.txt
#46 
unicycler -1 /home/cwm47/pig_pa/trim2/27_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/27_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/27_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pao1_46_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_46_PAO1_assembly -t 16 --mode conservative --keep 2 
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_46_PAO1_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_46_PAO1_assembly/stat_46_assembly.txt
#100
unicycler -1 /home/cwm47/pig_pa/trim2/35_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/35_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/35_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pa14_100_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_100_PA14_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_100_PA14_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_100_PA14_assembly/stat_100_assembly.txt
#126
unicycler -1 /home/cwm47/pig_pa/trim2/61_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/61_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/61_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pao1_126_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_126_PAO1_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_126_PAO1_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_126_PAO1_assembly/stat_126_assembly.txt
#146
unicycler -1 /home/cwm47/pig_pa/trim2/81_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/81_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/81_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/all_reads_cat/pao1_146_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_146_PAO1_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_146_PAO1_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_146_PAO1_assembly/stat_146_assembly.txt
# S54485
unicycler -1 /home/cwm47/pig_pa/trim/94_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim/94_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim/94_forward_unpaired.fq.gz  -l /home/cwm47/pig_pa/nanopore/all_reads_cat/PAS54485_12_final.fastq -o /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_S54485_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_S54485_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_S54485_assembly/stat_S54485_assembly.txt
# b23 assembly
unicycler -1 /home/cwm47/pig_pa/trim/91_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim/91_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim/91_forward_unpaired.fq.gz -l /home/cwm47/pig_pa/nanopore/06042019/b23_reads/b23_nanopore_reads.fastq -o /home/cwm47/pig_pa/nanopore/06042019/b23_unicycler_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/06042019/b23_unicycler_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/06042019/b23_unicycler_assembly/stat_b23_assembly.txt
# 130 assembly
unicycler -1 /home/cwm47/pig_pa/trim2/65_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/65_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/65_forward_unpaired.fq.gz -l /home/dmux/nanopore/190604-albacore-calls/workspace/pass/barcode12/*.fastq -o /home/cwm47/pig_pa/nanopore/06042019/130_unicycler_assembly -t 16 --mode conservative --keep 2
#stats
ruby /home/cwm47/build/gscripts/fasta_contig_len_distribution.rb /home/cwm47/pig_pa/nanopore/06042019/130_unicycler_assembly/assembly.fasta > /home/cwm47/pig_pa/nanopore/06042019/130_unicycler_assembly/stat_130_assembly.txt


#Prokka v1.13 annotation of genomes
for i in 2 15 16 31 39 46 100 126 146 S54485;do prokka --centre X  --force --locustag npore_"$i" --outdir /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_"$i"_*_assembly/"$i"_prokka --prefix "$i"_npore --gffver 3 --cpus 8 /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/sample_"$i"_*_assembly/assembly.fasta; done
#b23
prokka --centre X  --force --locustag b23_npore --outdir /home/cwm47/pig_pa/nanopore/06042019/b23_prokka --prefix b23_npore --gffver 3 --cpus 8 /home/cwm47/pig_pa/nanopore/06042019/b23_unicycler_assembly/assembly.fasta
#130
prokka --centre X  --force --locustag 130_npore --outdir /home/cwm47/pig_pa/nanopore/06042019/130_prokka --prefix 130_npore --gffver 3 --cpus 8 /home/cwm47/pig_pa/nanopore/06042019/130_unicycler_assembly/assembly.fasta


##############################
## PHASTER ##
##############################

#Predict prophage sequences in assemblies
for i in 2 05 15 16 31 39 46 100 126 146 S54485; do wget --post-file="/Users/chrismarshall/Documents/Pitt/Cooper_Lab/Misc_projects/Pig_wound/nanopore/final_nanopore_multirun/assemblies/"$i"_assembly.fasta" "http://phaster.ca/phaster_api?contigs=1" -O "$i"_phaster;done

for i in 2 05 15 16 31 39 46 100 126 146 S54485; do url=$(cut -d "\"" -f4 "$i"_phaster) && wget phaster.ca/phaster_api?acc=$url && mv phaster_api\?acc\=$url "$i"_phaster.txt;done

#in R
library(jsonlite)
for(i in c("2","05", "14", "15", "16", "31", "39", "46", "100", "126", "146", "S54485")) {
  p_i <- fromJSON(paste(i,"phaster.txt",sep="_"))$summary
  cat(p_i)
}


##############################
## Pangenomics ##
##############################

# assemblies for Anvio and MGEfinder


module load miniconda/miniconda-3
module load blast/blast-2.9.0
module load prokka/prokka-1.14.0
module load tbl2asn/tbl2asn-27.0
#redo spades on all PAO1 Illumina samples to keep consistent numbers
#5
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/5_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/5_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/5_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/02_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/02_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/02_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/5_spades
#11
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/11_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/11_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/11_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/03_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/03_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/03_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/11_spades
#15
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/15_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/15_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/15_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/05_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/05_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/05_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/15_spades
#18
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/18_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/18_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/18_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/06_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/06_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/06_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/18_spades
#19
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/19_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/19_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/19_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/07_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/07_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/07_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/19_spades
#21
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/21_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/21_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/21_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/08_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/08_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/08_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/21_spades
#22
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/22_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/22_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/22_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/09_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/09_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/09_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/22_spades
#25
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/25_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/25_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/25_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/12_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/12_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/12_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/25_spades
#26
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/26_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/26_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/26_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/13_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/13_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/13_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/26_spades
#30
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/30_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/30_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/30_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/17_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/17_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/17_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/30_spades
#31
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/31_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/31_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/31_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/18_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/18_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/18_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/31_spades
#32
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/32_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/32_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/32_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/19_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/19_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/19_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/32_spades
#33
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/21_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/21_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/21_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/33_spades
#34
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/34_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/34_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/34_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/22_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/22_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/22_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/34_spades
#35
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/35_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/35_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/35_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/23_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/23_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/23_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/35_spades
#39
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/39_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/39_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/39_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/24_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/24_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/24_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/39_spades
#41
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/41_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/41_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/41_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/25_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/25_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/25_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/41_spades
#44
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/44_forward_paired.fq.gz --pe1-2 /home/cwm47/pig_pa/trim/44_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/44_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/26_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/26_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/26_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/44_spades
#46
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 --pe1-1 /home/cwm47/pig_pa/trim/46_forward_paired.fq --pe1-2 /home/cwm47/pig_pa/trim/46_reverse_paired.fq.gz --pe1-s /home/cwm47/pig_pa/trim/46_forward_unpaired.fq.gz --pe2-1 /home/cwm47/pig_pa/trim2/27_forward_paired.fq.gz --pe2-2 /home/cwm47/pig_pa/trim2/27_reverse_paired.fq.gz --pe2-s /home/cwm47/pig_pa/trim2/27_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/46_spades
#106
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/41_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/41_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/41_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/106_spades
#110
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/45_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/45_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/45_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/110_spades
#111
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/46_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/46_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/46_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/111_spades
#117
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/52_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/52_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/52_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/117_spades
#118
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/53_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/53_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/53_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/118_spades
#119
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/54_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/54_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/54_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/119_spades
#120
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/55_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/55_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/55_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/120_spades
#121
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/56_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/56_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/56_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/121_spades
#122
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/57_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/57_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/57_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/122_spades
#123
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/58_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/58_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/58_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/123_spades
#124
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/59_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/59_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/59_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/124_spades
#125
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/60_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/60_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/60_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/125_spades
#126
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/61_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/61_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/61_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/126_spades
#130
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/65_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/65_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/65_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/130_spades
#131
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/66_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/66_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/66_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/131_spades
#134
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/69_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/69_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/69_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/134_spades
#136
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/71_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/71_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/71_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/136_spades
#137
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/72_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/72_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/72_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/137_spades
#138
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/73_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/73_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/73_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/138_spades
#139
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/74_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/74_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/74_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/139_spades
#140
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/75_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/75_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/75_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/140_spades
#141
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/76_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/76_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/76_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/141_spades
#142
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/77_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/77_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/77_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/142_spades
#143
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/78_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/78_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/78_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/143_spades
#144
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/79_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/79_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/79_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/144_spades
#145
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/80_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/80_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/80_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/145_spades
#146
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/81_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/81_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/81_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/146_spades
#149
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim2/84_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim2/84_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim2/84_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/149_spades
#152
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim3/77_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim3/77_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim3/77_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/152_spades
#153
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim3/78_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim3/78_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim3/78_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/153_spades
#154
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim3/79_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim3/79_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim3/79_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/154_spades
#155
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim3/80_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim3/80_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim3/80_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/155_spades
#156
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim3/81_forward_paired.fq.gz -2 /home/cwm47/pig_pa/trim3/81_reverse_paired.fq.gz -s /home/cwm47/pig_pa/trim3/81_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/156_spades
#89 - P14
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim/89_forward_paired.fq.gz -2 89_reverse_paired.fq.gz -s 89_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/PA14_spades
#90 - PAO1
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim/90_forward_paired.fq.gz -2 90_reverse_paired.fq.gz -s 90_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/PAO1_spades
#91 - B23-2
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim/91_forward_paired.fq.gz -2 91_reverse_paired.fq.gz -s 91_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/B23_spades
#92 - CF18-1
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim/92_forward_paired.fq.gz -2 92_reverse_paired.fq.gz -s 92_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/CF18_spades
#93 - MSH10-2
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim/93_forward_paired.fq.gz -2 93_reverse_paired.fq.gz -s 93_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/MSH10_spades
#94 - S54485-1
/home/cwm47/build/SPAdes-3.14.0-Linux/bin/spades.py --isolate -t 24 -k 33,55,77,91,111 -1 /home/cwm47/pig_pa/trim/94_forward_paired.fq.gz -2 94_reverse_paired.fq.gz -s 94_forward_unpaired.fq.gz -o /home/cwm47/pig_pa/PAO1_spades_update_nums/S54485_spades

#prokka
for i in 5 11 15 18 19 21 22 25 26 30 31 32 33 34 35 39 41 44 46 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156;do prokka --centre X --force --locustag "$i" --outdir /home/cwm47/pig_pa/PAO1_spades_update_nums/"$i"_prokka --prefix "$i" --gffver 3 --cpus 8 /home/cwm47/pig_pa/PAO1_spades_update_nums/"$i"_spades/contigs.fasta;done

#blast to custom database of all 6 ancestral strains
#concatenated all NCBI refseq genomes + B23 assembly into 6PAstrains.ffn 
makeblastdb -in /home/cwm47/database/blastdb/6PAstrains.ffn -dbtype nucl -parse_seqids
#assign locus tags to all genes based on closest blast hit (helps to find regions where insertions are from another strain)
for i in  5 11 15 18 19 21 22 25 26 30 31 32 33 34 35 39 41 44 46 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156; do blastn -task blastn -db 6PAstrains.ffn -query /home/cwm47/pig_pa/PAO1_spades_update_nums/"$i"_prokka/"$i".ffn -num_threads 8 -outfmt "6 qseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sseqid" -out /home/cwm47/pig_pa/PAO1_spades_update_nums/"$i"_blastresult.ffn -num_alignments 1;done


##############################
## MGE Finder ##
##############################

source activate mgefinder
#Current version of snakemake: 3.13.3
#Expected version of snakemake: 3.13.3
#Current version of einverted: EMBOSS:6.6.0.0
#Expected version of einverted: EMBOSS:6.6.0.0
#Current version of bowtie2: 2.3.5
#Expected version of bowtie2: 2.3.5
#Current version of samtools: 1.9
#Expected version of samtools: 1.9
#Current version of cd-hit: 4.8.1
#Expected version of cd-hit: 4.8.1
#index reference genome - only need to do this the first time
bwa index /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna
module load miniconda/miniconda-3
source activate mgefinder
cd /home/cwm47/pig_pa/MGEfinder/00.bam
#bwa version 0.7.17-r1188
#5
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/5_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/5_reverse_paired.fq.gz > 5_map.sam
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/02_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/02_reverse_paired.fq.gz > 5_map2.sam 
samtools sort -o 5_out1.bam -@ 16 5_map.sam
samtools sort -o 5_out2.bam -@ 16 5_map2.sam
samtools merge 5.pao1.bam 5_out1.bam 5_out2.bam
samtools index -b -@ 16 5.pao1.bam
#11
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/11_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/11_reverse_paired.fq.gz > 11_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/03_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/03_reverse_paired.fq.gz > 11_map2.sam 
samtools sort -o 11_out1.bam -@ 16 11_map.sam
samtools sort -o 11_out2.bam -@ 16 11_map2.sam
samtools merge 11.pao1.bam 11_out1.bam 11_out2.bam
samtools index -b -@ 16 11.pao1.bam
#15
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/15_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/15_reverse_paired.fq.gz > 15_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/05_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/05_reverse_paired.fq.gz > 15_map2.sam 
samtools sort -o 15_out1.bam -@ 16 15_map.sam
samtools sort -o 15_out2.bam -@ 16 15_map2.sam
samtools merge 15.pao1.bam 15_out1.bam 15_out2.bam
samtools index -b -@ 16 15.pao1.bam
#18
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/18_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/18_reverse_paired.fq.gz > 18_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/06_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/06_reverse_paired.fq.gz > 18_map2.sam 
samtools sort -o 18_out1.bam -@ 16 18_map.sam
samtools sort -o 18_out2.bam -@ 16 18_map2.sam
samtools merge 18.pao1.bam 18_out1.bam 18_out2.bam
samtools index -b -@ 16 18.pao1.bam
#19
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/19_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/19_reverse_paired.fq.gz > 19_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/07_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/07_reverse_paired.fq.gz > 19_map2.sam 
samtools sort -o 19_out1.bam -@ 16 19_map.sam
samtools sort -o 19_out2.bam -@ 16 19_map2.sam
samtools merge 19.pao1.bam 19_out1.bam 19_out2.bam
samtools index -b -@ 16 19.pao1.bam
#21
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/21_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/21_reverse_paired.fq.gz > 21_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/08_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/08_reverse_paired.fq.gz > 21_map2.sam 
samtools sort -o 21_out1.bam -@ 16 21_map.sam
samtools sort -o 21_out2.bam -@ 16 21_map2.sam
samtools merge 21.pao1.bam 21_out1.bam 21_out2.bam
samtools index -b -@ 16 21.pao1.bam
#22
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/22_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/22_reverse_paired.fq.gz > 22_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/09_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/09_reverse_paired.fq.gz > 22_map2.sam 
samtools sort -o 22_out1.bam -@ 16 22_map.sam
samtools sort -o 22_out2.bam -@ 16 22_map2.sam
samtools merge 22.pao1.bam 22_out1.bam 22_out2.bam
samtools index -b -@ 16 22.pao1.bam
#25
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/25_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/25_reverse_paired.fq.gz > 25_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/12_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/12_reverse_paired.fq.gz > 25_map2.sam
samtools sort -o 25_out1.bam -@ 16 25_map.sam
samtools sort -o 25_out2.bam -@ 16 25_map2.sam
samtools merge 25.pao1.bam 25_out1.bam 25_out2.bam
samtools index -b -@ 16 25.pao1.bam
#26
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/26_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/26_reverse_paired.fq.gz > 26_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/13_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/13_reverse_paired.fq.gz > 26_map2.sam
samtools sort -o 26_out1.bam -@ 16 26_map.sam
samtools sort -o 26_out2.bam -@ 16 26_map2.sam
samtools merge 26.pao1.bam 26_out1.bam 26_out2.bam
samtools index -b -@ 16 26.pao1.bam
#30
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/30_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/30_reverse_paired.fq.gz > 30_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/17_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/17_reverse_paired.fq.gz > 30_map2.sam 
samtools sort -o 30_out1.bam -@ 16 30_map.sam
samtools sort -o 30_out2.bam -@ 16 30_map2.sam
samtools merge 30.pao1.bam 30_out1.bam 30_out2.bam
samtools index -b -@ 16 30.pao1.bam
#31
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/31_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/31_reverse_paired.fq.gz > 31_map.sam
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/18_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/18_reverse_paired.fq.gz > 31_map2.sam
samtools sort -o 31_out1.bam -@ 16 31_map.sam
samtools sort -o 31_out2.bam -@ 16 31_map2.sam
samtools merge 31.pao1.bam 31_out1.bam 31_out2.bam
samtools index -b -@ 16 31.pao1.bam
#32
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/32_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/32_reverse_paired.fq.gz > 32_map.sam
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/19_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/19_reverse_paired.fq.gz > 32_map2.sam
samtools sort -o 32_out1.bam -@ 16 32_map.sam
samtools sort -o 32_out2.bam -@ 16 32_map2.sam
samtools merge 32.pao1.bam 32_out1.bam 32_out2.bam
samtools index -b -@ 16 32.pao1.bam
#33
bwa mem -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/21_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/21_reverse_paired.fq.gz > 33_map.sam
mgefinder formatbam 33_map.sam 33_out.bam
#34
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/34_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/34_reverse_paired.fq.gz > 34_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/22_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/22_reverse_paired.fq.gz > 34_map2.sam 
samtools sort -o 34_out1.bam -@ 16 34_map.sam
samtools sort -o 34_out2.bam -@ 16 34_map2.sam
samtools merge 34.pao1.bam 34_out1.bam 34_out2.bam
samtools index -b -@ 16 34.pao1.bam
#35
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/35_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/35_reverse_paired.fq.gz > 35_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/23_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/23_reverse_paired.fq.gz > 35_map2.sam
samtools sort -o 35_out1.bam -@ 16 35_map.sam
samtools sort -o 35_out2.bam -@ 16 35_map2.sam
samtools merge 35.pao1.bam 35_out1.bam 35_out2.bam
samtools index -b -@ 16 35.pao1.bam
#39
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/39_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/39_reverse_paired.fq.gz > 39_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/24_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/24_reverse_paired.fq.gz > 39_map2.sam
samtools sort -o 39_out1.bam -@ 16 39_map.sam
samtools sort -o 39_out2.bam -@ 16 39_map2.sam
samtools merge 39.pao1.bam 39_out1.bam 39_out2.bam
samtools index -b -@ 16 39.pao1.bam
#41
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/41_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/41_reverse_paired.fq.gz > 41_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/25_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/25_reverse_paired.fq.gz > 41_map2.sam 
samtools sort -o 41_out1.bam -@ 16 41_map.sam
samtools sort -o 41_out2.bam -@ 16 41_map2.sam
samtools merge 41.pao1.bam 41_out1.bam 41_out2.bam
samtools index -b -@ 16 41.pao1.bam
#44
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/44_forward_paired.fq.gz  /home/cwm47/pig_pa/trim/44_reverse_paired.fq.gz > 44_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/26_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/26_reverse_paired.fq.gz > 44_map2.sam
samtools sort -o 44_out1.bam -@ 16 44_map.sam
samtools sort -o 44_out2.bam -@ 16 44_map2.sam
samtools merge 44.pao1.bam 44_out1.bam 44_out2.bam
samtools index -b -@ 16 44.pao1.bam
#46
bwa mem  -t 24  /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim/46_forward_paired.fq  /home/cwm47/pig_pa/trim/46_reverse_paired.fq.gz > 46_map.sam 
bwa mem -t 24 ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/27_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/27_reverse_paired.fq.gz > 46_map2.sam 
samtools sort -o 46_out1.bam -@ 16 46_map.sam
samtools sort -o 46_out2.bam -@ 16 46_map2.sam
samtools merge 46.pao1.bam 46_out1.bam 46_out2.bam
samtools index -b -@ 16 46.pao1.bam
#106
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/41_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/41_reverse_paired.fq.gz > 106_map.sam
mgefinder formatbam 106_map.sam 106_out.bam
#110
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/45_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/45_reverse_paired.fq.gz > 110_map.sam
mgefinder formatbam 110_map.sam 110_out.bam
#111
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/46_forward_paired.fq.gz  /home/cwm47/pig_pa/trim2/46_reverse_paired.fq.gz > 111_map.sam
mgefinder formatbam 111_map.sam 111_out.bam
#117
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/52_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/52_reverse_paired.fq.gz > 117_map.sam
mgefinder formatbam 117_map.sam 117_out.bam
#118
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/53_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/53_reverse_paired.fq.gz > 118_map.sam
mgefinder formatbam 118_map.sam 118_out.bam
#119
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/54_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/54_reverse_paired.fq.gz > 119_map.sam
mgefinder formatbam 119_map.sam 119_out.bam
#120
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/55_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/55_reverse_paired.fq.gz > 120_map.sam
mgefinder formatbam 120_map.sam 120_out.bam
#121
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/56_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/56_reverse_paired.fq.gz > 121_map.sam
mgefinder formatbam 121_map.sam 121_out.bam
#122
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/57_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/57_reverse_paired.fq.gz > 122_map.sam
mgefinder formatbam 122_map.sam 122_out.bam
#123
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/58_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/58_reverse_paired.fq.gz > 123_map.sam
mgefinder formatbam 123_map.sam 123_out.bam
#124
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/59_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/59_reverse_paired.fq.gz > 124_map.sam
mgefinder formatbam 124_map.sam 124_out.bam
#125
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/60_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/60_reverse_paired.fq.gz > 125_map.sam
mgefinder formatbam 125_map.sam 125_out.bam
#126
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/61_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/61_reverse_paired.fq.gz > 126_map.sam
mgefinder formatbam 5_map.sam 126_out.bam
#130
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/65_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/65_reverse_paired.fq.gz > 130_map.sam
mgefinder formatbam 130_map.sam 130_out.bam
#131
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/66_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/66_reverse_paired.fq.gz > 131_map.sam
mgefinder formatbam 131_map.sam 131_out.bam
#134
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/69_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/69_reverse_paired.fq.gz > 134_map.sam
mgefinder formatbam 134_map.sam 134_out.bam
#136
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/71_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/71_reverse_paired.fq.gz > 136_map.sam
mgefinder formatbam 136_map.sam  136_out.bam
#137
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/72_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/72_reverse_paired.fq.gz > 137_map.sam
mgefinder formatbam 137_map.sam 137_out.bam
#138
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/73_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/73_reverse_paired.fq.gz > 138_map.sam
mgefinder formatbam 138_map.sam 138_out.bam
#139
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/74_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/74_reverse_paired.fq.gz > 139_map.sam
mgefinder formatbam 139_map.sam 139_out.bam
#140
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/75_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/75_reverse_paired.fq.gz > 140_map.sam
mgefinder formatbam 140_map.sam 140_out.bam
#141
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/76_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/76_reverse_paired.fq.gz > 141_map.sam
mgefinder formatbam 141_map.sam 141_out.bam
#142
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/77_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/77_reverse_paired.fq.gz > 142_map.sam
mgefinder formatbam 142_map.sam 142_out.bam
#143
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/78_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/78_reverse_paired.fq.gz > 143_map.sam
mgefinder formatbam 143_map.sam 143_out.bam
#144
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/79_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/79_reverse_paired.fq.gz > 144_map.sam
mgefinder formatbam 144_map.sam  144_out.bam
#145
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/80_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/80_reverse_paired.fq.gz > 145_map.sam
mgefinder formatbam 145_map.sam 145_out.bam
#146
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/81_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/81_reverse_paired.fq.gz > 146_map.sam
mgefinder formatbam 146_map.sam 146_out.bam
#149
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim2/84_forward_paired.fq.gz /home/cwm47/pig_pa/trim2/84_reverse_paired.fq.gz > 149_map.sam
mgefinder formatbam 149_map.sam 149_out.bam
#152
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim3/77_forward_paired.fq.gz /home/cwm47/pig_pa/trim3/77_reverse_paired.fq.gz > 152_map.sam
mgefinder formatbam 152_map.sam 152_out.bam
#153
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim3/78_forward_paired.fq.gz /home/cwm47/pig_pa/trim3/78_reverse_paired.fq.gz > 153_map.sam
mgefinder formatbam 153_map.sam 153_out.bam
#154
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim3/79_forward_paired.fq.gz  /home/cwm47/pig_pa/trim3/79_reverse_paired.fq.gz > 154_map.sam
mgefinder formatbam 154_map.sam 154_out.bam
#155
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim3/80_forward_paired.fq.gz /home/cwm47/pig_pa/trim3/80_reverse_paired.fq.gz > 155_map.sam
mgefinder formatbam 155_map.sam 155_out.bam
#156
bwa mem  -t 24 /home/cwm47/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna /home/cwm47/pig_pa/trim3/81_forward_paired.fq.gz /home/cwm47/pig_pa/trim3/81_reverse_paired.fq.gz > 156_map.sam
mgefinder formatbam 156_map.sam 156_out.bam

#make sure file structure for MGE is correct
for i in 5 11 15 18 19 21 22 25 26 30 31 32 33 34 35 39 41 44 46 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156;do cp /home/cwm47/pig_pa/PAO1_spades_update_nums/"$i"_spades/contigs.fasta /home/cwm47/pig_pa/MGEfinder/00.assembly/"$i".fna;done

for i in 33 46 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156;do mv "$i"_out.bam "$i".pao1.bam;done

for i in 33 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156;do mv "$i"_out.bam.bai "$i".pao1.bam.bai;done

srun -J mge mgefinder workflow --cores 24 MGEfinder/

##############################
## Anvi'o ##
##############################

source activate anvio-6.2

#fix fasta headers
anvi-script-reformat-fasta /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/assemblies/S54485_assembly.fasta -o S54485_contigs-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta /home/cwm47/pig_pa/nanopore/unicycler_assemblies_9_2018/b23_unicycler_assembly/assembly.fasta -o B23_contigs-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta ~/ref_genomes/Paeruginosa/PA14_GCF_000014625.1_ASM1462v1_genomic.fna -o PA14_contigs-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta ~/ref_genomes/Paeruginosa/PAO1_GCF_000006765.1_ASM676v1_genomic.fna -o PAO1_contigs-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta ~/ref_genomes/Paeruginosa/GCF_000407905.1_Pseu_aeru_MSH-10_V1_genomic.fna -o MSH10_contigs-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta ~/ref_genomes/Paeruginosa/GCF_000481925.1_Pseu_aeru_CF18_V1_genomic.fna -o CF18_contigs-fixed.fa -l 0 --simplify-names

#create contig databases of references
#Prodigal (v2.6.3)
anvi-gen-contigs-database -f S54485_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/S54485_contigs.db
anvi-gen-contigs-database -f B23_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/B23_contigs.db
anvi-gen-contigs-database -f PA14_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/PA14_contigs.db
anvi-gen-contigs-database -f PAO1_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/PAO1_contigs.db
anvi-gen-contigs-database -f MSH10_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/MSH10_contigs.db
anvi-gen-contigs-database -f CF18_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/CF18_contigs.db

#create contig databases of references
#Prodigal (v2.6.3)
anvi-gen-contigs-database -f S54485_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/S54485_contigs.db
anvi-gen-contigs-database -f B23_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/B23_contigs.db
anvi-gen-contigs-database -f PA14_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/PA14_contigs.db
anvi-gen-contigs-database -f PAO1_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/PAO1_contigs.db
anvi-gen-contigs-database -f MSH10_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/MSH10_contigs.db
anvi-gen-contigs-database -f CF18_contigs-fixed.fa -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/CF18_contigs.db

#Create contigs databses from fasta files
for i in 5 11 15 18 19 21 22 25 26 30 31 32 33 34 35 39 41 44 46 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156 PAO1 PA14 B23 S54485 MSH10 CF18; do anvi-gen-contigs-database -f /home/cwm47/pig_pa/MGEfinder/00.assembly/"$i".fna -n 'PAO1 Isolates from Pig Wounds' -o /home/cwm47/pig_pa/anvio_final/contigs_db/"$i"_contigs.db;done

for i in 5 11 15 18 19 21 22 25 26 30 31 32 33 34 35 39 41 44 46 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156 PAO1 PA14 B23 S54485 MSH10 CF18; do anvi-run-hmms --num-threads 24 --also-scan-trnas -c /home/cwm47/pig_pa/anvio_final/contigs_db/"$i"_contigs.db; done

for i in 5 11 15 18 19 21 22 25 26 30 31 32 33 34 35 39 41 44 46 106 110 111 117 118 119 120 121 122 123 124 125 126 130 131 134 136 137 138 139 140 141 142 143 144 145 146 149 152 153 154 155 156 PAO1 PA14 B23 S54485 MSH10 CF18; do anvi-run-ncbi-cogs -c /home/cwm47/pig_pa/anvio_final/contigs_db/"$i"_contigs.db --num-threads 24 --cog-data-dir /opt/miniconda/miniconda3/envs/anvio3/lib/python3.5/site-packages/anvio/data/misc/COG/;done #this searches with DIAMOND, could set to blast with --search-with blastp

#store contigs in database
anvi-gen-genomes-storage -e /home/cwm47/pig_pa/anvio_final/genomes_storage_nocontam.txt -o /home/cwm47/pig_pa/anvio_final/all_final_nocontam-GENOMES.db 

#pangenome analysis with diamond annotation (cite Diamond)
srun -J anvi anvi-pan-genome -g /home/cwm47/pig_pa/anvio_final/all_final_nocontam-GENOMES.db -n anvi_final_nocontam_pangenome --num-threads 24 --exclude-partial-gene-calls 
#####
#display
sudo anvi-display-pan -p /Users/chrismarshall/Documents/Pitt/Cooper_Lab/Misc_projects/Pig_wound/anvio_pig/anvio_final/anvi_final_nocontam_pangenome/anvi_final_nocontam_pangenome-PAN.db -g /Users/chrismarshall/Documents/Pitt/Cooper_Lab/Misc_projects/Pig_wound/anvio_pig/anvio_final/anvi_final_nocontam_pangenome/all_final_nocontam-GENOMES.db
