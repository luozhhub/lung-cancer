#!/usr/bin/python3

import os
import sys,re
from subprocess import *
import pandas as pd

class quality_control():

    def __init__(self):
        self.fastq_dir = "/home/zhluo/Project/lung_cancer/data"
        self.Trimmomatic = "/home/nazhang/luozhihui/software/Trimmomatic-0.38/trimmomatic-0.38.jar"
        self.outputDir = "/home/nazhang/luozhihui/project/CRC/trimmoResult"
        self.adaptor = "/home/nazhang/luozhihui/project/CRC/TruSeq2-PE.fa"
        self.ref = "/home/zhluo/Project/CRC/data_nazhang/refernece_genome/bowtie_index/GRCm38.primary_assembly.genome"

    def run(self, cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p

    def fastp(self, fastq1=None, fastq2=None, outputDir=None):
        sample = os.path.basename(fastq1).split(".")[0]
        output1 = sample + ".fastp.fq.gz"
        output1 = os.path.join(outputDir, output1)
        sample = os.path.basename(fastq2).split(".")[0]
        output2 = sample + ".fastp.fq.gz"
        output2 = os.path.join(outputDir, output2)
        sample = sample.split("_")[0]
        htmlPath = os.path.join(outputDir, sample + ".html")
        jsonPath = os.path.join(outputDir, sample + ".json")
        cmd = "fastp -z 4 -i %s -I %s -o %s -O %s -h %s -j %s" % (fastq1, fastq2, output1, output2, htmlPath, jsonPath)
        return cmd

    def extract(self, zipFile=None):
        """
        unzip the fastq file
        """

    def trimAdapterByTrimmomatic(self, fastq1=None, fastq2=None, outputDir=None):
        """
        fastq eg. "823-RNA_L8_1.fq.gz"
        ILLUMINACLIP:%s:2:30:10:
        %s is the adaptor file
        2 is bigest mismatch
        30 is palindrome method threshold value
        
        :param fastq1: string
        :param fastq2: string
        :return:
        """
        fastq1_name = os.path.basename(fastq1)
        fastq2_name = os.path.basename(fastq2)
        forward_paired = os.path.join(outputDir, fastq1_name.replace("_R1.fastp.fq.gz", "_1_paired.fq.gz"))
        forward_unpaired = os.path.join(outputDir, fastq1_name.replace("_R1.fastp.fq.gz", "_1_unpaired.fq.gz"))
        reverse_paired = os.path.join(outputDir, fastq2_name.replace("_R2.fastp.fq.gz", "_2_paired.fq.gz"))
        reverse_unpaired = os.path.join(outputDir, fastq2_name.replace("_R2.fastp.fq.gz", "_2_unpaired.fq.gz"))
        logFile = os.path.join(outputDir, fastq1_name.replace("_R1.fastp.fq.gz", ".log"))
        cmd = "java -jar %s PE -threads 4 -phred33 -trimlog %s\
         %s  %s %s %s %s %s \
         ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50" % \
              (self.Trimmomatic, logFile, fastq1, fastq2, forward_paired, forward_unpaired, reverse_paired, reverse_unpaired, self.adaptor)
        return cmd

    def fastqc(self, fastq=None, outputDir=None):
        #dirName = os.path.dirname(fastq)
        cmd = "/home/nazhang/source/FastQC/fastqc -t 6 -o %s %s" % (outputDir, fastq )
        return cmd

    def star(self, fastq1=None, fastq2=None, outputDir=None):
        fastq1_name = os.path.basename(fastq1)
        prefix = os.path.join(outputDir, fastq1_name.replace("_1_paired.fq.gz", ""))
        cmd = "STAR --genomeDir star_genome --readFilesIn %s %s \
        --readFilesCommand zcat --outSAMstrandField intronMotif --runThreadN 8 --outFileNamePrefix %s" % \
              (fastq1, fastq2, prefix)
        return cmd

    def run_bwa(self, core=6, ref=None, fastq_1=None, fastq_2=None, outsam=None, rg=None):
        cmd = "%s mem -R '@RG\\tID:1\\tPL:ILLUMINA\\tSM:%s' -t %s %s %s %s > %s"%("bwa", rg, core, ref, fastq_1, fastq_2, outsam)
        return cmd
        
    def run_bowtie(self, core=10, ref="", fastq_1=None, fastq_2=None, outbam=None):
        cmd = "bowtie2 -p %s -x %s -1 %s -2 %s | samtools view -Sb -q 10 - | samtools sort -O bam -@ %s -o - > %s" % (core, ref, fastq_1, fastq_2, core, outbam)
        return(cmd)
        
    def sort_bed(self, inputFile=None, outputFile=None):
        cmd = "sort -k1,1 -k2,2n %s > %s" %(inputFile, outputFile)
        self.run(cmd=cmd)
        
    def merge_bed(self, inputFile=None, outputFile=None):
        cmd = "bedtools merge -i %s >%s" %(inputFile, outputFile)
        self.run(cmd=cmd)  
        
    def merge_all_peaks(self, select_peak_files=None, state_total_file=None, zcat=True):
        files_str = " ".join(select_peak_files)
        if zcat is True:
            cmd = "zcat %s > %s" % (files_str, state_total_file)
        else:
            cmd = "cat %s > %s" % (files_str, state_total_file)
        self.run(cmd)
        
        sort_file = state_total_file + ".sort"
        self.sort_bed(inputFile = state_total_file, outputFile = sort_file)

        merged_file = sort_file + ".merged"
        self.merge_bed(inputFile = sort_file, outputFile = merged_file)
 

    def run_multiBamSummary(self):
        #calculate read count
        #bam files in /home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles
        bamDir = "/home/zhluo/Project/lung_cancer/step2_bam/"
        output_dir = "/home/zhluo/Project/lung_cancer/step5_read_count"
        sample_list = ["/home/zhluo/Project/lung_cancer/step2_bam/AM1_TAAGGCGA.bam",  "/home/zhluo/Project/lung_cancer/step2_bam/AM2_AGGCAGAA.bam",
         "/home/zhluo/Project/lung_cancer/step2_bam/AM3_TCCTGAGC.bam", "/home/zhluo/Project/lung_cancer/step2_bam/AM4_GGACTCCT.bam", "/home/zhluo/Project/lung_cancer/step2_bam/AM5_CTCTCTAC.bam",
         "/home/zhluo/Project/lung_cancer/step2_bam/AM6_CAGAGAGG.bam", "/home/zhluo/Project/lung_cancer/step2_bam/AM7_GCTACGCT.bam", "/home/zhluo/Project/lung_cancer/step2_bam/IM1_CGTACTAG.bam",
         "/home/zhluo/Project/lung_cancer/step2_bam/IM2_TAGGCATG.bam", "/home/zhluo/Project/lung_cancer/step2_bam/IM4_CGAGGCTG.bam", "/home/zhluo/Project/lung_cancer/step2_bam/IM5_AAGAGGCA.bam",
         "/home/zhluo/Project/lung_cancer/step2_bam/AM6_CAGAGAGG.bam"]
        lable_list = ["AM1", "AM2", "AM3", "AM4", "AM5", "AM6", "AM7", "IM1", "IM2", "IM4", "IM5", "IM6"] 
         
        
        bamfiles = " ".join(sample_list)
        lables = " ".join(lable_list)
        cmd = "multiBamSummary BED-file --BED /home/zhluo/Project/lung_cancer/step4_mergedbed/master.bed.sort.merged --numberOfProcessors 25 --bamfiles %s  --minMappingQuality 30 --labels %s -out %s --outRawCounts %s" %(bamfiles, \
        lables, os.path.join(output_dir, "readCounts.npz"), os.path.join(output_dir, "readCounts.tab"))
        pbs_handle = open("/home/zhluo/Project/lung_cancer/step5_read_count/multibamsummary.pbs", "w")
        pbs_handle.write(cmd)
        pbs_handle.close()
        
        
    def modify_read_count_file(self):
        #remove the "#" first
        df = pd.read_csv("/home/zhluo/Project/lung_cancer/step5_read_count/readCounts.tab", header=0, sep="\t", quotechar="'")
    
        df_sub = df.drop(df.columns[[0,1, 2]], axis=1)
        #print(df_sub[0:5])
        peaks = ["peak%s" % index for index, row in df_sub.iterrows()]
        #print(peaks)
        df_sub.index = peaks
        #df_sub.insert(0, 'peak', peaks)
        print(df_sub[0:5])
    
        df_sub.to_csv("/home/zhluo/Project/lung_cancer/step5_read_count/deseq2_read_count.txt", header=True, index=True)
        
    def modify_merged_bed(self):
        df = pd.read_csv("/home/zhluo/Project/lung_cancer/step5_read_count/readCounts.tab", header=0, sep="\t", quotechar="'")
    
        df_sub = df.iloc[:,[0,1, 2]]
        #print(df_sub[0:5])
        peaks = ["peak%s" % index for index, row in df_sub.iterrows()]
        #print(peaks)
        df_sub.index = peaks
        #df_sub.insert(0, 'peak', peaks)
        print(df_sub[0:5])
        
        df_sub.columns = ["chr", "start", "end"]
        df_sub["peak_ID"] = peaks
        df_sub["score"] = [1000 for i in range(0, len(df_sub))]
        df_sub["strand"] = ["+" for i in range(0, len(df_sub))]
        print(df_sub[0:5])
        df_sub.to_csv("/home/zhluo/Project/lung_cancer/step4_mergedbed/pre_annot_peak.bed", header=False, index=False, sep="\t")
        
        
if __name__ == "__main__":
    qc = quality_control()
    step = 6
    #step 1, generate fastqc file
    if step < 1:
        fastqFiles = os.listdir(qc.fastq_dir)
        for fi in fastqFiles:
            fastq = os.path.join(qc.fastq_dir, fi)
            cmd = qc.fastqc(fastq=fastq, outputDir="/home/zhluo/Project/lung_cancer/step1_qc")
            qc.run(cmd=cmd)
            
            
    if step < 2:
        fastq_dir = "/home/zhluo/Project/lung_cancer/fastq_ln"
        fastqFiles = os.listdir(fastq_dir)
        if len(fastqFiles) == 0:
            exit(1)

        sample_num = []
        for file in fastqFiles:
            if re.search(".log", file):
                continue
            regx = re.compile("(.*)_R(.*).fastq.gz")
            result = regx.search(file)
            if not result:
                continue
            number = result.groups()[0]
            if number not in sample_num:
                sample_num.append(number)

        for sample in sample_num:
            fastq1 = "%s_R1.fastq.gz" % (sample)
            fastq1 = os.path.join(fastq_dir, fastq1)
            fastq2 = "%s_R2.fastq.gz" % (sample)
            fastq2 = os.path.join(fastq_dir, fastq2)
            bam = os.path.join("/home/zhluo/Project/lung_cancer/step2_bam", sample + ".bam")
            print (sample)
            cmd = qc.run_bowtie(core=10, ref=qc.ref, fastq_1=fastq1, fastq_2=fastq2, outbam=bam)
            
            OP = open("/home/zhluo/Project/lung_cancer/pbs/bowtie_" + sample + ".pbs", "w")
            OP.write("cd /home/zhluo/Project/lung_cancer/pbs;" + cmd + "\n")
            OP.close()
    

    if step < 3:
        bam_dir = "/home/zhluo/Project/lung_cancer/step2_bam/"
        """
        nohup macs2 callpeak --treatment AM1_TAAGGCGA.bam AM2_AGGCAGAA.bam AM3_TCCTGAGC.bam AM4_GGACTCCT.bam --name AM_uninfect -g mm  -q 0.05 --keep-dup all -f BAMPE -B --nomodel --outdir /home/zhluo/Project/lung_cancer/step3_macs2/AM_uninfect &
        nohup macs2 callpeak --treatment AM5_CTCTCTAC.bam AM6_CAGAGAGG.bam AM7_GCTACGCT.bam --name AM_infect -g mm  -q 0.05 --keep-dup all -f BAMPE -B --nomodel --outdir /home/zhluo/Project/lung_cancer/step3_macs2/AM_infect &
        nohup macs2 callpeak --treatment IM1_CGTACTAG.bam IM2_TAGGCATG.bam --name IM_uninfect -g mm  -q 0.05 --keep-dup all -f BAMPE -B --nomodel --outdir /home/zhluo/Project/lung_cancer/step3_macs2/IM_uninfect &
        nohup macs2 callpeak --treatment IM4_CGAGGCTG.bam IM5_AAGAGGCA.bam IM6_GTAGAGGA.bam --name IM_uninfect -g mm  -q 0.05 --keep-dup all -f BAMPE -B --nomodel --outdir /home/zhluo/Project/lung_cancer/step3_macs2/IM_infect &
        """
        
    if step < 4:
        bed_files = ["/home/zhluo/Project/lung_cancer/step3_macs2/AM_infect/AM_infect_peaks.narrowPeak", "/home/zhluo/Project/lung_cancer/step3_macs2/AM_uninfect/AM_uninfect_peaks.narrowPeak", "/home/zhluo/Project/lung_cancer/step3_macs2/IM_infect/IM_uninfect_peaks.narrowPeak", "/home/zhluo/Project/lung_cancer/step3_macs2/IM_uninfect/IM_uninfect_peaks.narrowPeak" ]
        qc.merge_all_peaks(select_peak_files=bed_files, state_total_file="/home/zhluo/Project/lung_cancer/step4_mergedbed/master.bed", zcat=False)
        """
        nohup annotatePeaks.pl /home/zhluo/Project/lung_cancer/step4_mergedbed/master.bed.sort.merged hg19 > /home/zhluo/Project/lung_cancer/step4_mergedbed/master.bed.sort.merged.anno &
        """
        
    if step < 5:
        qc.run_multiBamSummary()
        
     
    if step < 6:
        qc.modify_read_count_file()
        
        
    if step < 7:
        qc.modify_merged_bed()
        """
        nohup annotatePeaks.pl /home/zhluo/Project/lung_cancer/step4_mergedbed/pre_annot_peak.bed mm10 > /home/zhluo/Project/lung_cancer/step4_mergedbed/aft_annot_peak.bed &
        """
    exit(1)
        
        
        
        