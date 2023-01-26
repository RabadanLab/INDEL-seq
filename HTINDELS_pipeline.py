#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
from multiprocessing import Pool

fq1 = sys.argv[1]
fq2 = sys.argv[2]
barcodes = sys.argv[3]

##preprocess
os.makedirs('preprocess/',exist_ok=True)
#cut the first 6 bp NNNNNN
subprocess.call('cutadapt -u 6 -U 6 -j 0 -o preprocess/' + os.path.basename(fq1).split('.')[0] + '.trimmed.fastq.gz -p preprocess/' + os.path.basename(fq2).split('.')[0] + '.trimmed.fastq.gz ' + fq1 + ' ' + fq2,shell=True)

#demultx
subprocess.call('/opt/fastq-multx/fastq-multx -b -B ' + barcodes + ' preprocess/' + os.path.basename(fq1).split('.')[0] + '.trimmed.fastq.gz preprocess/' + os.path.basename(fq2).split('.')[0] + '.trimmed.fastq.gz -o preprocess/%_R1.fq.gz preprocess/%_R2.fq.gz',shell=True)

#merge forward and reverse paired-end reads
def merge(ID):
    subprocess.call('/opt/pear-0.9.11-linux-x86_64/bin/pear -f preprocess/' + ID + '_R1.fq.gz -r preprocess/' + ID + '_R2.fq.gz -o preprocess/' + ID,shell=True)
N = os.cpu_count()
IDs = list(pd.read_csv(barcodes,sep="\t",header=None)[0].astype(str).values)
with Pool(N) as p:
    p.map(merge, IDs)

#alignment
os.makedirs('results/',exist_ok=True)

def run_needle(ID):
    subprocess.call('/opt/EMBOSS-6.6.0/emboss/needle TRIM37.fa preprocess/' + ID + '.assembled.fastq -outfile results/' + ID + '.sam -aformat sam -gapopen 10 -gapextend 0.5',shell=True)
with Pool(N) as p:
    p.map(run_needle, IDs)

def run_gzip(ID):
    subprocess.call('gzip results/' + ID + '.sam',shell=True)
with Pool(N) as p:
    p.map(run_gzip, IDs)

