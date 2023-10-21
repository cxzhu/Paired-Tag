#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='Add barcodes and umi from Paired-Tag readname to field CB and UMI.')
parser.add_argument('-i', type=str, dest="input", help='input Paired-Tag bam file')
parser.add_argument('-p', type=int, dest="cpu", default = 1, help='number of process to parallel. Default 1.')
parser.add_argument('--bc_len', type=int, dest="bc_len", default = "11", help='Length of barcodes is 8 or 11? Default 11')
parser.add_argument('--barcode_tag', type=str, dest="bc_fild", default = "CB", help='barcode tag assigned. Default CB')
parser.add_argument('--umi_tag', type=str, dest="umi_fild", default = "UB", help='umi tag assigned. Default UB')

args = parser.parse_args()

import os
import sys
import pysam
import multiprocessing
import re
from time import perf_counter as pc

def run():
    ''' Paired-Tag format to 10X format '''
    inbam = args.input
    cpu = args.cpu
    bc_len = args.bc_len
    bc_fild = args.bc_fild
    umi_fild = args.umi_fild
    start_time = pc()
    print("Add barcodes to field for file", inbam)
    fsam = generate_bams(inbam, bc_len, bc_fild, umi_fild, cpu)
    obam = "".join((inbam, ".10X.bam"))
    merge_command = [obam] + fsam
    pysam.merge(*merge_command)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_bams(inbam, bc_len, bc_fild, umi_fild, cpu):
    if bc_len != 11 and bc_len != 8:
        sys.exit("Currently barcode length should be 8 or 11. Exit")
    if cpu > 1:
        print("multithreading with", cpu, "threads...")
    pool = multiprocessing.Pool(cpu)
    ### parallel by chrom name
    inbamF = pysam.AlignmentFile(inbam, "rb")
    if not inbamF.has_index():
        pysam.index(inbam)
        inbamF.close()
    inbamF = pysam.AlignmentFile(inbam, "rb")
    idxstats = inbamF.get_index_statistics()
    inbamF.close()
    contigs = []
    fsam = []
    for i in idxstats:
        if i.mapped > 0:
            contigs.append(i.contig)
            fname = ".".join((inbam, i.contig, "tmp.bam"))
            fsam.append(fname)
    [pool.apply_async(generate_bams_worker, (inbam, bc_len, bc_fild, umi_fild, chrom, )) for chrom in contigs]
    pool.close()
    pool.join()
    return fsam

def generate_bams_worker(inbam, bc_len, bc_fild, umi_fild, chrom):
    tmpBam = ".".join((inbam, chrom, "tmp.bam"))
    inbamF = pysam.AlignmentFile(inbam, "rb")
    outb = pysam.AlignmentFile(tmpBam, "wb", template = inbamF)
    umi_regex = ":(\w+$)"
    for read in inbamF.fetch(chrom):
        ### determine barcode
        if bc_len == 11:
            bc_regex = "(..:..:..:..):\w+$"
        elif bc_len == 8:
            bc_regex = "(..:..:..):\w+$"
        bc = re.search(bc_regex, read.qname).group(1)
        umi = re.search(umi_regex, read.qname).group(1)
        read.tags = read.tags + [(bc_fild, bc)] + [(umi_fild, umi)]
        outb.write(read)
    inb.close()
    outb.close()

if __name__ == "__main__":
    run()


