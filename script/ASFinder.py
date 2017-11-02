#! /usr/bin/env python

import sys, os
import pybedtools
from collections import defaultdict

def load_xcount(bam, bed):
    print('loading xcount....')
    a = pybedtools.BedTool(bam)
    b = pybedtools.BedTool(bed)
    c = a.intersect(b, bed=True, wa=True, wb=True, split=True, f=1.0)
    xcount = defaultdict(int)
    for x in c:
        xcount['|'.join(x[12:15])] += 1
    return xcount

def load_xcount_RI(bam, bed):
    print('loading xcount for RI....')
    a = pybedtools.BedTool(bam)
    b = pybedtools.BedTool(bed)
    c = a.intersect(b, bed=True, wa=True, wb=True, split=True)
    xcount = defaultdict(int)
    for x in c:
        xcount['|'.join(x[12:15])] += 1
    return xcount

def load_jcount(junc):
    print('loading jcount....')
    jcount = {}
    for line in open(junc):
        chrom, left, right, count = line.strip().split('\t')
        jcount['%s|%s|%s' % (chrom, left, right)] = count
    return jcount 


def ASFinder(bam, bed, bed_RI, junc, out_dir):
    xcount = load_xcount(bam, bed)
    jcount = load_jcount(junc)
    xcount_RI = load_xcount_RI(bam, bed_RI)
    AS_dir = out_dir
    print('searching AS....')
    
    in_file = open(bed)
    SE_out = open(AS_dir+'/SE.txt', 'w')
    RI_out = open(AS_dir+'/RI.txt', 'w')
    A5SS_out = open(AS_dir+'/A5SS.txt', 'w')
    A3SS_out = open(AS_dir+'/A3SS.txt', 'w')
    AFE_out = open(AS_dir+'/AFE.txt', 'w')
    ALE_out = open(AS_dir+'/ALE.txt', 'w')
    MXE_out = open(AS_dir+'/MXE.txt', 'w')

    SE_out.write('%s\t%s\t%s\t%s\t%s\n' % ('AS_event', 'exon', 'junc12', 'junc34', 'junc14'))
    RI_out.write('%s\t%s\t%s\n' % ('AS_event', 'intron', 'junc23'))
    A5SS_out.write('%s\t%s\t%s\t%s\n' % ('AS_event', 'exon', 'junc_short', 'junc_long'))
    A3SS_out.write('%s\t%s\t%s\t%s\n' % ('AS_event', 'exon', 'junc_short', 'junc_long'))
    AFE_out.write('%s\t%s\t%s\t%s\t%s\n' % ('AS_event', 'exon_short', 'junc_short', 'exon_long', 'junc_long'))
    ALE_out.write('%s\t%s\t%s\t%s\t%s\n' % ('AS_event', 'exon_short', 'junc_short', 'exon_long', 'junc_long'))
    MXE_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('AS_event', 'exon23', 'junc12', 'junc36', 'exon45', 'junc14', 'junc56'))
    
    events = {}
    for line in in_file:
        lst = line.strip('\n').split('\t')
        events[lst[-2]+'\t'+lst[-1]] = '%s|%s|%s' % (lst[0], lst[1], lst[2])
    
    for event in events:
        AS, AS_type = event.split('\t')
    
        if AS_type == 'SE':
            chrom, s1, s2, s3, s4, strand, gene = AS.split('|')
            exon = events[event]
            junc12 = '%s|%s|%s' % (chrom, s1, s2)
            junc34 = '%s|%s|%s' % (chrom, s3, s4)
            junc14 = '%s|%s|%s' % (chrom, s1, s4)
            exon_count = xcount.get(exon, '0')
            junc12_count = jcount.get(junc12, '0')
            junc34_count = jcount.get(junc34, '0')
            junc14_count = jcount.get(junc14, '0')
            SE_out.write('%s\t%s\t%s\t%s\t%s\n' % (AS, exon_count, junc12_count, junc34_count, junc14_count))
        
        elif AS_type == 'RI':
            chrom, s1, s2, s3, s4, strand, gene = AS.split('|')
            intron = events[event]                
            junc23 = '%s|%s|%s' % (chrom, s2, s3)
            intron_count = xcount_RI.get(intron, '0')        
            junc23_count = jcount.get(junc23, '0')
            RI_out.write('%s\t%s\t%s\n' % (AS, intron_count, junc23_count))
    
        elif AS_type == 'A5SS':
            chrom, s1, s2, s3, s4, strand, gene = AS.split('|')
            exon = events[event]        
            exon_count = xcount.get(exon, '0')
            if strand == '+':    
                junc24 = '%s|%s|%s' % (chrom, s2, s4)
                junc34 = '%s|%s|%s' % (chrom, s3, s4)
                junc24_count = jcount.get(junc24, '0')
                junc34_count = jcount.get(junc34, '0')
                A5SS_out.write('%s\t%s\t%s\t%s\n' % (AS, exon_count, junc34_count, junc24_count))
            else:
                junc12 = '%s|%s|%s' % (chrom, s1, s2)
                junc13 = '%s|%s|%s' % (chrom, s1, s3)
                junc12_count = jcount.get(junc12, '0')
                junc13_count = jcount.get(junc13, '0')
                A5SS_out.write('%s\t%s\t%s\t%s\n' % (AS, exon_count, junc12_count, junc13_count))
        
        elif AS_type == 'A3SS':
            chrom, s1, s2, s3, s4, strand, gene = AS.split('|')
            exon = events[event]        
            exon_count = xcount.get(exon, '0')
            if strand == '-':    
                junc24 = '%s|%s|%s' % (chrom, s2, s4)
                junc34 = '%s|%s|%s' % (chrom, s3, s4)
                junc24_count = jcount.get(junc24, '0')
                junc34_count = jcount.get(junc34, '0')
                A3SS_out.write('%s\t%s\t%s\t%s\n' % (AS, exon_count, junc34_count, junc24_count))
            else:
                junc12 = '%s|%s|%s' % (chrom, s1, s2)
                junc13 = '%s|%s|%s' % (chrom, s1, s3)
                junc12_count = jcount.get(junc12, '0')
                junc13_count = jcount.get(junc13, '0')
                A3SS_out.write('%s\t%s\t%s\t%s\n' % (AS, exon_count, junc12_count, junc13_count))
        
        elif AS_type == 'AFE':
            chrom, s1, s2, s3, s4, s5, strand, gene = AS.split('|')
            if strand == '+':    
                junc25 = '%s|%s|%s' % (chrom, s2, s5)
                junc45 = '%s|%s|%s' % (chrom, s4, s5)
                junc25_count = jcount.get(junc25, '0')
                junc45_count = jcount.get(junc45, '0')
                exon12 = '%s|%s|%s' % (chrom, s1, s2)
                exon34 = '%s|%s|%s' % (chrom, s3, s4)
                exon12_count = xcount.get(exon12, '0')
                exon34_count = xcount.get(exon34, '0')
                AFE_out.write('%s\t%s\t%s\t%s\t%s\n' % (AS, exon34_count, junc45_count, exon12_count, junc25_count))
            else:
                junc12 = '%s|%s|%s' % (chrom, s1, s2)
                junc14 = '%s|%s|%s' % (chrom, s1, s4)
                junc12_count = jcount.get(junc12, '0')
                junc14_count = jcount.get(junc14, '0')
                exon23 = '%s|%s|%s' % (chrom, s2, s3)
                exon45 = '%s|%s|%s' % (chrom, s4, s5)
                exon23_count = xcount.get(exon23, '0')
                exon45_count = xcount.get(exon45, '0')
                AFE_out.write('%s\t%s\t%s\t%s\t%s\n' % (AS, exon23_count, junc12_count, exon45_count, junc14_count))
    
        elif AS_type == 'ALE':
            chrom, s1, s2, s3, s4, s5, strand, gene = AS.split('|')
            if strand == '-':    
                junc25 = '%s|%s|%s' % (chrom, s2, s5)
                junc45 = '%s|%s|%s' % (chrom, s4, s5)
                junc25_count = jcount.get(junc25, '0')
                junc45_count = jcount.get(junc45, '0')
                exon12 = '%s|%s|%s' % (chrom, s1, s2)
                exon34 = '%s|%s|%s' % (chrom, s3, s4)
                exon12_count = xcount.get(exon12, '0')
                exon34_count = xcount.get(exon34, '0')
                ALE_out.write('%s\t%s\t%s\t%s\t%s\n' % (AS, exon34_count, junc45_count, exon12_count, junc25_count))
            else:
                junc12 = '%s|%s|%s' % (chrom, s1, s2)
                junc14 = '%s|%s|%s' % (chrom, s1, s4)
                junc12_count = jcount.get(junc12, '0')
                junc14_count = jcount.get(junc14, '0')
                exon23 = '%s|%s|%s' % (chrom, s2, s3)
                exon45 = '%s|%s|%s' % (chrom, s4, s5)
                exon23_count = xcount.get(exon23, '0')
                exon45_count = xcount.get(exon45, '0')
                ALE_out.write('%s\t%s\t%s\t%s\t%s\n' % (AS, exon23_count, junc12_count, exon45_count, junc14_count))
    
        elif AS_type == 'MXE':
            chrom, s1, s2, s3, s4, s5, s6, strand, gene = event.split('|')
            junc12 = '%s|%s|%s' % (chrom, s1, s2)
            junc36 = '%s|%s|%s' % (chrom, s3, s6)
            junc14 = '%s|%s|%s' % (chrom, s1, s4)
            junc56 = '%s|%s|%s' % (chrom, s5, s6)
            junc12_count = jcount.get(junc12, '0')
            junc36_count = jcount.get(junc36, '0')
            junc14_count = jcount.get(junc14, '0')
            junc56_count = jcount.get(junc56, '0')
            exon23 = '%s|%s|%s' % (chrom, s2, s3)
            exon45 = '%s|%s|%s' % (chrom, s4, s5)
            exon23_count = xcount.get(exon23, '0')
            exon45_count = xcount.get(exon45, '0')
            MXE_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (AS, exon23_count, junc12_count, junc36_count, exon45_count, junc14_count, junc56_count))


if __name__ == '__main__':
    
    bam = sys.argv[1]
    junc = sys.argv[2]
    bed = sys.argv[3]
    bed_RI = sys.argv[4]
    out_dir = sys.argv[5]
    tmp_dir = out_dir+'/tmp'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    pybedtools.set_tempdir(tmp_dir)
    ASFinder(bam, bed, bed_RI, junc, out_dir)
    os.system('rm -r %s' % tmp_dir)
    
