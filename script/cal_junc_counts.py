#! /usr/bin/env python

import sys

read2junc = {}
read2chrom = {}

junc2counts = {}

f = open(sys.argv[1])
for line in f:
	chrom, start, end, read, n, strand = line.strip().split('\t')
	start = int(start) + 1
	end = int(end)
	if read not in read2junc:
		read2junc[read] = [start, end]
		read2chrom[read] = chrom
	else:
		read2junc[read].extend([start, end])

for i in read2junc:
	junc = [str(x) for x in sorted(read2junc[i])]
	if len(junc) == 4:
		j12 = '%s\t%s\t%s' % (read2chrom[i], junc[1], junc[2])
		if j12 in junc2counts:
			junc2counts[j12] += 1
		else:
			junc2counts[j12] = 1
	else:
		j12 = '%s\t%s\t%s' % (read2chrom[i], junc[1], junc[2])
		j34 = '%s\t%s\t%s' % (read2chrom[i], junc[3], junc[4])
		if j12 in junc2counts:
			junc2counts[j12] += 0.5
		else:
			junc2counts[j12] = 0.5
		if j34 in junc2counts:
			junc2counts[j34] += 0.5
		else:
			junc2counts[j34] = 0.5

f.close()

for i in sorted(junc2counts.keys()):
	print('%s\t%s' % (i, junc2counts[i]))
