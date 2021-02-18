#! /usr/bin/env python3
from optparse import OptionParser
import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)
import argparse
import subprocess


class annotation_record:
	def __init__(self, seqname, source, feature, start, end, score, strand, frame, group):
		self.seqname = seqname
		self.source = source
		self.feature = feature
		self.start = start
		self.end = end
		self.score = score
		self.strand = strand
		self.frame = frame
		self.group = self.group2dict(group)
	def group2dict(self, group):
		retdict = {}
		for each_item in group.split(";"):
			key, value = each_item.split("=")
			retdict[key] = value
		return retdict
	def __str__(self):
		retstr = ""
		retstr += (f"seqname:{self.seqname}\t")
		retstr += (f"source:{self.source}\t")
		retstr += (f"feature:{self.feature}\t")
		retstr += (f"start:{self.start}\t")
		retstr += (f"end:{self.end}\t")
		retstr += (f"score:{self.score}\t")
		retstr += (f"strand:{self.strand}\t")
		retstr += (f"frame:{self.frame}\t")
		group_str = ""
		for each_item in self.group:
			group_str += each_item + ":" + self.group[each_item] + "\t"
		retstr += (f"{group_str}\t")
		return retstr


def gff3_parser(filename):
	return_array = []
	with open(filename, "r") as f:
		for each_line in f:
			if each_line.startswith("#"):
				continue
			seqname, source, feature, start, end, score, strand, frame, group = each_line.split()
			return_array.append(annotation_record(seqname, source, feature, start, end, score, strand, frame, group))
	return return_array

def gff3toBed(gff3_array):
	bed4_array = []
	for line in gff3_array:
		if int(line.end) - int(line.start) >= 100:
			bed4_array.append("\t".join([line.seqname, str(line.start), str(line.end), ":".join([line.seqname, line.feature, line.group["gene_name"], str(line.start), str(line.end)])]))
	return bed4_array



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "parsing GFF3 format gene annotation data")
	parser.add_argument("GFF3", metavar = "GFF3", type = str, help = "input  GFF3  file name")
	parser.add_argument("Hg19", metavar = "Hg19", type = str, help = "input  FASTA file name which contains all chromosome")
	parser.add_argument("BED4", metavar = "BED4", type = str, help = "output BED4  file name")
	parser.add_argument("FASTA", metavar = "FASTA", type = str, help = "output FASTA file name")
	args = parser.parse_args()
	gff3_filename = args.GFF3
	bed4_filename = args.BED4
	hg19_filename = args.Hg19
	fasta_filename = args.FASTA
	print(args, file = sys.stderr)
	allrecord = gff3_parser(gff3_filename)
	bed4_array = gff3toBed(allrecord)
	with open(bed4_filename, "w") as f:
		for each_annotation in bed4_array:
			print(each_annotation, file = f)
	seqkit_cmd = ["seqkit", "subseq", "--bed", bed4_filename, hg19_filename, "-o", fasta_filename]
	print(seqkit_cmd)
	subprocess.call(seqkit_cmd)
	sys.exit()
