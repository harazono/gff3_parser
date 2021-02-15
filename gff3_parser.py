#! /usr/bin/env python3
from optparse import OptionParser
import sys
import pprint
pp = pprint.PrettyPrinter(indent=4)



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
		retstr += (f"seqname: {self.seqname}\t")
		retstr += (f"source: {self.source}\t")
		retstr += (f"feature: {self.feature}\t")
		retstr += (f"start: {self.start}\t")
		retstr += (f"end: {self.end}\t")
		retstr += (f"score: {self.score}\t")
		retstr += (f"strand: {self.strand}\t")
		retstr += (f"frame: {self.frame}\t")
		group_str = ""
		for each_item in self.group:
			group_str += each_item + ":" + self.group[each_item] + "\t"
		retstr += (f"group: {group_str}\t")
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

"""
def gff3_info_parser(inputdict_array):
	for each_dict in inputdict_array:
		info = each_dict["info"]
		each_item_of_info = info.split(";")
		tmpdict = {}
		for each in each_item_of_info:
			key, value = each.split("=")
			tmpdict[key] = value
		each_dict["info"] = tmpdict
"""
def extract_gene_from_file(filename):
	grch37_dict = gff3_parser(filename)
	gff3_info_parser(grch37_dict)
	#print("chrom, start, end, gene_name, gene_status")
	for each_gene in grch37_dict:
		if each_gene['info']['gene_type'] == "protein_coding" and int(each_gene['start']) < int(each_gene['end']):
			#print(f"{each_gene['chrom'].replace('chr', '')}, {each_gene['start']}, {each_gene['end']}, {each_gene['info']['gene_name']}, {each_gene['info']['gene_status']}")
			print(f"{each_gene['chrom']}\t{each_gene['start']}\t{each_gene['end']}\t{each_gene['info']['gene_name']}\t\t")

if __name__ == "__main__":
	sys.exit() if len(sys.argv) < 2
	print(f"{sys.argv[1]} was given", file = sys.stderr)
	allrecord = gff3_parser(sys.argv[1])
	for each_annotation in allrecord:
		print(each_annotation)

	sys.exit()
	#genes = extract_gene_from_file("gencode.v19.gene_annotation.gff3")



"""
def gff3_dict_modifier(inputdict_array):
	ret_dict = {}
	for each_dict in inputdict_array:
		gene = each_dict["info"]["gene_name"]
		chrom = each_dict["chrom"]
		start = each_dict["start"]
		end = each_dict["end"]
		if not (chrom, gene) in ret_dict.keys():
			ret_dict[(chrom, gene)] = [(chrom, start, end)]
		else:
			ret_dict[(chrom, gene)].append((chrom, start, end))
			#print(f"gene name {gene} at chrom {chrom} is already enrolled in a dictionary.\n{ret_dict[(chrom, gene)][0]}: {ret_dict[(chrom, gene)][1]}-{ret_dict[(chrom, gene)][2]}\n{chrom}: {start}-{end}\n", file = sys.stderr)
	return ret_dict

def get_union(dict1, dict2):
	retarray = []
	for gene in dict1:
		if gene in dict2:
			retarray.append(gene)
"""



"""
def main():
	if len(sys.argv) < 2:
		sys.exit()

	grch37_annotations_gff3_filepath = sys.argv[1]
	grch38_annotations_gff3_filepath = sys.argv[2]

	grch37_dict = gff3_parser(grch37_annotations_gff3_filepath)
	grch38_dict = gff3_parser(grch38_annotations_gff3_filepath)
	gff3_info_parser(grch37_dict)
	gff3_info_parser(grch38_dict)
	grch37_modified_dict = gff3_dict_modifier(grch37_dict)
	grch38_modified_dict = gff3_dict_modifier(grch38_dict)

	gene_in_both = grch37_modified_dict.keys() & grch38_modified_dict.keys()
	gene_in_only37 = grch37_modified_dict.keys() - grch38_modified_dict.keys()
	gene_in_only38 = grch38_modified_dict.keys() - grch37_modified_dict.keys()

	print(f"GRCh37(total)\t{len(grch37_dict)}")
	print(f"GRCh38(total)\t{len(grch38_dict)}")
	print(f"GRCh37(unique)\t{len(grch37_modified_dict)}")
	print(f"GRCh38(unique)\t{len(grch38_modified_dict)}")
	print(f"both reference\t{len(gene_in_both)}")
	print(f"only GRCh37\t{len(gene_in_only37)}")
	print(f"only GRCh38\t{len(gene_in_only38)}")
	multiple_loci_37 = [x for x in grch37_modified_dict if len(grch37_modified_dict[x]) > 1]
	multiple_loci_38 = [x for x in grch38_modified_dict if len(grch38_modified_dict[x]) > 1]
	print(f"registered in multiple loci(GRCh37)\t{len(multiple_loci_37)}")
	print(f"registered in multiple loci(GRCh38)\t{len(multiple_loci_38)}")
	print(f"registered in multiple loci(both)\t{len(set(multiple_loci_37) & set(multiple_loci_38))}")
	#print(f"registered in multiple loci(GRCh38): {len([x for x in grch37_modified_dict if len(grch38_modified_dict[x]) > 1])}")

	for gene in gene_in_both:
		if len(grch37_modified_dict[gene]) > 1 and len(grch38_modified_dict[gene]) > 1:
			output_str_37 = "\t".join([f"{x[0]}:{x[1]}-{x[2]}" for x in grch37_modified_dict[gene]])
			output_str_38 = "\t".join([f"{x[0]}:{x[1]}-{x[2]}" for x in grch38_modified_dict[gene]])
			print(f"{gene}\n37: {len(grch37_modified_dict[gene])}箇所\t{output_str_37}\n38: {len(grch38_modified_dict[gene])}箇所\t{output_str_38}\n")
"""

if __name__ == "__main__":
	extract_gene_from_file("gencode.v19.annotation.gff3")