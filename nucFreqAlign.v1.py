import sys
import math
import argparse
from Bio import AlignIO
#This takes in an alignment in fasta format and calculates the nucleotide frequencies at each position for each seq or group of seqs for a given category
#Further processing of data is required if indels are present in the reference genome that serves as the coordinate system
#The fasta file identifier is required to be in the following format, where "Country' can be any grouping variable:
#>UniqueName_Country
#ie, Country must be the first item following the underscore

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--infile", required = True, help = "Specify path to input alignment file in fasta format")
ap.add_argument("-t", "--type", required = True, help = "Specify file type: multiple sequences per group/country? '-t m', or single sequence per group? '-t s'")
ap.add_argument("-r", "--redundant", required = True, help = "Specify whether redundant bases are to be incorporated in analysis (T) or ignored (F) or '")
ap.add_argument("-o", "--outfile", required = True, help = "Specify output filename")
args = ap.parse_args()

infile =  (args.infile)
intype = (args.type).lower()
redundancy = (args.redundant).upper()
if intype != 'm' and intype != 's':
	print ('multiple sequences per country > -t m, single sequence > -t s')
	sys.exit(1)
if redundancy != 'T' and redundancy != 'F':
	print ('Redundancy analysis should be analyzed (T) or not (F)')
	sys.exit(1)
outfile = open((args.outfile),'w')

def baseCount(countryOrg,cpos,s):
	translate = {'A':0,'C':1,'G':2,'T':3,'N':4,'-':5}
	RedundantBases= {'R':['A','G'],'Y':['C','T'], 'K':['G','T'],'M':['A','C'],'S':['C','G'],'W':['A','T'],'B':['G','C', 'T'],'D':['G','A', 'T'],'H':['A','C', 'T'],'V':['G','C', 'A']}
	s = s.upper()
	if s in RedundantBases:
		if redundancy == 'T':
			for r in RedundantBases[s]:
				r_trans = translate[r]
				countryOrg[cpos][r_trans] = countryOrg[cpos][r_trans] + 1
		else:
			countryOrg[cpos][4] = countryOrg[cpos][4] + 1
	else:
		s_trans = translate[s]
		countryOrg[cpos][s_trans] = countryOrg[cpos][s_trans] + 1



countryOrg = {}
nucs = [0,0,0,0,0,0]
seqs=[]
countryCounts = {}

if intype == 'm':
	outfile.write('Country,CountryTotal,Pos,A,C,G,T,N,Gap'+'\n')
else:
	outfile.write('Country,Pos,A,C,G,T,N,Gap'+'\n')
for alignment in AlignIO.parse(infile, "fasta"):
	for record in alignment:
		if intype == 'm':
			try:
				country = record.id.split('_')[1]
				if country in countryCounts:
					countryCounts[country] += 1
				else:
					countryCounts[country] = 1
			except:
				print ('Country or grouping variable should be underscore delimited (ex: ">uid_country")')
				sys.exit(1)
		else:
			country = record.id
		seqLen = len(record.id)
		count = 0			
		for s in record.seq:
			count += 1
			cpos = country+','+str(count)
			if not cpos in countryOrg:
				countryOrg[cpos] = list(nucs)
			baseCount(countryOrg,cpos,s)
for p in countryOrg:
	modnucs = map(str, countryOrg[p])
	if intype == 'm':
		psp = p.split(',')
		countryTotal = str(countryCounts[psp[0]])
		outfile.write(psp[0]+','+countryTotal+','+psp[1]+','+','.join(modnucs)+'\n')
	else:
		outfile.write(p+','+','.join(modnucs)+'\n')