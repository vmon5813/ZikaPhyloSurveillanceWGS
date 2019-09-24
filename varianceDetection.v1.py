import sys
import pandas as pd
import argparse
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--infile", required = True, help = "Specify path to input csv file")
ap.add_argument("-t", "--type", required = True, help = "Specify file type: multiple sequences per group/country? '-t m', or single sequence per group? '-t s'")

args = ap.parse_args()

infile =  (args.infile)
intype = (args.type).lower()
if intype != 'm' and intype != 's':
	print ('multiple sequences per country > -t m, single sequence > -t s')
	sys.exit(1)
elif intype == 'm':
	type_flag = 'MultipleSeqs'
else:
	type_flag = 'SingleSeq'

df = pd.read_csv(infile)
shortName = infile.split('.csv')[0]

countryTotals = {}
countryCounts = {}
posInfo = {}
bases = {}
nucOut = open(shortName+'.'+type_flag+'.NucleotideVariance.csv','w')
nucOut.write('Country,Position,A,C,G,T,N,Gap,InterCountryVariance'+'\n')
countryOut = open(shortName+'.'+type_flag+'.CountryVariance.csv','w')
countryOut.write('Position,Base,CountryCount,CountryList'+'\n')
trans = {0:'A',1:'C',2:'G',3:'T',4:'N',5:'-'}
count = 0
def baseDictFormat(bases):
	translate = {'A':0,'C':1,'G':2,'T':3,'N':4,'-':5}
	counts = [0,0,0,0,0,0]
	for base in bases:
		c_trans = translate[base]
		counts[c_trans] = str(bases[base])
	return counts
def dictSetup(country,pos):
	if not country in posInfo:
		posInfo[country] = dict()
	if not pos in posInfo[country]:
		posInfo[country][pos]=dict({'A':0,'G':0,'C':0,'T':0,'N':0,'-':0})
def varCalc(d):
	varPos = []
	for p in d:
		ps = p.split('-')
		pos = ps[0]
		if len(set(d[p])) == 1:
			varPos.append(str(pos)+'Var')
	return varPos
def countryCounter(l):
	c = dict(Counter(l))
	counts = []
	for i in c:
		counts.append(i+':'+str(round(float(c[i])/float(countryTotals[i]),2)))
	return counts

print ('Collecting variant counts per group...')
for index, row in df.iterrows():
	pos = str(row['Pos'])
	country = row['Country']
	dictSetup(country,pos)
	if intype == 'm':
		countryTotals[row['Country']] = row['CountryTotal']
	nucs = [row['A'],row['C'],row['G'],row['T'],row['N'],row['Gap']]
	count = 0
	variants = []
	for n in nucs:
		if str(n) != '0':
			base = trans[int(count)]
			posInfo[country][pos][base] += 1
			uid = str(pos)+'-'+base
			for x in range(0,int(n)):
				countryCounts.setdefault(uid,[]).append(country)
			bases.setdefault(pos,[]).append(base)
		count += 1

#nucOut format:		
#Country	Pos	A	C	G	T	N	Gap	Variance
#countryOut format:
#Pos,Nuc,CountryCount,CountryList

print ('Completed data processing, writing files...')
posCheck = []
for c in posInfo:
	for pos in posInfo[c]:
		counts = baseDictFormat(posInfo[c][pos])
		if len(set(bases[pos])) > 1:
			nucOut.write(c+','+pos+','+','.join(counts)+',Variance'+'\n')
		else:
			nucOut.write(c+','+pos+','+','.join(counts)+',NoVariance'+'\n')
		for b in posInfo[c][pos]:
			if posInfo[c][pos][b] > 0:
				uidCheck = pos+'-'+b+'-'+country
				uid = pos+'-'+b
				if not uidCheck in posCheck:
					if intype == 'm':
						countryOut.write(pos+','+b+','+str(len(set(countryCounts[uid])))+','+','.join(countryCounter(countryCounts[uid]))+'\n')
					else:
						countryOut.write(pos+','+b+','+str(len(set(countryCounts[uid])))+','+','.join(set(countryCounts[uid]))+'\n')
					posCheck.append(uidCheck)
