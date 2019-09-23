import sys
import numpy as np
from sys import *
import os
import time
import dendropy
from dendropy.calculate import treemeasure
import argparse



#Ensure all of the required files and settings were input correctly
#The fasta file identifier is required to be in this format, where "Country' can be any grouping variable:
#>UniqueName_Country
#Country must be the first item after first underscore
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--infile", required = True, help = "Specify path to input newick file")
ap.add_argument("-o1", "--out1", required = True, help = "Specify total patristic distance output filename")
ap.add_argument("-o2", "--out2", required = True, help = "Specify organized patristic distance spectrum output filename")
args = ap.parse_args()

infile =  (args.infile)
PatDistOutput = (args.out1)
outfile = (args.out2)


#Generate pairwise pat distances for all variants below the specified threshold 
def dendro(infile, dendroOut):
    PatDistOutput = open(dendroOut,'w')
    tree = dendropy.Tree.get_from_path(infile, "newick")
    pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree)
    for i, t1 in enumerate(tree.taxon_namespace):
        for t2 in tree.taxon_namespace[i+1:]:
            PatDistOutput.write("%s,%s,%s,%s" % (t1.label.replace(' ','_'), t2.label.replace(' ','_'), pdm(t1, t2),str(pdm.mrca(t1, t2).label) )+'\n')

#Store pat dist in a dictionary for subsequent percentile calculations
def PatDistAppend(line):
    linesplit = line.split(",")
    node = linesplit[0]
    mate = linesplit[1]
    patdist = linesplit[2].rstrip()
    nodeshort = node.split('_')[1]
    mateshort = mate.split('_')[1]
    comb = nodeshort+"__"+mateshort
    comb2 = mateshort+"__"+nodeshort
    shortdouble = nodeshort+"__"+mateshort
    if nodeshort != mateshort:      
        if not comb2 in dictInterPatDist:
            if not comb in dictInterPatDist:
                dictInterPatDist[comb]=[]
            dictInterPatDist[comb].append(float(patdist))
        else:
            dictInterPatDist[comb2].append(float(patdist))
    else:
        if not shortdouble in dictInterPatDist:
            dictInterPatDist[shortdouble]=[]
        dictInterPatDist[shortdouble].append(float(patdist))
    return dictInterPatDist
 
#Calculate pat dist percentiles   
def Spectrum():
    Percentiles={}
    for k in dictInterPatDist:
        dictInterPatDist[k].sort()
        a = np.asarray(dictInterPatDist[k])
        Q0 = np.percentile(a, 0)
        Q1 = np.percentile(a, 1)
        Q5 = np.percentile(a, 5)
        Q10 = np.percentile(a, 10)
        Q20 = np.percentile(a, 20)
        Q25 = np.percentile(a, 25)
        Q30 = np.percentile(a, 30)
        Q35 = np.percentile(a, 35)
        Q40 = np.percentile(a, 40)
        Q45 = np.percentile(a, 45)
        Q50 = np.percentile(a, 50)
        Q75 = np.percentile(a, 75)
        Q90 = np.percentile(a, 90)
        Q99 = np.percentile(a, 99)
        Q100 = np.percentile(a, 100)
        if not k in Percentiles:
            Percentiles[k]=[]
        Percentiles[k].extend([Q0,Q1,Q5,Q10,Q20,Q25,Q30,Q35,Q40,Q45,Q50,Q75,Q90,Q99,Q100])
    return Percentiles
                
if __name__ == '__main__':
    
    if not os.path.isfile(PatDistOutput):
        if os.path.isfile(infile):
            dendro(infile,PatDistOutput)
        else:
            "WARNING: Cannot locate input tree file"
            sys.exit(1)
    else:
        print ("Pairwise patristic distances extracted already for %s, calculating percentiles..." %(PatDistOutput))
    
    dictInterPatDist={}
    outfile=open(outfile,'w')
    outfile.write("Samples,0,1,5,10,20,25,30,35,40,45,50,75,90,99,100"+"\n")

    with open(PatDistOutput) as f:
        for line in f:
            PatDistAppend(line.rstrip("\n"))
        PDspectrum = Spectrum()
        for k in PDspectrum:
            outfile.write(k +","+str(PDspectrum[k]).replace("[","").replace("]","")+"\n")