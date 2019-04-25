#! /usr/bin/env python2.7

import sys
import os, shutil
import subprocess as sp

#########################################################################################################################

samtools = "/Users/nk9/Programs/samtools-1.3.1/samtools"
bcftools = "/Users/nk9/Programs/bcftools-1.3.1/bcftools"
mge = "scripts/MGE_rem_staph15.bed"
ref = "scripts/reference.fa"

#########################################################################################################################
######################### Usage : mec_cov_analysis.py downloadFloder summaryFile.txt outputfile #########################
#########################################################################################################################


in1 = sys.argv[1] ### folder where the files have been downloaded
#in2 = sys.argv[2] ### summary file downloaded for the run
out1 = sys.argv[2] ### Name for the analysis output file

#########################################################################################################################


#summary_db = {base.split("\t")[0]:base.strip().split() for base in open(in2).readlines()}

bamFiles = [base for base in os.listdir(in1) if base.endswith(".bam") and not base.startswith(".")]

mecFiles = [base for base in os.listdir(in1) if base.endswith("report.tsv") and not base.startswith(".")]

report = {}
for each in mecFiles:
	mec_data = [base.strip().split("\t") for base in open(in1+"/"+each).readlines()[1:]][0]
	name = each.split("_")[0]
	if len(mec_data) >= 1:
		t1 = mec_data
		cov = "{0:.2f}".format(int(t1[8])/float(t1[7]))
		inf1 = [t1[0], t1[5], cov, t1[9], t1[12]]
		report[name] = inf1
	else:
		inf1 = "- - - - - -".split(" ")
		report[name] += inf1

	

os.mkdir("results")

header = "strainID hetSnps homSnps hethomRatio refCov20x geneName reads geneCovered identity mecAcovDepth genomeCovDepth stdDev"

mec_report = ["\t".join(header.split(" "))]

for each1 in bamFiles:

	each = each1.split(".")[0]

	os.mkdir("results/"+each)

	os.chdir("results/"+each)

	bamcopy, err = (sp.Popen("cp ../../%s/%s ./" % (in1, each1), stdout=sp.PIPE, shell=True)).communicate()

	#bamcopy, err = (sp.Popen("cp ../../%s/%s/ariba/mec.report.tsv ./" % (in1, each), stdout=sp.PIPE, shell=True)).communicate()
	

	faIndex, err = (sp.Popen("%s faidx ../../%s" % (samtools, ref), stdout=sp.PIPE, shell=True)).communicate()

	samSort, err = (sp.Popen("%s sort -T temp -O BAM -o samtools_sorted.bam %s" % (samtools,each1), stdout=sp.PIPE, shell=True)).communicate()

	samMpileup, err = (sp.Popen("%s mpileup -d 8000 -t AD,DP,SP -Bug -f ../../%s -o samtools.bcf samtools_sorted.bam" % (samtools, ref),stdout=sp.PIPE, shell=True)).communicate()

	#print err

	samDeth, err = (sp.Popen("%s depth -aa -Q 30 samtools_sorted.bam | awk -F '\t' '{ if ( $3 >= 20 ) print $0}' | wc -l > cov" % samtools, stdout=sp.PIPE, shell=True)).communicate()
		
	samInd, err = (sp.Popen("%s index samtools_sorted.bam" % samtools, stdout=sp.PIPE, shell=True)).communicate()

	covstd, err = (sp.Popen("%s depth -aa -Q 30 -b ../../%s samtools_sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}'" % (samtools, mge), stdout=sp.PIPE, shell=True)).communicate()
	
	varCalls, err = (sp.Popen("%s call -P 0.01 -mv -o samtools.vcf -O v samtools.bcf" % bcftools, stdout=sp.PIPE, shell=True)).communicate()

	remClust, err = (sp.Popen("python2.7 ../../scripts/filtervcf_v4.py samtools.vcf %s ../../scripts/MGE_staph.bed cov 50" % each, stdout=sp.PIPE, shell=True)).communicate()
		
	stats = (open(each+"_50stats").readlines()[1]).strip().split("\t")
	
	if each in report:
		z1 = [each] + stats+ report[each]+ covstd.strip().split(" ")
	else:
		z1 = [each] + stats+ "- - - - -".split(" ")+ covstd.strip().split(" ")
	print (covstd.strip().split(" "))
	mec_report.append("\t".join(z1))
#		

	#os.remove("mapping.bam")
	os.remove("samtools_sorted.bam")
	os.remove("samtools.bcf")
	os.remove("cov")
	os.remove("samtools_sorted.bam.bai")
	os.chdir("../../")
	
open(out1, "a").writelines("\n".join(mec_report)+"\n")	

shutil.rmtree("results")