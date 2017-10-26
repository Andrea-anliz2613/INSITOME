# Quality Control Check 
# wd: INSITOME/pipeline

# identify individuals whose reported sex in PED files does not match estimated sex 
# 1 = male, 2 = female, 0 = not determined 
echo "Identify individuals with inconsistent reported sexes"
../bin/plink --bfile ../data/hgdp/00_data/hgdp  --check-sex --out ../data/hgdp/01_QC/hgdp
grep PROBLEM ../data/hgdp/01_QC/hgdp.sexcheck > ../data/hgdp/01_QC/hgdp_problem

# Identification of individuals with elevated missing data rates or outlying heterozygosity rate 
# find missing SNP; output = hgdp.imiss and hgdp.lmiss file
echo "Identify individuals with elevated missing data rates"

../bin/plink --bfile ../data/hgdp/00_data/hgdp --missing --out ../data/hgdp/01_QC/hgdp
# find heterozygous genotypes; output = hgdp.het file
echo "Identify individuals with elevated outlying heterozygosity rates"
../bin/plink --bfile ../data/hgdp/00_data/hgdp --het --out ../data/hgdp/01_QC/hgdp

# calculate heterozygosity rate per individual and save as pdf 
# investigate plot to see if any individuals need to be removed 
# -3 sd away from mean: 0.2278569
# +3 sd away from mean: 0.3365223
# mean = 0.2821896
# sd = 0.02716635
# add the family and individual ids of those who do not meet these thresholds: 
# one: exclude if genotype failure rate ≥ 0.03
# two: exclude if is not in +- 3 SD away from mean heterozygosity rate 
# 77 samples will be excluded: see fail-imisshet-qc.txt
echo "Saving heterozygosity rate as pdf"
echo "Outputting fail-imisshet-qc.txt"
R CMD BATCH ../src/qc/imiss-vs-het.Rscript ../data/hgdp/01_QC/hgdp

# identify duplicated or replicated individuals 
# skipping to keep high LD 
# plink --file  --exclude 00_data/high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out 

# generate pair-wise IBS identity by state - detect individuals who looked more different than 
# they would in a random, homogenous sample 
# did not include extract argument because skipped removing high LD 
echo "Identifying duplicated individuals"
../bin/plink --bfile ../data/hgdp/00_data/hgdp -genome --out ../data/hgdp/01_QC/hgdp
# nohup plink --noweb --bfile 00_data/hgdp -genome --out 01_QC/03_duplication/hgdp &


# identify all pairs of individuals with IBD > 0.185 to remove individual with lowest call rate 
echo "Identifying pairs of individuals with IBD > 0.185"
perl ../src/qc/run-IBD-QC.pl ../data/hgdp/01_QC/hgdp


# error: extracting all variants so nothing left in bim file????? SKIPPING STEPS 14-19 OF PAPER
# identify individuals of divergent ancestry 
# merge study group to 1000 genome phase 3 
# paper used hapmap phase 3 instead 
# trying to stay consistent with reference panel for imputation  
# make a new bed file that excludes SNPs not featured in 1000 genome phase 3 panel
# plink --bfile 00_data/hgdp -extract 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt --make-bed --out 01_QC/hgdp.1000-snps
# # error - no variants remaining after extract 
# plink --bfile 00_data/hgdp -extract 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt --make-bed --allow-no-vars --out 01_QC/hgdp.1000-snps_novar
# # merge and extract new target to reference panel 
# plink --bfile 01_QC/hgdp.1000-snps_novar --bmerge 1000.bed 1000.bim 1000.fam --make-bed --out 01_QC/hgdp_1000


# remove individuals failing QC 
echo "Combining failed files"
cat ../data/hgdp/01_QC/fail-* | sort -k1 | uniq > ../data/hgdp/01_QC/fail-qc-inds.txt
echo "Removing individuals failing QC"
../bin/plink --bfile ../data/hgdp/00_data/hgdp --remove ../data/hgdp/01_QC/fail-qc-inds.txt --make-bed --out ../data/hgdp/00_data/clean-inds-hgdp-data

# identify markers with excessive missing data rate 
echo "Identifying markers with excessiving missing data rate"
../bin/plink --bfile ../data/hgdp/00_data/clean-inds-hgdp-data --missing --out ../data/hgdp/01_QC/clean-inds-hgdp-data

# Plot a histogram of the missing genotype rate to identify a threshold for extreme
# genotype failure rate
# paper used call-rate threshold of 3% 
echo "Plotting histogram of missing genotype rate"
R CMD BATCH ../src/qc/lmiss-hist.Rscript ./data/hgdp/01_QC/hgdp-R-miss-geno

# SKIPPED: Test markers for different genotype call rates between cases and contols

# Remove all markers failing QC
# geno = filters out samples exceeding missing genotype rate of 3% 
# maf = filters out below minor allele frequency threshold 
# hwe = filters out below Hardy-Weinberg equilibrium exact test p-value
plink --bed ./data/hgdp/00_data/clean-inds-hgdp-data.bed --bim ./data/hgdp/00_data/clean-inds-hgdp-data.bim --fam ./data/hgdp/00_data/clean-inds-hgdp-data.fam -maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ./data/hgdp/00_data/clean-hgdp-data


