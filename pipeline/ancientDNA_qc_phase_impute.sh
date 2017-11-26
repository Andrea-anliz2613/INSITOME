
ROOT_DIR=./
INSITOME_DATA_DIR=00_data/insitomedata/
QC_DIR=${ROOT_DIR}01_QC/ancientDNA_qc/
PHASE_DIR=00_data/ancientDNA_impute/phased/
IMPUTE_DIR=00_data/ancientDNA_impute/imputeresults/
ALIGNMENTS_DIR=00_data/ancientDNA_impute/alignments/

 
# executable
# PLINK_EXEC=${ROOT_DIR}bin/plink
# SHAPEIT_EXEC=${ROOT_DIR}bin/shapeit
 
# reference data files
GENMAP_FILE=${ROOT_DIR}00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt
 
# insitome files 
list_files=("Europe" "FuQ0.9" "FuQ2.2M" "Haak" "Mathfull" "ScythianSarmatian" "Sko" "StMarys") 


# remove chromosome 90 information from Europe DNA


for i in "${list_files[@]}"; do

	# declare bed, bim, fam names 
	GWASDATA_BED=${INSITOME_DATA_DIR}$i.bed
	GWASDATA_BIM=${INSITOME_DATA_DIR}$i.bim
	GWASDATA_FAM=${INSITOME_DATA_DIR}$i.fam

	# echo $GWASDATA_BED
	# echo $GWASDATA_BIM
	# echo $GWASDATA_FAM

	# rename chromosomes x, y, xy, and mitochondrial to 23, 24, 25, and 26 
	sed -i '' -- 's/X/23/g' ${GWASDATA_BIM}
	sed -i '' -- 's/Y/24/g' ${GWASDATA_BIM}
	sed -i '' -- 's/XY/25/g' ${GWASDATA_BIM}
	sed -i '' -- 's/90/26/g' ${GWASDATA_BIM}

	# echo ${GWASDATA_BIM}

	# STEP 1: QUALITY CONTROL

	# discordant sex information 
	plink --bed ${GWASDATA_BED} --bim ${GWASDATA_BIM} --fam ${GWASDATA_FAM} --check-sex --out ${QC_DIR}$i

	## if [ -f ${QC_DIR}$i.sexcheck ]; then 
	## 	grep "PROBLEM" ${QC_DIR}$i.sexcheck > ${QC_DIR}$i_problem
	## fi
	## Haak and FuQ0.9only has one gender, so was excluded 
	## StMarys didn't have gender issues 

	grep "PROBLEM" 01_QC/ancientDNA_qc/Europe.sexcheck > 01_QC/ancientDNA_qc/Europe_problem
	grep "PROBLEM" 01_QC/ancientDNA_qc/FuQ2.2M.sexcheck > 01_QC/ancientDNA_qc/FuQ2.2M_problem
	grep "PROBLEM" 01_QC/ancientDNA_qc/Mathfull.sexcheck > 01_QC/ancientDNA_qc/Mathfull_problem
	grep "PROBLEM" 01_QC/ancientDNA_qc/ScythianSarmatian.sexcheck > 01_QC/ancientDNA_qc/ScythianSarmatian_problem
	grep "PROBLEM" 01_QC/ancientDNA_qc/Sko.sexcheck > 01_QC/ancientDNA_qc/Sko_problem


	## Identification of individuals with elevated missing data rates or outlying heterozygosity rate 
	## find missing SNP; output = hgdp.imiss and hgdp.lmiss file
	plink --bfile ${INSITOME_DATA_DIR}$i --missing --out ${QC_DIR}$i
	# find heterozygous genotypes; output = hgdp.het file
	plink --bfile ${INSITOME_DATA_DIR}$i --het --out ${QC_DIR}$i
	# calculate heterozygosity rate per individual and save as pdf 
	# investigate plot to see if any individuals need to be removed 
	# did not loop in R file 
	R CMD BATCH script/imiss-vs-het.Rscript ${QC_DIR}$i

	# identify duplicated or replicated individuals 
	# skipping to keep high LD 
	plink --file  --exclude 00_data/high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out 

	## generate pair-wise IBS identity by state - detect individuals who looked more different than 
	## they would in a random, homogenous sample 
	## did not include extract argument because skipped removing high LD 
	plink --bfile ${INSITOME_DATA_DIR}$i -genome --out ${QC_DIR}$i
	# identify all pairs of individuals with IBD > 0.185 to remove individual with lowest call rate 
	perl script/run-IBD-QC-ancientDNA.pl ${QC_DIR}$i

	# skip steps 14-19 of paper

	# note: couldnt figure out how to rename the files appropriately in perl so must rename now for next steps to work
	
	# remove individuals failing QC 
	cat ${QC_DIR}fail-*-$i.txt | sort -k1 | uniq > ${QC_DIR}fail-$i-qc-inds.txt

	# ISSUE: Error: Line 1 of --remove file has fewer tokens than expected.
	# proceeding with only Sko bed format files 

	plink --bfile --bed ${GWASDATA_BED} --bim ${GWASDATA_BIM} --fam ${GWASDATA_FAM} \
	--no-fid --no-parents --no-sex --no-pheno \
	--remove ${QC_DIR}fail-$i-qc-inds.txt \
	--make-bed \
	--out ${QC_DIR}clean-inds-$i-data


done


# Sko bed format files 
# identify markers with excessive missing data rate 
plink --bfile ${QC_DIR}clean-inds-Sko-data --missing --out ${QC_DIR}clean-inds-Sko-data-excess-missing

# Plot a histogram of the missing genotype rate to identify a threshold for extreme
# genotype failure rate
# paper used call-rate threshold of 3% 
R CMD BATCH script/lmiss-hist.Rscript ${QC_DIR}Sko-R-miss-geno


# SKIPPED: Test markers for different genotype call rates between cases and contols

# Remove all markers failing QC
# geno = filters out samples exceeding missing genotype rate of 3% 
# maf = filters out below minor allele frequency threshold 
# hwe = filters out below Hardy-Weinberg equilibrium exact test p-value

#Error: All variants excluded due to missing genotype data (--geno).
# ERROR FOUND IN Mathfull, FUQ2.2M and FUQ0.9
# 17 Nov 25: Will phase the raw files for those three above 


plink --bed ${QC_DIR}clean-inds-Sko-data.bed --bim ${QC_DIR}clean-inds-Sko-data.bim --fam ${QC_DIR}clean-inds-Sko-data.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ${QC_DIR}clean-Sko-data

plink --bed ${INSITOME_DATA_DIR}StMarys.bed --bim ${INSITOME_DATA_DIR}StMarys.bim --fam ${INSITOME_DATA_DIR}StMarys.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ${QC_DIR}clean_QC_marker_StMarys

plink --bed ${INSITOME_DATA_DIR}Europe.bed --bim ${INSITOME_DATA_DIR}Europe.bim --fam ${INSITOME_DATA_DIR}Europe.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ${QC_DIR}clean_QC_marker_Europe

plink --bed ${INSITOME_DATA_DIR}FuQ0.9.bed --bim ${INSITOME_DATA_DIR}FuQ0.9.bim --fam ${INSITOME_DATA_DIR}FuQ0.9.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ${QC_DIR}clean_QC_marker_FuQ0.9

plink --bed ${INSITOME_DATA_DIR}FuQ2.2M.bed --bim ${INSITOME_DATA_DIR}FuQ2.2M.bim --fam ${INSITOME_DATA_DIR}FuQ2.2M.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ${QC_DIR}clean_QC_marker_FuQ2.2M

plink --bed ${INSITOME_DATA_DIR}ScythianSarmatian.bed --bim ${INSITOME_DATA_DIR}ScythianSarmatian.bim --fam ${INSITOME_DATA_DIR}ScythianSarmatian.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ${QC_DIR}clean_QC_marker_ScythianSarmatian

plink --bed ${INSITOME_DATA_DIR}Haak.bed --bim ${INSITOME_DATA_DIR}Haak.bim --fam ${INSITOME_DATA_DIR}Haak.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed --out ${QC_DIR}clean_QC_marker_Haak

plink --bed ${INSITOME_DATA_DIR}Mathfull.bed --bim ${INSITOME_DATA_DIR}Mathfull.bim --fam ${INSITOME_DATA_DIR}Mathfull.fam \
-maf 0.01 --geno 0.03 --hwe 0.00001 --make-bed \
--out ${QC_DIR}clean_QC_marker_Mathfull


# no lifting over necessary for sko, stmarys, FUQ, haak, and europe
# unknown versions for ScythianSarmatian and mathfull

# phasing 

# extract chromosome 1 data 
# make ped and map files including only chr 1 
plink --bfile ${QC_DIR}clean_QC_marker_Europe  --chr 1 --recode --out ${INSITOME_DATA_DIR}clean_QC_marker_Europe_chr1

# make bed, bim, fam files including only chr 1
plink --file ${INSITOME_DATA_DIR}clean_QC_marker_Europe_chr1 --make-bed --out ${INSITOME_DATA_DIR}clean_QC_marker_Europe_chr1

plink --bfile ${QC_DIR}clean-Sko-data --chr 1 --recode --out ${INSITOME_DATA_DIR}clean-Sko-data_chr1
plink --file ${INSITOME_DATA_DIR}clean-Sko-data_chr1 --make-bed --out ${INSITOME_DATA_DIR}clean-Sko-data_chr1

plink --bfile ${QC_DIR}clean_QC_marker_Haak --chr 1 --recode --out ${INSITOME_DATA_DIR}clean_QC_marker_Haak_chr1
plink --file ${INSITOME_DATA_DIR}clean_QC_marker_Haak_chr1 --make-bed --out ${INSITOME_DATA_DIR}clean_QC_marker_Haak_chr1

plink --bfile ${QC_DIR}clean_QC_marker_ScythianSarmatian --chr 1 --recode --out ${INSITOME_DATA_DIR}clean_QC_marker_ScythianSarmatian_chr1
plink --file ${INSITOME_DATA_DIR}clean_QC_marker_ScythianSarmatian_chr1 --make-bed --out ${INSITOME_DATA_DIR}clean_QC_marker_ScythianSarmatian_chr1

plink --bfile ${QC_DIR}clean_QC_marker_StMarys --chr 1 --recode --out ${INSITOME_DATA_DIR}clean_QC_marker_StMarys_chr1
plink --file ${INSITOME_DATA_DIR}clean_QC_marker_StMarys_chr1 --make-bed --out ${INSITOME_DATA_DIR}clean_QC_marker_StMarys_chr1

# check strand alignment before pre-phasing 
# -P = ped format 
shapeit \
-check \
-P ${INSITOME_DATA_DIR}clean_QC_marker_Europe_chr1 \
--input-ref 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz  00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log ${ALIGNMENTS_DIR}unphased-clean_Europe_chr1.alignments 

shapeit \
-check \
-P ${INSITOME_DATA_DIR}clean_QC_marker_Haak_chr1 \
--input-ref 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz  00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log ${ALIGNMENTS_DIR}unphased-clean_Haak_chr1.alignments 

shapeit \
-check \
-P ${INSITOME_DATA_DIR}clean_QC_marker_ScythianSarmatian_chr1 \
--input-ref 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz  00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log ${ALIGNMENTS_DIR}unphased-clean_ScythianSarmatian.alignments 

shapeit \
-check \
-P ${INSITOME_DATA_DIR}clean_QC_marker_StMarys_chr1 \
--input-ref 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz  00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log ${ALIGNMENTS_DIR}unphased-clean_StMarys_chr1.alignments 

shapeit \
-check \
-P ${INSITOME_DATA_DIR}clean-Sko-data_chr1 \
--input-ref 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz  00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log ${ALIGNMENTS_DIR}unphased-clean_Sko_chr1.alignments 

# "Europe" 
# new_list_files=("Haak" "ScythianSarmatian" "Sko" "StMarys") 
# for i in "${new_list_files[@]}"; do

# 	# declare bed, bim, fam names 
# 	GWASDATA_BED=${INSITOME_DATA_DIR}$i.bed
# 	GWASDATA_BIM=${INSITOME_DATA_DIR}$i.bim
# 	GWASDATA_FAM=${INSITOME_DATA_DIR}$i.fam

# 	# echo $GWASDATA_BED
# 	# echo $GWASDATA_BIM
# 	# echo $GWASDATA_FAM

# 	shapeit \
# 	-check \
# 	-P ${INSITOME_DATA_DIR}clean*$i*_chr1 \
# 	--input-ref 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz  00_data/1000GP_Phase3/1000GP_Phase3.sample \
# 	--output-log ${ALIGNMENTS_DIR}unphased-clean_$i_chr1.alignments 

# done




## make a list of sites to flip and exclude 
# grab the main id position in the 3rd field 
# remove duplicates from final list 
cat ${ALIGNMENTS_DIR}unphased-clean_Europe_chr1.alignments.snp.strand | grep "Strand" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_Europe_chr_1_strand_flip.txt
cat ${ALIGNMENTS_DIR}unphased-clean_Europe_chr1.alignments.snp.strand | grep "Missing" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_Europe_chr_1_missing.txt

cat ${ALIGNMENTS_DIR}unphased-clean_ScythianSarmatian_chr1.alignments.snp.strand | grep "Strand" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_ScythianSarmatian_chr_1_strand_flip.txt
cat ${ALIGNMENTS_DIR}unphased-clean_ScythianSarmatian_chr1.alignments.snp.strand | grep "Missing" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_ScythianSarmatian_chr_1_missing.txt

cat ${ALIGNMENTS_DIR}unphased-clean_Haak_chr1.alignments.snp.strand | grep "Strand" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_Haak_chr_1_strand_flip.txt
cat ${ALIGNMENTS_DIR}unphased-clean_Haak_chr1.alignments.snp.strand | grep "Missing" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_Haak_chr_1_missing.txt

cat ${ALIGNMENTS_DIR}unphased-clean_StMarys_chr1.alignments.snp.strand | grep "Strand" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_StMarys_chr_1_strand_flip.txt
cat ${ALIGNMENTS_DIR}unphased-clean_StMarys_chr1.alignments.snp.strand | grep "Missing" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_StMarys_chr_1_missing.txt

cat ${ALIGNMENTS_DIR}unphased-clean_Sko_chr1.alignments.snp.strand | grep "Strand" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_Sko_chr_1_strand_flip.txt
cat ${ALIGNMENTS_DIR}unphased-clean_Sko_chr1.alignments.snp.strand | grep "Missing" | cut -f 4 | uniq > ${ALIGNMENTS_DIR}unphased-clean_Sko_chr_1_missing.txt

## flip strand and exclude SNPs missing in reference panel using plink
## make new bed files 

## Europe: ERROR 
## Duplicate id error in map file: Error: Duplicate ID 'Aff23-50009726'
## tried removing duplicate; ERROR: Error: Line 1 of --remove file has fewer tokens than expected.
plink \
--file ${INSITOME_DATA_DIR}clean_QC_marker_Europe_chr1 \
--remove ${ALIGNMENTS_DIR}remove_indiv_Europe.txt \
--recode \
--out ${ALIGNMENTS_DIR}no_dup_unphased_Europe_chr_1 

plink \
--file ${INSITOME_DATA_DIR}clean_QC_marker_Europe_chr1 \
--flip ${ALIGNMENTS_DIR}unphased-clean_Europe_chr_1_strand_flip.txt \
--exclude ${ALIGNMENTS_DIR}unphased-clean_Europe_chr_1_missing.txt \
--recode \
--out ${ALIGNMENTS_DIR}aligned_unphased_Europe_chr_1 

# Scythian 
plink \
--file ${INSITOME_DATA_DIR}clean_QC_marker_ScythianSarmatian_chr1 \
--exclude ${ALIGNMENTS_DIR}unphased-clean_ScythianSarmatian_chr_1_missing.txt \
--recode \
--out ${ALIGNMENTS_DIR}aligned_unphased_ScythianSarmatian_chr_1 

# Haak 
# plink \
# --file ${INSITOME_DATA_DIR}clean_QC_marker_Haak_chr1 \
# --flip ${ALIGNMENTS_DIR}unphased-clean_Haak_chr_1_strand_flip.txt \
# --exclude ${ALIGNMENTS_DIR}unphased-clean_Haak_chr_1_missing.txt \
# --recode \
# --out ${ALIGNMENTS_DIR}aligned_unphased_Haak_chr_1 
plink \
--file ${INSITOME_DATA_DIR}clean_QC_marker_Haak_chr1 \
--exclude ${ALIGNMENTS_DIR}combined_unphased-clean_Haak_chr_1_missing.txt \
--recode \
--out ${ALIGNMENTS_DIR}aligned_unphased_Haak_chr_1 

# StMarys ERROR 
# Same error as Europe 
# Error: Line 1 of --remove file has fewer tokens than expected.
plink \
--file ${INSITOME_DATA_DIR}clean_QC_marker_StMarys_chr1 \
--remove ${ALIGNMENTS_DIR}remove_indiv_StMarys.txt \
--recode \
--out ${ALIGNMENTS_DIR}no_dup_unphased_StMarys_chr_1 
# Error: Duplicate ID 'rs10737260'.
plink \
--file ${INSITOME_DATA_DIR}clean_QC_marker_StMarys_chr1 \
--flip ${ALIGNMENTS_DIR}unphased-clean_StMarys_chr_1_strand_flip.txt \
--exclude ${ALIGNMENTS_DIR}unphased-clean_StMarys_chr_1_missing.txt \
--recode \
--out ${ALIGNMENTS_DIR}aligned_unphased_StMarys_chr_1 

# Sko 
# plink \
# --file ${INSITOME_DATA_DIR}clean-Sko-data_chr1 \
# --flip ${ALIGNMENTS_DIR}unphased-clean_Sko_chr_1_strand_flip.txt \
# --exclude ${ALIGNMENTS_DIR}unphased-clean_Sko_chr_1_missing.txt \
# --recode \
# --out ${ALIGNMENTS_DIR}aligned_unphased_Sko_chr_1 
plink \
--file ${INSITOME_DATA_DIR}clean-Sko-data_chr1 \
--flip ${ALIGNMENTS_DIR}unphased-clean_Sko_chr_1_strand_flip.txt \
--exclude ${ALIGNMENTS_DIR}combined_unphased-clean_Sko_chr_1_missing.txt \
--recode \
--out ${ALIGNMENTS_DIR}aligned_unphased_Sko_chr_1 


# phase using SHAPEIT
## ScythianSarmatian
shapeit \
--input-ped ${ALIGNMENTS_DIR}aligned_unphased_ScythianSarmatian_chr_1.ped ${ALIGNMENTS_DIR}aligned_unphased_ScythianSarmatian_chr_1.map \
--input-map 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-R 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz 00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-max ${PHASE_DIR}aligned_ScythianSarmatian_chr_1.phased \
--thread 4 \
--output-log ${PHASE_DIR}aligned_ScythianSarmatian_chr_1_LOG.phased

## Haak 
shapeit \
--input-ped ${ALIGNMENTS_DIR}aligned_unphased_Haak_chr_1.ped ${ALIGNMENTS_DIR}aligned_unphased_Haak_chr_1.map \
--input-map 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-R 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz 00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-max ${PHASE_DIR}aligned_Haak_chr_1.phased \
--thread 4 \
--output-log ${PHASE_DIR}aligned_Haak_chr_1_LOG.phased

# ERROR
# ERROR: Reference and Main panels are not well aligned:
#   * #Missing sites in reference panel = 0
#   * #Misaligned sites between panels = 13
#   * #Multiple alignments between panels = 0
# will include the flipped strands into the mssing file and remove; file combined_unphased-clean_Haak_chr_1_missing.txt
# Aff23-15487538
# Aff23-10183760
# Aff23-10986945
# Aff23-6298600
# Aff23-5505842
# Aff23-9585960
# Aff23-13063949
# Aff23-16273477
# Aff23-50008367
# Aff23-5274155
# Aff23-5372333
# Aff23-5527576
# Aff23-5796562
# Aff23-6575511
# Aff23-7319400
# Aff23-7322716
# Aff23-7406976
# Aff23-7507973
# Aff23-50020381


## Sko
## ERROR: Reference and Main panels are not well aligned: combined_unphased-clean_Sko_chr_1_missing.txt
  # * #Missing sites in reference panel = 0
  # * #Misaligned sites between panels = 9
  # * #Multiple alignments between panels = 0
  # ADDED the following to new missing file to remove 
#   rs476527
# rs4655836
# rs7532826
# rs10912657
# rs12031354
# rs1795244
# rs10798410
# rs2609351
# rs12405469
# rs2502408
# rs637280
# rs570661

shapeit \
--input-ped ${ALIGNMENTS_DIR}aligned_unphased_Sko_chr_1.ped ${ALIGNMENTS_DIR}aligned_unphased_Sko_chr_1.map \
--input-map 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-R 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz 00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-max ${PHASE_DIR}aligned_Sko_chr_1.phased \
--thread 4 \
--output-log ${PHASE_DIR}aligned_Sko_chr_1_LOG.phased


# Imputation

## Sko
impute2 \
-use_prephased_g \
-known_haps_g ${PHASE_DIR}aligned_Sko_chr_1.phased.haps \
-m 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-h 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
-l 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
-int 168549811 170549811\
-Ne 20000 \
-allow_large_regions \
-seed 367946 \
-phase \
-o ${IMPUTE_DIR}Sko_imputed.gen

## Haak
impute2 \
-use_prephased_g \
-known_haps_g ${PHASE_DIR}aligned_Haak_chr_1.phased.haps \
-m 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-h 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
-l 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
-int 168549811 170549811\
-Ne 20000 \
-allow_large_regions \
-seed 367946 \
-phase \
-o ${IMPUTE_DIR}Haak_imputed.gen

## ScythianSarmatian
impute2 \
-use_prephased_g \
-known_haps_g ${PHASE_DIR}aligned_ScythianSarmatian_chr_1.phased.haps \
-m 00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-h 00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
-l 00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
-int 168549811 170549811\
-Ne 20000 \
-allow_large_regions \
-seed 367946 \
-phase \
-o ${IMPUTE_DIR}ScythianSarmatian_imputed.gen


