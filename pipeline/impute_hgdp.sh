# extract chr 1 data, phase, and impute 

# extract chromosome 1 data from lifted over HGDP data 
# make ped and map files including only chr 1 
../bin/plink --file ../data/hgdp/00_data/clean-hgdp-data.HG19 --chr 1 --recode --out ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1
# make bed, bim, fam files including only chr 1
../bin/plink --file ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1 --make-bed --out ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1

# check strand alignment before pre-phasing 
../bin/shapeit -check \
-P ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1 \
--input-ref ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz  ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log ../data/hgdp/00_data/unphased-clean-hgdp-data_HG19_chr_1.alignments

## make a list of sites to flip and exclude 
# grab the main id position in the 4th field, which has rs numbers 
# remove duplicates from final list 
cat ../data/hgdp/00_data/unphased-clean-hgdp-data_HG19_chr_1.alignments.snp.strand | grep "Strand" | cut -f 4 | uniq > ../data/hgdp/00_data/hgdp-data_HG19_chr_1_strand_flip.txt 
cat ../data/hgdp/00_data/unphased-clean-hgdp-data_HG19_chr_1.alignments.snp.strand | grep "Missing" | cut -f 4 | uniq > ../data/hgdp/00_data/hgdp-data_HG19_chr_1_missing.txt

## flip strand and exclude SNPs missing in reference panel using plink
## make new bed files 
plink \
--file ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1 \
--flip ../data/hgdp/00_data/hgdp-data_HG19_chr_1_strand_flip.txt  \
--exclude ../data/hgdp/00_data/00_data/hgdp-data_HG19_chr_1_missing.txt \
--recode \
--make-bed \
--out ../data/hgdp/00_data/aligned_unphased_hgdp-data_HG19_chr_1 



# extract Northern European Individuals 
cat ../data/hgdp/00_data/aligned_unphased_hgdp-data_HG19_chr_1.txt | grep -E "French_Basque|Sardinian|Russian|Sardinian|French|North_Italian|Tuscan" |cut -d ' ' -f 1-2 | uniq > ../data/hgdp/00_data/NE_keep.txt

plink --file ../data/hgdp/00_data/aligned_unphased_hgdp-data_HG19_chr_1 --keep ../data/hgdp/00_data/NE_keep.txt --recode --out ../data/hgdp/00_data/NE_unphased_hgdp_data_chr_1
plink --file ../data/hgdp/00_data/NE_unphased_hgdp_data_chr_1 --make-bed --out ../data/hgdp/00_data/NE_unphased_hgdp_data_chr_1
# phase using North European data using Shapeit
../bin/shapeit --input-ped ../data/hgdp/00_data/NE_unphased_hgdp_data_chr_1.ped ../data/hgdp/00_data/NE_unphased_hgdp_data_chr_1.map \
--input-map ../data/hgdp/00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
--output-max ../data/hgdp/00_data/NE_hgdp_data_chr_1.phased \
--thread 4 \
--output-log ../data/hgdp/00_data/NE_hgdp_data_chr_1.phased

# phase full HGDP using SHAPEIT
../bin/shapeit --input-ped ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1.ped ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1.map \
--input-map ../data/hgdp/00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
--output-max ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1.phased \
--thread 4 \
--output-log ../data/hgdp/00_data/log_clean-hgdp-data_HG19_chr_1.phased

# imputation with impute2: NE HGDP
impute2 \
-use_prephased_g \
-known_haps_g ../data/hgdp/00_data/22Nov2017/NE_hgdp_data_chr_1.phased.haps \
-m ../data/hgdp/00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-h ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
-l ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
-int 168549811 170549811\
-Ne 20000 \
-allow_large_regions \
-seed 367946 \
-phase \
-o ../data/hgdp/00_data/NE_imputed_hgdp.gen

# imputation with impute2: full HGDP
impute2 \
-use_prephased_g \
-known_haps_g ../data/hgdp/00_data/aligned_clean-hgdp-data_HG19_chr_1.phased.haps \
-m ../data/hgdp/00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-h ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
-l ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
-int 168549811 170549811\
-Ne 20000 \
-allow_large_regions \
-seed 367946 \
-phase \
-o 00_data/22Nov2017/mputed_hgdp.gen

