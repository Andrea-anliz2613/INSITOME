# extract chr 1 data, phase, and impute 

# extract chromosome 1 data from lifted over HGDP data 
# make ped and map files including only chr 1 
../bin/plink --file ../data/hgdp/00_data/clean-hgdp-data.HG19 --chr 1 --recode --out ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1
# make bed, bim, fam files including only chr 1
../bin/plink --file ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1 --make-bed --out ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1

# check strand alignment before pre-phasing 
# ../bin/shapeit -check -P chr${chr}.unphased --input-ref $REF_DIR/$REF_PREFIX$chr$REF_SUFFIX.hap.gz $REF_DIR$REF_PREFIX$chr$REF_SUFFIX.legend.gz  $REF_DIR$REF_PREFIX.sample --output-log $OUT_PREFIXchr${chr}.alignments

# phase using SHAPEIT
../bin/shapeit --input-ped ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1.ped ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1.map \
--input-map ../data/hgdp/00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
--output-max ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1.phased \
--thread 8 \
--output-log ../data/hgdp/00_data/log_clean-hgdp-data_HG19_chr_1.phased

# imputation with impute2 

impute2 \
-use_prephased_g \
-known_haps_g ../data/hgdp/00_data/clean-hgdp-data_HG19_chr_1.phased.haps \
-m ../data/hgdp/00_data/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
-h ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
-l ../data/hgdp/00_data/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
-int 168549811 170549811\
-Ne 20000 \
-allow_large_regions \
-seed 367946 \
-o ../data/hgdp/00_data/imputed_hgdp.impute2