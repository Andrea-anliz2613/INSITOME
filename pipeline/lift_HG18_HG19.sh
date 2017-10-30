# LiftOver 
# convert files between build versions
# wd: INSITOME/pipeline 

# download chain files if you have not 
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz -O ../liftover/data
# gunzip ../liftover/hg18ToHg19.over.chain.gz 

# create map and ped files 
plink --bfile ../data/hgdp/00_data/clean-hgdp-data --recode --out ../data/hgdp/00_data/clean-hgdp-data

# create BED file based on map file 
gawk '{print "chr"$1, $4, $4+1, $2}' OFS="\t" ../data/hgdp/00_data/clean-hgdp-data.map > ../data/hgdp/00_data/clean-hgdp-data_HG18.BED

# perform liftover 
# left out bedPlus = 4 argument bc would hve treated first 4 cols as BED; unsure if appropriate  
./liftOver ../data/hgdp/00_data/clean-hgdp-data_HG18.BED ../liftover/hg18ToHg19.over.chain ../data/hgdp/00_data/clean-hgdp-data.HG19.BED ../data/hgdp/00_data/clean-hgdp-data_unmapped.txt

# create list of unmapped SNPs 
gawk '/^[^#]/ {print $4}' ../data/hgdp/00_data/clean-hgdp-data_unmapped.txt > ../data/hgdp/00_data/clean-hgdp-data_unmappedSNPs.txt

# create map file using new bed file 
gawk '{print $4, $2}' OFS="\t" ../data/hgdp/00_data/clean-hgdp-data.HG19.BED > ../data/hgdp/00_data/clean-hgdp-data.HG19.mapping.txt

# remove unmapped SNP from target data set 
# create new ped and map files 
plink --file ../data/hgdp/00_data/clean-hgdp-data --exclude ../data/hgdp/00_data/clean-hgdp-data_unmappedSNPs.txt --update-map ../data/hgdp/00_data/clean-hgdp-data.HG19.mapping.txt --make-bed --out ../data/hgdp/00_data/clean-hgdp-data.HG19.temp
plink --bfile ../data/hgdp/00_data/clean-hgdp-data.HG19.temp --recode --out ../data/hgdp/00_data/clean-hgdp-data.HG19

# create new snp list for data set 
gawk '{print $2}' ../data/hgdp/00_data/clean-hgdp-data.HG19.map > ../data/hgdp/00_data/clean-hgdp-data.HG19.snplist