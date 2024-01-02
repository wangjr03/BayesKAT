# Dataset in ped format can be converted to a binary ped file (*.bed) using this simple line of code using PLINK software in https://zzz.bwh.harvard.edu/plink/data.shtml#bed : 
# plink --file mydata --make-bed
# This line of code creates three files:
#     plink.bed      ( binary file, genotype information )
#     plink.fam      ( first six columns of mydata.ped ) 
#     plink.bim      ( extended MAP file: two extra cols = allele names)
#
# R code to convert these files into a genotype matrix data 
#

library("genio")
data1=read_bim(file="plink.bim",verbose=TRUE)
data2=read_fam(file="plink.fam",verbose=TRUE)
geno_data=read_bed("plink",names_loci = data1$id,
       names_ind = data2$id,
       ext = "bed",
       verbose = TRUE
)

# The function fn_Z can be used for creating input genotype matrix Z 
#for BayesKAT for a given set of SNPs.
# @param geno_data is the genotype matrix where rows are SNPs and columns are individuals
# @param SNP_set is the array of rsids of SNPs of interest.

fn_Z<-function(geno_data, SNP_set){
    snp_id=rownames(geno_data)
    ind_id=colnames(geno_data)
    Z=data[which(snp_id %in% SNP_set),]
    return(Z)
  }

