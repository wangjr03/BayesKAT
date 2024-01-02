# Dataset in ped format can be converted to a binary ped file (*.bed) using this simple line of code using PLINK software in https://zzz.bwh.harvard.edu/plink/data.shtml#bed : 
# plink --file mydata --make-bed
# This line of code creates three files:
#     mydata.bed      ( binary file, genotype information )
#     mydata.fam      ( first six columns of mydata.ped ) 
#     mydata.bim      ( extended MAP file: two extra cols = allele names)
#
# Quality control step using PLINK:
# Output are afterQC_mydata.bed, afterQC_mydata.fam, afterQC_mydata.bim
# plink --bfile mydata --geno 0.1 --mind 0.1 --maf 0.05 --hwe 0.0000000001 --allow-no-sex --nonfounders --make-bed --out afterQC_mydata
#
#
# R code to convert these files into a genotype matrix data 
library("genio")
data1=read_bim(file="afterQC_mydata.bim",verbose=TRUE)
data2=read_fam(file="afterQC_mydata.fam",verbose=TRUE)
geno_data=read_bed("afterQC_mydata",names_loci = data1$id,
       names_ind = data2$id,
       ext = "bed",
       verbose = TRUE
)
#
# If geno_data matrix contains missing values, those can be imputed using this code
#
fn_impute<-function(x){
	pos=which(is.na(x)==1)
	prob=table(x[-pos])/length(x[-pos])
	if(length(pos)!=0){
		for(i in pos){
			x[i]=which(rmultinom(1,1,prob)==1)-1		
		}
	}
	return(x)
}
#set.seed(2570)
geno_data=t(apply(geno_data,1,fn_impute))
#
#
# The function fn_Z creates input genotype matrix Z for BayesKAT for a given set of SNPs.
# @param geno_data is the genotype matrix where rows are SNPs and columns are individuals
# @param SNP_set is the array of rsids of SNPs of interest.
# @return a genotype matrix Z with rows as individuals and columns as SNPs.

fn_Z<-function(geno_data, SNP_set){
    snp_id=rownames(geno_data)
    ind_id=colnames(geno_data)
    Z=data[which(snp_id %in% SNP_set),]
    return(t(Z))
  }
input_Z=fn_Z(geno_data=geno_data, SNP_set=SNP_set)


