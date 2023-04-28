# About ADNI data

We have used the individual-level genotype, phenotype, and demographic covariates data from The ADNI project available on this website https://adni.loni.usc.edu/ to conduct some real data analysis. More information about the ADNI datasets and the data generation details are described in paper by  Saykin AJ et al. The phenotype and covariate data are also available in the ADNIMERGE R package, which can be installed from loni website. We choose the whole brain volume variable as trait of interest and Age, gender, and occupation data as covariates. We use Individual level genotype data to test if any given group of SNPs are associated with the trait. ADNI website has information on over 800 individuals, but we chose 755 individuals for our study as they have information on all the required variables.  

We choose different biologically meaningful and functionally similar SNP groups to perform genome-wide testing. We conduct 
(1)Gene-wise analysis for 18999 protein-coding genes. We simultaneously test each SNP group located within 5KB downstream or upstream of each gene.
(2)Pathway-wise analysis considering all 352 KEGG pathways. For each pathway, we find the genes and test each SNP group located within 5KB downstream or upstream of the genes.We perform this testing simultaneously for each KEGG pathway. 

For each of these analysis, we have made available the final results. The Ranked_kegg_pathways.txt file and Ranked_genes.txt file contains the ranked group and their corresponding posterior probability. Higher posterior probability indicates higher chances for the corresponding SNP group to be associated with the trait. 

### reference:
1. Saykin AJ, Shen L, Foroud TM, Potkin SG, Swaminathan S, Kim S, Risacher SL, Nho K, Huentelman MJ, Craig DW, Thompson PM, Stein JL, Moore JH, Farrer LA, Green RC, Bertram L, Jack CR Jr, Weiner MW. (2010)
Alzheimer's Disease Neuroimaging Initiative biomarkers as quantitative phenotypes: Genetics core aims, progress, and plans.

# Remark:
Note that, all the results for KEGE pathways and genes and described here are based on ADNI dataset only. There might be other SNPs(and genes or pathways) which might be strongly related to the trait but not genotyped in the ADNI dataset. 
