#This is the codes for two-sample Mendelian randomization and colocalization analyses for 
#the causal effects of T2D-related pathways, insulin secretion and insulin sensitivity on birth weight.

library(TwoSampleMR)
library(LDlinkR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(coloc)

#load instrumental variables for exposures
T2DGGI_SNPs = read.csv("./T2DGGI_SNPs.csv", header = TRUE)  
Insulin_secretion_SNPs = read.csv("./Insulin_secretion_SNPs.csv", header = TRUE)  
Insulin_sensitivity_SNPs = read.csv("./Insulin_sensitivity_SNPs.csv", header = TRUE)  

#load outcome summary data
NG2019_BW_MA_EF_EUR<-fread("./Maternal_Effect_European_meta_NG2019.txt.gz", sep=' ',header = TRUE) #maternal genetically-determined birth weight
NG2019_BW_FE_EF_EUR<-fread("./Fetal_Effect_European_meta_NG2019.txt.gz", sep=' ',header = TRUE) #fetal genetically-determined birth weight

result_MR<-data.frame()
n=1
exposure_data<-subset(T2DGGI_SNPs,cluster=='obesity') #extract T2D pathway related to obesity
exposure_data<-exposure_data[,c('SNP','risk','Other','Beta_EUR','SE_EUR','Pvalue_EUR','EAF_EUR')] # extract necessary columns
colnames(exposure_data)<-c('SNP','effect_allele','other_allele','beta','se','pval','eaf') #rename column names
missing<-subset(exposure_data,!SNP %in% NG2019_BW_MA_EF_EUR$RSID) #missing SNPs in outcome summary data
m<-nrow(exposure_data)

#find proxy SNPs (LD R2>0.8) in outcome summary data for missing SNPs 
for (i in missing$SNP){
  b<-LDproxy(
    snp=i,
    pop = c("EUR"),
    r2d = "r2",
    token ='' , #fill in your token
    genome_build = "grch37",
  )
  if(nrow(b)==1){
    next()
  }
  data3<-merge(b,NG2019_BW_MA_EF_EUR,by.x="RS_Number",by.y='RSID',all=FALSE)
  proxy1<-subset(data3,R2>0.8)
  if(nrow(proxy1)==0){
    next()
  }
  proxy<-subset(proxy1,Distance==min(abs(proxy1$Distance)) | Distance==-min(abs(proxy1$Distance))) #find the closest SNP with highest R2 as the proxy SNP
  #add the proxy SNP into the exposure data
  proxy_snp<-proxy[1,'RS_Number']
  exposure_data[m+1,'SNP']=proxy_snp
  exposure_data[m+1,'effect_allele']=toupper(subset(NG2019_BW_MA_FE_EUR,RSID==proxy_snp)$ea)
  exposure_data[m+1,'other_allele']=toupper(subset(NG2019_BW_MA_FE_EUR,RSID==proxy_snp)$nea)
  exposure_data[m+1,'beta']=subset(exposure_data,SNP==i)$beta
  exposure_data[m+1,'se']=subset(exposure_data,SNP==i)$se
  exposure_data[m+1,'pval']=subset(exposure_data,SNP==i)$pval
  exposure_data[m+1,'eaf']=subset(exposure_data,SNP==i)$eaf
  m<-m+1
}

outcome_data<-dplyr::filter(NG2019_BW_MA_EF_EUR,RSID %in% exposure_data$SNP)
outcome_data<-outcome_data[,c('RSID','ea','nea','beta','se','p','eaf')]
colnames(outcome_data)<-c('SNP','effect_allele','other_allele','beta','se','pval','eaf') #change the column names of outcome data
outcome_data$effect_allele<-toupper(outcome_data$effect_allele)
outcome_data$other_allele<-toupper(outcome_data$other_allele)

exposure_data<-dplyr::filter(exposure_data,SNP %in% outcome_data$SNP)
exposure_data<-data.frame(exposure_data)

random_exp_dat <- format_data(exposure_data, type = "exposure")
random_out_dat <- format_data(outcome_data, type = "outcome")

dat<- harmonise_data(
  exposure_dat = random_exp_dat,
  outcome_dat = random_out_dat
)

dat$r2.exposure<-(get_r_from_bsen(dat$beta.exposure,dat$se.exposure,751754))^2
dat$F.exposure<-(751754-2)*dat$r2.exposure/(1-dat$r2.exposure)
dat<-subset(dat, F.exposure>10) # filter out weak instrumental variables

  result_MR[n,'exposure']<-'obesity'
  result_MR[n,'outcome']<-'NG2019_BW_MA_EF_EUR'
  res <- mr(dat) #MR results
  res<-generate_odds_ratios(res)
  #MR-Egger
  result_MR[n,'Egger_b']<-res[1,7]
  result_MR[n,'Egger_se']<-res[1,8]
  result_MR[n,'Egger_p']<-res[1,9]
  result_MR[n,'Egger_low']<-res[1,10]
  result_MR[n,'Egger_up']<-res[1,11]
  # MR-median-weighted
  result_MR[n,'WMedian_b']<-res[2,7]
  result_MR[n,'WMedian_se']<-res[2,8]
  result_MR[n,'WMedian_p']<-res[2,9]
  result_MR[n,'WMedian_low']<-res[2,10]
  result_MR[n,'WMedian_up']<-res[2,11]
  #MR-IVW
  result_MR[n,'IVW_b']<-res[3,7]
  result_MR[n,'IVW_se']<-res[3,8]
  result_MR[n,'IVW_p']<-res[3,9]
  result_MR[n,'IVW_low']<-res[3,10]
  result_MR[n,'IVW_up']<-res[3,11]
  #MR-weighted-mode
  result_MR[n,'WMode_b']<-res[5,7]
  result_MR[n,'WMode_se']<-res[5,8]
  result_MR[n,'WMode_p']<-res[5,9]
  result_MR[n,'WMode_low']<-res[5,10]
  result_MR[n,'WMode_up']<-res[5,11]
  #MR-PRESSO
  res_presso<-run_mr_presso(dat)
  result_MR[n,'presso_b']<-res_presso[[1]]$'Main MR results'[1,3]
  result_MR[n,'presso_se']<-res_presso[[1]]$'Main MR results'[1,4]
  result_MR[n,'presso_p']<-res_presso[[1]]$'Main MR results'[1,6]
  result_MR[n,'presso_low']<-res_presso[[1]]$'Main MR results'[1,3]-1.96*res_presso[[1]]$'Main MR results'[1,4]
  result_MR[n,'presso_up']<-res_presso[[1]]$'Main MR results'[1,3]+1.96*res_presso[[1]]$'Main MR results'[1,4]

  #heterogeneity test
  het<-mr_heterogeneity(dat)
  result_MR[n,'Egger_HET_Q']<-het[1,6]
  result_MR[n,'Egger_HET_Q_df']<-het[1,7]
  result_MR[n,'Egger_HET_pvalue']<-het[1,8]

  result_MR[n,'IVW_HET_Q']<-het[2,6]
  result_MR[n,'IVW_HET_Q_df']<-het[2,7]
  result_MR[n,'IVW_HET_pvalue']<-het[2,8]

  # horizontal pleiotropy test
  hplei<-mr_pleiotropy_test(dat)
  result_MR[n,'Egger_intercept']<-hplei[1,5]
  result_MR[n,'Egger_intercept_se']<-hplei[1,6]
  result_MR[n,'Egger_intercept_pvalue']<-hplei[1,7]

  result_MR[n,'presso_global_RSSobs']<-res_presso[[1]][["MR-PRESSO results"]][["Global Test"]][["RSSobs"]]
  result_MR[n,'presso_global_p']<-res_presso[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]

  #leave-one-out test
  res_loo <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(res_loo)

  #scatterplot
  res <- mr(dat)
  mr_scatter_plot(res, dat)

#colocalization analysis
T2DGGI_EUR<-fread("./EUR_Metal_LDSC-CORR_Neff.v2.txt", sep='\t', header = TRUE) 
CIR_EUR<-fread("./meta1cir.filthetmafn.rsid.selectedcolumns", sep='\t',header = TRUE)
CIR_EUR<-separate_wider_delim(CIR_EUR, cols = MarkerName, delim = ":", names = c("CHR", "BP"),cols_remove=F)
ISI_adj_EUR<-fread("./MAGIC_postchallengeIR_ISI_adjBMI_EUR.tsv.gz", sep='\t',header = TRUE)

SNPs<-read.table("clipboard",sep='\t',header=TRUE)
result_coloc<-data.frame()
for (i in rownames(SNPs)){
data1<-T2DGGI_EUR%>% filter(Chromsome==SNPs[i,'CHR'],Position >=SNPs[i,'BP']-500000,Position <= SNPs[i,'BP']+500000)
data2<-NG2019_BW_FE_EF_EUR%>% filter(chr==SNPs[i,'CHR'],pos >=SNPs[i,'BP']-500000,pos <= SNPs[i,'BP']+500000)
data1$chr_pos<-paste(data1$Chromsome,data1$Position,sep=':')
data = merge(data1,data2,by.x='chr_pos',by.y="chr_pos")
data = data[!duplicated(data$chr_pos),]
data$ea<-toupper(data$ea)
data$nea<-toupper(data$nea)
data = data %>% filter((EffectAllele==ea & NonEffectAllele==nea)|(EffectAllele==nea & NonEffectAllele==ea)) 
data=data%>%mutate(beta.y=ifelse(EffectAllele==ea,beta,-beta))
data$VAR.x = data$SE^2
data$VAR.y = data$se^2
data = data[data$VAR.x!=0 & data$VAR.y!=0 ,]
data$MAF.y <- ifelse(data$eaf<0.5,data$eaf,1-data$eaf)
data1 = data[,c("Effect","VAR.x","chr_pos")]
data2 = data[,c("beta.y","VAR.y","chr_pos",'MAF.y','n_offBW')]
colnames(data1)=c("beta","varbeta","snp")
colnames(data2)=c("beta","varbeta","snp",'MAF','N')
data1 = as.list(data1)
data2 = as.list(data2)
data1$type = "cc"
data2$type = "quant"
res = coloc.abf(data1,data2,p1=1e-4,p2=1e-4,p12=1e-5)
result_coloc[i,'SNP']<-SNPs[i,'SNP']
result_coloc[i,'PP.H0']<-res$summary[2]
result_coloc[i,'PP.H1']<-res$summary[3]
result_coloc[i,'PP.H2']<-res$summary[4]
result_coloc[i,'PP.H3']<-res$summary[5]
result_coloc[i,'PP.H4']<-res$summary[6]
}

  



