library("readxl")

##################################################################################################################
################### Upload the dataset: sheet 8 is volume-Batch normalized imputed dataset provided by the Metabolon
################### Sheet5 is normalized non-imputed dataset
##################################################################################################################
meta <- read_excel("/home/sahmad/metabolomics/gut_liver_brain/DUKE-0304-19ML+_CLIENT DATA TABLE.XLSX",sheet = 8)
metadata<-as.data.frame(meta)
metadata<-metadata[,-c(1,2,4)] 
meta2 <- read_excel("/home/sahmad/metabolomics/gut_liver_brain/DUKE-0304-19ML+_CLIENT DATA TABLE.XLSX",sheet = 5)
metadata2<-as.data.frame(meta2)
metadata2<-metadata2[,-c(1,2,4)]
### placing NA values back into the normalized data
for(i in 1:dim(metadata2)[1]){
for(j in 1:dim(metadata2)[2]){
metadata[i,j][which(is.na(metadata2[i,j]))]<-NA
}
}
rownames(metadata)<-metadata$CLIENT_IDENTIFIER
metadata<-metadata[,-1]

################### save(metadata,file="Metabolomics_data_after_NAs_in.Rdata")
################### load("Metabolomics_data_after_NAs_in.Rdata")
##################################################################################################################
################### Extract Rotterdam study sample ids: subjects have ids (starting from 0-9, so extracted ids only), N=1082, metabolites = 1387
##################################################################################################################
all_ids<-rownames(metadata)
participants<-grep("^[1-9]",all_ids,value=TRUE)
data<-metadata[participants,]


################################################################################################################## STEP 1
################# STEP1.1 ==> Remove samples with missingness more than 5xSD of the mean  
################# 			  mean(missingInd*100) [1] 19.06627 > sd(missingInd*100) [1] 7.073257
#################             Missing 5 times SD 54.43%, 0.55
################# STEP1.2==> Save the metabolites with more than 30% missing as separate dataframe with binomial coding
################# STEP1.3==> Remove missing metabolites more than 70% percent to prepare dataset for next step
##################################################################################################################

############ STEP 1.1 [ 14 individuals are removed, so now N= 1068 ]
missingInd<-rowMeans(is.na(data))
good_ind<-missingInd[which(missingInd<=0.55)]
ind_ids<-rownames(as.data.frame(good_ind))
data_participants<-data[ind_ids,]

########### STEP 1.2 [ 299 metabolites are included in the extra data ]
missing_metabolites_1.2<-colMeans(is.na(data_participants))
list_missing_30<-missing_metabolites_1.2[which(missing_metabolites_1.2>0.30)]
metaname_30<-rownames(as.data.frame(list_missing_30))
data_participants_for_binary_data<-data_participants[,metaname_30]
data_participants_for_binary_data[!is.na(data_participants_for_binary_data)]<-1
data_participants_for_binary_data[is.na(data_participants_for_binary_data)]<-0
save(data_participants_for_binary_data,file="Metabolomics_data_greater_than_30percent_missing_binary.Rdata")

############ STEP 1.3 [ 140 metabolites are removed so metabolites = 1247 out of 1387 ]
missing_metabolites<-colMeans(is.na(data_participants))
good_metab<-missing_metabolites[which(missing_metabolites<=0.70)]
metaname<-rownames(as.data.frame(good_metab))
data_participants_V2<-data_participants[,metaname]



################################################################################################################## STEP 2
################# STEP2.1 ==> Extract NIST samples            
################# STEP2.2 ==> Missingness in NIST samples (if extreme outliers remove)- only 1 is removed due to 5X SD of mean missingness
#################			  mean(t$missing_nist_samples), #[1] 0.1569822, sd(t$missing_nist_samples) [1] 0.05152471
#################             0.05152471*5 = [1] 0.2576236  -- mean + 5XSD= [1] 0.4146058
################# STEP2.3 ==> Missingness in metabolites, remove greater than 30%
################# STEP2.4 ==> calculate the CV in NIST samples and flag metabolites
################# STEP2.5 ==> Remove metabolites with >30CV from Rotterdam participant data generated in 1.3 step 
##################################################################################################################

########## STEP 2.1
nistsamples<-grep("NIST",all_ids,value=TRUE)
nistdata<-metadata[nistsamples,]

######### STEP 2.2
missing_nist_samples<-rowMeans(is.na(nistdata))
df_missing_nist<-as.data.frame(missing_nist_samples)
good_NIST_samples<-missing_nist_samples[which(missing_nist_samples<=0.4146058)]
good_NIST_ids<-rownames(as.data.frame(good_NIST_samples))
nistdata2<-nistdata[good_NIST_ids,]

######## STEP 2.3
missing_nist<-colMeans(is.na(nistdata2))
good_metab_nist<-missing_nist[which(missing_nist<=0.30)]
good_metab_nist<-rownames(as.data.frame(good_metab_nist))
nistdata_qc<-nistdata2[,good_metab_nist]

####### STEP 2.4 
mean_per_metabolite <- apply(nistdata_qc, 2, mean, na.rm = TRUE)
      sd_per_metabolite <- apply(nistdata_qc, 2, sd, na.rm = TRUE)
      cv_per_metabolite <- data.frame(mean = mean_per_metabolite,
                                sd = sd_per_metabolite,
                                cv = sd_per_metabolite/mean_per_metabolite)



cv_per_metabolite$flag<-NA
cv_per_metabolite$flag[which(cv_per_metabolite$cv>0.3)]<-"flagged"
cv_per_metabolite$flag[which(cv_per_metabolite$cv<=0.3)]<-"ok"

###### STEP 2.5 (136 metabolites removed so left 1111)
flagged_metabolites<-rownames(cv_per_metabolite)[which(cv_per_metabolite$flag=="flagged")]
'%ni%'<-Negate('%in%')
data_participants_V3<-data_participants_V2[,which(colnames(data_participants_V2)%ni%flagged_metabolites)]


################################################################################################################## STEP 3
################# STEP3.1 ==> Log2 transformation of the participants data generated in step 2.5
##################################################################################################################

names_meta<-colnames(data_participants_V3)
for(i in 1:length(names_meta)){
data_participants_V3[,names_meta[i]]<-log2(data_participants_V3[,names_meta[i]])
}


################################################################################################################## STEP 4
################# STEP4.1 ==> remove metabolites with missing greater than 30 percent in data generated STEP 3.1 
##################################################################################################################

missing_metabolites<-colMeans(is.na(data_participants_V3))
good_metab<-missing_metabolites[which(missing_metabolites<=0.30)]
metaname<-rownames(as.data.frame(good_metab))
data_participants_V4<-data_participants_V3[,metaname]

## This step we leave 991 metabolites for imputation

################################################################################################################## STEP 4
################# Imputation function from Gabi
################# Function to impute data points in metabolomics dataset using 
################# K nearest neighbours (KNN) per observation with selected metabolites 
################# (i.e., considering both KNN in sample as well as in metabolite dimensions)
################# implemented by Simone Wahl and Gabi Kastenmüller
##################################################################################################################
impute.knn.obs.sel <- function(dat, K=10) {
  
  results <- list()
  cor.cutoff <- 0.2     # use only variables with cor>0.2 for distance computation
  
  da1 <- dat 
  da1list <- da2list <- rep(list(dat),length(K)) 
  
  incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
  incom.obs <- which(apply(da1,1,function(x) any(is.na(x))))
  
  Cor <- cor(da1,use="p")
  
  D2list <- lapply(incom.vars, function(j) {
    varsel <- which(abs(Cor[j,])>cor.cutoff)  
    if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
    if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
    D2 <- as.matrix(dist(scale(da1[,varsel])),upper=T,diag=T) 
    if(any(is.na(D2))) {
      D2a <- as.matrix(dist(scale(da1)),upper=T,diag=T)*sqrt(length(varsel)/ncol(da1)) 
      D2[is.na(D2)] <- D2a[is.na(D2)] 
    }
    diag(D2) <- NA
    D2})
  names(D2list) <- incom.vars
  
  for (i in incom.obs){
    comvars <-  complete.cases(as.numeric(da1[i,]))
    for (j in which(!comvars)) {
      D2 <- D2list[[as.character(j)]]                                 
      if(any(!is.na(D2[i,]))) {
        KNNids <- order(D2[i,],na.last=NA)
        KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(da1[ii,j])))] 
      } else {
        KNNids  <- NULL
      }
      da1list <- lapply(1:length(da1list),function(ii) {
        k <- K[ii] 
        da <-  da1list[[ii]]
        if(!is.null(KNNids)) {
          KNNids_sel <- intersect(KNNids[1:min(k,length(KNNids))],KNNids_naomit)
        }
        if(length(KNNids_sel)<1) {
          KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))]
        } else if (length(which(sapply(KNNids_sel,function(ii) !is.na(da1[ii,j])))) < floor(k/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))] 
        if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
          da_sel <- da[KNNids_sel,j]
          da[i,j] <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) }
        da}) 
    }
  }
  da1list <- lapply(da1list, function(da) {
    da <- apply(da,2, function(x) {
      if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
      x}) 
    da})
  
  results <- c(results,list(da1list))
  names(results)[length(results)] <- "knn.sample.euc.sel" 
  rm(da1list,da2list)  
  return(results$knn.sample.euc.sel[[1]])
}


#################################################
################################################# Run imputations and write results
imputeddata<-impute.knn.obs.sel(data_participants_V4)
save(imputeddata,file="/home/sahmad/metabolomics/gut_liver_brain/quality_control_2021/KNN_imputed_metabolon_data_3rdMarch.RData")

