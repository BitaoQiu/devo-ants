# Detect canalizatin related genes.

source('shared_functions.R')
library(plyr)
# sampleTable.filter.used includes 2nd to 3rd instar larvae, pre-pupae, early and old pupae and imagoes (adults).
# The code here is for M.pharaonis.

exp_data.full = quick_subset(subset.age = age_levels[c(22:27)],all.sampletable = sampleTable.filter.used,
  deg_test = F,ab.type = 'vst.contrast',deg_contrast = 'caste*age')
# We use variance-stabilizing transformation to obtain a robust estimation of variance.
# Other normalization methods give similar results.
                             
exp_data.full$abundance.filtered = exp_data.full$abundance

cal_canalized = function(x, sampleInfo,pseudo_count = .1){
  x1_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Gyne')]
  x2_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Worker')]
  x1 = x[x1_id]
  x2 = x[x2_id]
  n1 = length(x1)
  n2 = length(x2)
  exp_diff = (median(x1)- median(x2))
  c_score =sign(exp_diff)*(abs(exp_diff)+pseudo_count)/(pseudo_count+sqrt((var(x1) + var(x2))))
  if ((var(x1) == 0 & var(x2) == 0)| abs(exp_diff) < pseudo_count ){c_score = 0}  
  return(c_score)
}

cal_canalized_length = function(x, sampleInfo,pseudo_count = .1,size = 'body_length'){
  x1_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Gyne' & !is.na(sampleInfo[names(x),size]))]
  x2_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Worker' & !is.na(sampleInfo[names(x),size]))]
  x1 = x[x1_id]
  x2 = x[x2_id]
  coef_diff = rlm(x[c(x1_id,x2_id)]~ log(sampleInfo[c(x1_id,x2_id),size]) + sampleInfo[c(x1_id,x2_id),'caste'])$coefficients[3]
  x1_res = resid(rlm(x1 ~ log(sampleInfo[x1_id,size])))
  x2_res = resid(rlm(x2 ~ log(sampleInfo[x2_id,size])))
  c_score = sign(coef_diff)*(abs(coef_diff)+pseudo_count)/(pseudo_count+sqrt((var(x1_res)+ var(x2_res))))
  if ((var(x1_res) == 0 & var(x2_res) == 0)| abs(coef_diff) < pseudo_count ){c_score = 0}
  return(c_score)
}

cal_canalized_score = function(x, sampleInfo,pseudo_count = .1, quantile.q = .75){
  # Using upper quantile to obstain a robust estimation. Mostly for adult samples.
  x1_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Gyne')]
  x2_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Worker')]
  x1 = x[x1_id]
  x2 = x[x2_id]
  n1 = length(x1)
  n2 = length(x2)
  exp_diff = (quantile(x1,quantile.q)- quantile(x2,quantile.q))
  c_score =sign(exp_diff)*(abs(exp_diff)+pseudo_count)/(pseudo_count+sqrt((var(x1) + var(x2))))
  if ((var(x1) == 0 & var(x2) == 0)| abs(exp_diff) < pseudo_count ){c_score = 0}
  return(c((var(x1) + var(x2)),exp_diff))
}

cal_score = function(exp_diff, var_sum,pseudo_count = .1){
  c_score =sign(exp_diff)*(abs(exp_diff)+pseudo_count)/(pseudo_count+sqrt(var_sum))
  if ((var_sum == 0)| abs(exp_diff) < pseudo_count ){c_score = 0}
  return(c_score)
}

canalized_score = data.frame(row.names = rownames(exp_data.full$abundance.filtered))
# We first calculate the ratio of between-caste expression difference and within-caste expression variation in each stage.
for(i in age_levels[c(22:27)]){
  if (i %in% c('2nd','3rd')){
  # For 2nd and 3rd instar larvae, we account for the variation due to body length difference.
    canalized_score[,i] = future_apply(exp_data.full$abundance.filtered[,which(exp_data.full$sampleInfo$age %in% i)],
                                       1, FUN = function(x) cal_canalized_length(x,exp_data.full$sampleInfo))}
  else{
    canalized_score[,i] = future_apply(exp_data.full$abundance.filtered[,which(exp_data.full$sampleInfo$age %in% i)], 
                                       1, FUN = function(x) cal_canalized(x,exp_data.full$sampleInfo))}
}



c.score_old_pupa = t(future_apply(exp_data.full$abundance.filtered[,which(exp_data.full$sampleInfo$age %in% "Pupa.Old")], 
                                   1, FUN = function(x) cal_canalized_score(x,exp_data.full$sampleInfo)))
c.score_imago = t(future_apply(exp_data.full$abundance.filtered[,which(exp_data.full$sampleInfo$age %in% "Imago")], 
                                1, FUN = function(x) cal_canalized_score(x,exp_data.full$sampleInfo)))

# Because of the huge variation in imagos (tissue lysis), the ratio in imagos are calculated by 
# combining the 3rd quantile expression difference in imagos and the variation in old pupae.

canalized_score$Imago = future_sapply(rownames(canalized_score),FUN = function(x){cal_score(c.score_imago[x,2],c.score_old_pupa[x,1])})

canalized_score[apply(canalized_score[,c(1:6)],1,anyNA),]

canalized_score[,c(1:6)][is.na(canalized_score[,c(1:6)])] = 0
examined_age = c(1:6) 
# We then quantify the trend of caste-expression-difference ratio across developmental stages.
canalized_score[,c('pvalue','cor')] = t(apply(canalized_score[,examined_age],1,FUN = function(x){
  if(anyNA(x)){return(c(NA,NA))}
  else{
    cor.tmp = cor.test(abs(as.numeric(x)), examined_age,method = 'p',alternative = 'g')
    return(c(cor.tmp$p.value,cor.tmp$estimate))}
  }))

canalized_score$same_direction = ((canalized_score$Imago*canalized_score$Pupa.Old > 0))
canalized_score$c.score = -log10(canalized_score$pvalue)

# Calculation of canalization score by combining (1) the trend of increasing caste-difference and (2) the caste-expression difference at the end stage:
canalized_score$combined = canalized_score$c.score*canalized_score$Pupa.Old # We use caste-expression difference in old pupae because of the poor sample quality in imagos.
canalized_score$combined[which(canalized_score$same_direction == F)] = 0
canalized_score = canalized_score[order(abs(canalized_score$combined),decreasing = T),]

canalized_score$p.same = apply(canalized_score[,c(1:6)],1,FUN = function(x)wilcox.test(x)$p.value)
canalized_score[,c('blastp','gene.name','NCBI','dmel','coverage','bitScore','Evalue','ortholog',"Aech_Ortholog")] = 
  mpha.info[match(rownames(canalized_score),mpha.info$geneID),c('blastp','gene.name','NCBI','dmel','coverage','bitScore','Evalue','ortholog','aech_ortholog')]

# Plotting one canalized gene as an example.
plot_single_candidate('LOC105837931',exp_data.full, abundance = 'abundance.filtered')
