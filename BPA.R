# Backward progressive algorithm that predicts caste identities of a target stage based on the principle of developmental continuity. 
source('shared_functions.R')
# sampleTable.used is the sample information.
# Using prediction (validation) of 2nd instar larvae in M. pharaonis as an example.
train_age = '2nd'
test_age = '3rd'

# This step reads in the gene expression matrix and sample information into a list, exp_data. See shared_functions.R for detail.
exp_data = quick_subset(c(train_age,test_age),all.sampletable = sampleTable.used,deg_test = F,ab.type = 'count.upper')

# Step 0: Identificaton of caste DEGs at 3rd instar as features (candidate genes). Note: One can use all genes as features, but the signal to noise ratio will be very low.
exp_deg = quick_subset(test_age,all.sampletable = sampleTable.used,deg_test = T,deg_contrast = 'caste',LRT = F)
exp_deg$res = results(exp_deg$dds,contrast = c('caste','Gyne','Worker'),alpha = 5e-2,lfcThreshold = 0)
summary(exp_deg$res)

candidate_genes = rownames(exp_deg$res)[which(exp_deg$res$padj < 5e-2)]
candidate_genes = candidate_genes[which(candidate_genes %in% rownames(exp_data$abundance))]

train_data = which(exp_data$sampleInfo$age %in% train_age)
test_data = which(exp_data$sampleInfo$age %in% test_age)

# Code for BPA:
### Beginning of code:
reverse_pca.exp = function(train_data.exp, test_data.exp, train_age, test_age, sample.info,
                           correct_dev = T, method_d = 'combat',n_gene = 512){
  filter_table = sample.info[c(colnames(train_data.exp),colnames(test_data.exp)),]
  batch = droplevels(filter_table$age)
  modcombat = model.matrix(~1, data= filter_table)
  shared_genes = intersect(rownames(train_data.exp),rownames(test_data.exp))
  shared_genes = rownames(train_data.exp[shared_genes,])[order(apply(train_data.exp[shared_genes,],1,var),decreasing = T)[c(1:n_gene)]]
  exp.filter = cbind(train_data.exp[shared_genes,],test_data.exp[shared_genes,])
  filter_genes = shared_genes
  for(i in c(train_age, test_age)){
    filter_genes = intersect(filter_genes,rownames(exp.filter)[which(apply(exp.filter[,which(filter_table$age %in% i)],1,var) > 0)])
  }
  exp.filter = exp.filter[filter_genes,]
  if (correct_dev == T){
    if (method_d == 'combat'){
      combat_edata = ComBat(dat=exp.filter, batch=batch, mod=modcombat,mean.only = F,
                            par.prior=TRUE,  prior.plots=FALSE)}
    if (method_d == 'scale'){
      combat_edata = exp.filter
      for(i in c(train_age, test_age)){
        combat_edata[,which(filter_table$age %in% i)] =  t(scale(t(exp.filter[,which(filter_table$age %in% i)])))
      }
    }
    if (method_d == 'mean'){
      combat_edata = exp.filter
      for(i in c(train_age, test_age)){
        combat_edata[,which(filter_table$age %in% i)] = exp.filter[,which(filter_table$age %in% i)] - rowMeans(exp.filter[,which(filter_table$age %in% i)])}
    }
  }
  if (correct_dev != T)
  {combat_edata = exp.filter}
  trainning_data = combat_edata[,which(colnames(combat_edata) %in% rownames(filter_table)[which(filter_table$age %in% train_age)])]
  trainning_data.centre = trainning_data - rowMeans(trainning_data)
  n = 10
  tmp_1 = svd(trainning_data.centre,nu = n,nv = n)
  testing_data = combat_edata[,which(filter_table$age %in% c(test_age,train_age))]
  testing_data.centre = testing_data - rowMeans(testing_data)
  pc_data_test = data.frame(t(t(tmp_1$u) %*% testing_data.centre))
  pc_data_test = cbind(pc_data_test, filter_table[rownames(pc_data_test),])
  names(pc_data_test)[c(1:n)] = c(paste0('PC',c(1:n)))
  return(pc_data_test)
}
### End of code.

# Step 1: Construction of PCs from transcriptomes of 2nd instar larvae and 
# Step 2: Projection of the 3rd instar transcriptomes onto the constructed PCs.
pca_data = reverse_pca.exp(exp_data$abundance[candidate_genes,train_data],
                           exp_data$abundance[candidate_genes,test_data],
                           train_age,test_age,
                           sample.info = exp_data$sampleInfo,
                           n_gene = length(candidate_genes))

# Step 3.a: Visually examine the PC value similarity between 2nd instar and 3rd instar larvae. 
ggplot(pca_data,aes(x = PC1, y = PC2,size = log(body_length), col = caste))+
  geom_point()+
  scale_size_continuous(range = c(.5,3))+
  facet_wrap(~age)+
  theme_bw()

# Step 3.b: Stastitically identify PCs that associated with caste identities in 3rd instar larvae.  
summary(glm(PC1 ~ caste + log(body_length), data = droplevels(pca_data[which(pca_data$age %in% test_age),])))

# Step 3.c: Train a caste prediction model (with 3rd instar transcriptomic data), based on the identified caste PCs on Step 3.b.
pca.3rd =lda(caste ~ PC1 + PC2 , data =  droplevels(pca_data[which(pca_data$age %in% test_age),]),probability = TRUE)

# Step 4: Prediction of caste identities among 2nd instar larvae.
pca_data$predict_caste = predict(pca.3rd,pca_data)$class
pca_data$predict_prob = predict(pca.3rd,pca_data)$posterior[,'Gyne']
table(pca_data$predict_caste,pca_data$caste,pca_data$age)

pca_data$age_plot = pca_data$age
levels(pca_data$age_plot) = c("2nd instar",'3rd instar')
pca_data$caste = factor(pca_data$caste,levels = c("Gyne",'Worker'))

# Visualization of the prediction result.
ggplot(pca_data,aes(x =  PC1, y = PC2,size = body_length, fill = predict_prob, shape= paste(caste,age),col = paste(age,caste)),colour = NULL)+
  geom_point(stroke = .3)+
  scale_size_binned(range = c(0.3,2),n.breaks = 3,name = 'Body length')+
  scale_fill_gradient2(low = brewer.pal(n = 3, name = "Set1")[c(2)],high = brewer.pal(n = 3, name = "Set1")[c(1)],
                       mid = 'white',midpoint = .5,name = 'Probability\n of being gyne')+
  scale_color_manual(values = c("black",'black',brewer.pal(3,"Set1")[c(1,2)]))+
  facet_wrap(~age_plot)+
  scale_shape_manual(values = c(24,24,22,22),name = "Observed")+
  theme_bw()+
  xlab(label = 'PC1: 37.4% explained variance (2nd instar)')+ # This number is calculated in a separated analysis.
  ylab(label = 'PC2: 12% explained variance (2nd instar)')+ # This number is calculated in a separated analysis.
  theme(legend.position = "none",strip.text = element_text(size = 6,margin = margin(.02,0,.02,0, "in")),
        plot.margin = unit(c(0.01,0.01, 0,0),'in'),
        panel.spacing = unit(.5, "lines"),
        axis.title = element_text(size = 6),axis.text = element_text(size = 6,face = 'plain'))
