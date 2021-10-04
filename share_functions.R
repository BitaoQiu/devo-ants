# Prepare libraries
library(tximport)
library(DESeq2)
library(preprocessCore)
library(pheatmap)
library("RColorBrewer")
library("FactoMineR")
library('factoextra')
library(tidyverse)
library(EDASeq)
library("igraph")
library(vegan)
library(ppcor)
library(gage)
library(e1071)
library(sva)
library("dplyr")

#load('~/Documents//KU_data_analysis/evo_devo/Salmon/Mpha_NCBI_ver2/common/GO.Rdata') 
#load('~/Documents//KU_data_analysis/evo_devo/Salmon/Dmel_NCBI/common/Dmel_tissues.Rdata')
#source('~/Documents//KU_data_analysis/evo_devo/Salmon/Mpha_NCBI_ver2/common/prepare_dmel_tissues_adult.R')
caste_col = rev(colorRampPalette(brewer.pal(n = 3,name = "Set1"))(3))
caste_col[1] = 'grey'
caste_aech_col = rev(colorRampPalette(brewer.pal(n = 5,name = "Set1"))(5))
caste_aech_col[1] = 'grey'
caste_aech_col[5] = 'red'
age_col = (colorRampPalette(brewer.pal(9, "YlOrRd"),bias = .6)(16))

colors <- rev(colorRampPalette(brewer.pal(9, "Blues"),bias = 1)(255))
colors_order <- rev(colorRampPalette(brewer.pal(9, "Blues"),bias = 0.5,interpolate = 'linear')(255))
age_col = gray.colors(16, start = 1, end = 0,gamma = 2.4)

#################################################################
# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1.3/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.color[v]
  }
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          fg = vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("square")$clip,
          plot=mytriangle)
#######

get_abundance = function(x, gene_file, sample_list, pseudo_count = 1e-5){
  abundance = log2(x[which(rownames(x) %in% gene_file),sample_list] + pseudo_count)
  abundance.quantile = normalize.quantiles(abundance)
  rownames(abundance.quantile) = rownames(abundance)
  colnames(abundance.quantile) = colnames(abundance)
  abundance.quantile[which(abundance.quantile < 0)] = 0
  return(abundance.quantile)
}

quick_pca = function(tmp.abundance, tmp.sampletable, return.PCA = F, 
                     pc.x = 'Dim.1', pc.y = 'Dim.2', ica = F,
                     annotate_color = 'caste',annotate_shape = 'age', n_gene = dim(tmp.abundance)[1],n = min(n_gene,10),size = 2,scale_unit = F){
  tmp.abundance = tmp.abundance[order(apply(tmp.abundance,1,var),decreasing = T)[c(1:n_gene)],]
  #tmp.abundance = tmp.abundance[which(apply(tmp.abundance,1,var)>0),]
  tmp.sampletable = tmp.sampletable[colnames(tmp.abundance),]
  if (ica == F){
    tmp.pca <- PCA(t(tmp.abundance),scale.unit = scale_unit, ncp = n, graph = FALSE)
    # Take a look at the amount of variations explained by each PC.
    #fviz_eig(abundance.quantile.pca, addlabels = TRUE,main = 'Explained variance for each PC')
    tmp.pca.var = tmp.pca$eig #Top 6 PC
    rownames(tmp.pca.var) = paste('Dim',c(1:length(rownames(tmp.pca.var))), sep = '.')
    tmp.pca.data = cbind(tmp.pca$ind$coord[,c(1:n)],tmp.sampletable)
  }
  else
  {
    n = 2
    tmp.pca = data.frame(fastICA(t(tmp.abundance),n.comp = n,alg.typ = 'deflation',method = 'C')$S)
    colnames(tmp.pca) = paste('Dim',c(1:n), sep = '.')
    tmp.pca.data = cbind(tmp.pca,tmp.sampletable)
    p1 = ggplot(tmp.pca.data, aes_string(x = pc.x, y = pc.y, colour = annotate_color, shape = annotate_shape)) +
      geom_point(size=size) +
      #  xlab(paste0(pc.x,':', "Explained variance:", round(tmp.pca.var[pc.x,'percentage of variance'],2), '%'))+
      # ylab(paste0(pc.y,':',"Explained variance:", round(tmp.pca.var[pc.y,'percentage of variance'],2), '%'))+
      # coord_fixed()+
      theme_bw()
    return(p1)
  }
  if (return.PCA == TRUE){
    return(tmp.pca)}
  else
    ggplot(tmp.pca.data, aes_string(x = pc.x, y = pc.y, colour = annotate_color, shape = annotate_shape)) +
    geom_point(size=size)+
    xlab(paste0(pc.x,':', "Explained variance:", round(tmp.pca.var[pc.x,'percentage of variance'],1), '%'))+
    ylab(paste0(pc.y,':',"Explained variance:", round(tmp.pca.var[pc.y,'percentage of variance'],1), '%'))+
    #  coord_fixed()+
    theme_bw()
}

cal_dist = function(tmp.abundance, dist_fun){
  if (dist_fun %in% c('s','p','k')){
    tmp.Dists = as.dist(1- cor(tmp.abundance,method = dist_fun,use = 'p'))} 
  else{
    tmp.Dists = dist(t(tmp.abundance),method = dist_fun)
  }
  return(tmp.Dists)}


quick_heatmap = function(tmp.abundance, tmp.sampletable, sample_annotate= c("age",'caste'), dist_fun = 's',order = T){
  callback = function(hc, mat,var = sample_annotate){
    age_caste = factor(c(paste(tmp.sampletable[,var[1]],tmp.sampletable[,var[2]],sep = '.')))
    sample_order = order(tmp.sampletable[,var[1]],tmp.sampletable[,var[2]],decreasing = T)
    age_caste = factor(age_caste, 
                       levels = unique(paste(tmp.sampletable[sample_order,var[1]],tmp.sampletable[sample_order,var[2]],sep = '.')))
    dend = reorder(hc, wts = age_caste,agglo.FUN = 'mean')
  }
  tmp.sampletable= droplevels(tmp.sampletable[colnames(tmp.abundance),])
  new_order = c()
  var_levels = levels(tmp.sampletable[,sample_annotate[1]])
  for(i in c(1:length(levels(tmp.sampletable[,sample_annotate[1]])))){
    sub_samples = rownames(tmp.sampletable)[which(tmp.sampletable[,sample_annotate[1]] %in% var_levels[c(i-1,i,i+1)])]
    tmp_cluster = hclust(cal_dist(tmp.abundance[,sub_samples], dist_fun))
    tmp_cluster.order = reorder(tmp_cluster,wts = tmp.sampletable[sub_samples,sample_annotate[1]], agglo.FUN = 'mean')
    tmp_cluster.order.lables = tmp_cluster.order$labels[tmp_cluster.order$order]
    new_order = c(new_order,
                  tmp_cluster.order.lables[which(tmp.sampletable[tmp_cluster.order.lables,sample_annotate[1]] %in% var_levels[i])])
  }
  tmp.abundance = tmp.abundance[,new_order]
   
  tmp.Dists = cal_dist(tmp.abundance, dist_fun)
  tmp.Dists.plot = as.matrix(tmp.Dists)
  colnames(tmp.Dists.plot) = colnames(tmp.abundance)
  rownames(tmp.Dists.plot) = colnames(tmp.abundance)
  #
  pheatmap(tmp.Dists.plot,annotation_col = tmp.sampletable[,sample_annotate],legend = T, annotation_legend = T, 
           clustering_callback = callback,cluster_rows = order,cluster_cols = order,
           show_colnames = F,
           show_rownames =  F,
           clustering_distance_rows = tmp.Dists,border_color = NA,
           clustering_distance_col = tmp.Dists,
         #  annotation_colors = ann_colors,
           color = colors)
}


quick_subset = function(subset.age, all.files = files.used, all.sampletable = sampleTable.filter.2,
                        deg_test = F, deg_contrast = "caste", ab.type = 'count.quantile',
                        LRT = F, balance = F, lrt_reduce = '1',fitT = 'parametric', tx2gene = tx2gene_gene,txOut = F){
  if ('all' %in% subset.age){
    subset.names = rownames(all.sampletable)
  } else{
    subset.names = rownames(all.sampletable[which(all.sampletable$age %in% subset.age),])}
  subset.sampleTable = droplevels(all.sampletable[subset.names,])
  # Create pseudo samples for replacement
  if (balance == T){
    subset.names.dulplicate = subset.names[grep('.',subset.names,fixed = T)]
    all.files.duplicate = all.files[sapply(subset.names.dulplicate,FUN = function(x) return(strsplit(x,split = '.',fixed = T)[[1]][1]))]
    names(all.files.duplicate) = subset.names.dulplicate
    all.files = c(all.files, all.files.duplicate)
  }
  subset.salmon = tximport(all.files[subset.names], type = "salmon", tx2gene = tx2gene, #Read in counts for gene only, or all (tx2gene)
                           countsFromAbundance = 'no', txOut = txOut)
  if (deg_test == T){
    subset.ddsTxi <- DESeqDataSetFromTximport(subset.salmon,
                                              colData = subset.sampleTable,
                                              design = as.formula(paste('~',deg_contrast)))
    if (LRT == T){
      subset.ddsTxi = DESeq(subset.ddsTxi,test = 'LRT',reduced = as.formula(paste('~',lrt_reduce)), sfType = 'poscounts',parallel = T, quiet = F,minmu=1e-6)
      subset.abundance = assay(vst(subset.ddsTxi,blind = F))
    }
    else if (LRT == F){
      subset.ddsTxi = DESeq(subset.ddsTxi,useT = T, sfType = 'poscounts',parallel = T,quiet = F,minmu=1e-6)
      subset.abundance = assay(vst(subset.ddsTxi,blind = F))
      
    }
  }
  else if (ab.type == 'vst.contrast'){
    subset.ddsTxi <- DESeqDataSetFromTximport(subset.salmon,
                                              colData = subset.sampleTable,
                                              design = as.formula(paste('~',deg_contrast)))
    subset.ddsTxi = estimateSizeFactors(subset.ddsTxi,type = 'poscounts')
    subset.ddsTxi = estimateDispersions(subset.ddsTxi)
    subset.abundance = assay(vst(subset.ddsTxi,blind = F,fitType = fitT))
  }
  else if (ab.type == 'rlog.age'){
    subset.ddsTxi <- DESeqDataSetFromTximport(subset.salmon,
                                              colData = subset.sampleTable,
                                              design = ~age)
    subset.abundance = assay(rlog(subset.ddsTxi,blind = F))
  }
  else{subset.ddsTxi <- DESeqDataSetFromTximport(subset.salmon,
                                                 colData = subset.sampleTable,
                                                 design = ~1)
  if (ab.type == 'count.raw'){
    subset.abundance = subset.salmon$counts + 1}
  else if (ab.type == 'count.quantile'){
    subset.abundance = subset.salmon$counts
    subset.abundance =  betweenLaneNormalization(subset.abundance, offset=FALSE, round=F,which = 'full')
    subset.abundance = log2(subset.abundance +1)
    #tmp_data.exp.quantile = tmp_data$abundance
    subset.abundance = subset.abundance[apply(subset.abundance,1,var)>0,]}
  else if (ab.type == 'count.quantile.filtered'){
    subset.abundance = subset.salmon$counts
    subset.abundance = subset.abundance[which(apply(subset.salmon$abundance,1,min) > 5), ]
    subset.abundance =  betweenLaneNormalization(subset.abundance, offset=FALSE, round=F,which = 'full')
    subset.abundance = log2(subset.abundance)
    #tmp_data.exp.quantile = tmp_data$abundance
    subset.abundance = subset.abundance[apply(subset.abundance,1,var)>0,]}
  else if (ab.type == 'count.upper'){ #Upper quantile, common approach, but the mean-variance problem still exist
    subset.abundance = subset.salmon$counts
    subset.abundance =  betweenLaneNormalization(subset.abundance, offset=FALSE, round=F,which = 'upper')
    subset.abundance = log2(subset.abundance +1)
    #tmp_data.exp.quantile = tmp_data$abundance
    subset.abundance = subset.abundance[apply(subset.abundance,1,var) >0,]}
  else if (ab.type == 'vst.null'){
    subset.abundance = assay(vst(subset.ddsTxi,blind = T,fitType = fitT))}
  else if (ab.type == 'rlog.null'){
    subset.abundance = assay(rlog(subset.ddsTxi,blind = T))}
  else if (ab.type == 'abundance.log'){
    subset.abundance = log2(subset.salmon$abundance + 1)}
  else if (ab.type == 'abundance.raw'){
    subset.abundance = subset.salmon$abundance + 1} 
  else if (ab.type == 'abundance.quantile'){
    subset.abundance = subset.salmon$abundance
    subset.abundance =  betweenLaneNormalization(subset.abundance, offset=FALSE, round=F,which = 'full')
    subset.abundance = log2(subset.abundance +1)
    #tmp_data.exp.quantile = tmp_data$abundance
    subset.abundance = subset.abundance[apply(subset.abundance,1,var)>0,]}
  else if (ab.type == 'cosine'){
    subset.abundance = subset.salmon$abundance
    subset.abundance = apply(subset.abundance,2,function(x) return(x/sum(x)))
      }
  else if (ab.type == 'sd'){
    subset.abundance = t(scale(t(get_abundance(subset.salmon$abundance,tx2gene$V2,subset.names))))
  }
  }
  if (txOut == T){
    subset.abundance = subset.abundance[which(rownames(subset.abundance) %in% tx2gene$V1),]
  }
  return(list(abundance = subset.abundance, sampleInfo = subset.sampleTable, dds = subset.ddsTxi, counts =  subset.salmon$counts + 1))
}

#####


reverse_pca.exp = function(train_data.exp, test_data.exp, train_age, test_age, sample.info = sampleTable.filter.2, 
                           correct_dev = T, check_function = T, method_d = 'combat',n_gene = 512){
  #exp_data.both = quick_subset(c(test_age,train_age),all.sampletable = sample.info,deg_test = F)
  #exp_data.both$sampleInfo$caste = as.character(exp_data.both$sampleInfo$caste)
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
  
  if (check_function == T){
    pc_data_test$body_size = sqrt(filter_table[rownames(pc_data_test),'body_size'])
    pc_data_test$box = filter_table[rownames(pc_data_test),'box']
    pc_data_test$e.box = filter_table[rownames(pc_data_test),'e.box']
    pc_data_test$library = filter_table[rownames(pc_data_test),'library']
    pc_data_test$sex = filter_table[rownames(pc_data_test),'gender']
    pc_data_test$lifestage = filter_table[rownames(pc_data_test),'lifestage']
    pc_data_test$pseudo.age = 0#filter_table[rownames(pc_data_test),'pseudo.age']
    if(pc_data_test$age[1] %in% age_levels[c(1:9)]){
      pc_data_test$body_size = 1
      pc_data_test$pseudo.age = 1}
    for(i in levels(pc_data_test$age)){
      for(j in levels(factor(pc_data_test$caste))){
        
        tmp_sub = which(pc_data_test$age %in% i & pc_data_test$caste %in% j)
        pc_data_test[tmp_sub,'body_size'] = scale(pc_data_test[tmp_sub,'body_size'])
        pc_data_test[tmp_sub,'pseudo.age'] = scale(pc_data_test[tmp_sub,'pseudo.age'])
      }
    }
    pc_data_test$caste.score = filter_table$caste.score
    pc_data_test$age.true = filter_table$age.true
    pc_data_test$caste.true = filter_table$caste.true
  }
  pc_data_test$sampleID = filter_table$sampleID
  return(pc_data_test)
}

reverse_pca.exp.simple = function(train_data.exp, test_data.exp,sample.info = sampleTable.filter.2, 
                           correct_dev = T, center = T, method_d = 'combat',n_gene = 512){
  filter_table = sample.info[c(colnames(train_data.exp),colnames(test_data.exp)),]
  filter_table$batch_T = 'A'
  filter_table[colnames(train_data.exp),'batch_T'] = "A"
  filter_table[colnames(test_data.exp),'batch_T'] = "B"
  filter_table$batch_T = factor(filter_table$batch_T)
  batch = droplevels(filter_table$batch_T)
  modcombat = model.matrix(~1, data= filter_table)
  shared_genes = intersect(rownames(train_data.exp),rownames(test_data.exp))
  shared_genes = rownames(train_data.exp[shared_genes,])[order(apply(train_data.exp[shared_genes,],1,var),decreasing = T)[c(1:n_gene)]]
  exp.filter = cbind(train_data.exp[shared_genes,],test_data.exp[shared_genes,])
  filter_genes = shared_genes
  exp.filter = exp.filter[filter_genes,]
  if (correct_dev == T){
    if (method_d == 'combat'){
      combat_edata = ComBat(dat=exp.filter, batch=batch, mod=modcombat,mean.only = F,
                            par.prior=TRUE,  prior.plots=FALSE)} 
    
  }
  if (correct_dev != T)
  {combat_edata = exp.filter}
  trainning_data = combat_edata[,colnames(train_data.exp)]
  testing_data = combat_edata[,colnames(test_data.exp)]
  if (center == T){
    trainning_data.centre = trainning_data - rowMeans(trainning_data)
    testing_data.centre = testing_data - rowMeans(testing_data)}
  else{
    trainning_data.centre = trainning_data
    testing_data.centre = testing_data  }
  
  n = 10
  tmp_1 = svd(trainning_data.centre,nu = n,nv = n)
  pc_data_train = data.frame(t(t(tmp_1$u) %*% trainning_data.centre))
  pc_data_test = data.frame(t(t(tmp_1$u) %*% testing_data.centre))
  #print(pc_data_test)
  pc_data_out = rbind(pc_data_train,pc_data_test)
 # print(pc_data_out)
  pc_data_out = cbind(pc_data_out, filter_table[rownames(pc_data_out),])
  names(pc_data_out)[c(1:n)] = c(paste0('PC',c(1:n)))
  pc_data_out$sampleID = filter_table[,'sampleID']
  #print(rownames(pc_data_out))
  return(pc_data_out)
}

reverse_pca.exp_raw = function(train_data.exp, test_data.exp, train_age, test_age, sample.info,scale_unit = F){
  filter_table = sample.info[c(colnames(train_data.exp),colnames(test_data.exp)),]
  shared_genes = intersect(rownames(train_data.exp),rownames(test_data.exp))
  exp.filter = cbind(train_data.exp[shared_genes,],test_data.exp[shared_genes,])
  trainning_data = exp.filter[shared_genes,which(filter_table$age %in% train_age)]
  pc.train = PCA(t(trainning_data),scale.unit = scale_unit,graph = F,ncp = 10)
  testing_data = exp.filter[shared_genes,which(filter_table$age %in% c(test_age,train_age))]
  pc.test = predict.PCA(pc.train,newdata = t(testing_data),scale.unit = scale_unit)
  pc_data_test  = pc.test$coord
  pc_data_test = cbind(pc_data_test, filter_table)
  colnames(pc_data_test)[c(1:10)] = paste0("PC",c(1:10))
  return(pc_data_test)
}


format_pca = function(exp, sampleInfo,scale_unit = F){
  pc_data = quick_pca(exp,tmp.sampletable =  sampleInfo,return.PCA = T,n = 10,scale_unit = scale_unit)
  pc_data.data_frame  = data.frame(pc_data$ind$coord)
  pc_data.data_frame[,colnames(sampleInfo)] = sampleInfo[rownames(pc_data.data_frame),]
  pc_data.data_frame = droplevels(pc_data.data_frame)
  return(pc_data.data_frame)
}

canalization_aov = function(pca_data, co_var = '~ age*caste + batch'){
  res_data = pca_data
  for(i in c(1:10)){
    model = as.formula(paste("Dim.",i,co_var,sep = ''))
    res_data[,i] = resid(glm(formula = model,data = pca_data,na.action = na.exclude))
  }
  return(res_data)
}


quick_combat = function(exp_data, sampleInfo, co_var = 'age',method_d = 'combat',var_interest = '~1', ref.batch = NULL,mean.only = F){
  var_interest = as.formula(var_interest)
  filter_table = droplevels(sampleInfo[colnames(exp_data),])
  filter_table[,co_var] = factor(filter_table[,co_var])
  batch = droplevels(filter_table[,co_var])
  co_type = levels(filter_table[,co_var])
  modcombat = model.matrix(var_interest, data= filter_table)
  exp.filter = exp_data[which(apply(exp_data,1,var) > 0),]
  
  if (method_d == 'combat'){
    combat_edata = ComBat(dat=exp.filter, batch=batch, mod=modcombat,mean.only = mean.only,ref.batch = ref.batch, 
                          par.prior=TRUE,  prior.plots=FALSE)} 
  if (method_d == 'scale'){
    combat_edata = exp.filter
    for(i in co_type){
      combat_edata[,which(filter_table[,co_var] %in% i)] =  t(scale(t(exp.filter[,which(filter_table[,co_var] %in% i)]),scale = mean.only))
    }
  }
  return(combat_edata)
}


transcriptome_dist = function(check_age, candidates = tx2gene_gene$V2,sampleInfo = sampleTable.filter.used, 
                              dist_fun = 'eu',trait = c('age','caste')){
  exp_check_data = quick_subset(check_age,all.sampletable = sampleInfo,deg_test = F,ab.type = 'count.quantile')
  exp_check_data$abundance =  exp_check_data$abundance[which(rownames(exp_check_data$abundance) %in% candidates),]
  n_gene = dim(exp_check_data$abundance)[1]
  exp_check_data$abundance_center = data.frame(center = apply(exp_check_data$abundance,1, median),exp_check_data$abundance)
  trans_dist = data.frame(sample_dist = as.matrix(dist(t(exp_check_data$abundance_center),method = dist_fun))[1,-1]/n_gene)
  trans_dist[,trait] =  sampleInfo[colnames(exp_check_data$abundance),trait]
  return(trans_dist)
}


find_age_gene = function(target_age, exp_all, n_gene = 50, caste = T,target_pheno = 'Gyne',method = 'p',abundance ='abundance'){
  subset_samples = which(exp_all$sampleInfo$age %in% target_age & exp_all$sampleInfo$caste %in% target_pheno)
  subset_exp = exp_all[[abundance]][,subset_samples]
  subset_info = exp_all$sampleInfo[subset_samples,]
  if (caste == T){
    gene_age_cor = future_apply(subset_exp,1,FUN = function(x)pcor.test(x,as.numeric(subset_info$age),
                                                                        z = as.numeric(subset_info$caste),method = method )$p.value)
  }else{
    gene_age_cor = future_apply(subset_exp,1,FUN = function(x)cor.test(x,as.numeric(subset_info$age),method = method)$p.value)
  }
  age_gene = rownames(subset_exp)[order(gene_age_cor)[c(1:n_gene)]]
  return(age_gene)
}

find_age_gene_glm = function(target_age, exp_all, n_gene = 50, model_formula = 'as.numeric(age) + caste + library_2'){
  subset_samples = which(exp_all$sampleInfo$age %in% target_age)
  subset_exp = exp_all$abundance[,subset_samples]
  subset_info = exp_all$sampleInfo[subset_samples,]
  target_formula = as.formula(paste('exp ~', model_formula))
  gene_age.pvalue = future_apply(subset_exp,1,FUN = function(x){
    gene_data = data.frame(exp = x)
    gene_data[,colnames(subset_info)] = subset_info
    gene_glm = summary(glm(target_formula,data = gene_data))
    return(coef(gene_glm)[2,4])
  })
  gene_age.pvalue = p.adjust(gene_age.pvalue)
  print(length(which(gene_age.pvalue < 5e-2)))
  age_gene = names(gene_age.pvalue)[order(gene_age.pvalue)[c(1:n_gene)]]
  return(age_gene)
}

correct_age = function(exp_raw, m = 'age'){
  if(m == 'age'){
    f  = function(x) resid(rlm(formula = x ~ exp_raw$sampleInfo$pseudo.age + exp_raw$sampleInfo$age))
  } 
  if (m=='pseudo'){
    f  = function(x) resid(loess(formula = x ~ exp_raw$sampleInfo$pseudo.age))  #Using loess model.
  } 
  if (m == 'caste'){
    f  = function(x) resid(rlm(formula = x ~ exp_raw$sampleInfo$pseudo.age + 
                                 exp_raw$sampleInfo$age +  
                                 exp_raw$sampleInfo$caste))}
  exp_raw.res = t(future_apply(exp_raw$abundance,1,FUN = f))
  return(exp_raw.res)
}

plot_candidate = function(geneid, exp_data = exp_norm.train,
                          gene_info = mpha_annotate, x = 'age', col = 'caste',ncol = 1,
                          shape = 'caste',abundance = 'abundance',geneOut = T){
  geneid = geneid[which(geneid %in% rownames(exp_data[[abundance]]))]
  if (length(geneid) < 2){
    return(print(geneid))
  }
  gene.exp = t(exp_data[[abundance]][geneid,])
  gene.data = cbind(gene.exp,exp_data$sampleInfo)
  gene.data.long = reshape(gene.data,varying = colnames(gene.exp),direction = 'long',v.names = 'Expression.level',timevar =  'geneID')
  gene.data.long$geneID = colnames(gene.exp)[gene.data.long$geneID]
  #  gene.data.long$gene.name = paste(gene.data.long$geneID,gene_info$blastp[match(gene.data.long$geneID, gene_info$geneID)])
  if (geneOut == T){
    gene.data.long$gene.name = gsub('%2C transcript variant.*','',gene_info$blastp[match(gene.data.long$geneID, gene_info$geneID)])} else{
      gene.data.long$gene.name = gsub('%2C transcript variant ','_',gene_info$blastp[match(gene.data.long$geneID, gene_info$txID)])}
  gene.data.long$gene.name = factor(gene.data.long$gene.name,levels =unique(gene.data.long$gene.name[match(geneid,gene.data.long$geneID)])) 
  gene.data.long[,'smooth_x'] = as.numeric(gene.data.long[,x])
  
  ggplot(gene.data.long, aes_string(x = x, shape = shape, colour = col, y = 'Expression.level'))+
    geom_boxplot(position =position_dodge(0.8))+
    geom_jitter(position=position_dodge(0.8), size = 2, alpha = 0.5)+
    geom_smooth(data = gene.data.long, method = 'loess',
                aes_string(x = 'smooth_x',y = 'Expression.level', colour = col),se = F,show.legend = FALSE,linetype = "dashed")+
    theme_bw()+
  #  scale_color_manual(values = caste_col)+
    xlab(label = 'Developmental stage')+
    #ylab(label = 'Expression abundance: log2(TPM)')+
    facet_wrap(facets = ~geneID + gene.name,scales = 'free_y',ncol = ncol)
}

plot_candidate_size = function(geneid, exp_data = exp_norm.train,gene_info = mpha.info, x = 'log(body_length)', 
                               col = 'caste',ncol = 1,shape = 'caste',abundance = 'abundance',method_f = 'rlm',geneOut = T){
  geneid = geneid[which(geneid %in% rownames(exp_data[[abundance]]))]
  if (length(geneid) <= 1){
    return(print(geneid))
  }
  gene.exp = t(exp_data[[abundance]][geneid,])
  gene.data = cbind(gene.exp,exp_data$sampleInfo)
  gene.data.long = reshape(gene.data,varying = colnames(gene.exp),direction = 'long',v.names = 'Expression.level',timevar =  'geneID')
  gene.data.long$geneID = colnames(gene.exp)[gene.data.long$geneID]
  if (geneOut == T){
    gene.data.long$blastp = gsub('%2C transcript variant.*','',gene_info$blastp[match(gene.data.long$geneID, gene_info$geneID)])} else{
      gene.data.long$blastp = gsub('%2C transcript variant ','_',gene_info$blastp[match(gene.data.long$geneID, gene_info$txID)])}
  gene.data.long$blastp = factor(gene.data.long$blastp,levels =unique(gene.data.long$blastp[match(geneid,gene.data.long$geneID)])) 
  gene.data.long$gene.name = gene_info$gene.name[match(gene.data.long$geneID, gene_info$geneID)]
  gene.data.long$gene.name = factor(gene.data.long$gene.name,levels = unique(gene.data.long$gene.name[match(geneid,gene.data.long$geneID)]))
  #gene.data.long[,'smooth_x'] = as.numeric(gene.data.long[,x])
  print('a')
  ggplot(gene.data.long, aes_string(x = x, shape = shape, colour = col, y = 'Expression.level'))+
    geom_point(size = 1)+
    geom_smooth(data = gene.data.long, method = method_f,
                aes_string(x = x,y = 'Expression.level', col = col,shape = shape),se = T,show.legend = FALSE,linetype = "dashed")+
    theme_bw()+
    #scale_color_manual(values = caste_col)+
    #ylab(label = 'Expression abundance: log10(TPM)')+
    facet_wrap(facets = ~geneID + blastp,scales = 'free_y',ncol = ncol)
}

plot_single_candidate = function(geneid, exp_data = exp_norm.train,gene_info = mpha_annotate, shape = 'library_2',
                                 x = 'age', col = 'caste',ncol = 1,abundance = 'abundance', geneOut = T,size = 2){
  gene.data = data.frame(Expression.level = exp_data[[abundance]][geneid,],exp_data$sampleInfo, geneID = geneid)
  if (geneOut == T){
    gene.data$gene.name = gene_info$blastp[match(gene.data$geneID, gene_info$geneID)]
    gene.data$gene.name = gsub('%2C transcript variant.*','',gene.data$gene.name)} else{
      gene.data$gene.name = gene_info$blastp[match(gene.data$geneID, gene_info$txID)]}
  gene.data[,'smooth_x'] = as.numeric(gene.data[,x])
  ggplot(gene.data, aes_string(x = x, colour = col,shape = shape, y = 'Expression.level'))+
    geom_boxplot(position =position_dodge(0.8),outlier.size = NULL)+
    geom_jitter(position=position_dodge(0.8), size = size, alpha = 0.5)+
    #geom_smooth(data = gene.data, method = 'loess',aes_string(x = 'smooth_x',y = 'Expression.level', colour = col),se = F,show.legend = FALSE,linetype = "dashed")+
    theme_bw()+
    #ylim(c(0,max(gene.data$Expression.level + 1)))+
    
    #scale_color_manual(values = caste_col)+
    ggtitle(label = gene.data$gene.name[1])+
    #ylab(label = 'Expression abundance: log2(TPM)')+
    xlab(label = 'Developmental stage')
}


plot_single_candidate_size = function(geneid, exp_data = exp_norm.train,gene_info = mpha_annotate,
                                      x = 'log(body_length)', col = 'caste',shape = 'caste',
                                      abundance = 'abundance',method_f = 'rlm',geneOut = T,
                                      lm_formula = "Expression.level ~ age*(log(body_length) + caste)"){
  gene.data = data.frame(Expression.level = exp_data[[abundance]][geneid,],exp_data$sampleInfo, geneID = geneid)
  if (geneOut == T){
    gene.data$gene.name = gene_info$blastp[match(gene.data$geneID, gene_info$geneID)]
    gene.data$gene.name = gsub('%2C transcript variant.*','',gene.data$gene.name)} else{
      gene.data$gene.name = gene_info$blastp[match(gene.data$geneID, gene_info$txID)]}
  #gene.data.predict = droplevels(gene.data[!(is.na(gene.data$body_length)),])
 # mod<-rlm(as.formula(lm_formula),data = gene.data.predict)
 # rlm_mod <-cbind(gene.data.predict,predict(mod,interval="confidence"))
  ggplot(gene.data, aes_string(x = x, colour = col, shape = shape, y = 'Expression.level'))+
    geom_point(size = 1)+
    geom_smooth(data = gene.data, method = method_f,
                aes_string(x = x,y = 'Expression.level', col = col,shape = shape),se = T,show.legend = FALSE,linetype = "dashed")+
    theme_bw()+
  #  ylim(c(0,max(gene.data$Expression.level + 1)))+
   # geom_line(aes(y=fit),data = rlm_mod,linetype = 2)+
    #scale_color_manual(values = caste_col)+
    #ylab(label = 'Expression abundance: log2(TPM)')+
    ggtitle(label = gene.data$gene.name[1])
}

plot_single_candidate_size_simple = function(geneid, exp_data = exp_norm.train,gene_info = mpha_annotate,
                                             x = 'log(body_length)', abundance = 'abundance',geneOut = T){
  gene.data = data.frame(Expression.level = exp_data[[abundance]][geneid,],exp_data$sampleInfo, geneID = geneid)
  if (geneOut == T){
    gene.data$gene.name = gene_info$blastp[match(gene.data$geneID, gene_info$geneID)]
    gene.data$gene.name = gsub('%2C transcript variant.*','',gene.data$gene.name)} else{
      gene.data$gene.name = gene_info$blastp[match(gene.data$geneID, gene_info$txID)]}
  ggplot(gene.data, aes_string(x = x, y = 'Expression.level'))+
    geom_point(size = 1)+
    theme_bw()+
    ggtitle(label = gene.data$gene.name[1])
}

plot_cor_genes = function(exp, x, y, geneInfo = mpha_annotate, col = 'caste',shape = 'caste'){
  exp_plot = t(exp$abundance)
  exp_plot = cbind(exp_plot, exp$sampleInfo[rownames(exp_plot),])
  ggplot(exp_plot, aes_string(x = x, y = y, col = col,shape = shape))+
    geom_point()+
    geom_smooth(method = 'rlm')+
    xlab(label = geneInfo$blastp[which(geneInfo$geneID == x)])+
    ylab(label = geneInfo$blastp[which(geneInfo$geneID == y)])+
    theme_bw()
}

exp_dist = function(exp_data,sampleInfo = sampleTable.filter.used, dist_fun = 's',trait = c('age','caste')){
  n_gene = dim(exp_data)[1]
  if (dist_fun %in% c('s','p','k')){
    n_gene = 1}
  exp_data.abundance_center = data.frame(center = apply(exp_data,1, mean),exp_data)
  trans_dist = data.frame(sample_dist = as.matrix(cal_dist(exp_data.abundance_center,dist_fun = dist_fun))[1,-1]/n_gene)
  trans_dist[,trait] =  sampleInfo[colnames(exp_data),trait]
  return(trans_dist)
}

check_gene_contr = function(res.pca, mpha.info, pc, res.train, log2 = F, npc = 10,geneID = 'geneID', 
                            check_stat = c("geneID",'coverage','length','blastp','NCBI','gene.name','GO_term')){
  gene_contribution = res.pca$var$contrib[order(res.pca$var$contrib[,pc],decreasing = T),c(1:npc)]
  gene_cor = res.pca$var$cor[order(res.pca$var$contrib[,pc],decreasing = T),c(1:npc)]
  colnames(gene_cor) = paste('cor',colnames(gene_cor),sep = '_')
  pc_gene.info = data.frame(mpha.info[match(rownames(gene_contribution),mpha.info[[geneID]]),check_stat])
  pc_gene.info = cbind(pc_gene.info,gene_contribution,gene_cor)
  if(log2 == T){
    pc_gene.info$fold_change = res.train[pc_gene.info$geneID,'log2FoldChange']
    pc_gene.info$padj = res.train[pc_gene.info$geneID,'padj']}
  return(pc_gene.info)
}
quick_pc_function = function(exp, sampleInfo, target_pc, geneInfo = mpha.info, f_type = 'GO'){
  pc_data = quick_pca(exp,sampleInfo,return.PCA = T)
  pc_gene = check_gene_contr(pc_data,geneInfo,pc = target_pc,log2 = F)
  pc_gene_cor = pc_gene[,paste('cor_Dim',target_pc,sep = '.')]
  names(pc_gene_cor) = pc_gene$NCBI
  if (f_type == 'GO'){
    gsets = go.bp}
  if (f_type == "KEGG"){
    gsets = kegg.gs}
  fc.go.p <- gage(pc_gene_cor, gsets = gsets, ref = NULL, samp = NULL,same.dir = T)
  return(fc.go.p)
}

quick_enrichment = function(res_data,geneInfo = mpha.info, stat = 'stat',same.dir = F, f_type = 'GO',rank = F){
  if (f_type == 'GO'){
    gsets = go.bp}
  if (f_type == "KEGG"){
    gsets = kegg.gs}
  stat_value = res_data[,stat]
  names(stat_value) = geneInfo$NCBI[match(rownames(res_data),geneInfo$geneID)]
  fc.go.p <- gage(stat_value, gsets = gsets, ref = NULL, samp = NULL,same.dir = same.dir,set.size = c(3,500),rank.test = rank)
  return(fc.go.p)
}

quick_enrichment_paired = function(res_data, stat, ref, samp, geneInfo = mpha.info, same.dir = F, f_type = 'GO',compare = 'as.group',rank = F){
  if (f_type == 'GO'){
    gsets = go.bp}
  if (f_type == "KEGG"){
    gsets = kegg.gs}
  stat_value = as.matrix(res_data[,stat])
  rownames(stat_value) = geneInfo$NCBI[match(rownames(res_data),geneInfo$geneID)]
  fc.go.p <- gage(stat_value, gsets = gsets, ref = ref, samp =  samp,same.dir = same.dir,set.size = c(3,500),compare = compare,rank.test = rank)
  return(fc.go.p)
}

cal_canalized = function(x, sampleInfo, size = 'body_length'){
  # Calculate canalization score controlled with body size
  x1_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Gyne')]
  x2_id = names(x)[which(sampleInfo[names(x),'caste'] %in% 'Worker')]
  x1 = x[x1_id]
  x2 = x[x2_id]
  coef_diff = glm(x[c(x1_id,x2_id)]~ log(sampleInfo[c(x1_id,x2_id),'body_length']) + sampleInfo[c(x1_id,x2_id),'caste'])$coefficients[3]
  x1_res = resid(glm(x1 ~ log(sampleInfo[x1_id,'body_length'])))
  x2_res = resid(glm(x2 ~ log(sampleInfo[x2_id,'body_length'])))
  c_score = (coef_diff)/sqrt((var(x1_res) + var(x2_res)))
  #c_score = (mean(x1)- mean(x2))/sqrt((var(x1_res) + var(x2_res)))
  return(c_score)
}

summary_dist = function(exp,dist_fun, sampleInfo, target_age, ref_age, central = F, 
                        mean_fun = median,caste = 'caste', age = 'age',c1 = "Worker",c2 = "Gyne",n_gene = 1000){
  sampleInfo = sampleInfo[colnames(exp),]
  exp = exp[order(apply(exp[,which(sampleInfo[,age] == ref_age)],1,var),decreasing = T)[c(1:n_gene)],]
  if (central == F){
    samples_dist = as.matrix(cal_dist(exp,dist_fun))
    samples_dist.target = data.frame(samples_dist[which(sampleInfo[,age] %in% target_age),
                                                  which(sampleInfo[,age] %in% ref_age)])
    #samples_dist.target = samples_dist.target[order(apply(samples_dist.target,1,sum),decreasing = T)[-c(1,2)],]
    samples_dist.target[,c1] = apply(samples_dist.target[,sampleInfo[colnames(samples_dist.target),caste] %in% c1], 1,mean_fun)
    samples_dist.target[,c2] = apply(samples_dist.target[,sampleInfo[colnames(samples_dist.target),caste] %in% c2], 1,mean_fun)
    samples_dist.target = samples_dist.target[,c(c2,c1)]
    samples_dist.ref= samples_dist[which(sampleInfo[,age] %in% ref_age & sampleInfo[,caste] == c1),
                                   which(sampleInfo[,age] %in% ref_age & sampleInfo[,caste] == c2)]
    samples_dist.ref = samples_dist.ref[,order(apply(samples_dist.ref,2,sum),decreasing = T,na.last = T)[-c(1:2)]]
    samples_dist.ref = samples_dist.ref[order(apply(samples_dist.ref,1,sum),decreasing = T,na.last = T)[-c(1:2)],]
    samples_dist.ref_between = as.numeric(samples_dist.ref)
    ref_between = mean_fun(samples_dist.ref_between)
    samples_dist.target$ref_between = ref_between
  }
  else {
    exp_target = exp[,which(sampleInfo[,age] %in% target_age)]
    exp_ref = exp[,which(sampleInfo[,age] %in% ref_age)]
    exp_target.out = data.frame(exp_target,
                                Worker = apply(exp_ref[,which(sampleInfo[colnames(exp_ref),caste] %in% c1)],1, mean_fun))
    exp_target.out = data.frame(exp_target.out,
                                Gyne = apply(exp_ref[,which(sampleInfo[colnames(exp_ref),caste] %in% c2)],1, mean_fun))
    colnames(exp_target.out)[ncol(exp_target.out)-1] = c1
    colnames(exp_target.out)[ncol(exp_target.out)] = c2
    ref_between = cal_dist(exp_target.out[,c(c2,c1)],dist_fun)[1]
    samples_dist.target = data.frame(as.matrix(cal_dist(exp_target.out,dist_fun)))[colnames(exp_target),c(c1,c2)]
    samples_dist.target$ref_between = ref_between
  }
  samples_dist.target = cbind(samples_dist.target,
                              sampleInfo[rownames(samples_dist.target),])
  #samples_dist.target[,c('worker','gyne')] = samples_dist.target[,c('worker','gyne')]/rowSums(samples_dist.target[,c('worker','gyne')])
  return (samples_dist.target)
}


exp_dist_glm = function(exp_data,sampleInfo, model = '1', p = 1){
  n_gene = dim(exp_data)[1]
  gene_dist = t(future_apply(exp_data,1, FUN = function(x){
    tmp_matrix = data.frame(exp = x)
    tmp_matrix[,colnames(sampleInfo)] = sampleInfo[rownames(tmp_matrix),]
    return(resid(glm(formula = as.formula(paste('exp ~',model)),data =  tmp_matrix)))}))
  trans_dist = data.frame(sample_dist = future_apply(gene_dist,2,FUN = function(x) sum(abs(x)^p)/n_gene))
  #  print(trans_dist)
  trans_dist =  cbind(trans_dist, sampleInfo[colnames(exp_data),])
  return(trans_dist)
}


glm_contrast = function(geneID, exp_data, subsamples, model = 'exp ~ age*(caste_2 + log(body_size)) + library_2',
                        co_var = c('library_2B','library_2C'), body_length = 'body_length',lm_method = rlm, abundance = 'abundance'){
  tmp_data = data.frame(exp = exp_data[[abundance]][geneID, subsamples])
  tmp_data = cbind(tmp_data, exp_data$sampleInfo[subsamples,])
  if (!is.na(body_length)){
    tmp_data = tmp_data[which(!is.na(tmp_data[[body_length]])),]}
  tmp_data = droplevels(tmp_data)
  fit = lm_method(as.formula(model),data = tmp_data) # Using rlm to estimate the batch effect is more accurate is less prone to biases.
  tryCatch(
    {sum_fit = summary(fit)
    fit_coefficients = sum_fit$coefficients[co_var,]
    return(fit_coefficients)},
    error = function(e) return(rep(0,6)))
}

glm_contrast_predict = function(geneID, exp_data,subsamples, model = 'exp ~ age_2*(caste_2 + log(body_length)) + library_2', check_age = age_levels[c(22:28)],
                                abundance = 'abundance'){
  tmp_data = data.frame(exp = exp_data[[abundance]][geneID,subsamples])
  tmp_data = cbind(tmp_data, exp_data$sampleInfo[subsamples,])
  if (1 %in% grep('body_length',model) ){
    tmp_data = tmp_data[which(!is.na(tmp_data$body_length)),]} else{
      tmp_data = tmp_data
    }
  tmp_length = min(tmp_data[,'body_length'],na.rm = T)
  tmp_data = droplevels(tmp_data)
  fit = rlm(as.formula(model),data = tmp_data)
  fit_predict = predict(fit, data.frame(age_2 = rep(check_age, each = 2),age = rep(check_age, each = 2),
                                        caste_2 = c("Gyne",'Worker'),caste = c("Gyne",'Worker'), 
                                        body_length = tmp_length,library_2 = "A"))
  return((fit_predict))
}


draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 1, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

sim_GO = function(x, by = 'pvalue', cutoff = .8,decreasing = F,Cluster = 'Cluster', count_filter = 2, q_cutoff = 1e-1){
  # Cluster eq Caste means separating the output by caste
  x$Cluster_ID = paste(x[,Cluster],x$ID,sep = "_")
  # Clustering GO functions based on similarity of IDs. 
  x.filter = x[which(x$qvalue < q_cutoff & x$Count >= count_filter),] 
  retain_terms = x.filter[order(x.filter$ID,x.filter$Count,decreasing = T),] 
  # When there are multiple groups, retain the one that includes the most genes.
  retain_terms = retain_terms[!duplicated(retain_terms$Cluster_ID),c('Cluster_ID','Cluster','ID','geneID')]
  rownames(retain_terms) = retain_terms$Cluster_ID
  full_members = unique(unlist(sapply(retain_terms$geneID,FUN = strsplit,split ='/')))
  go_sim = t(apply(retain_terms,1,FUN = function(x) full_members %in% strsplit(x[4],split ='/')[[1]]))
  colnames(go_sim) = full_members
  
  go_sim[go_sim == FALSE] = 0
  go_sim[go_sim == TRUE] = 1
  # GO similarity clustering (separately for gyne and worker caste), using ID as vector \n
  # GO similarity based on shared memberships: intersect/union
  go_sim.cor = matrix(nrow = nrow(go_sim),ncol = nrow(go_sim),dimnames = list(rownames(go_sim),rownames(go_sim)))
   for(i in rownames(go_sim)){
    for(j in rownames(go_sim)){
      go_sim.cor[i,j] = sum(go_sim[i,]*go_sim[j,])/sum(sapply(go_sim[i,] + go_sim[j,],min,1))}
  }
  go_clusters = cutree(hclust(as.dist(1 - go_sim.cor),method = 'average'),h = cutoff) # Select the right cut-off
  plot(hclust(as.dist(1 - go_sim.cor)))
  abline(h = cutoff)
  x.filter$GO_cluster = go_clusters[match(x.filter$Cluster_ID,names(go_clusters))] 
  # Within each cluster, retain one term that with the lowest qvalue.
  go_filtered = tapply(x.filter$Cluster_ID,INDEX = x.filter$GO_cluster,function(x){ 
    tmp = x.filter[which(x.filter$Cluster_ID %in% x),]
    return(tmp$Cluster_ID[order(tmp[,by],decreasing = decreasing)][1])})
  return(go_filtered)
}
