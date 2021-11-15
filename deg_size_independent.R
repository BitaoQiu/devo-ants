# We detect 2nd instar body-size-independent differentially expressed genes (DEGs) as an example.
source('shared_functions.R')
# sampleTable.used includes all samples in M.pharaonis.

exp_data = quick_subset(age_levels[c(22)],all.sampletable = sampleTable.used,deg_test = F,ab.type = 'count.quantile') 

glm_contrast = function(geneID, exp_data, subsamples, model = 'exp ~ age*(caste_2 + log(body_size)) + library_2', co_var = c('library_2B','library_2C'), body_length = 'body_length',lm_method = rlm, abundance = 'abundance'){
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

target_samples = rownames(exp_data$sampleInfo)
# Using a single gene as an example:
glm_contrast('LOC105837931',exp_data,subsamples = target_samples,
             model = 'exp ~ caste + log(body_length)',body_length = 'body_length',
             co_var = 'casteGyne',lm_method = rlm,abundance = 'abundance')

# Examine whole transcriptome.
mpha_2nd_res = data.frame(t(future_sapply(rownames(exp_data$abundance), glm_contrast, 
                                          exp_data = exp_data, 
                                          subsamples = target_samples, 
                                          model = 'exp ~ caste + log(body_length)',
                                          body_length = 'body_length',co_var = 'casteGyne')))

# Calculate the p-values and adjusted pvalues:
mpha_2nd_res$pvalue = sapply(mpha_2nd_res$t.value,FUN = 
                               function(x) 2*pt(q = abs(x),df = length(mpha_samples) - 3,lower.tail = F))
mpha_2nd_res$padj = p.adjust(mpha_2nd_res$pvalue,method = 'fdr')

head(mpha_2nd_res[order(mpha_2nd_res$Value,decreasing = T),])
