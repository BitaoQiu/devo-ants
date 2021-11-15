# Identify the size threshold for each candidate genes within the metamorphosis gene module.
source('shared_functions.R'))
library(chngpt)
library(lmtest)

sampleTable.used = sampleTable.filter.2[which(sampleTable.filter.2$caste %in% c("Gyne",'Worker')),]

# Getting 2nd and 3rd instar samples in M. pharaonis.
exp_data = quick_subset(age_levels[c(22:23)], all.sampletable = sampleTable.used,deg_test = F,ab.type = 'abundance.log',txOut = F)
exp_data$abundance.filtered = exp_data$abundance[which(apply(exp_data$abundance,1, FUN = function(x) sort(x, decreasing = T)[10]) > 1),]
exp_data$abundance.filtered = betweenLaneNormalization(exp_data$abundance.filtered,which = 'full',round = F)

subset_samples = which(exp_data$sampleInfo$age %in% c('2nd','3rd'))

# Using LOC105832848 (E93) as an example:
tmp_data = data.frame(exp = exp_data$abundance.filtered['LOC105832848',subset_samples],
                      size = log(exp_data$sampleInfo[subset_samples,'body_length']),
                      caste = exp_data$sampleInfo[subset_samples,'caste'])

# We calculate the threshold for gyne and worker samples separately.
fit.gyne = chngptm(formula.1 = exp ~ 1, formula.2 = ~ size, data = tmp_data[which(tmp_data$caste %in% 'Gyne'),], type="M11", family="gaussian")
fit.worker = chngptm(formula.1 = exp ~ 1, formula.2 = ~ size, data = tmp_data[which(tmp_data$caste %in% 'Worker'),], type="M11", family="gaussian")

summary(fit.gyne)
summary(fit.worker)

# Test the significance.
# Null model (without threshold expression pattern)
fit.gyne.0=lm(exp ~ size, tmp_data[which(tmp_data$caste %in% 'Gyne'),]) 
fit.worker.0=lm(exp ~ size, tmp_data[which(tmp_data$caste %in% 'worker'),]) 

lrtest(fit.gyne, fit.gyne.0)
lrtest(fit.worker, fit.worker.0)

# Illustration:
plot_single_candidate_size('LOC105832848',exp_data, gene_info = mpha.info,abundance = 'abundance.filtered',method_f = 'loess')+
  facet_wrap(~life_stage,scales = 'free_x',nrow = 3)+
  geom_vline(xintercept = c(2.28,3.13),linetype = c(2,2),col = c('blue','red'))+
  scale_color_manual(values = c(brewer.pal(9, 'RdYlBu')[c(9,1)],'black'))

