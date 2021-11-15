# Using the between stage transition to reveal the canalization pattern (developmental potential).

source('shared_functions.R')

# sampleTable.filter.used includes 1st to 3rd instar larvae, pre-pupae, early and old pupae and imagoes (adults).
# Caste identities of 1st instar larvae were predicted by BPA.
# The code here is for M.pharaonis.

exp_data = quick_subset(subset.age = age_levels[c(21:27)], # We exlude embryo samples.
                        all.sampletable = sampleTable.filter.used[which(sampleTable.filter.used$caste %in% c('Gyne','Worker','Unknown')),], 
                        deg_test = F,ab.type = 'abundance.quantile')

exp_data$sampleInfo$caste = factor(exp_data$sampleInfo$caste, levels = c("Gyne",'Worker','Unknown'))

# Becasue sample categories are not balanced (e.g. 30 workers vs 20 gynes), we perform bootstrapping to obtain a balanced set of samples.
exp_data_balanced = exp_data
exp_data_balanced$abundance = exp_data_balanced$abundance[apply(exp_data$counts,1,function(x) sort(x,decreasing = T)[10]) > 5,]
subsamples = c()
for(i in age_levels[c(21:27)]){
  n_sample = 100
  for (j in c("Gyne",'Worker')){
    subsamples = c(subsamples,sample(rownames(exp_data_balanced$sampleInfo)[which(exp_data_balanced$sampleInfo$age %in% i & 
                                                                                  exp_data_balanced$sampleInfo$caste %in% j)],n_sample,replace = T))
  }
}
exp_data_balanced$abundance = exp_data_balanced$abundance[,subsamples]
exp_data_balanced$sampleInfo = exp_data_balanced$sampleInfo[subsamples,]

# The following step normalizes for the between-age expression difference, so that the mean expression level of each age group is the same.
# Note that it's important to have a balanced set of gynes and workers in each age group.
exp_combat = quick_combat(exp_data_balanced$abundance,exp_data_balanced$sampleInfo,co_var = 'age',mean.only = F)

# Calculation of developmental potential
dist_caste = data.frame()
for (i in c(21:26)){
  between_stage = c(i,i+1)
  age_id = age_levels[i]
  tmp_candidate = rownames(exp_combat)
  test_tmp = summary_dist(exp_combat[tmp_candidate,which(exp_data_balanced$sampleInfo$caste %in% c('Unknown','Gyne','Worker'))],'man',
                          exp_data_balanced$sampleInfo,
                          age_levels[between_stage[1]],age_levels[between_stage[2]],
                          central = T,mean_fun = median,caste = 'caste',n_gene = 2000)
  test_tmp[,'n_gene'] = length(tmp_candidate)
  dist_caste = rbind(dist_caste, test_tmp)
}
dist_caste = droplevels(dist_caste)
levels(dist_caste$age) = c(levels(dist_caste$age)[c(1:3)],'Pre-pupa','Early pupa','Late pupa')
levels(dist_caste$age)[c(1:3)] = paste(levels(dist_caste$age)[c(1:3)],'instar')

# A simple visualization.
ggplot(dist_caste,aes(x = (Worker - Gyne)/(ref_between), y = log(body_length) ,col = caste))+
  geom_point()+
  scale_color_manual(values = c(brewer.pal(n = 3, name = "Set1")[c(1,2)],'black'))+
  facet_wrap(vars(age), ncol =1, scales = 'free_y',strip.position = 'right')+
  ylab(label = 'log (body length) (mm)')+
  xlab(label = expression(paste("Developmental potential",' (',delta,')')))+
  geom_density_2d(contour_var = "ndensity",alpha = .4)+
  theme(strip.text = element_text(size = 15))+
  theme_bw()

# Generate Figure 3C.
library(extrafont)
p_list =list()
range_list = list(c(1,1.4),c(1.4,2.7),c(1.5,3.5), 
                   c(2.5,3.4),c(2.5,3.4),c(2.5,3.4))
names(range_list) = levels(dist_caste$age)
levels(dist_caste$age)[c(4:6)] = c("Pre-pupae",'Early pupae','Late pupae')
for(i in levels(dist_caste$age)[c(1:6)]){
  if (i == levels(dist_caste$age)[1]){plot_col = c(brewer.pal(n = 3, name = "Set1")[c(1,2)])}
  else{
    plot_col = c(brewer.pal(n = 3, name = "Set1")[c(1,2)])}
  p_list[[i]] = ggplot(dist_caste[which(dist_caste$age %in% i),],aes(x = (Worker - Gyne)/(ref_between), y = (body_length)/10 ,col = caste))+
    geom_point(size = 1)+
    scale_color_manual(values = plot_col)+
    facet_wrap(vars(age), ncol =1,strip.position = 'right')+
    ylab(label = '')+
    scale_y_log10(labels = scales::number_format(accuracy = 0.1),limits = exp(range_list[[i]])/10)+
    xlab(label = expression(paste("Developmental potential",' (',delta,')')))+
    geom_density_2d(contour_var = "ndensity",alpha = .4)+
    xlim(c(-.7,.7))+
    theme_bw()+
    theme(legend.position = "none",plot.margin = unit(c(0,0.01,0,0),'in'),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 6,face = 'plain'),
          strip.text = element_text(size = 6,margin = margin(0,0.02,0,0.03, "in")))
}
p_list[c(1:6)] = lapply(p_list[c(1:6)],  "+", theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()))
p_list[[1]] = p_list[[1]] + theme(plot.margin = unit(c(.01,0.01,0,0),'in'))

p_list
