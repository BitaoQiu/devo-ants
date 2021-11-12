# Construct the transcriptomic developmental trajectory of a target species.
source('shared_functions.R')
# abundance.quantile is the (quantile normalized) gene expression matrix (log2 transformed) of a target species (e.g., M. pharaonis) in all developmental stages. 
# Rows are genes, columns are samples.
# sampleInfo.Table is the sample information of a target species.

devo_matrix = cor(abundance.quantile,method = 's')
devo_matrix[which(devo_matrix < 0.81)] = 0 # Two samples are connected if their transcriptomic correlation coefficient > 0.81.
filter_sample = which(apply(devo_matrix,1,function(x) sum(x > .9)) < 5) # Remove outlier samples.
devo_matrix = devo_matrix[-filter_sample, -filter_sample]

exp_data_info = sampleInfo.Table 
exp_data_info = exp_data_info[-filter_sample,]
exp_data_info = droplevels(exp_data_info)
exp_data_info$caste = factor(exp_data_info$caste,levels = c('Gyne','Worker','Unknown'))

devo_net = graph_from_adjacency_matrix(adjmatrix = devo_matrix, weighted = T,mode = 'undirected',diag = F) # Transform correlation coefficient matrix into network.
V(devo_net)$label <- NA
V(devo_net)$Age = exp_data_info$age
V(devo_net)$Caste = exp_data_info$caste

age_col_net = c(colorRampPalette(brewer.pal(9, "Blues"),bias = .6)(9)[c(1:9)], # For embryos
                (colorRampPalette(brewer.pal(9, "YlOrRd"),bias = 1)(7))) # For larvae and imago

age_text =  levels(exp_data_info$age)
age_text[c(1:9)] = age_text[c(1:9)]
age_text[c(10:12)] = paste(age_text[c(10:12)],'instar')
age_text[c(13:15)] = c("Pre-pupae",'Early pupae','Late pupae')
age_text[16] = "Adults"
ggraph(devo_net,layout = 'igraph', algorithm = 'fr')+
  geom_edge_link(n = 2,alpha = .8, aes(col = weight^1.5, width = (weight^2/5)), show.legend = FALSE) +
  scale_edge_color_gradient(low = 'white',high = weigth_col[16])+
  scale_edge_width_continuous(range = c(0,.4))+
  geom_node_point(size = 1,stroke = .1,aes(shape = factor(Caste),fill = factor(Age))) +
  scale_fill_manual(values = age_col_net,labels = age_text)+
  scale_shape_manual(values = c(24,22,21),labels = levels(exp_data_info$caste))+
  theme_graph(background = NA, foreground = NA, fg_text_colour = NA)
