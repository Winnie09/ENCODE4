rm(list=ls())
ddir <- '/home/whou10/data/zji/encode/analysis/rnaatac_cor/promoter/sc/res/'
pdir <- '/home/whou10/scratch4/whou10/encode4/ca_ge_cor_sc_promoter/plot/'
rdir <- '/home/whou10/scratch4/whou10/encode4/ca_ge_cor_sc_promoter/result/'
dir.create(pdir, recursive = T)
dir.create(rdir, recursive = T)
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
theme_set(.new_theme)


pd <- readRDS(paste0(ddir, 'acrosssample.rds'))
str(pd)
stu = sub('_.*','',pd[,6])

tb <- read.csv('/home/whou10/data/zji/encode/data/samplemeta/list.csv')
str(tb)
stu2 <- paste0(stu, tb[match(stu, tb[,1]), 'Tissue'])
pd <- data.frame(pd, experiment = stu, exp_tissue = stu2, 
                 tissue = tb[match(stu, tb[,1]), 'Tissue'],
                 species = tb[match(stu, tb[,1]), 'Species'])



i = 2
for (i in c(2,3)){
  br <- seq(min(pd[,i]), max(pd[,i]), length.out = 6)
  br[2:(length(br)-1)] <- round(br[2:(length(br)-1)], 3)
  c <- cut(pd[,i], breaks= br)
  if (i == 2) pd <- data.frame(pd, cut.rnazp = c)
  if (i == 3) pd <- data.frame(pd, cut.ataczp = c)  
}
str(br)
pd <- pd[complete.cases(pd), ]
dim(pd)

library(ggplot2)
library(gridExtra)
assign('var1', 'cut.rnazp')
assign('var2', 'cor')
assign('var3', 'cut.ataczp')
assign('var4', 'experiment')
str(pd)

pdf(paste0(pdir, 'pcc_by_sd_acrosssample.pdf'), width = 2.2, height = 4.2)
plist <- list()
plist[[1]] <- ggplot()+
  geom_jitter(data = pd, aes_string(x = var1, y = var2, color = var1), 
              size = 0.1, stroke = 0, alpha = 0.6, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var1, y = var2, fill = var1), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('snRNA-seq non-zero proportion') + 
  ylab('Cross-sample Pearson correlation')

plist[[2]] <- ggplot()+
  geom_jitter(data = pd, aes_string(x = var3, y = var2, color = var3), 
              size = 0.1, stroke = 0, alpha = 0.6, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var3, y = var2, fill = var3), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('snATAC-seq non-zero proportion') + 
  ylab('Cross-sample Pearson correlation')
grid.arrange(grobs=plist,nrow=2)
dev.off()


pdf(paste0(pdir, 'pcc_by_experiment_acrosssample.pdf'), width = 8.3, height = 2)
ggplot()+
  geom_jitter(data = pd, aes_string(x = var4, y = var2, color = var4), 
              size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var4, y = var2, fill = var4), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  #scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  #scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('Experiment') + 
  ylab('Cross-sample Pearson correlation')
dev.off()

assign('var5', 'tissue')
pdf(paste0(pdir, 'pcc_by_tissue_acrosssample.pdf'), width = 9, height = 2.5)
ggplot()+
  #geom_jitter(data = pd[sample(1:nrow(pd), 1e5),], aes_string(x = var5, y = var2, color = var4), 
  #            size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  +
  geom_violin(data = pd, aes_string(x = var4, y = var2, fill = var5), 
              scale = 'width', alpha = 0.9,
              draw_quantiles = 0.5) +
  #scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  #scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('Tissue') + 
  ylab('Cross-sample Pearson correlation')+
  guides(colour = "none")
dev.off()


assign('var5', 'species')
pdf(paste0(pdir, 'pcc_by_species_acrosssample.pdf'), width = 9, height = 2.5)
ggplot()+
  #geom_jitter(data = pd[sample(1:nrow(pd), 1e5),], aes_string(x = var5, y = var2, color = var4), 
  #            size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 0.6))  +
  geom_violin(data = pd, aes_string(x = var4, y = var2, fill = var5), 
              scale = 'width', alpha = 0.5,
              draw_quantiles = 0.5) +
  #scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  #scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('Species') + 
  ylab('Cross-gene Pearson correlation') +
  guides(colour = "none")
dev.off()


#########################################################################
pdf(paste0(pdir, 'hist_crosssample.pdf'), 
    width = 2, height = 1.2)
ggplot(data = pd) + 
  geom_histogram(aes(x = cor), 
                 color = 'black', fill = 'skyblue', alpha = 0.2,
                 binwidth = 0.01) +
  xlab('Cross-sample Pearson correlation') + 
  ylab('Count')+
  geom_vline(aes(xintercept=mean(cor)),
             color="blue", linetype="dashed", size=1)
dev.off()


