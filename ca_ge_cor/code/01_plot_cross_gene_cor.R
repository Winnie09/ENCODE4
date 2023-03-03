rm(list=ls())
ddir <- '/home/whou10/data/zji/encode/analysis/promoter_gene_cor/sc/res/'
pdir <- '/home/whou10/scratch4/whou10/encode4/ca_ge_cor/plot/'
rdir <- '/home/whou10/scratch4/whou10/encode4/ca_ge_cor/result/'
pd <- readRDS(paste0(ddir, 'acrossgene.rds'))
stu = sub('_.*','',pd[,4])
pd <- data.frame(pd, experiment = stu)

for (i in c(2,3)){
  br <- seq(min(pd[,i]), max(pd[,i]), length.out = 6)
  br[2:(length(br)-1)] <- round(br[2:(length(br)-1)], 3)
  c <- cut(pd[,i], breaks= br)
  if (i == 2) pd <- data.frame(pd, cut.rnasd = c)
  if (i == 3) pd <- data.frame(pd, cut.atacsd = c)  
}
pd <- pd[complete.cases(pd), ]
dim(pd)

library(gridExtra)
pdf(paste0(pdir, 'pcc_by_sd_acrossgene.pdf'), 
    width = 2.8, height = 6.5)
assign('var1', 'cut.rnasd')
assign('var2', 'cor')
assign('var3', 'cut.atacsd')
assign('var4', 'experiment')
str(pd)
plist <- list()
plist[[1]] <- ggplot()+
  geom_jitter(data = pd, aes_string(x = var1, y = var2, color = var1), 
              size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var1, y = var2, fill = var1), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('snRNA-seq signals sd') + 
  ylab('Cross-gene Pearson correlation')

plist[[2]] <- ggplot()+
  geom_jitter(data = pd, aes_string(x = var3, y = var2, color = var3), 
              size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var3, y = var2, fill = var3), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('snATAC-seq signals sd') + 
  ylab('Cross-gene Pearson correlation')

plist[[3]] <- ggplot()+
  geom_jitter(data = pd, aes_string(x = var4, y = var2, color = var4), 
              size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var4, y = var2, fill = var4), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('Experiment') + 
  ylab('Cross-gene Pearson correlation')
grid.arrange(grobs=plist,nrow=3)
dev.off()

#########################################################################
pdf(paste0(pdir, 'hist_crossgene.pdf'), 
    width = 2, height = 1.2)
ggplot(data = pd) + 
  geom_histogram(aes(x = cor), 
                 color = 'black', fill = 'skyblue', alpha = 0.2,
                 binwidth = 0.005) +
  xlab('Cross-gene Pearson correlation') + 
  ylab('Count')+
  geom_vline(aes(xintercept=mean(cor)),
             color="blue", linetype="dashed", size=1)
dev.off()
