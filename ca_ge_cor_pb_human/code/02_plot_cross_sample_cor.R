rm(list=ls())
ddir <- '/home/whou10/data/zji/encode/analysis/promoter_gene_cor/pb/res/'
pdir <- '/home/whou10/scratch4/whou10/encode4/ca_ge_cor_pb_human/plot/'
rdir <- '/home/whou10/scratch4/whou10/encode4/ca_ge_cor_pb_human/result/'
pd <- readRDS(paste0(ddir, 'human_acrosssample.rds'))
str(pd) ## 58721 obs. of  3 
head(pd)

for (i in c(2,3)){
  br <- seq(min(pd[,i]), max(pd[,i]), length.out = 6)
  br[2:(length(br)-1)] <- round(br[2:(length(br)-1)], 3)
  c <- cut(pd[,i], breaks= br)
  if (i == 2) pd <- data.frame(pd, cut.rnasd = c)
  if (i == 3) pd <- data.frame(pd, cut.atacsd = c)  
}
pd <- pd[complete.cases(pd), ]
dim(pd) ##  51148     5
set.seed(123)
pd.tmp <- pd[sample(1:nrow(pd), 1e4, replace = FALSE), ]
str(pd.tmp)
library(gridExtra)
pdf(paste0(pdir, 'pcc_by_sd_acrosssample.pdf'), 
    width = 2.8, height = 4.5)
assign('var1', 'cut.rnasd')
assign('var2', 'cor')
assign('var3', 'cut.atacsd')
str(pd)
plist <- list()
plist[[1]] <- ggplot()+
  geom_jitter(data = pd.tmp, aes_string(x = var1, y = var2, color = var1), 
              size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var1, y = var2, fill = var1), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('snRNA-seq signals sd') + 
  ylab('Cross-gene Pearson correlation')

plist[[2]] <- ggplot()+
  geom_jitter(data = pd.tmp, aes_string(x = var3, y = var2, color = var3), 
              size = 0.2, stroke = 0, alpha = 0.8, width = 0.15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4), legend.position = 'none')  +
  geom_violin(data = pd, aes_string(x = var3, y = var2, fill = var3), scale = 'width', alpha = 0.1,
              draw_quantiles = 0.5) +
  scale_fill_brewer(palette = 'Set1', direction = -1)   + 
  scale_color_brewer(palette = 'Set1', direction = -1)  + 
  xlab('snATAC-seq signals sd') + 
  ylab('Cross-sample Pearson correlation')
grid.arrange(grobs=plist,nrow=2)
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

