library(tidyverse)
library(ggpubr)
# combine ATAC and RNA-seq data -------------------------------------------
# read data 
SE <- read_csv("results/data/atac/annotated/SE.csv") %>% 
  select(Geneid = "Gene Name", tss_dist = "Distance to TSS", l2fc, N, atac_state) %>%
  mutate(type = "SE")
TE <- read_csv("results/data/atac/annotated/TE.csv") %>%
  select(Geneid = "Gene Name", tss_dist = "Distance to TSS", l2fc, N, atac_state) %>% 
  mutate(type = "TE")
fano_counts <- read_csv("results/data/rna_preprocessing/fano_mean_counts.csv")
markers <- read_csv("results/data/rna_markers/activated_markers.csv")

dir.create("results/plots/ATAC_RNA")

# remove TE in SE
TE <- TE %>% filter(!Geneid %in% SE$Geneid)

# remove infinite L2FC
TE <- TE %>% filter(is.finite(l2fc))

# extract only one TE/SE by distance to TSS
TE <- TE %>% arrange(abs(tss_dist)) %>% group_by(Geneid) %>% slice_head(n = 1) %>% ungroup()
SE <- SE %>% arrange(abs(tss_dist)) %>% group_by(Geneid) %>% slice_head(n = 1) %>% ungroup()

SE_TE <- bind_rows(SE, TE)

# RNA l2fc across SE and TE conditions
inner_join(markers, SE_TE, by = c("gene" = "Geneid")) %>% select(type, atac_state)%>% table()
  
inner_join(markers, SE_TE, by = c("gene" = "Geneid")) %>%
  mutate(atac_state = factor(atac_state, levels = c("gained", "unchanged", "lost"))) %>%
  ggplot(aes(x = atac_state, y = avg_log2FC, fill = atac_state)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_grid(cols = vars(type), scales="fixed")+ cowplot::theme_cowplot(30) + cowplot::panel_border(color = "black")+
  coord_cartesian(ylim= c(-1,1))+
  scale_fill_manual(values = c("#FCBBB6",  "white", "#80DFE2")) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave("results/plots/ATAC_RNA/SE_TE_state_RNA_lfc.svg", height = 7, width = 10)

# fano factor plot
fano_counts <- fano_counts %>%
  pivot_wider(names_from = sample, values_from = c(fano, mean)) %>%
  mutate(fano_ratio = fano_2/fano_1)

# scatter plot
inner_join(markers, SE_TE, by = c("gene" = "Geneid")) %>% 
  inner_join(fano_counts, by = c("gene" = "Geneid")) %>%
  ggplot(aes(x = l2fc, y = log2(fano_ratio), color = l2fc)) + 
  geom_point(size = 3, alpha = 1)  + 
  stat_smooth(method = "lm", col = "black") +
  stat_cor(method = "spearman", size = 10) +
  scale_color_gradientn(colours=colorRampPalette(c("#386cb0", "#ef3b2c"))(100))+
  geom_hline(yintercept = 0)+
  guides(color = F)+
  cowplot::theme_cowplot(font_size = 40, line_size = 1)+
  facet_grid(cols = vars(type))

ggsave("results/plots/ATAC_RNA/SE_TE_fano_scatter.svg", height = 7, width = 10)

# fano factor boxplot between SE and TE
inner_join(markers, SE_TE, by = c("gene" = "Geneid")) %>% 
  inner_join(fano_counts, by = c("gene" = "Geneid")) %>% 
  pivot_longer(cols = c("fano_1", "fano_2")) %>%
  ggplot(aes(x = type, y = value, fill = type)) + 
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(0,0.7)) +
  cowplot::theme_cowplot(30) +
  cowplot::panel_border(color = "black") +
  scale_fill_manual(values = c("#FCBBB6", "#80DFE2")) + 
  facet_grid(cols = vars(name)) +
  theme(legend.position = "none", strip.background = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave("results/plots/ATAC_RNA/SE_TE_fano_boxplot.svg", height = 7, width = 10)





