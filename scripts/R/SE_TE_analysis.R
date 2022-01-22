
# analysis of ATAC + RNA --------------------------------------------------
library(tidyverse)

SE <- read_tsv("results/data/atac/annotated/merged_superEnhancers_annotated_intensity.txt")
TE <- read_tsv("results/data/atac/annotated/merged_typicalEnhancers_annotated_intensity.txt")

SE <- SE %>% rename(neg = 20, pos = 21)  
TE <- TE %>% rename(neg = 20, pos = 21)  

# calculate log2 fold change
SE <- SE %>% mutate(l2fc = log2((pos/neg)))
TE <- TE %>% mutate(l2fc = log2((pos/neg)))

# add rank
SE <- SE %>% arrange(desc(l2fc)) %>% mutate(N = seq.int(nrow(.)))
TE <- TE %>% arrange(desc(l2fc)) %>% mutate(N = seq.int(nrow(.)))

# gained and lost SE and TE graphs
dir.create("results/plots/ATAC_gained_lost")
ggplot(SE,aes(x = N , y = l2fc, fill = -N))+
  geom_bar(stat = "identity", width = 1)+
  geom_hline(yintercept = c(SE %>% pull(l2fc) %>% quantile(0.25),SE %>% pull(l2fc) %>% quantile(0.75)), linetype="dashed", size = 2)+
  coord_flip()+
  scale_x_reverse()+
  guides(fill = F)+
  cowplot::theme_cowplot(40)+
  cowplot::panel_border("black", size = 2)+
  scale_fill_viridis_c()+
  ylab("") +
  xlab("") +
  ylab(element_blank())+
  xlab(element_blank())

ggsave("results/plots/ATAC_gained_lost/SE.png", height = 7, width = 7)

ggplot(TE,aes(x = N , y = l2fc, fill = -N))+
  geom_bar(stat = "identity", width = 1)+
  geom_hline(yintercept = c(TE %>% pull(l2fc) %>% quantile(0.25),TE %>% pull(l2fc) %>% quantile(0.75)), linetype="dashed", size = 2)+
  coord_flip()+
  scale_x_reverse()+
  guides(fill = F)+
  cowplot::theme_cowplot(40)+
  cowplot::panel_border("black", size = 2)+
  scale_fill_viridis_c()+
  ylab(element_blank())+
  xlab(element_blank())

ggsave("results/plots/ATAC_gained_lost/TE.png", height = 7, width = 7)


# calculate TE SE gain and lost
SE <- SE %>% 
  mutate(atac_state = case_when(l2fc > l2fc %>% quantile(0.75) ~ "gained",
                                l2fc < l2fc %>% quantile(0.25) ~ "lost",
                                l2fc <= l2fc %>% quantile(0.75) & l2fc >= l2fc %>% quantile(0.25) ~ "unchanged"
  ))

TE <- TE %>% 
  mutate(atac_state = case_when(l2fc > l2fc %>% quantile(0.75) ~ "gained",
                                l2fc < l2fc %>% quantile(0.25) ~ "lost",
                                l2fc <= l2fc %>% quantile(0.75) & l2fc >= l2fc %>% quantile(0.25) ~ "unchanged"
  ))

write_csv(SE, "results/data/atac/annotated/SE.csv")
write_csv(TE, "results/data/atac/annotated/TE.csv")

SE$atac_state %>% table
TE$atac_state %>% table

