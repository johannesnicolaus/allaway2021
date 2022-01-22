library(tidyverse)


# calculate fano factor ---------------------------------------------------
counts <- read_csv("results/data/rna_preprocessing/baynorm_normalized.csv")

# sample cells randomly
cells <- colnames(counts)[-1]

# check the number of cells per condition
cells %>% str_extract(".$") %>% table()

# downsample to 3000 cells
set.seed(1)
cells <- c(cells %>% str_subset(".*1") %>% sample(size = 3000),
  cells %>% str_subset(".*2") %>% sample(size = 3000)
  )

counts <- counts[, c("Geneid", cells)]

# make long 
counts <- counts %>% pivot_longer(!Geneid)

# calculate mean and fano factor
counts <- counts %>% mutate(sample = str_extract(name, ".$"))

counts_calc <- counts %>% group_by(Geneid, sample) %>% summarise(fano = var(value)/mean(value), mean = mean(value))

counts_calc %>% write_csv("results/data/rna_preprocessing/fano_mean_counts.csv")


counts_calc



