suppressPackageStartupMessages(library(tidyverse))

AAs <- str_split("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

d <- tibble(pos = 1:21,
       WT = c(AAs, "*")) %>%
    group_by(pos, WT) %>%
    do(tibble(var = c(AAs, "-"))) %>%
    ungroup() %>%
    mutate(score = runif(n(), -3, 1) %>% signif(4),
           extra = str_c("a", sample(n())))

write_tsv(d, "data.txt")

d %>%
    select(extra, var, WT, score, pos) %>%
    write_tsv("data_order.txt")

d %>%
    select(extra, Mutant = var, wt = WT, fitness = score, Position = pos) %>%
    write_tsv("data_colnames.txt")

d %>%
    filter(WT != "A") %>%
    write_tsv("data_incomplete.txt")

d %>%
    slice(1:200) %>%
    write_tsv("data/data1.txt")

d %>%
    slice(201:n()) %>%
    write_tsv("data/data2.txt")
