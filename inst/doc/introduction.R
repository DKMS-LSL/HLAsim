## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE, collapse = TRUE, comment = "#>")
library(HLAsim)
library(dplyr)
library(ggplot2)

## ----load1---------------------------------------------------------------
data(list = c("drb1", "dqb1", "dpb1"), package = "HLAsim")
drb1

## ----load2---------------------------------------------------------------
data(list = c("drb1_eag1412", "dqb1_eag1412", "dpb1_eag1412"), package = "HLAsim")
drb1_eag1412

## ----make_samplers, cache=TRUE, cache.lazy=FALSE, cache.comments=FALSE----
sample_drb_de <- make_genotype_sampler(drb1[provenance == "DE"], drb1_eag1412)
sample_drb_pl <- make_genotype_sampler(drb1[provenance == "PL"], drb1_eag1412)
sample_drb_uk <- make_genotype_sampler(drb1[provenance == "UK"], drb1_eag1412)

sample_dqb_de <- make_genotype_sampler(dqb1[provenance == "DE"], dqb1_eag1412)
sample_dqb_pl <- make_genotype_sampler(dqb1[provenance == "PL"], dqb1_eag1412)
sample_dqb_uk <- make_genotype_sampler(dqb1[provenance == "UK"], dqb1_eag1412)

sample_dpb_de <- make_genotype_sampler(dpb1[provenance == "DE"], dpb1_eag1412)
sample_dpb_pl <- make_genotype_sampler(dpb1[provenance == "PL"], dpb1_eag1412)
sample_dpb_uk <- make_genotype_sampler(dpb1[provenance == "UK"], dpb1_eag1412)

## ----show_sampler--------------------------------------------------------
sample_drb_de

## ----load3---------------------------------------------------------------
data("concentration", package = "HLAsim")
concentration

## ----load4, fig.width=6, fig.height=4, fig.align='center'----------------
data("cdfA", package = "HLAsim")
plot(cdfA)

## ----settings------------------------------------------------------------
# Total number of random samples
n <- 1e6
# Binwidth for DNA concentrations
bin_width <- 2.5

## ----sample, cache=TRUE--------------------------------------------------
## DRB1
gtbl_drb1_de <- sample_drb_de(concentration[provenance == "DE"], n, bin_width)
gtbl_drb1_pl <- sample_drb_pl(concentration[provenance == "PL"], n, bin_width)
gtbl_drb1_uk <- sample_drb_uk(concentration[provenance == "UK"], n, bin_width)

## DQB1
gtbl_dqb1_de <- sample_dqb_de(concentration[provenance == "DE"], n, bin_width)
gtbl_dqb1_pl <- sample_dqb_pl(concentration[provenance == "PL"], n, bin_width)
gtbl_dqb1_uk <- sample_dqb_uk(concentration[provenance == "UK"], n, bin_width)

## DPB1
gtbl_dpb1_de <- sample_dpb_de(concentration[provenance == "DE"], n, bin_width)
gtbl_dpb1_pl <- sample_dpb_pl(concentration[provenance == "PL"], n, bin_width)
gtbl_dpb1_uk <- sample_dpb_uk(concentration[provenance == "UK"], n, bin_width)

gtbl_drb1_de

## ----inject_error, cache=TRUE--------------------------------------------
## DRB1
err_tbl_drb1_de <- inject_errors(gtbl_drb1_de, cdfA, perr = 0.01, odds = 0.25)
err_tbl_drb1_pl <- inject_errors(gtbl_drb1_pl, cdfA, perr = 0.01, odds = 0.25)
err_tbl_drb1_uk <- inject_errors(gtbl_drb1_uk, cdfA, perr = 0.01, odds = 0.25)

## DQB1
err_tbl_dqb1_de <- inject_errors(gtbl_dqb1_de, cdfA, perr = 0.01, odds = 1)
err_tbl_dqb1_pl <- inject_errors(gtbl_dqb1_pl, cdfA, perr = 0.01, odds = 1)
err_tbl_dqb1_uk <- inject_errors(gtbl_dqb1_uk, cdfA, perr = 0.01, odds = 1)

## DPB1
err_tbl_dpb1_de <- inject_errors(gtbl_dpb1_de, cdfA, perr = 0.01, odds = 0.25)
err_tbl_dpb1_pl <- inject_errors(gtbl_dpb1_pl, cdfA, perr = 0.01, odds = 0.25)
err_tbl_dpb1_uk <- inject_errors(gtbl_dpb1_uk, cdfA, perr = 0.01, odds = 0.25)

err_tbl_drb1_de

## ----analysis, fig.width=7.5, fig.height=4-------------------------------
library(foreach)
library(iterators)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

summarise_gtbl <- function(x) {
  tbl <- summary(x) %>%
    dplyr::select(-pHet, -nerr) %>%
    tidyr::gather(type, proportion, pHom:pNA) %>%
    dplyr::tbl_df()
  levels(tbl$type) <- c("erroneous HM", "unchanged HT", "erroneous HT", "incomp. exons")
  tbl
}

tbls <- iter(ls(pattern = "^err_tbl_.+_.+"))
tbls.summary <- tbl_dt(foreach(tbl = tbls, .combine = "rbind") %do% {
  gene <- toupper(strsplit(tbl, "_")[[1]][3])
  provenance <- toupper(strsplit(tbl, "_")[[1]][4])
  summarise_gtbl(get(tbl)) %>%
    mutate(provenance = provenance, gene = gene)
})

ggplot(tbls.summary, aes(x = type, y = proportion, fill = provenance)) +
  geom_bar(position = "dodge", stat = "identity", colour = "white") +
  facet_grid(. ~ gene) +
  geom_text(aes(label = ifelse(proportion > 0.01, round(proportion, 2), ""),
                hjust = ifelse(proportion > 0.2, 1.2, -0.2)),
            angle = 90, colour = "black", position = position_dodge(0.9), size = 2.5) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(x = "Genotyping error category", y = "Proportion", fill = "Provenance") +
  theme(axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 1),
        axis.text.x = element_text(size = 8, angle = 30, vjust = 1, hjust = 1)) +
  theme_hc() + 
  theme(legend.position = "top")

## ------------------------------------------------------------------------
gtf_table(sample_drb_de)

## ------------------------------------------------------------------------
drb1.01.de <- map_genotype_frequencies(err_tbl_drb1_de) %>%
  dplyr::mutate(provenance = "DE")
drb1.01.de %>% 
  group_by(Class) %>% 
  sample_n(5)

