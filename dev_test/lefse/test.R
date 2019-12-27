var <- "SampleType"
taxa_rank <- "OTU"
mbset <- MicrobiomeAnalystR:::.get.mbSetObj(mbSet)
ps <- mbset$dataSet$norm.phyobj
class <- sample_data(ps)[[var]]
otus <- otu_table(ps)

if (taxa_are_rows(otus)) {
  otus <- t(otus)
}

otus <- as_tibble(otus@.Data, rownames = NA)
kw_p <- map_dbl(otus, ~ kruskal.test(.x, class)$p.value)
kw_fdr <- p.adjust(kw_p, method = "fdr")

# lda analysis
lda_res <- MASS::lda(class~., data = otus)
lda_mean <- t(lda_res$means) %>%
  as_tibble(rownames = "otu")
lda_mean <- rowwise(lda_mean) %>%
  mutate(max = max(c(gut, `left palm`, `right palm`, tongue)),
    min = min(c(gut, `left palm`, `right palm`, tongue)))
lda_mean <- mutate(lda_mean,
  lda_score = signif(log10(1 + abs(max - min)/2), digits = 3))
lda_mean$p_value <- kw_p
lda_mean$fdr <- kw_fdr


class_unique <-  c("gut", "left palm", "right palm", "tongue")
sig_feature <-
  mutate(sig_feature,
    enrich = class_unique[which.max(c(gut, `left palm`, `right palm`, tongue))])


if (taxa_are_rows(otus)) {
  otus <- t(otus)
}
otus <- as_tibble(otus@.Data, rownames = NA)

## kw rank sum test
set.seed(56290)
# kw_p <- apply(otus, 2,
#   function(x) kruskal.test(x, class)$p.value
# )
kw_p <- purrr::map_dbl(otus, ~ kruskal.test(.x, class)$p.value)
# ord_indx <- order(kw_p)
# kw_p <- kw_p[ord_indx]
# a <- sort(kw_p)
kw_fdr <- p.adjust(kw_p, method = "fdr")
# otus <- otus[ord_indx, ]

# otus_tb <- as_tibble(otus, rownames = NA)

# wilcoxon rank sum test is not preformed if there is no subclass

# lda analysis
lda_res <- MASS::lda(class~., data = otus)
lda_mean <- t(lda_res$means) %>%
  tibble::as_tibble(rownames = "otu")
lda_mean <- rowwise(lda_mean) %>%
  mutate(max = max(c(gut, `left palm`, `right palm`, tongue)),
    min = min(c(gut, `left palm`, `right palm`, tongue)))
lda_mean <- mutate(lda_mean,
  lda_score = signif(log10(1 + abs(max - min)/2), digits = 3))
lda_mean$p_value <- kw_p
lda_mean$fdr <- kw_fdr

# by default, as_tibble add ` around invalid names, and many otu is repersentd as
# a nubmer
lda_mean <- mutate(lda_mean, otu = gsub("`", "", otu))

# significant feature
sig_feature <- filter(lda_mean, fdr < 0.05, lda_score > 2)
enrich <-

  if (nrow(sig_feature) == 0) {
    warning("No significant featrue were identified")
  }

# order result by p value, fdr and lda_score
class_unique <-  c("gut", "left palm", "right palm", "tongue")
sig_feature <- arrange(sig_feature, desc(lda_score))
sig_feature <-
  mutate(sig_feature,
    enrich = class_unique[which.max(c(gut, `left palm`, `right palm`, tongue))])

top_sig_feature <- slice(sig_feature, 1:15) %>%
  arrange(desc(lda_score))

# 两种方法修改离散坐标的顺序
# bar plot
ggplot(a[1:15, ], aes(otu, lda_score, fill = enrich_group)) +
  geom_col() +
  labs(x = "Features", y = "LDA score", fill = NULL) +
  scale_x_discrete(limits = rev(a[1:15, ]$otu)) +
  coord_flip() +
  theme_bw()

ggplot(top_sig_feature, aes(fct_reorder(otu, lda_score), lda_score, fill = enrich)) +
  geom_col() +
  labs(x = "Features", y = "LDA score", fill = NULL) +
  coord_flip() +
  theme_bw()

# online
library(MicrobiomeAnalystR)
mbSet<-Init.mbSetObj()
mbSet<-SetModuleType(mbSet, "mdp")
mbSet<-Read16SAbundData(mbSet, "data/treebiom/otu_table_mc2_w_tax_no_pynast_failures.biom","biom","GreengenesID","F");
mbSet<-ReadSampleTable(mbSet, "data/treebiom/map.txt");
mbSet<-ReadTreeFile(mbSet, "data/treebiom/rep_set.tre");
mbSet<-SanityCheckData(mbSet, "biom");
mbSet<-PlotLibSizeView(mbSet, "norm_libsizes_0","png");
mbSet<-CreatePhyloseqObj(mbSet, "biom","GreengenesID","F")
mbSet<-ApplyAbundanceFilter(mbSet, "prevalence", 4, 0.2);
mbSet<-ApplyVarianceFilter(mbSet, "iqr", 0.1);
mbSet<-PerformNormalization(mbSet, "none", "colsum", "none");
mbSet<-PerformLefseAnal(mbSet, 0.1, "fdr", 2.0, "SampleType","F","NA","OTU");
mbSet<-PlotLEfSeSummary(mbSet, 15, "dot", "bar_graph_0","png");
mbSet<-PlotLEfSeSummary(mbSet, 15, "bar", "bar_graph_1","png");
# mbSet<-PerformLefseAnal(mbSet, 0.05, "fdr", 2.0, "SampleType","F","NA","OTU");
# mbSet<-PlotLEfSeSummary(mbSet, 15, "bar", "bar_graph_2","png")
#
# export to lefse
otus_lefse <- t(otus)
lefse_in <- rbind(as.character(class) %>% gsub(" ", "", .),
  otus_lefse
)
row.names(lefse_in) <- c("otu", row.names(otus_lefse))
write.table(lefse_in, file = "data/lefse_in.tsv",
  col.names = FALSE, quote = FALSE, sep = "\t")

lefse_in <- as_tibble(lefse_in, rownames = "otu")
write_tsv(lefse_in, path = "data/lefse_in.tsv", col_names = FALSE)

# p 值计算为`_`??, 会带来错误
lefse_res <- read_tsv("data/lefse_lda.res", col_names = FALSE)
lefse_res <- mutate(lefse_res, X3 = gsub(" ", "", X3))
lefse_res <- filter(lefse_res, X5 != "-")
write_tsv(lefse_res, "data/lefse_lda2.res", col_names = FALSE)

lefse_out <- lefse(caporaso_phyloseq, class = "SampleType")

lefse_out_arranged <- arrange(lefse_out,enrich_group, desc(lda_score)) %>%
  mutate(otu = factor(otu, levels = otu))

ggplot(lefse_out_arranged, aes(otu, lda_score, fill = enrich_group)) +
  geom_col() +
  labs(x = "Features", y = "LDA score", fill = NULL) +
  scale_x_discrete(limits = rev(lefse_out_arranged$otu)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "LDA SCORE (log10)") +
  coord_flip() +
  theme_bw()

microbiomeViz::lefse.lda.plot("lefse/lefse_lda2.res", negate.class = "gut", lda.threshold = 2)

# lefse cladogram
library(microbiomeViz)
data("SRS014459_Stool_profile")
tr <- parseMetaphlanTSV(SRS014459_Stool_profile)

df <- read.table("http://bailab.genetics.ac.cn/markdown/R/microbiomeViz/merged_abundance_table.txt", head=TRUE, stringsAsFactors = FALSE)
dat <- data.frame(V1=df[,1], V2=rowMeans(df[,-1]), stringsAsFactors = FALSE)
tr <- parseMetaphlanTSV(dat, node.size.offset=2, node.size.scale=0.2)
p <- tree.backbone(tr, size=0.5)

lefse_lists = data.frame(node=c('s__Haemophilus_parainfluenzae','p__Proteobacteria',
  'f__Veillonellaceae','o__Selenomonadales',
  'c__Negativicutes', 's__Streptococcus_parasanguinis',

  'p__Firmicutes','f__Streptococcaceae',
  'g__Streptococcus','o__Lactobacillales',
  'c__Bacilli','s__Streptococcus_mitis'),
  color=c(rep('darkgreen',6), rep('red','6')),
  stringsAsFactors = FALSE
)
clade.anno(p, lefse_lists, alpha=0.3)

tr <- parsePhyloseq(caporaso_phyloseq, use_abundance = FALSE)
p <- tree.backbone(tr, size=0.5, )


# phyloseq2lefse ----------------------------------------------------------
otus <- otu_table(caporaso_phyloseq) %>% .@.Data
otu_id <- row.names(otus) %>% tibble(otu_id = .)
otus <- as_tibble(otus) %>%
  bind_cols(otu_id, .)
ps_meta <- sample_data(caporaso_phyloseq)
cls <- ps_meta$SampleType
all(row.names(ps_meta) == names(otus)[-1])
names(otus) <- c("otu_id", as.character(cls))
ps_otu_anno <-
# otus_lefse <- t(otus)
#lefse_in <- rbind(as.character(class) %>% gsub(" ", "", .),
  otus_lefse
)
row.names(lefse_in) <- c("otu", row.names(otus_lefse))
write.table(lefse_in, file = "data/lefse_in.tsv",
  col.names = FALSE, quote = FALSE, sep = "\t")

lefse_in <- as_tibble(lefse_in, rownames = "otu")
write_tsv(lefse_in, path = "data/lefse_in.tsv", col_names = FALSE)

