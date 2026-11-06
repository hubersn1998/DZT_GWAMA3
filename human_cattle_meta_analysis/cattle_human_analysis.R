### Import libraries ----
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ashr)
library(ggrepel)
library(hexbin)  



### Define functions ----
read_magma <- function(path, tag) {
  # Read a MAGMA .genes.out file and suffix ZSTAT / P with the dataset tag
  fread(path, data.table = FALSE) %>%
    as_tibble() %>%
    select(GENE, CHR, START, STOP, ZSTAT, P) %>%
    rename_with(~ paste0(sub("ZSTAT", "Z", .x), "_", tag), c(ZSTAT, P))
}

calibrate_z <- function(z, name) {
  # MAD-scale Z-scores to adjust for genomic inflation before Bayesian shrinkage
  med_z <- median(z)
  mad_z <- mad(z, center = med_z, constant = 1)
  
  if (!is.finite(mad_z) || mad_z <= 0) {
    stop(sprintf("Invalid MAD for %s: %s", name, mad_z))
  }
  
  null_sd <- mad_z / qnorm(0.75)
  z_cal <- (z - med_z) / null_sd
  lambda_gc <- median(z_cal^2) / qchisq(0.5, df = 1)
  
  cat(sprintf(
    "\n=== %s ===\nn genes: %d\nmedian(Z): %.4f\nMAD: %.4f\nestimated null SD: %.4f\nmean calibrated Z: %.4f\n",
    name, length(z), med_z, mad_z, null_sd, lambda_gc, mean(z_cal)
  ))
  
  list(z = z_cal, median = med_z, null_sd = null_sd, lambda_gc = lambda_gc)
}

run_ash <- function(z, s = 1) {
  # Adaptive shrinkage of Z-scores (non-finite values are returned as NA)
  ok <- is.finite(z)
  z_shrunk <- rep(NA_real_, length(z))
  res <- ash(betahat = z[ok], sebetahat = rep(s, sum(ok)), mixcompdist = "normal")
  z_shrunk[ok] <- get_pm(res)
  z_shrunk
}

plot_pair <- function(z1, z2, xlab, ylab, title = NULL) {
  df <- data.frame(Z1 = z1, Z2 = z2)
  df <- df[complete.cases(df), ]
  
  cor_test <- cor.test(df$Z1, df$Z2, method = "pearson")
  
  ggplot(df, aes(Z1, Z2)) +
    geom_hex(bins = 60) +
    scale_fill_viridis_c(trans = "log10", name = "count") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "red", linewidth = 0.5) +
    labs(
      title = title,
      subtitle = sprintf("Pearson r = %.3f, p = %.3g",
                         unname(cor_test$estimate), cor_test$p.value),
      x = xlab, y = ylab
    ) +
    coord_equal() +
    theme_minimal()
}

stouffer_meta <- function(df, w_irish_bulls, w_german_dams) {
  # Stage 1: weighted Stouffer across the two cattle datasets -> Z_meta_cattle
  # Stage 2: unweighted Stouffer of Z_meta_cattle and Z_human_meta -> Z_meta
  # Only genes with a human Z and at least one cattle Z are returned.
  valid_ib <- !is.na(df$Z_irish_bulls)
  valid_gd <- !is.na(df$Z_german_dams)
  
  w_ib <- ifelse(valid_ib, w_irish_bulls, 0)
  w_gd <- ifelse(valid_gd, w_german_dams, 0)
  z_ib <- ifelse(valid_ib, df$Z_irish_bulls, 0)
  z_gd <- ifelse(valid_gd, df$Z_german_dams, 0)
  
  df$Z_meta_cattle <- (w_ib * z_ib + w_gd * z_gd) / sqrt(w_ib^2 + w_gd^2)
  df$Z_meta_cattle[!valid_ib & !valid_gd] <- NA_real_
  df$p_meta_cattle <- pnorm(df$Z_meta_cattle, lower.tail = FALSE)
  
  df %>%
    filter(!is.na(Z_human_meta), !is.na(Z_meta_cattle)) %>%
    mutate(
      Z_meta = (Z_meta_cattle + Z_human_meta) / sqrt(2),
      p_meta = pnorm(Z_meta, lower.tail = FALSE)
    )
}



### Load files ----
human_meta <- read_magma("human_DZT.genes.out", "human_meta") %>%
  select(GENE, ends_with("human_meta"), CHR, START, STOP)

german_dams <- read_magma("german_dams_twinning.genes.out", "german_dams") %>%
  select(GENE, ends_with("german_dams"),
         cattle_CHR = CHR, cattle_START = START, cattle_STOP = STOP)

irish_bulls <- read_magma("irish_bulls_twinning.genes.out", "irish_bulls") %>%
  select(GENE, ends_with("irish_bulls"))

orthologs <- fread(
  "Gene_orthologs.txt",
  select = c("Gene stable ID", "Cattle gene stable ID", "Cattle gene name",
             "Cattle orthology confidence [0 low, 1 high]"),
  data.table = FALSE
) %>%
  as_tibble() %>%
  rename(
    human_gene = `Gene stable ID`,
    cattle_gene = `Cattle gene stable ID`,
    cattle_symbol = `Cattle gene name`,
    orthology_conf = `Cattle orthology confidence [0 low, 1 high]`
  ) %>%
  distinct() %>%
  filter(orthology_conf == 1)

### Merge cattle and human gene-level stats on orthologs ----
Z_stats <- orthologs %>%
  left_join(irish_bulls, by = c("cattle_gene" = "GENE")) %>%
  left_join(german_dams, by = c("cattle_gene" = "GENE")) %>%
  left_join(human_meta,  by = c("human_gene" = "GENE")) %>%
  mutate(CHR = suppressWarnings(as.numeric(ifelse(CHR == "X", "23", CHR)))) %>%
  filter(!is.na(CHR),
         !(is.na(Z_irish_bulls) & is.na(Z_german_dams)))

### Calibrate and shrink Z-scores for correlation plots ----
# Genes with finite Z in all three datasets
Z_stats_cc <- Z_stats %>%
  filter(is.finite(Z_irish_bulls), is.finite(Z_german_dams), is.finite(Z_human_meta))

dataset_labels <- c(
  irish_bulls = "Irish Bulls",
  german_dams = "German Dams",
  human_meta  = "Human meta-analysis"
)

# MAD calibration, then Bayesian adaptive shrinkage
# MAD calibration is only done so data matches expected shrinkage distribution
shrunk <- lapply(names(dataset_labels), function(d) {
  cal <- calibrate_z(Z_stats_cc[[paste0("Z_", d)]], paste0("Z_", d))
  run_ash(cal$z, s = 1)
})
names(shrunk) <- names(dataset_labels)

# Pairwise correlations and plots
pair_defs <- list(
  c("irish_bulls", "german_dams"),
  c("irish_bulls", "human_meta"),
  c("german_dams", "human_meta")
)

pair_cors <- sapply(pair_defs, function(p) {
  cor(shrunk[[p[1]]], shrunk[[p[2]]], method = "pearson",
      use = "pairwise.complete.obs")
})
names(pair_cors) <- sapply(pair_defs, paste, collapse = "_vs_")
print(pair_cors)

pair_plots <- lapply(pair_defs, function(p) {
  plot_pair(shrunk[[p[1]]], shrunk[[p[2]]],
            dataset_labels[[p[1]]], dataset_labels[[p[2]]])
})

p_z <- wrap_plots(pair_plots, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
print(p_z)

### Calculate attenuation ratio ----
# ratio of mean p-values between cattle sumstats (variants with p<5e-8)
# Needs access to full sumstats
# Rename columns internally
# gwas1 <- gwas1[, .(
#   SNP = get('SNP'),
#   Z1 = get('Z'),
#   P1 = get('P')
# )]
# 
# gwas2 <- gwas2[, .(
#   SNP = get('SNP'),
#   Z2 = get('Z'),
#   P2 = get('P')
# )]
# 
# # Merge
# dat <- merge(gwas1, gwas2, by="SNP")
# alpha <- 5e-8
# threshold <- alpha 
# hits <- dat[
#   P1 < threshold & P2 < threshold
# ]
# hits[, attenuation := abs(Z2) / abs(Z1)]


### Meta-analyze ----
meta_analysis <- stouffer_meta(Z_stats, weight_irish_bulls, weight_german_dams)
write.table(meta_analysis, 'human_cattle_meta_analysis.txt', sep = "\t", row.names = FALSE, quote = FALSE)

### Manhattan plot ----
manhattan_df <- meta_analysis %>%
  filter(is.finite(p_meta), !is.na(START), !is.na(STOP)) %>%
  mutate(
    pos = (START + STOP) / 2,
    logP = -log10(p_meta),
    SYMBOL = cattle_symbol,
    fdr = p.adjust(p_meta, method = "BH")
  ) %>%
  arrange(CHR, pos)

# Largest p-value passing BH FDR < 0.05 (NA if nothing passes)
sig <- manhattan_df$fdr < 0.05
threshold <- if (any(sig)) max(manhattan_df$p_meta[sig]) else NA_real_

# Cumulative genomic position
chr_offsets <- manhattan_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = cumsum(chr_len) - chr_len)

manhattan_df <- manhattan_df %>%
  left_join(chr_offsets, by = "CHR") %>%
  mutate(pos_cum = pos + offset)

axis_set <- manhattan_df %>%
  group_by(CHR) %>%
  summarise(center = mean(pos_cum), .groups = "drop")

# Split points once, handling the case where no threshold is defined
if (!is.na(threshold)) {
  sig_df    <- filter(manhattan_df, p_meta <= threshold)
  nonsig_df <- filter(manhattan_df, p_meta >  threshold)
} else {
  sig_df    <- manhattan_df[0, ]
  nonsig_df <- manhattan_df
}

manhattan <- ggplot() +
  # Non-significant points (alternating chromosome colors)
  geom_point(
    data = nonsig_df,
    aes(x = pos_cum, y = logP, color = as.factor(CHR %% 2)),
    size = 0.8
  ) +
  # Significant points (green)
  geom_point(
    data = sig_df,
    aes(x = pos_cum, y = logP),
    color = "green3",
    size = 0.8
  ) +
  scale_color_manual(values = c("blue3", "cornflowerblue")) +
  scale_x_continuous(
    breaks = axis_set$center,
    labels = ifelse(axis_set$CHR == 23, "X", axis_set$CHR)
  ) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

# Threshold line and gene labels only if a threshold exists
if (!is.na(threshold)) {
  manhattan <- manhattan +
    geom_hline(yintercept = -log10(threshold), color = "red", linewidth = 0.6) +
    geom_text_repel(
      data = sig_df,
      aes(x = pos_cum, y = logP, label = SYMBOL),
      color = "black", size = 3, nudge_y = 0.3, max.overlaps = 50
    )
}

print(manhattan)
ggsave(manhattan_file, manhattan, width = 12, height = 5, dpi = 300)

