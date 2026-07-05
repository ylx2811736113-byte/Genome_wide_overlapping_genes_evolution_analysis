############################################################
# expression_analysis.R
#
# Description:
# This script performs comparative expression analyses
# between overlapping and non-overlapping genes,
# including:
#   1. Expression matrix preprocessing
#   2. Expression profile extraction
#   3. Differential expression comparison
#   4. Spearman correlation analysis
#   5. Linear regression analysis

############################################################

### Example: Human gene expression analysis
############################################################
# Load required R packages
############################################################

library(readr)
library(data.table)
library(limma)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(colorspace)

############################################################
## Read and normalize gene expression matrix
############################################################

# Input file
input_file <- "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct"

# Read the GCT file
expr <- read.table(
  input_file,
  sep = "\t",
  header = TRUE,
  skip = 2,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Preview the first few columns
head(expr[,1:5])

# Remove Ensembl version numbers
expr$Name <- sub("\\..*$", "", expr$Name)

# Identify expression matrix columns (starting from the third column)
expr_cols <- 3:ncol(expr)

# log2(TPM + 1)
expr[, expr_cols] <- log2(expr[, expr_cols] + 1)

# Save the results
write.csv(
  expr,
  file = "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/human_data_normalized.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("Finished!\n")
cat("Output file: human_data_normalized.csv\n")


############################################################
## Extract expression profiles for overlapping and non-overlapping genes
############################################################

### Extract expression profiles for overlapping genes
data_normalized <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/human_data_normalized.csv")
RNA_seq <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/overlapping_gene_pairs.csv")
data_normalized_df <- as.data.frame(data_normalized)
merged_data <- merge(RNA_seq, data_normalized_df, by = "Gene")
write.csv(merged_data,file="D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/merged_data_overlap_1.csv")

### Extract expression profiles for non-overlapping genes
## Random non-overlapping gene control group
# Extract non-overlapping genes
gene_info <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/Homo_sapiens_overlap_genes.csv")
all_gene_info <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/Homo_sapiens_extracted_genes.csv")
gene1_id <- gene_info["Gene1_ID"]
gene2_id <- gene_info["Gene2_ID"]
logical_vector <- all_gene_info$GeneID %in% gene1_id$Gene1_ID
all_gene_info <- all_gene_info[!logical_vector, ]
logical_vector <- all_gene_info$GeneID %in% gene2_id$Gene2_ID
all_gene_info <- all_gene_info[!logical_vector, ]
write.csv(all_gene_info,file="D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/all_gene_info.csv")

# Retrieve expression profiles for non-overlapping genes
RNA_seq <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/all_gene_info.csv")
data_normalized <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/human_data_normalized.csv")
data_normalized_df <- as.data.frame(data_normalized)
merged_data <- merge(RNA_seq, data_normalized_df, by = "Gene")
write.csv(merged_data,file="D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/merged_data_non_overlap.csv")

# Randomly sample the same number of non-overlapping genes as the overlapping gene set
tra_gene_info <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/merged_data_non_overlap.csv")
num_pairs = 1890
set.seed(123)
non_overlap_genes <- tra_gene_info %>%
  sample_n(num_pairs, replace = FALSE)

# Save the sampled genes to a new file
write.csv(non_overlap_genes, file = "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/non_overlap_genes.csv", row.names = FALSE)


## Adjacent non-overlapping gene control group
# Retrieve expression profiles for non-overlapping genes
data_normalized <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/human_data_normalized.csv")
RNA_seq <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/adjacent_gene_pairs.csv")
data_normalized_df <- as.data.frame(data_normalized)
merged_data <- merge(RNA_seq, data_normalized_df, by = "Gene")
write.csv(merged_data,file="D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/merged_data_ad_2.csv")

# Randomly sample the same number of non-overlapping genes as the overlapping gene set
tra_gene_info <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/control.csv")
set.seed(123)
sampled_data <- sample_n(tra_gene_info, 945)

# Save the sampled genes to a new file
write.table(sampled_data, file = "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/ad_control.csv", sep = ",", row.names = FALSE, col.names = TRUE)

###########################Comparison of expression levels between overlapping and non-overlapping gene pairs
group1 <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/RNA_seq3.csv")
group2 <- read_csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/non_overlap_genes.csv")
group1_mean_values <- mean(group1$Average)
print(group1_mean_values)
group2_mean_values <- mean(group2$Average)
print(group2_mean_values)

## Mann–Whitney U test
# Perform two-tailed Mann-Whitney U test
result <- wilcox.test(group1$Average, group2$Average)
print(result)

# Perform one-tailed Mann-Whitney U test with the hypothesis that group 1 is significantly smaller than group 2
result <- wilcox.test(group1$Average, group2$Average, alternative = "less")
print(result)

# Perform one-tailed Mann-Whitney U test with the hypothesis that group 1 is significantly larger than group 2
result <- wilcox.test(group1$Average, group2$Average, alternative = "greater")
print(result)

# Adjust p-values using the Benjamini–Hochberg method
set.seed(123)
# Perform pairwise Wilcoxon tests
p_values <- c(wilcox.test(group1$Average, group2$Average)$p.value)
# Print unadjusted p-values
print(p_values)
# Multiple testing correction
adjusted_p_values <- p.adjust(p_values, method = "BH")
# Print adjusted p-values
print(adjusted_p_values)

## Kruskal–Wallis test
sub_analyse_1 <- read.table("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/sub_overlap_convergent.csv", sep = ",", header = TRUE)
sub_analyse_2 <- read.table("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/sub_overlap_co.oriented.csv", sep = ",", header = TRUE)
sub_analyse_3 <- read.table("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/sub_overlap_divergent.csv", sep = ",", header = TRUE)
sub_analyse_4 <- read.table("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/sub_overlap_nested.csv", sep = ",", header = TRUE)
group1 <- sub_analyse_1$Convergent
group2 <- sub_analyse_2$Co.oriented
group3 <- sub_analyse_3$Divergent
group4 <- sub_analyse_4$Nested
data <- c(group1, group2, group3, group4)

# Create a grouping factor (group sizes may differ)
group <- factor(c(rep("group1", length(group1)), 
                  rep("group2", length(group2)), 
                  rep("group3", length(group3)),
                  rep("group4", length(group4))))

# Perform Kruskal-Wallis test
result <- kruskal.test(data ~ group)
# Display the results
print(result)
# Perform pairwise comparisons
result <- wilcox.test(group3, group4)
print(result)

# Adjust p-values using the Benjamini–Hochberg method
# Perform pairwise Wilcoxon tests
p_values <- c(
  wilcox.test(group1, group2)$p.value,
  wilcox.test(group1, group3)$p.value,
  wilcox.test(group1, group4)$p.value,
  wilcox.test(group2, group3)$p.value,
  wilcox.test(group2, group4)$p.value,
  wilcox.test(group3, group4)$p.value
)

# Print unadjusted p-values
print(p_values)

# Multiple testing correction
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Print adjusted p-values
print(adjusted_p_values)


############################################################
## Spearman correlation and linear regression analyses
############################################################

### Spearman correlation analysis
# Read overlapping gene pair data
data_overlap <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/overlap_gene_expression_update.csv", header = TRUE, sep = ",")

# Read adjacent and random control gene pair datasets
data_ad_control <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/control.csv", header = TRUE, sep = ",")
data_random_control <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/non_overlap_genes_1.csv", header = TRUE, sep = ",")

# Add category labels for control groups
data_ad_control$Category <- "Neighboring_genes_control"
data_random_control$Category <- "Random_genes_control"

# Merge the datasets
data <- rbind(data_overlap, data_ad_control, data_random_control)

# Automatically identify columns ending with '_g1' and '_g2'
g1_cols <- grep("_g1$", colnames(data), value = TRUE)
g2_cols <- grep("_g2$", colnames(data), value = TRUE)

# Calculate Spearman correlation coefficients
spearman_results <- t(sapply(1:nrow(data), function(i) {
  g1_expr <- as.numeric(data[i, g1_cols])
  g2_expr <- as.numeric(data[i, g2_cols])
  
  if (all(is.na(g1_expr)) | all(is.na(g2_expr))) {
    return(c(correlation = NA, p_value = NA))
  }
  
  test_result <- cor.test(g1_expr, g2_expr, method = "spearman")
  return(c(correlation = test_result$estimate, p_value = test_result$p.value))
}))

# Convert results to a data frame and append gene pair information
spearman_df <- as.data.frame(spearman_results)
spearman_df$Gene1 <- data$Gene1
spearman_df$Gene2 <- data$Gene2
spearman_df$Category <- data$Category

# Save the results
write.csv(spearman_df, "spearman_gene_correlation.csv", row.names = FALSE)

# View the first few rows
head(spearman_df)

# Read Spearman correlation results
spearman_df <- read.csv("spearman_gene_correlation.csv")

# Ensure correlation values are numeric
spearman_df$correlation <- as.numeric(spearman_df$correlation)

# Remove missing values
spearman_df <- na.omit(spearman_df)

# Rename categories
spearman_df$Category <- recode(spearman_df$Category,
                               "Neighboring_genes_control" = "Neighboring genes control",
                               "Random_genes_control" = "Random genes control",
                               "Co-oriented" = "Co-oriented",
                               "Convergent" = "Convergent",
                               "Divergent" = "Divergent",
                               "Same_direction_nested" = "Same direction nested",
                               "Reverse_direction_nested" = "Reverse direction nested"
)

# Reset the category order to ensure that "Control" appears first
spearman_df$Category <- factor(spearman_df$Category, levels = c("Neighboring genes control","Random genes control", "Co-oriented", "Convergent", "Divergent", "Same direction nested", "Reverse direction nested"))

# draw a boxplot
custom_colors <- c(
  "Neighboring genes control" = "#BEBADA",
  "Random genes control" = "#80B1D3",
  "Co-oriented" = "#8DD3C7",
  "Convergent" = "#FDB462",
  "Divergent" = "#FFED6F",
  "Same direction nested" = "#FCCDE5",
  "Reverse direction nested" =  "#D9D9D9"
)

deeper_colors <- darken(custom_colors, amount = 0.3)

p <- ggplot(spearman_df, aes(x = Category, y = correlation, fill = Category)) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.7,
    outlier.shape = NA,
    color = NA
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2, aes(color = Category)) +
  geom_boxplot(
    width = 0.5,
    alpha = 0,
    aes(color = Category), 
    outlier.shape = NA
  ) +
  stat_summary(
    fun.data = function(y) {
      stats <- fivenum(y)
      outliers <- y[y < stats[2] - 1.5 * IQR(y) | y > stats[4] + 1.5 * IQR(y)]
      data.frame(y = outliers)
    },
    geom = "point",
    aes(color = Category),
    shape = 21,
    size = 1.5,
    alpha = 0.7
  ) +
  scale_fill_manual(values = deeper_colors) +
  scale_color_manual(values = deeper_colors) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

# Add significance annotations
p + stat_compare_means(
  comparisons = list(
    c("Neighboring genes control", "Random genes control"),
    c("Neighboring genes control", "Co-oriented"),
    c("Neighboring genes control", "Convergent"),
    c("Neighboring genes control", "Divergent"),
    c("Neighboring genes control", "Same direction nested"),
    c("Neighboring genes control", "Reverse direction nested")
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  hide.ns = FALSE,
  bracket.size = 1
)

ggsave("spearman.png", width = 8, height = 10, dpi = 300)

### Linear regression analysis
# Read input datasets
df_overlap <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/overlap_gene_expression_update.csv", header = TRUE)
df_ad <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/control.csv", header = TRUE)
df_random <- read.csv("D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/non_overlap_genes_1.csv", header = TRUE)

df_ad$Category <- "Neighboring_genes_control"
df_random$Category <- "Random_genes_control"

# Merge the datasets
df <- rbind(df_overlap, df_ad, df_random)

# Extract Gene1 and Gene2 expression columns
gene1_cols <- grep("_g1$", colnames(df), value = TRUE)
gene2_cols <- grep("_g2$", colnames(df), value = TRUE)

# Create a data frame to store regression results
results <- data.frame(Gene1 = character(),
                      Gene2 = character(),
                      Category = character(),
                      Slope = numeric(),
                      Intercept = numeric(),
                      R_squared = numeric(),
                      P_value = numeric(),
                      stringsAsFactors = FALSE)

# Fit a linear regression model for each gene pair
for (i in 1:nrow(df)) {
  gene1_expr <- as.numeric(df[i, gene1_cols])
  gene2_expr <- as.numeric(df[i, gene2_cols])
  
  if (length(gene1_expr) > 0 & length(gene2_expr) > 0) {
    skip <- FALSE
    tryCatch({
      model <- lm(gene2_expr ~ gene1_expr)
      summary_model <- summary(model)
      results <- rbind(results, data.frame(
        Gene1 = df$Gene1[i],
        Gene2 = df$Gene2[i],
        Category = df$Category[i],
        Slope = coef(model)[2],
        Intercept = coef(model)[1],
        R_squared = summary_model$r.squared,
        P_value = summary_model$coefficients[2, 4]
      ))
    }, error = function(e) {
      warning(paste("Model fitting failed at the", i, "row"))
      skip <- TRUE
    })
    if (skip) {
      next
    }
  }
}

# Save regression results
write.csv(results, "linear_regression_results.csv", row.names = FALSE)

# Read regression results
results <- read.csv("linear_regression_results.csv", header = TRUE)

# Rename categories
results$Category <- recode(results$Category,
                           "Neighboring_genes_control" = "Neighboring genes control",
                           "Random_genes_control" = "Random genes control",
                           "Co-oriented" = "Co-oriented",
                           "Convergent" = "Convergent",
                           "Divergent" = "Divergent",
                           "Same_direction_nested" = "Same direction nested",
                           "Reverse_direction_nested" = "Reverse direction nested"
)

results$Category <- factor(results$Category, 
                          levels = c("Neighboring genes control","Random genes control", "Co-oriented", "Convergent", "Divergent", "Same direction nested", "Reverse direction nested"))

# draw a boxplot
custom_colors <- c(
  "Neighboring genes control" = "#BEBADA",
  "Random genes control" = "#80B1D3",
  "Co-oriented" = "#8DD3C7",
  "Convergent" = "#FDB462",
  "Divergent" = "#FFED6F",
  "Same direction nested" = "#FCCDE5",
  "Reverse direction nested" =  "#D9D9D9"
)

deeper_colors <- darken(custom_colors, amount = 0.3)

p <- ggplot(results, aes(x = Category, y = R_squared, fill = Category)) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.7,
    outlier.shape = NA,
    color = NA
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2, aes(color = Category)) +
  geom_boxplot(
    width = 0.5,
    alpha = 0,
    aes(color = Category),
    outlier.shape = NA
  ) +
  stat_summary(
    fun.data = function(y) {
      stats <- fivenum(y)
      outliers <- y[y < stats[2] - 1.5 * IQR(y) | y > stats[4] + 1.5 * IQR(y)]
      data.frame(y = outliers)
    },
    geom = "point",
    aes(color = Category),
    shape = 21,
    size = 1.5,
    alpha = 0.7
  ) +
  scale_fill_manual(values = deeper_colors) +
  scale_color_manual(values = deeper_colors) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

# Add significance annotations
p + stat_compare_means(
  comparisons = list(
    c("Neighboring genes control", "Random genes control"),
    c("Neighboring genes control", "Co-oriented"),
    c("Neighboring genes control", "Convergent"),
    c("Neighboring genes control", "Divergent"),
    c("Neighboring genes control", "Same direction nested"),
    c("Neighboring genes control", "Reverse direction nested")
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  hide.ns = FALSE,
  bracket.size = 1
)

ggsave("linear_regression.png", width = 8, height = 10, dpi = 300)


















