# draw scatter plot of predicted indel freq against the experimentally measured indel freq

library(openxlsx)
library(dplyr)
library(ggplot2)

df <- read.xlsx('Evalulation_Brien2023.xlsx')

# ================ Brien2023 =================
df.comp.brien <- select(df, c(pred_Brien2023, indel_freq))

p.brien <- ggplot(data = df.comp.brien, aes(x = pred_Brien2023, y = indel_freq)) +
  geom_point()

# Add a linear regression line to the plot
p.brien <- p.brien + geom_smooth(method = "lm", se = FALSE, color = "red")

# Calculate Spearman correlation and p-value
spearman.corr<- cor.test(df.comp.brien$pred_Brien2023, df.comp.brien$indel_freq, method = "spearman")

spearman.brien <- spearman.corr$estimate
pvalue.brien <- spearman.corr$p.value


# Add the correlation and p-value to the plot
p.brien <- p.brien + labs(title = paste0("Spearman's r = ", spearman.brien, ", p = ", pvalue.brien))

# Display the plot
pdf("scatter_brien2023.pdf")
print(p.brien)
dev.off()