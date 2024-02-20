setwd("C:/Users/davey/OneDrive/Desktop/UCSF Work/OneDrive - UCSF/modified_data/biomarker_corr")
library(corrplot)

df = read.csv("correlation_change4wk_primary.csv", row.names=1)
m = cor(df)
corrplot(m, order="hclust", tl.col="black", tl.cex = 0.5, addrect=10,
rect.col="black", title="Primary Group: 4 Week Biomarker Correlation")

df = read.csv("correlation_change20wk_primary.csv", row.names=1)
m = cor(df)
corrplot(m, order="hclust", tl.col="black", tl.cex = 0.5, addrect=10,
rect.col="black", title="Primary Group: 20 Week Biomarker Correlation")

df = read.csv("correlation_change4wk_both.csv", row.names=1)
m = cor(df)
corrplot(m, order="hclust", tl.col="black", tl.cex = 0.5, addrect=10,
rect.col="black", title="Primary and Secondary Groups: 4 Week Biomarker Correlation")

df = read.csv("correlation_change20wk_both.csv", row.names=1)
m = cor(df)
corrplot(m, order="hclust", tl.col="black", tl.cex = 0.5, addrect=10,
rect.col="black", title="Primary and Secondary Groups: 20 Week Biomarker Correlation")
