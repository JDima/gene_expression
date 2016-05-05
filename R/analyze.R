setwd("~/PycharmProjects/gene_expression/R")
source("common.R")

setwd("/Users/JDima/PycharmProjects/gene_expression/src")

data <- readData("model_output/stats/means.txt", "means_head.txt")


plot.ts(getEveryKniNRow(data, 25), y = 1:40 * 25)
plot.ts(getEveryRNANRow(data, 25), y = 1:40 * 25)

plotSliceEveryTfNRow(data, 250, "kni")
plotSliceEveryRNANRow(data, 250)


df <- as.numeric(getKniNRow(data, 500))
plot(df, type = "l")

old.par <- par(mfrow=c(2, 2))
df <- getEveryKniNRow(data, 252)
for(i in 1:nrow(df)) {
  row <- as.numeric(df[i,])
  plot(row, type = "l", ylab = "Сoncentration", xlab = "Core")
}
par(old.par)

