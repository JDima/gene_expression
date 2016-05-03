source("common.R")

setwd("/Users/JDima/PycharmProjects/gene_expression/src")

data <- readData("means_head.txt", "model_output/stats/means.txt")

plot.ts(getEveryKniNRow(data, 25))

plotSliceEveryTfNRow(data, 250, "kni")


df <- as.numeric(getKniNRow(data, 500))
plot(df, type = "l")

old.par <- par(mfrow=c(2, 2))
df <- getEveryKniNRow(data, 252)
for(i in 1:nrow(df)) {
  row <- as.numeric(df[i,])
  plot(row, type = "l", ylab = "Сoncentration", xlab = "Core")
}
par(old.par)

