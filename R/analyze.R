setwd("/Users/JDima/PycharmProjects/gene_expression/src")
cols <- read.table("means_head.txt", sep = "\t", head = FALSE)
data <- read.table("model_output/stats/means.txt", sep = "\t", head = FALSE)
data <- data[, 1:length(cols)]
colnames(data) <- unlist(cols)

getNRow <- function(df, nrow)
{
  return(df[nrow, ])
}

getKniNRow <- function(df, nrow)
{
  return(df[nrow, grep("kni", colnames(df), perl = TRUE )])
}

getKNICount <- function(df)
{
  return(df[ , grep("kni", colnames(df), perl = TRUE )])
}

getEveryNRow <- function(df, n)
{
  maxT <- length(df[,1])
  nrows <- seq(1, maxT, by = n)
  return(df[nrows,])
}

getEveryKniNRow<- function(df, n)
{
  kniCounts <- getKNICount(df)
  return(getEveryNRow(kniCounts, n))
}

plot.ts(getEveryKniNRow(data, 25))

df <- as.numeric(getKniNRow(data, 500))
plot(df, type = "l")

old.par <- par(mfrow=c(2, 2))

df <- getEveryKniNRow(data, 252)

for(i in 1:nrow(df)) {
  row <- as.numeric(df[i,])
  plot(row, type = "l", ylab = "Сoncentration", xlab = "Core")
}
par(old.par)




df <- getEveryKniNRow(data, 250)
row <- as.numeric(df[1,])
plot(row, type = "l", ylab = "Сoncentration", xlab = "Core", xlim = c(1, length(df[1,])), ylim = c(0, max(df)))
cl <- rainbow(length(df[1,]))
for(i in 2:nrow(df)) {
  row <- as.numeric(df[i,])
  lines(row, type = "l", ylab = "Сoncentration", xlab = "Core",col = cl[i])
}








getRNANRow <- function(df, nrow)
{
  return(df[nrow, grep("mRNA", colnames(df), perl = TRUE )])
}

getRNACount <- function(df)
{
  return(df[ , grep("mRNA", colnames(df), perl = TRUE )])
}

getEveryNRow <- function(df, n)
{
  maxT <- length(df[,1])
  nrows <- seq(1, maxT, by = n)
  return(df[nrows,])
}

getEveryRNANRow<- function(df, n)
{
  kniCounts <- getRNACount(df)
  return(getEveryNRow(kniCounts, n))
}

df <- getEveryRNANRow(data, 250)
row <- as.numeric(df[1,])
plot(row, type = "l", ylab = "RNA", xlab = "Core", xlim = c(1, length(df[1,])), ylim = c(0, max(df)))
cl <- rainbow(length(df[1,]))

for(i in 2:nrow(df)) {
  row <- as.numeric(df[i,])
  lines(row, type = "l", ylab = "RNA", xlab = "Core",col = cl[i])
}

