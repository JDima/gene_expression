readData <- function(input, header)
{
  cols <- read.table(header, sep = "\t", head = FALSE)
  data <- read.table(input, sep = "\t", head = FALSE)
  data <- data[, 1:length(cols)]
  colnames(data) <- unlist(cols)
  return(data)
}

getNRow <- function(df, nrow)
{
  return(df[nrow, ])
}

getTFNRow <- function(df, nrow, tf)
{
  tmp <- df[nrow, grep(tf, colnames(df), perl = TRUE )]
  return(tmp[nrow, grep("free", colnames(tmp), perl = TRUE )])
}

getTfCount <- function(df, tf)
{
  tmp <- df[, grep(tf, colnames(df), perl = TRUE )]
  return(tmp[ , grep("free", colnames(tmp), perl = TRUE )])
}

getEveryTfNRow<- function(df, n, tf)
{
  kniCounts <- getTfCount(df, tf)
  return(getEveryNRow(kniCounts, n))
}

getEveryNRow <- function(df, n)
{
  maxT <- length(df[,1])
  nrows <- seq(1, maxT, by = n)
  return(df[nrows,])
}

plotSliceEveryTfNRow <- function(data, N, tf)
{
  df <- getEveryTfNRow(data, 250, tf)
  row <- as.numeric(df[1,])
  plot(row, type = "l", ylab = "Сoncentration", xlab = "Core", xlim = c(1, length(df[1,])), ylim = c(0, max(df)), main = tf)
  cl <- rainbow(length(df[1,]))
  for(i in 2:nrow(df)) {
    row <- as.numeric(df[i,])
    lines(row, type = "l", ylab = "Сoncentration", xlab = "Core",col = cl[i])
  }
}

# START KNI

getKniNRow <- function(df, nrow)
{
  return(getTFNRow(df, nrow, "kni"))
}

getKNICount <- function(df)
{
  return(getTfCount(df, "kni"))
}

getEveryKniNRow<- function(df, n)
{
  kniCounts <- getKNICount(df)
  return(getEveryNRow(kniCounts, n))
}

# END KNI

# START mRNA

getRNANRow <- function(df, nrow)
{
  return(df[nrow, grep("mRNA", colnames(df), perl = TRUE )])
}

getRNACount <- function(df)
{
  return(df[ , grep("mRNA", colnames(df), perl = TRUE )])
}

getEveryRNANRow<- function(df, n)
{
  rnaCounts <- getRNACount(df)
  return(getEveryNRow(rnaCounts, n))
}

plotSliceEveryRNANRow <- function(data, N)
{
  df <- getEveryRNANRow(data, 250)
  row <- as.numeric(df[1,])
  plot(row, type = "l", ylab = "Сoncentration", xlab = "Core", xlim = c(1, length(df[1,])), ylim = c(0, max(df)), main = "mRNA")
  cl <- rainbow(length(df[1,]))
  for(i in 2:nrow(df)) {
    row <- as.numeric(df[i,])
    lines(row, type = "l", ylab = "Сoncentration", xlab = "Core",col = cl[i])
  }
}
# END mRNA