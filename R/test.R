source("common.R")

setwd("/Users/JDima/PycharmProjects/gene_expression/test/result")

files <- list.files()
out.files <- files[grep("output", files, fixed=T)]
curdir <- getwd()
for (f in out.files)
{
  setwd(f)
  
  data <- readData("stats/means.txt", "../../means_head.txt")
  print(getwd())
  
  plotSliceEveryTfNRow(data, 250, "kni")
  plotSliceEveryTfNRow(data, 250, "mRNA")
  
  setwd(curdir)
}


data <- readData("means_head.txt", "model_output/stats/means.txt")

plotSliceEveryTfNRow(data, 250, "kni")
plotSliceEveryTfNRow(data, 250, "mRNA")
