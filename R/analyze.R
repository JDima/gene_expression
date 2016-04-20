setwd("../src")
cols <- read.table("means_head.txt", sep = "\t", head = FALSE)
data <- read.table("model_output/stats/means.txt", sep = "\t", head = FALSE)
data <- data[, 1:length(cols)]
colnames(data) <- unlist(cols)
