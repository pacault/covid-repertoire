###############################################################################
#
# script: covab_ms.R
#
# Manuscript: iScience manuscript ISCIENCE-D-22-03671R1
#
# objective : extract covabad human sequences from input file and get occurences
#
# input: covabdab_141121 CoV2 human patients and vaccinnees.xlsx
#
# output: table for each VJ pair its frequency (occurence) in input file
#
# 2023 D. Rossille
#
###############################################################################

# emptying envt
rm(list=ls())

# libraries
library(openxlsx)


# global vars
wd <- getwd()
wd
dirData <- paste(wd, "DATA", sep="/")
dirRes <- paste(wd, "RESULTS", sep="/")


###############################################################################
# covabdab db

covab <- read.xlsx (paste(dirData, "covabdab_141121 CoV2 human patients and vaccinnees.xlsx", sep="/"),
                    sheet="Sheet1", startRow = 2)
names(covab)

# restrict to useful cols
cols <- c("Name", "Heavy.V.Gene", "Heavy.J.Gene", "CDRH3")
covab <- covab[, which(names(covab) %in% cols)]


# restrict to human sequences
covab <- covab[grep("Human", covab$Heavy.V.Gene), ]
covab <- covab[grep("Human", covab$Heavy.J.Gene), ]
covab$Heavy.V.Gene <- sapply(covab$Heavy.V.Gene, function(x) unlist(strsplit(x, split=" "))[1])
covab$Heavy.J.Gene <- sapply(covab$Heavy.J.Gene, function(x) unlist(strsplit(x, split=" "))[1])


# concatenate to form VJ
for (i in 1:nrow(covab))
  covab$VJ[i] <- paste(covab$Heavy.V.Gene[i], covab$Heavy.J.Gene[i], sep="_")
covab <- covab[order(covab$VJ), ]
length(unique(covab$VJ))

# distribution VJ pairs
tab <- table(covab$VJ)
tab <- as.data.frame(tab)
tab <- tab[order(tab$Freq, decreasing=TRUE), ]
n <- tab$Var1
tab$Var1 <- factor(tab$Var1, levels=n)


# saving
write.xlsx(tab, "distrib_covab.xlsx")


