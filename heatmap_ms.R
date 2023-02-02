###############################################################################
#
# script: heatmap_ms.R
#
# Manuscript: iScience manuscript ISCIENCE-D-22-03671R1
#
# Objective: VJ heatmap for ARDS/HV cohorts, ADN/ARN 
#
# input: clonotypes_totaux_covabdab_fix_<ADN|ARN>.tsv
#
# method
#   - choose ADN or ARN (beware: nonARN are not the same pts as nonADN!)
#   - heatmap with library(ComplexHeatmap) 
#       -> row and col annotations + gap between clusters
#   - DE analysis using library(edgeR)
#   - heatmap
#       pre-filtering using Wald DE genes betw ARDS J0+J7 vs Control 
#         and heatmap on DE VJ for all data
#
# output 
#   Repertoire_<ADN|ARN>_<cond1_vs_cond2>_<allData>_heatmap.pdf
#   Repertoire_<ADN|ARN>_<cond1_vs_cond2>_<allData>_Top.xlsx
#
# 2023 D Rossille
#
###############################################################################

# emptying envt
rm(list=ls())

# libraries
library(openxlsx)
library(edgeR)
library(ComplexHeatmap)
library(ggplot2)


# global vars
wd <- getwd()
wd
dirData <- paste(wd, "DATA", sep="/")
dirRes <- paste(wd, "RESULTS", sep="/")




###############################################################################
#
# Parameters to be selected by user
# 
###############################################################################

# to be selected by User: heatmap on ARN or ADN data?
type <- "ARN"
type <- "ADN"




###############################################################################
#
# Script's Parameters
# 
###############################################################################


# dataset to be imported according to $type
d.name <- paste("clonotypes_totaux_covabdab_fix_", type, ".tsv", sep="")


# define non ARDS patients for ADN data
ADN_nonARDS <- c("DC_02", "KM_74", "PE_80", "SJ_50", "TOSUN", "YM_65") 


# define colors for groups
group.cols <- c("red",  "#FFDB00",  "#49FF00",  "#00FF92",  "#0092FF",  "#4900FF",  "#FF00DB" )
names(group.cols) <- c("ARDSneg_D0", "ARDSneg_D7", "ARDSneg_M4", 
                       "ARDSpos_D0", "ARDSpos_D7", "ARDSpos_M4", "Control")

group3.cols <- c("#FFDB00",  "purple",  "red" )
names(group3.cols) <- c("without ARDS", "with ARDS", "Control")

day.cols <- c("red", "green", "chartreuse4", "blue", "deepskyblue2")
names(day.cols) <- c("Control", "D0", "D7", "M4", "D14")

inCovabdab.cols <- c("red", "blue")
names(inCovabdab.cols) <- c("0", "1") 




#######################################################################################
#
# function
#
#######################################################################################

funct_heatmap <- function(DonneesFiltrees, name0=NULL, col_ha, row_ha, fontsize, clusters.nrow, clusters.ncol, col_fun,
                               repertoire, width, height, out=FALSE)
{
  # parameters
  ## 
  ## repertoire : name of the output file (pdf)
  # library(ComplexHeatmap)
  
  if (out==TRUE)  jpeg(repertoire, width=width, height = height, units="in", res=300)
  p <- Heatmap(DonneesFiltrees, 
          name = "scale", #title of legend
          row_names_gp = gpar(fontsize = fontsize), # Text size for row names
          show_column_names = FALSE,
          clustering_distance_rows = "euclidean",
          clustering_method_rows = "ward.D2",
          clustering_distance_columns = "pearson",
          clustering_method_columns = "ward.D2",
          row_dend_reorder = FALSE,
          column_dend_reorder = FALSE,
          col = col_fun, 
          top_annotation = col_ha,
          right_annotation = row_ha,
          row_split=clusters.nrow, column_split = clusters.ncol, # split rows and cols
          row_title = NULL, column_title = NULL  # no labels for the split rows and cols
  )
  print(p)

  if (out==TRUE) 
  {
    dev.off()
    cat(paste("\nSaving into", repertoire))
  }
    
}







###############################################################################
#
# Preparation of the input file
#
###############################################################################

# large table 
T <- read.delim(paste(dirData, d.name, sep="/"))
# table(T$type) # ARN or ADN

colnames(T)

T$TotalReads <- T$reads/T$freq
T$VJ <- paste(T$V, T$J, sep = "_")
T$PatientID <- paste(T$X.samplename, T$barcode, T$group, T$day, sep = "_")

#### construction of the VJ table for edgeR
Data <- T[T$type == type,] # choisir
T_VJ_clono <- Data[,c(1:7,12,27,28,29)]
colnames(T_VJ_clono)
VJ_Patient <- paste(T_VJ_clono$PatientID, T_VJ_clono$VJ, sep = ".") # associe le VJ e PaientID
T_VJ_clono$VJ_Patient <- VJ_Patient 

# identify all VJ combinations for each patient
T_VJ_tot <- T_VJ_clono[duplicated(VJ_Patient, incomparables = FALSE)==FALSE,]
T_VJ_tot <- T_VJ_tot[order(T_VJ_tot$VJ_Patient),]

# count all reads per VJ and per patient
VJ_count <- aggregate(T_VJ_clono[,8], by = list(VJ_Patient) , sum)
VJ_count <- VJ_count[order(VJ_count[,1]),]

# count VJ
all(VJ_count[,1] == T_VJ_tot$VJ_Patient)                          
T_VJ_tot$VJ_count <- VJ_count[,2]

# Patient metadata table
Patient <-T_VJ_tot[duplicated(T_VJ_tot$PatientID, incomparables = FALSE)==FALSE,c(1:5,11)] 

if (type=="ARN")
{
  Patient$ARDS <- grepl("R", Patient$barcode) 
  Patient$Groupe <- paste(Patient$day, Patient$ARDS, sep = "_") 
  Patient$Groupe[Patient$Groupe == "na_FALSE"] <- "Control" 
  Patient$Groupe[grep("J0_TRUE", Patient$Groupe)] <- "ARDSpos_D0"
  Patient$Groupe[grep("J7_TRUE", Patient$Groupe)] <- "ARDSpos_D7"
  Patient$Groupe[grep("J14_TRUE", Patient$Groupe)] <- "ARDSpos_D14"
  Patient$Groupe[grep("M4_TRUE", Patient$Groupe)] <- "ARDSpos_M4"
  Patient$Groupe[grep("J0_FALSE", Patient$Groupe)] <- "ARDSneg_D0"
  Patient$Groupe[grep("J7_FALSE", Patient$Groupe)] <- "ARDSneg_D7"
  Patient$Groupe[grep("J14_FALSE", Patient$Groupe)] <- "ARDSneg_D14"
  Patient$Groupe[grep("M4_FALSE", Patient$Groupe)] <- "ARDSneg_M4"
  Patient$Cov <- Patient$day  
  Patient$Cov[Patient$Cov == "na"] <- "Control"
  Patient$Cov[grep("J0", Patient$Cov)] <- "D0"
  Patient$Cov[grep("J7", Patient$Cov)] <- "D7"
  Patient$Cov[grep("M4", Patient$Cov)] <- "M4"
  Patient$Cov <- factor(Patient$Cov, levels=c("Control", "D0", "D7", "M4"))
  Patient$group3 <- sapply(Patient$Groupe, function(x) unlist(strsplit(x, split="_"))[1])  
  Patient$group3 <- factor(Patient$group3)
  levels(Patient$group3)[grep("neg", levels(Patient$group3))] <- "without ARDS"
  levels(Patient$group3)[grep("pos", levels(Patient$group3))] <- "with ARDS"
  Patient$group3 <- factor(Patient$group3, levels=c("Control", "without ARDS", "with ARDS"))
  levels(Patient$group3)
  Patient$id <- NA
  for (i in 1:nrow(Patient))
  {
    if (Patient$group[i]=="Control") Patient$id[i] <- Patient$barcode[i]
    if (Patient$group[i]=="Patient") Patient$id[i] <- paste(Patient$barcode[i], Patient$day[i], sep="_")
  }
  Patient$id <- gsub("_J", "_D", Patient$id)
  
}

if (type=="ADN")
{
  Patient$group3 <- "with ARDS"
  Patient$group3[which(Patient$group == "Control")] <- "Control"
  Patient$group3[which(Patient$X.samplename %in% ADN_nonARDS)] <- "without ARDS"
  table(Patient$group3)
  Patient$Cov <- Patient$day 
  Patient$Cov[Patient$Cov == "na"] <- "Control"
  Patient$Cov[grep("J0", Patient$Cov)] <- "D0"
  Patient$Cov[grep("J7", Patient$Cov)] <- "D7"
  Patient$Cov[grep("J14", Patient$Cov)] <- "D14"
  Patient$Cov <- factor(Patient$Cov, levels=c("Control", "D0", "D7", "D14", "M4"))
  Patient$Groupe <- "Control"
  Patient$Groupe[grep("ARDS", Patient$group3)] <- paste(Patient$group3[grep("ARDS", Patient$group3)],
                                                        Patient$Cov[grep("ARDS", Patient$group3)], sep="_")
  Patient$id <- NA
  for (i in 1:nrow(Patient))
  {
    if (Patient$group[i]=="Control") Patient$id[i] <- Patient$barcode[i]
    if (Patient$group[i]=="Patient") Patient$id[i] <- paste(Patient$barcode[i], Patient$day[i], sep="_")
  }
  Patient$id <- gsub("_J", "_D", Patient$id)
}

# nb subjects per group
length(unique(Patient$X.samplename[grep("Control", Patient$group3)]))
length(unique(Patient$X.samplename[grep("without", Patient$group3)]))
length(unique(Patient$X.samplename[grep("^with ", Patient$group3)]))

# list of VJ
VJ_liste <- T_VJ_tot$VJ[duplicated(T_VJ_tot$VJ)==FALSE]

# table of VJ and patient, used for heatmap display
Donnees <- matrix(nrow = length(VJ_liste), ncol = nrow(Patient), data = 0)
rownames(Donnees) <- VJ_liste
colnames(Donnees) <- Patient$id

for (i in 1:nrow(Donnees)) 
{
  for (j in 1:ncol(Donnees))
  {
    #  i<- 2; j<-3
    A <- as.numeric(which(T_VJ_tot$VJ == rownames(Donnees)[i]))
    B <- as.numeric(which(T_VJ_tot$PatientID == Patient$PatientID[which(Patient$id == colnames(Donnees)[j])] ))
    C <- A[A %in% B]
    if (length(C)>0) {Donnees[i,j] <- T_VJ_tot$VJ_count[C]}
  }
}



###############################################################################
#
# Pre-filtering using Wald DE genes betw ARDS J0+J7 vs Control 
#         and heatmap on DE VJ for all data
#
###############################################################################
# heatmap only on DE VJ (DESeq2 Wald test)

# name of heatmap output file
elt <- "w-wo_ARDS_D0_D7_vs_Control"
repertoire <- paste("Repertoire", type, elt, "allData", "heatmap.jpeg", sep="_")
name0 <- paste(type, elt, "allData", sep="_")

# --- identify DE VJ for comparison ARDS J0+J7 vs Control 
Patient$group2 <- NA
Patient$group2[grep("Control", Patient$Groupe)] <- "Control"
Patient$group2[grep("D0|D7", Patient$Groupe)] <- "ARDS_D0_D7"
group <- as.factor(Patient$group2) 
# edgeR on subtable of patients only ARDS J0 or J7 or Control
temp <- Donnees[, which(colnames(Donnees) %in% Patient$id[grep("ARDS|Control", Patient$group2)])]
# create DGEList object
y <- DGEList(counts = temp, group = na.omit(group), genes = rownames(temp)) 
# filter the data
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
# normalisation 
y <- calcNormFactors(y)
# estimation de la dispersion 
y <- estimateDisp(y, robust = TRUE)


levels(group)
et <- exactTest(y, pair = c("Control","ARDS_D0_D7"))

DE_et <- topTags(et, n=30000, p.value = 0.1, adjust.method = "BH") 
Top <- (DE_et$table)
# write.xlsx(Top, gsub("heatmap.pdf", "DE_VJ.xlsx", repertoire), overwrite = TRUE)

# --- filter all data on DE VJ
group <- as.factor(Patient$Groupe) 
y <- DGEList(counts = Donnees, group = group, genes = rownames(Donnees)) 
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, robust = TRUE)

logcpm <- cpm(y, log = TRUE, prior.count = 3)
DonneesFiltrees <- logcpm[rownames(logcpm) %in% Top[,1],]
DonneesFiltrees <- t(scale(t(DonneesFiltrees))) # Centers and scales data.



# --- heatmap
# --- colors for heatmap annotations
cols <- list(Group = group3.cols, Day = day.cols)
col_ha <- HeatmapAnnotation(Group=as.factor(Patient$group3), 
                            Day=as.factor(Patient$Cov), col=cols, 
                            annotation_name_gp= gpar(fontsize = 8, font=2))

#inCovabad annotations
covab <- read.xlsx("distrib_covab.xlsx")
names(covab)
inter <- intersect(covab$Var1, rownames(DonneesFiltrees))
cov <- data.frame(VJ=rownames(DonneesFiltrees))
cov$inCovabdab <- sapply(cov$VJ, function(x, y) if (length(grep(x, y))>0) 1 else 0, y=inter)
cov$Freq <- NA
for (i in 1:nrow(cov))
  if (cov$inCovabdab[i]==1) cov$Freq[i] <- covab$Freq[which(covab$Var1 == cov$VJ[i])]

# covabdab as a Freq with scale colors
q1 <- round(as.integer(summary(cov$Freq)["1st Qu."]))
q3 <- round(as.integer(summary(cov$Freq)["3rd Qu."]))
col_fun = colorRamp2(c(0, 1, q1, q3, 150), c("white", rainbow(20)[c(15,17,19)], "red"))
row_ha <- rowAnnotation(covabdab=cov$Freq, na_col = "white", col=list(covabdab=col_fun),
                        annotation_name_gp= gpar(fontsize = 10, font=2))


# create heatmap
funct_heatmap(DonneesFiltrees, name0=name0, col_ha, row_ha, fontsize=8, 
              clusters.nrow=4, clusters.ncol=5, 
              col_fun=NULL, repertoire=repertoire, width=10, height=8, out=TRUE)





