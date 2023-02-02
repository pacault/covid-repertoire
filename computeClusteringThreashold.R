#!/usr/bin/env Rscript
library(optparse)
library(alakazam)
library(shazam)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="change-o database file", metavar="character"),
  make_option(c("-m", "--model"), type="character", default="aa", 
              help="model to use for clusterize (aa:distance from aminoacid sequence, ham:distance for nucleotidic sequence)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=".", 
              help="Output directory for result files", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  # stop("No change-o database file provided", call.=FALSE)
}

database_file = opt$file
model = opt$model
outdir = opt$outdir

# database_file = "4068KL_db-pass_prod.tsv"
model = "aa"

# print(database_file)
# print(outdir)
# print(model)

if(!file.exists(database_file)){
  stop("The input file doesn't exists")
}

packageVersion('shazam')
packageVersion('alakazam')

db <- readChangeoDb(database_file)

dist_ham <- distToNearest(db, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="aa", normalize="len", nproc=1)
print(dist_ham)
# 
# p1 <- ggplot(subset(dist_ham, !is.na(dist_nearest)),
#              aes(x=dist_nearest)) +
#   theme_bw() +
#   xlab("Hamming distance") +
#   ylab("Count") +
#   ggtitle("Density Method p1 aa") +
#   scale_x_continuous(breaks=seq(0, 1, 0.1)) +
#   geom_histogram(color="white", binwidth=0.02) +
#   geom_vline(xintercept=0.12, color="firebrick", linetype=2)
# plot(p1)

# Find threshold using density method
threshold <- NA
tryCatch(
    expr = {
        output <- findThreshold(dist_ham$dist_nearest, method="density")
        threshold <- output@threshold
        # Plot distance histogram, density estimate and optimum threshold
        pdf(paste(outdir, "/rplot.pdf", sep=""))
        plot(output, title=paste("Density (Method", model, ")", sep=" "))
        dev.off()
    },
    error = function(e){ 
        threshold <- NA
    }
)
#threshold <- output@threshold

if (is.na(threshold)) {
  threshold <- 1 / 30
  cat("\n")
  cat(paste("WARNING : No bimodal distribution, impossible to compute threshold. Default threshold defined =", threshold, sep= " "))
  cat("\n\n")
} 

# fileConn<-file("threshold.txt")
# writeLines(cat(threshold), fileConn)
# close(fileConn)
cat(as.character(threshold),file=paste(outdir, "/threshold.txt", sep=""),sep="\n")


