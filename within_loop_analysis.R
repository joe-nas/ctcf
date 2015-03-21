## anova analysis
library(gdata)
library(plyr)
library(lattice)

# load chip-seq data
load("chipseq.RDat")

# load interaction data set
ctcf_mediated_interactome <- "http://www.nature.com/ng/journal/v43/n7/extref/ng.857-S4.xls"
ctcf_loops <- read.xls(xls = ctcf_mediated_interactome, sheet = 1, 
                       header = T, stringsAsFactors = FALSE, skip = 1)

loops <- with(ctcf_loops, 
              GRanges(Rle(chrom), 
              IRanges(start = start,
                      end = end.1),
              hcc = Histone.cluster.category))
# split loops by histone cluster category
loops_split <- split(loops, loops$hcc)
# remove loops without category
loops_split[[1]] <- NULL


library(plyr)
library(doMC)
registerDoMC(cores = 4)

# overlap loops with chip-seq data
loops_chip_intersection <- llply(GRangesList(ChIP_seq), 
                                 function(x){
                                   llply(loops_split, 
                                         .parallel  = T,
                                         function(y){
                                           countOverlaps(y,x)})
                                   })
# merge per loop intersections count into a data.frame
lci_tab <- Reduce(cbind,Map(loops_chip_intersection, 
                            f = function(x){
                              as.data.frame(as.numeric(unlist(x))))
                  })
colnames(lci_tab) <- names(loops_chip_intersection)

# add factor variable histone cluster class: y
df <- ldply(loops_split, length)
y <- Reduce(c,Map(f = function(x) rep(df[x,1],df[x,2]) , 1:nrow(df)))
lci_tab$y <- as.factor(y)

# perform anova and compute "Compute Tukey Honest Significant Differences"
within_loop_anova <- alply(lci_tab[,1:ncol(lci_tab)-1], 2,
                           function(x){
                             df <- as.data.frame(TukeyHSD(aov(x[[1]]~lci_tab[,"y"]), 
                                                 conf.level = 0.999)[[1]])
                             df <- df[order(df$`p adj`,decreasing = T),]
                             df <- cbind(df, y = factor(rownames(df), 
                                                        levels = rownames(df)))
                             df <- cbind(df, series = factor(colnames(x)))
                           })
names(within_loop_anova) <- colnames(lci_tab)[1:(ncol(lci_tab)-1)]


# plot function for resultsof tukeyHSK
tukey.plot <- function(data, title = ""){
  xrang <- range(data[,1:4])
  xrang <- xrang + sqrt(abs(xrang)) * sign(xrang)
  stripplot(y ~ diff, data = data, 
            scales = list(tck = c(1,0)),
            xlim = xrang, pch = 16,
            main = list(label = paste(title, "conf.level: .999", sep = "; "),
                        cex = .8),
            panel = function(x,y,...){
              panel.abline(v = 0, col = "grey", lwd = 2)
              panel.segments(data$lwr,data$y,data$upr,data$y, lwd = 1)
              panel.stripplot(data$lwr,data$y,cex = .5,...)
              panel.stripplot(data$diff,data$y,cex = .5,...)
              panel.stripplot(data$upr,data$y,cex = .5,...)
            },
            xlab = list(label = "differences in mean levels of class I-V loops", 
                        cex = .65))
}
# generate the actual plots
tukey_plots <- llply(names(within_loop_anova), function(x){
  tukey.plot(within_loop_anova[[paste(x)]],x)
})

# merge plots
library(gridExtra)
do.call(grid.arrange, c(lapply(tukey_plots[1:12], update), list(ncol = 3, nrow = 4)))
do.call(grid.arrange, c(lapply(tukey_plots[13:20], update), list(ncol = 3, nrow = 4)))
