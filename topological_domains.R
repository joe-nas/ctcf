### topological domains
library(rtracklayer)
library(gdata)
# load topological domain data
topological_domains <- "http://www.nature.com/nature/journal/v485/n7398/extref/nature11082-s2.xls"

# define domains of domain barriers as domain start-1: start and domain end: end+1
topod_mm9 = with(read.xls(xls = topological_domains, sheet = 4, 
                       header = F, stringsAsFactors = FALSE),
              GRanges(Rle(rep(V1,2)), IRanges(start = c(V2, V2-1), 
                                       end = c(V2 +1, V2))))


# import the mm9 to mm8 chain file
ch = import.chain("mm9ToMm8.over.chain")
# lift over mm9 to mm8
topod_mm8 <- unlist(liftOver(topod_mm9, ch))


# load head and tail anchors
ctcf_mediated_interactome <- "http://www.nature.com/ng/journal/v43/n7/extref/ng.857-S4.xls"
ctcf_loops <- read.xls(xls = ctcf_mediated_interactome, sheet = 1, 
                       header = T, stringsAsFactors = FALSE, skip = 1)

head_anchors <- with(ctcf_loops, 
                     GRanges(Rle(chrom),
                             IRanges(start = start, end = end),
                             hcc = Histone.cluster.category, 
                             distance = end - start, 
                             loopID = rownames(ctcf_loops),
                             peakID = as.integer(NA),
                             height = as.integer(NA))
                     )


tail_anchors <- with(ctcf_loops, 
                     GRanges(Rle(chrom.1),
                             IRanges(start = start.1, end = end.1), 
                             hcc = Histone.cluster.category, 
                             distance = end.1 - start.1, 
                             loopID = rownames(ctcf_loops),
                             peakID = as.integer(NA),
                             height = as.integer(NA))
                     )
                            

                             
# load ctcf binding sites
ctcf_bs_url <- "http://www.nature.com/ng/journal/v43/n7/extref/ng.857-S2.xls"
ctcf_bs_r <- read.xls(xls = ctcf_bs_url, sheet = 1, 
                       header = T, stringsAsFactors = F, skip = 1)
ctcf_bs_df <- cbind(do.call('rbind', strsplit(as.character(ctcf_bs_r$Binding.loci),'[:-]')), 
                    ctcf_bs_r[,2:3], stringsAsFactors = F)
colnames(ctcf_bs_df) = c("chrom", "start", "end", "height", "interaction")
ctcf_bs_gr <- with(ctcf_bs_df,
                   GRanges(Rle(chrom),
                           IRanges(start = as.numeric(start), 
                                   end = as.numeric(end)),
                           height = height, 
                           interaction = interaction,
                           peakID = rownames(ctcf_bs_df)))
                           

# if possible, annotate head/tail anchors with peak height and peakID
ov_ha <- findOverlaps(head_anchors, ctcf_bs_gr)
head_anchors[queryHits(ov_ha)]$height <- ctcf_bs_gr[subjectHits(ov_ha)]$height
head_anchors[queryHits(ov_ha)]$peakID <- ctcf_bs_gr[subjectHits(ov_ha)]$peakID

ov_ta <- findOverlaps(tail_anchors, ctcf_bs_gr)
tail_anchors[queryHits(ov_ta)]$height <- ctcf_bs_gr[subjectHits(ov_ta)]$height
tail_anchors[queryHits(ov_ta)]$peakID <- ctcf_bs_gr[subjectHits(ov_ta)]$peakID

# merge head and tail anchor sets in order to for distanc calculation later in the process
ht_anchors_df <- merge(as.data.frame(head_anchors), 
                       as.data.frame(tail_anchors), 
                       by = "loopID", all = T, sort = F)

loop_props <- with(ht_anchors_df[complete.cases(ht_anchors_df),], data.frame(
  distance  = end.y - start.x, 
  head_height = height.x, tail_height = height.y, 
  loopID = loopID
))

# remove loops larger than 1MB
loop_props <- loop_props[which(loop_props$distance <= 10^6),]

# create a GRanges object from head and tail anchors which corresponds to model-loops for artificial loops generation
obs_loops_gr <- with(ht_anchors_df[ht_anchors_df$loopID %in% loop_props$loopID,],
                     GRanges(Rle(seqnames.x),
                             IRanges(start = start.x,
                                     end = end.y),
                             head_height = height.x,
                             tail_height = height.y))

# generate function to generate random loops
ctcf_split <- split(ctcf_bs_gr, ctcf_bs_gr$interaction)

randomLoop <- function(gr, n){
# takes GRanges object and an integer as input
# samples 2 times n of gr and combines them to a new gr given a number of constraints
# returns a gr containing hypothetical loops
  n_gr <- length(gr)
  ha_idx <- sample(1: n_gr, replace = T, size = n)
  ta_idx <- sample(1: n_gr, replace = T, size = n)
  
  ov <- findOverlaps(resize(gr[ha_idx], width = 10^6, fix = "start"), gr[ta_idx])
  ha <- gr[ha_idx][queryHits(ov)]
  ta <- gr[ta_idx][subjectHits(ov)]
  loops <- GRanges(Rle(seqnames(ha)), 
                   IRanges(start = start(ha), 
                           end = end(ta)),
                   head_height = ha$height,
                   tail_height = ta$height)
  loops <- loops[which(width(loops) >= 10^4)]
  return(list(gr = loops, 
              sim_loops = data.frame(distance = width(loops),
                                     head_height = loops$head_height,
                                     tail_height = loops$tail_height)))
}

matchLoops <- function(obs_loops, sim_loops){
# takes 2 matrices as input, calculates distances between first and the second matrix
# returns a vector of indeces comprising indeces of the second matrix with min dist to first matrix 
  require(StatMatch)
  require(plyr)
  idx <- c()
  mahal_dist <- mahalanobis.dist(obs_loops, sim_loops)
  for(i in 1:nrow(mahal_dist)){
    tmp <- which(min(mahal_dist[i,]) == mahal_dist[i,])
    if(length(tmp) > 1){
      idx <- c(idx, sample(tmp, 1))
    }else{
      idx <- c(idx, tmp)
    }
  }
  return(idx)
}


# done on server
# library(doMC)
# library(plyr)
# registerDoMC(cores = 6)
# simulated_loops <- llply(1:100, function(x){
#   rl <- randomLoop(ctcf_split$no, 5000)
#   return(rl$gr[matchLoops(loop_props[,1:3], rl$sim_loops)])
# }, .parallel = T)


library(curl)
simulated_loops_url <- 
  "https://github.com/joe-nas/ctcf/raw/master/simulated_loops.Rdata"
load(curl(url = simulated_loops_url, open = "r"))

# count overlaps of simulated loops as well as observed loops with topological domain barriers
expected <- countOverlaps(GRangesList(simulated_loops), topod_mm8)
observed <- countOverlaps(GRangesList(obs_loops_gr), topod_mm8)
## visualizing observed and expected loop-topological domain overlap
library(lattice)
histogram(expected, xlim = range(c(expected+20,observed-20)), 
          breaks = 30, scales = list(tck = c(1,0), cex = 1.5), 
          auto.key=list(text=c("Observed", "Expected"), cex = 2,
                        columns = 2, lines = F, rectangles = T),
          par.settings = simpleTheme(col = c("red", "#08306B"), 
                                     border = "transparent"),
          panel = function(x,y, ...){
            panel.histogram(x, ..., col = "#08306B")
            panel.abline(v = observed, col = "red" , lwd = 4)
          }, 
          xlab = list(label = "Topological domain overlaps", cex = 2),
          ylab = list(cex = 2))

