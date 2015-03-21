# function to load genome info table from ucsc
getGenomeTable <- function(genome = "mm8"){
  require(plyr)
  gt_url <- sub(pattern = "GENOME",
                replacement = genome,
                x = "http://hgdownload.cse.ucsc.edu/goldenPath/GENOME/database/chromInfo.txt.gz",
                ignore.case = F)
  gt <- readLines(gzcon(url(gt_url, "rb")))
  gt_df <- ldply(gt, function(x){
    tmp <- unlist(strsplit(x , "\t"))
    as.data.frame(cbind(chromosome = tmp[1], 
                        size = as.numeric(tmp[2])), 
                  stringsAsFactors = F)
  })
  gt_df
}


# granges shuffle function, similar to bedtools shuffle
grShuffle <- function(gr_in = NULL, genome_info = NULL, n_sets = 1, method = "chromosome"){
  require(plyr)
  gi <- genome_info
  # remove chromosomes which are not in gr
  gi <- gi[gi$chromosome %in% levels(seqnames(gr_in)),]
  gi[,"size"] <- as.numeric(gi[,"size"])
  # chromosome ends in genomic length cooridinates
  chromosome_ends <- laply(1:nrow(gi), function(idx){
    sum(as.numeric(gi[1:idx,"size"]))})
  # adding chromosome ends in genome coordinates to genome info
  gi <- cbind(gi, chromosome_ends)
  
  gr_out <- GRanges()
  if(method == "chromosome"){
    grl <- split(gr_in, seqnames(gr_in))
    gr_out <- llply(1:nrow(gi), function(idx){
      n <- length(grl[[gi[idx,"chromosome"]]])
      chromosome <- gi[idx,"chromosome"]
      new_starts <- runif(n, 1, as.integer(gi[idx,"size"]))
      in_width <- width(grl[[gi[idx,"chromosome"]]])
      c(gr_out, GRanges(Rle(rep(chromosome, n)),
                        IRanges(start = new_starts, 
                                end = new_starts + in_width)))
    })
    gr_out <- unlist(GRangesList(gr_out))
    gr_out
  }else if(method == "genome"){
    genome_length <- as.numeric(sum(gi[,"size"]))
    # width of input genomic regions
    in_width <- width(gr_in)
    # new start positions in genomic coordinates
    new_starts <- runif(length(gr_in), 1, genome_length)
    # mapping genomic coordinates to chromosomal coordinates
    # and return GRanges object
    for(idx in 1:nrow(gi)){
      chromosome <- gi[idx,"chromosome"]
      if(idx == 1){
        starts <- which(new_starts > 0 & new_starts <= gi[idx,"chromosome_ends"])
      }else{
        starts <- which(new_starts > gi[idx-1,"chromosome_ends"] & new_starts <= gi[idx,"chromosome_ends"])
      }
      start_pos <- gi[idx,3] - new_starts[starts]
      end_pos <- gi[idx,3] - new_starts[starts] + in_width[starts]
      tmp_gr <- GRanges(Rle(rep(gi[idx,"chromosome"], length(starts))),
                        IRanges(start_pos, end_pos))
      gr_out <- c(gr_out, tmp_gr)
    }
    gr_out
  }
}

# non parametric p vslues
nonParP <- function(observed, expected){
  p.out <- sum(abs(expected) > observed)
  p.at <- sum(abs(expected) == observed)
  (p.out + p.at /2) / (length(expected) +1)
}

# outlier
outlier <- function(x, ...){
  idx <- which(x %in% boxplot.stats(x, ...)$out)
  idx
}



