### ----Function definitions----
check_config <- function(config_file) {
  # maximum of number of platforms we can support
  MAX_N_PLATFORM <- 7
  old_config <- jsonlite::fromJSON(config_file)
  if (! "project" %in% names(old_config)) {
    stop("please provide project description")
  }
  if (! "data" %in% names(old_config)) {
    stop("please provide project data information")
  }
  if (! "top_num" %in% names(old_config)) {
    old_config$top_num <- 50
  }
  config <- old_config$data
  n_platform <- nrow(config)
  if (n_platform > MAX_N_PLATFORM) {
    stop(paste0("currently supports up to ", MAX_N_PLATFORM, " platforms"))
  } else if (n_platform < 1) {
    stop("must provide at least data from 1 platform")
  }
  # check if data file valid
  for (i in seq_len(n_platform)) {
    cur_platform <- config[i, "platform"]
    if(!file.exists(config[i, "gmt_file"])) {
      stop(paste0("gmt_file does not exist for platform ", cur_platform))
    } else if (!file.exists(config[i, "score_file"])){
      stop(paste0("score_file does not exist for platform ", cur_platform))
    }
  }
  # TODO: check platform abbr is unique
  return(old_config)
}

# convert list of gene names to gmt file
# each list component has a gene set name
# the value of each list item is an array of gene names/ids
list2gmt <- function(rlist, out_file) {
  m <- do.call(rbind, lapply(seq_along(rlist), function(i) paste(c(names(rlist)[[i]],"", rlist[[i]]), collapse="\t")))
  write.table(m,file=out_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# for now, for my own use
extractScore <- function(rlist, out_file) {
  m <- do.call(rbind, lapply(seq_along(rlist), function(i) paste(c(names(rlist)[[i]], rlist[[i]][["signP"]]), collapse="\t")))
  write.table(m,file=out_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
}

gmt2list <- function(gmt_file){
  if (!file.exists(gmt_file)) {
    stop("There is no such gmt file.")
  }
  x <- scan(gmt_file, what="character", sep="\n", quiet=TRUE, strip.white = T) # added strip.white=T to remove the trailing tabs in my gmt files
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  gmt_list <- lapply(y, `[`, c(-1,-2))
}

# for cytoscape, set node size
computeNodeDim <- function(maxVal, minVal, val){
  maxCyNodeSize<-200;
  minCyNodeSize<-50;
  return((val-minVal)/(maxVal-minVal)*(maxCyNodeSize-minCyNodeSize)+minCyNodeSize)
}

# set color: min: blue, middle: white, max: red
computeNodeColor <- function(maxVal, minVal, val){
  # if all positive, map to [0.5, 1]
  if(minVal > 0){
    maxColor <- 1;
    minColor <- 0.5;
  } else if(maxVal < 0) {  # all negative
    maxColor <- 0.5;
    minColor <- 0;
  } else {
    larger = max(maxVal, abs(minVal));
    maxVal <- larger;
    minVal <- -1.0*larger;
    maxColor <- 1;
    minColor <- 0;
  }
  return((val-minVal)/(maxVal-minVal)*(maxColor-minColor)+minColor)
}
