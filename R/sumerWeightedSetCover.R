sumerWeightedSetCover <- function(config_file, output_dir, n_threads = 4) {
  full_config <- check_config(config_file)

  ## Make sure that the output directory exists and is empty
  if(dir.exists(output_dir)) {
    overwrite <- readline(prompt = "Output directory not empty, do you want to overwrite existing output? (enter Y or N)")
    if(toupper(overwrite) != "Y"){
      return(0)
    } else {
      unlink(output_dir, recursive = T)
      dir.create(output_dir)
    }
  } else {
    dir.create(output_dir)
    }

  config <- full_config$data
  n_platform <- nrow(config)

  all_platform_abbr <- c()
  all_data <- data.frame(platfrom=character(), name=character(), size=integer(), score=numeric())
  work_dirs <- list()

  # iterates over platforms; prepares data and saves as genesets.RData; runs topGeneSets and topSCGenesets
  # topGeneSets: saves topsets.txt with pathway:score:genes info for up to max_num pathways; plots the score and saves as topsets.{pdf,png}
  # topSCGeneSets:
  for(i in seq_len(n_platform)){
    work_dir <- uuid.gen()
    cur_output_dir <- file.path(output_dir, work_dir)
    if(!dir.exists(cur_output_dir)) {
      dir.create(cur_output_dir)
    }
    cur_config <- config[i,] # set config for current platform
    cur_config <- cbind(cur_config, top_num=full_config$top_num) # add max number of sets to current conf
    cur_platform <- cur_config[["platform_name"]] # set current platform name
    cur_platform_abbr <- cur_config[["platform_abbr"]] # set current platform abbreviation
    all_platform_abbr <- c(all_platform_abbr, cur_platform_abbr) # read all abbreviations
    work_dirs[[cur_platform]] <- work_dir # create dictionary platform:folder_name
    print(paste0("platform: ", cur_platform)) # print what platform is being processed
    geneset_ids <- gmt2list(cur_config[["gmt_file"]]) # create list pathway_id:genes for current platform
    score_list <- scan(cur_config[["score_file"]], what="", sep="\n", quiet=TRUE) # parse score file
    score_list <- strsplit(score_list, "\t") # split on tab
    names(score_list) <- sapply(score_list, `[[`, 1) # set names
    geneset_info <- lapply(score_list, `[`, c(-1)) # drop the names elements
    geneset_info <- lapply(geneset_info, "as.double") # now we have pathway:score dictionary
    save(geneset_info, geneset_ids, file=file.path(cur_output_dir, "genesets.RData")) # save pathway:score and pathway:genes dicts for each platform in separate folder
    top_ret <- topGeneSets(cur_config, geneset_ids, geneset_info, cur_output_dir) # select the best scoring pathways up to max allowed??; top_ret seems not to be used anywhere
    sc_ret <- topSCGeneSets(cur_config, geneset_ids, geneset_info, cur_output_dir, n_threads) # perform set cover??; sc_ret seems not to be used anywhere
    # order by name
    geneset_ids <- geneset_ids[order(names(geneset_ids))]
    geneset_info <- geneset_info[order(names(geneset_info))]
    set_platfrom <- rep(cur_platform_abbr, length(geneset_ids))
    set_name <- names(geneset_ids)
    set_size <-  unname(sapply(geneset_ids, length)) # get the number of pathways returned by set cover
    set_score <- unname(sapply(geneset_info, "as.double")) # !! need to repair this; should be equal to the number of genes in a set, but now every list holding the genes contains also empty elements, which makes them appear equal length; is created by gmt2list function
    all_data <- rbind(all_data, data.frame(platform=set_platfrom, name=set_name, size=set_size, score=set_score, stringsAsFactors=FALSE))
  }
}

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

gmt2list <- function(gmt_file){
  if (!file.exists(gmt_file)) {
    stop("There is no such gmt file.")
  }
  x <- scan(gmt_file, what="character", sep="\n", quiet=TRUE, strip.white = T) # added strip.white=T to remove the trailing tabs in my gmt files
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  gmt_list <- lapply(y, `[`, c(-1,-2))
}
