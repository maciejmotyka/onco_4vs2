require("jsonlite")
require("proxy")
require("apcluster")
require("uuid")
require("dplyr")

sumerWeightedSetCover <- function(config_file, full_gmts_dir, output_dir, n_threads = 4, project) {

  ## Check if config file makes sense
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

  # Retrieve some config data
  config <- full_config$data
  n_platform <- nrow(config)

  # Create empty variables
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

  ### ----Affinity Propagation----
  ## 1. Create input data for Affinity Propagation

  ap_data <- list()
  ap_data$output.dir <- output_dir # output_dir is an argument of the sumer function
  ap_data$md5val <- "sumer" # used when naming some of the output files

  # Iterate over platforms, read the results from Weighted Set Cover, substitute AGS genes with genes from the full GMT files
  genesetIds <- list()
  genesetInfo <- list()
  for(i in seq_len(n_platform)){ # iterates over platforms
    cur_config <- config[i,]
    cur_platform <- cur_config[["platform_name"]]
    cur_platform_abbr <- cur_config[["platform_abbr"]]
    cur_sc_output_dir <- work_dirs[[cur_platform]] # reads which folder belongs to which platform

    # reads the sctopsets.txt (pathway:score:genes) file from a platform specific folder into a dataframe
    sc_gs <- read.table(file.path(output_dir, cur_sc_output_dir, "sctopsets.txt"),
                        sep="\t", stringsAsFactors=FALSE, header=TRUE)

    # Determine which full GMT file to use
    if(file.exists(file.path(full_gmts_dir, paste0(cur_platform, "-", project, ".gmt")))){
      full_gmt_file <- file.path(full_gmts_dir, paste0(cur_platform, "-", project, ".gmt"))
    } else {
      full_gmt_file <- file.path(full_gmts_dir, paste0(project, "-", cur_platform, ".gmt"))
    }

    # Parse GMT file with all of the genes that belong te each pathway
    full_gmt <- read.table(full_gmt_file, sep = "\t", stringsAsFactors = F,
                           col.names = c("pathway_name", "pathway_id", "all_genes"), quote = "")

    # Substitute AGS genes interacting with a pathway with all of the genes belonging to that pathway
    sc_gs2 <- dplyr::left_join(sc_gs, full_gmt, by = c("geneset" = "pathway_name")) %>%
      dplyr::select(geneset, score, all_genes)

    # Populate genesetIds and genesetInfo with relevant pathway data from sc_gs2
    for (row in 1:nrow(sc_gs2)) {
      tmpname <- paste(cur_platform_abbr, sc_gs2[row, "geneset"], sep="_")
      genesetIds[[tmpname]] <- unlist(strsplit(sc_gs2[row, "all_genes"], ","))
      genesetInfo[[tmpname]] <- list()
      genesetInfo[[tmpname]][["signP"]] <- sc_gs2[row, "score" ] # store  score as signP
      genesetInfo[[tmpname]][["leNum"]] <- length(genesetIds[[tmpname]])
      # not used
      genesetInfo[[tmpname]][["PValue"]] <-  NULL
      genesetInfo[[tmpname]][["NES"]] <-  sc_gs2[row, "score" ] # store score as NES
      # not used
      genesetInfo[[tmpname]][["FDR"]] <- NULL
      genesetInfo[[tmpname]][["platform"]] <- cur_platform
    }
  }

  ap_data$genesetIds <- genesetIds # add genesetIds to ap_data
  ap_data$genesetInfo <- genesetInfo # add genesetInfo to ap_data

  ## 2. Perform the Affinity Propagation
  ap_results <- affinityPropagation(ap_data) # performs the AP
}
