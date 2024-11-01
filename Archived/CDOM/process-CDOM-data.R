## process DBS and DIT EEM's projects to get indices

DBS_file <- "~/1_Research/0_Misc/BSLE_Degradation_EEMs/RC3_DBS_CDOM"
DIT_file <- "~/1_Research/0_Misc/BSLE_Degradation_EEMs/RC3_DIT_CDOM"

library(fewsdom)
#process DIT files
run_eems(DIT_file, meta_name="RC3_DIT_CDOM_metadata.csv",
         get_doc=F, rayleigh_width="manual", rayleigh_mask = c(20, 10, 12, 10),
         raman_width="manual", type=2) ## changed the interpolation type for the raman line


run_eems(DIT_file, meta_name="RC3_DIT_CDOM_metadata.csv",
         get_doc=F)


run_eems(DBS_file, meta_name="RC3_DBS_CDOM_metadata.csv",
         get_doc=F)

##### Set workflow manually to figure out bugs
  #set the file and metadata name
    prjpath <- "~/1_Research/0_Misc/BSLE_Degradation_EEMs/RC3_DIT_CDOM"
    meta_name="RC3_DIT_CDOM_metadata.csv"
    meta_file <- paste(prjpath,"/", meta_name, sep="")

  #make clean files and rename samples
    meta <- clean_files(prjpath=prjpath, meta_file=meta_file,
                        meta_sheet = meta_sheet)

  #convert absorbance files from .dat to .csv  (gives warning of clipping rows, this should be fine though)
    abs_preprocess(prjpath=prjpath, "mixed", meta)

  #Load Data in R
    data<- load_eems(prjpath = prjpath)
    X <- data[[1]]
    X_blk <- data[[2]]
    Sabs <- data[[3]]

  #Check data with metadata, remove samples that don't have data
    test <- check_samps(meta, X, X_blk, Sabs)
    meta <- test[[1]]
    X <- test[[2]]
    X_blk <- test[[3]]
    Sabs <- test[[4]]

  ## Process the EEM's
    data_process <- eem_proccess(prjpath=prjpath, eemlist=X, blanklist=X_blk, abs=Sabs,
                                 process_file=F, meta=meta,
                                 rayleigh_width="manual", rayleigh_mask = c(20, 10, 12, 10),
                                 raman_width="manual")
    X_clean <- data_process[[1]]
    abs_clean <- data_process[[2]]


  #plot the EEM's
    plot_eems(prjpath = prjpath, meta=meta, eem=X_clean, doc_norm=F)


  #save indices
    get_indices(X_clean, abs_clean, meta, prjpath=prjpath, doc_norm="both",
                sampsascol=F, waves=NULL)

  #save raw files (eems, abs)
    save_eems(X_clean, abs_clean, meta, prjpath=prjpath)

