## script used to process EEM's data with DOC data and generate exploratory plots 
  #written by Katie Wampler on 2024-11-05 
#### only partially updated to run DBS data ######

#section 0: load functions and libraries ------- 
  library(fewsdom) #to install run this code, remotes::install_github("katiewampler/fewsdom")
  library(ggplot2)
  library(stringr)
  library(readr)
  library(data.table)
  library(readxl)

#replace NA's with N/A or -9999 and round numeric to less decimals
fix_na <- function(col){
  if(class(col) %in% c("character")){
    col[is.na(col) | col== ""] <- "N/A"
  }
  if(class(col) %in% c("numeric", "logical")){
    col <- round(col, 3)
    col[is.na(col) | col== ""] <- -9999
  }
  return(col)
}

#section 1: load data ----- 
  setwd("C:/Users/russ143/PNNL/Core Richland and Sequim Lab-Field Team - Data Generation and Files/RC3/02_Processed_data_by_activity_and_analysis/DBS/")
  save_loc <- "EEM"
  meta <- read_excel("DBS_Metadata.xlsx")
  doc <- read.csv("NPOC_TN/DBS_NPOC_TN_Merged_on_2024-10-02_by_gara009.csv")

  #load raw EEM's data (we'll put processed data in the save_loc file)  
  data_loc <- "~/1_Research/0_Misc/BSLE_Degradation_EEMs/RC3_DBS_CDOM"

  #metadata for processing EEM's 
  EEM_meta <- read.csv(file.path(data_loc, "RC3_DBS_CDOM_metadata.csv"))
  
#section 2: add DOC data to metadata  ------ 
  #remove flags from DOC data 
    flagged <- grep("[|]", doc$NPOC_mg_C_per_L)
    doc$NPOC_mg_C_per_L[flagged] <- gsub("_ppm_Final_Corrected", "", str_split_i(doc$NPOC_mg_C_per_L[flagged], "[|]",i=3))

  #make numeric 
    doc$NPOC_mg_C_per_L <- as.numeric(doc$NPOC_mg_C_per_L)

  #remove underscores to match EEM's data 
    doc$Sample_ID <- gsub("_", "", doc$Sample_ID)
  
  #remove integration time in metadata description 
    EEM_meta$description <- gsub(" [1-9]s", "", EEM_meta$description)
  
  #add to metadata 
    EEM_meta <- EEM_meta %>% 
      left_join((doc %>% select(Sample_ID, NPOC_mg_C_per_L)), by = c("description" = "Sample_ID"))
    EEM_meta$DOC_mg_L <- EEM_meta$NPOC_mg_C_per_L
    EEM_meta <- EEM_meta %>% select(-NPOC_mg_C_per_L)
    
  #save meta file
    write_csv(EEM_meta, file.path(data_loc, "RC3_DBS_CDOM_metadata.csv"))
        
#section 3: process data with DOC data ------ 
  #since there's a lot of samples this step can take a few minutes
  run_eems(data_loc, meta_name="RC3_DBS_CDOM_metadata.csv",
             rayleigh_width="manual", rayleigh_mask = c(20, 10, 12, 10),
             sum_plot=F) #sum_plot removes the summary plot, which is useless with this many samples 
    
  #run some tidying things for data package 
    #add wavelength to column header for EEM's, rename and save in processed data folder  
    files <- list.files(file.path(data_loc, "5_Processed"), pattern= "s.csv")
    for(x in files){
      df <- fread(file.path(data_loc, "5_Processed", x), header =T)
      colnames(df)[1] <- "wavelength"
      new_name <- paste0(substr(x, 1,3), "_", substr(x, 4,5),"_",substr(x, 6,8), "_", substr(x,10,11), "_DilCorr_IFE_RamNorm.csv")
      fwrite(df, file.path("C:/Users/russ143/PNNL/Core Richland and Sequim Lab-Field Team - Data Generation and Files/RC3/02_Processed_data_by_activity_and_analysis/DBS", save_loc, "Fluorescence", new_name), na="NA")
    }  

    #save absorbance data individually in processed data folder
    absorb <- read.csv(file.path(data_loc, "5_Processed/Processed_Absorbance.csv"))
    wave <- absorb$wavelength
    abs_name <- gsub("[.]", " ", colnames(absorb)[-1])
    new_name <- paste0(substr(abs_name, 1,3), "_", substr(abs_name, 4,5),"_",substr(abs_name, 6,8), "_", substr(abs_name,10,11), "_DilCorr_Abs.csv")
    
    absorb <- absorb[,-1]
    for(y in 1:ncol(absorb)){
        df <- data.frame(wavelength=wave, absorbance_dil_corr=as.numeric(absorb[,y]))
        file <- new_name[y]
        write_csv(df, file.path("C:/Users/russ143/PNNL/Core Richland and Sequim Lab-Field Team - Data Generation and Files/RC3/02_Processed_data_by_activity_and_analysis/DBS", save_loc, "Absorbance", file))
      }
    
  #tidy indices file 
    fluor_index <- read_excel(file.path(data_loc, "5_Processed/SpectralIndices_RC3_DBS_CDOM.xlsx"), 
                              sheet = "fluor_indices_DOC")
    
    abs_index <- read_excel(file.path(data_loc, "5_Processed/SpectralIndices_RC3_DBS_CDOM.xlsx"), 
                              sheet = "abs_indices")

    #replace any blanks with -9999 
    fluor_index <- as.data.frame(sapply(fluor_index,fix_na))
    abs_index <- as.data.frame(sapply(abs_index,fix_na))
    
    #remove DOC data from abs
    abs_index <- abs_index %>% select(-DOC_mgL)
    
    #rename samples to match files 
    fluor_index$sample <- paste0(substr(fluor_index$sample, 1,3), "_", substr(fluor_index$sample, 4,5),"_",substr(fluor_index$sample, 6,8), "_", substr(fluor_index$sample,10,11))
    abs_index$sample <- paste0(substr(abs_index$sample, 1,3), "_", substr(abs_index$sample, 4,5),"_",substr(abs_index$sample, 6,8), "_", substr(abs_index$sample,10,11))
    
    #merge together for one file 
    indices <- fluor_index %>% left_join(abs_index, by="sample")
    
    #save in processed folder
    write_csv(indices, "EEM/DBS_CDOM_Indices_DOCnorm.csv")
    
#section 4: make plots of indices ------ 
    #read index file 
      indices <- read_csv("EEM/DBS_CDOM_Indices_DOCnorm.csv")
    
    #merge with metadata 
      #get sample 
        indices <- indices %>% mutate(int_time = substr(sample, 12,12),
                                      sample = substr(sample, 1,10)) %>% 
          left_join(meta, by=c("sample" = "Sample_ID"))
        
    #make treatment a factor to organize nicely 
      indices$Treatment <- factor(indices$Treatment, 
                                  levels=c("Saltwater_control","Mixture_control","Freshwater_control",
                                           "Salwater_leachate", "Mixture_leachate", "Freshwater_leachate"),
                                  ordered = T)
    #make plots and save
      plot_dit <- function(index){
        plot_data <- indices %>% select(any_of(c("Treatment","Day", index))) %>%
          rename("metric" = index) %>% group_by(Treatment, Day) %>% 
          summarise(mean = mean(metric),
                    sd = sd(metric))
        
        p1 <- ggplot(plot_data,aes(x=Day, y=mean)) + geom_point() + 
          geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd)), width=0) + 
          facet_wrap(~Treatment) + 
          labs(x="Day", y=index) 
       
        png(paste0("~/1_Research/0_Misc/BSLE_Degradation_EEMs/rcsfa-RC3-BSLE-DIT-degradation/plots/DIT_",index,".png"), res=300, height = 15, width=20, units="cm")
        print(p1)
        dev.off()
      }
      metrics <- colnames(indices)[!(colnames(indices) %in% c("sample", "pE",
                                                              "int_time","Treatment","Day",
                                                              "Replicate_number"))]
      sapply(metrics, plot_dit)
      
