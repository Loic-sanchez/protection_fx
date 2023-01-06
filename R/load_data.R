load_rdata = function() { 
  
  files <- list.files(here::here("data", "raw_data"), 
                      pattern = ".RData|Rdata", 
                      full.names = TRUE)

  load_path <- lapply(files, 
                      load, 
                      .GlobalEnv)  
  
  }

