## set correct path to results
library(readr)
mypath="Results/"
# setwd(mypath)

## function to combine data
combine_dat <- function(ptrn="SIM_infprior_n300_POS_sd005_iter*",   ## ptrn should be the common prefix INCLUDING _iter, followed by asterisk
                        orderby="iter_no"){ ## name of a column by which to order the results (ie iteration number)
  ## set to NA to not reorder
  
  ## get list of files
  txt_files_ls = list.files(path=mypath, pattern=ptrn)
  
  ## load files
  txt_files_df <- lapply(paste0(mypath,txt_files_ls), function(x) {read_csv(file = x)})
  
  ## combine
  combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame))
  
  ## save
  filename=substr(ptrn,1,nchar(ptrn)-6) ## remove the "_iter*"
  if(!is.na(orderby)){ ## reorder by a columnname
    combined_df <- combined_df[order(combined_df[,orderby]),]
  }
  write.csv(combined_df,file=paste0(mypath,"Full/",filename,".csv"),row.names=F)## Added "Full/" to put files in a separate folder, making things a little tidier
  
}


# combine_dat("bkmr_infprior_n300_POS_sd005_iter*")
# combine_dat("SIM_infprior_n300_POS_sd005_iter*")
# combine_dat("SIMpos_infprior_n300_POS_sd005_iter*")
# combine_dat("SIMrank_infprior_n300_POS_sd005_iter*")
# combine_dat("SIMdir_infprior_n300_POS_sd005_iter*")
# combine_dat("SIMdirequal_infprior_n300_POS_sd005_iter*")
#
# combine_dat("bkmr_infprior_n300_POSwrong_sd005_iter*")
# combine_dat("SIM_infprior_n300_POSwrong_sd005_iter*")
# combine_dat("SIMpos_infprior_n300_POSwrong_sd005_iter*")
# combine_dat("SIMrank_infprior_n300_POSwrong_sd005_iter*")
# combine_dat("SIMdir_infprior_n300_POSwrong_sd005_iter*")
# combine_dat("SIMdirequal_infprior_n300_POSwrong_sd005_iter*")
#
# combine_dat("bkmr_infprior_n300_NEG_sd005_iter*")
# combine_dat("SIM_infprior_n300_NEG_sd005_iter*")
# combine_dat("SIMpos_infprior_n300_NEG_sd005_iter*")
# combine_dat("SIMrank_infprior_n300_NEG_sd005_iter*")
# combine_dat("SIMdir_infprior_n300_NEG_sd005_iter*")
# combine_dat("SIMdirequal_infprior_n300_NEG_sd005_iter*")
#
# combine_dat("bkmr_infprior_n300_ZERO_sd005_iter*")
# combine_dat("SIM_infprior_n300_ZERO_sd005_iter*")
# combine_dat("SIMpos_infprior_n300_ZERO_sd005_iter*")
# combine_dat("SIMrank_infprior_n300_ZERO_sd005_iter*")
# combine_dat("SIMdir_infprior_n300_ZERO_sd005_iter*")
# combine_dat("SIMdirequal_infprior_n300_ZERO_sd005_iter*")

