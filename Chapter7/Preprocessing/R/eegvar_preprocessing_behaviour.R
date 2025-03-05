######################################################################
######################### Preprocessing script #######################
#1. Get data from OSF format to dataframe
#2. Descriptives and transformations
######################################################################


setwd("~\\EEGvar\\Behaviour") # select the path 


# Load RawData & Specify Factors ------------------------------------------

# Code the Colum names according to the order of variables in DiplayTrial function
# of the Matlab functions for the Sternberg task
data_colnames <- c("Subject", "Session","Block", "TrialNum", "TaskDescription", "Probe",
                   "Match", "SetSize", "Position",
                   "Memset1", "Memset2", "Memset3", "Memset4", "Memset5",
                   "Fix", "ISI1", "ISI2", "ISI3", "ISI4", "ISI5", "CUE", "ITI", # presentation times of the fixation cross, the inter-stimulus-intervals, the questionmark as indicator for the following probe, and the inter-trial-intervall
                   "Accuracy", "RT","Response", "CorrResp")
data_colclasses <- c("integer", "integer", "integer","integer", "character","integer",
                     "integer", "integer", "integer",
                     "integer", "integer", "integer","integer", "integer",
                     "numeric","numeric","numeric", "numeric","numeric","numeric","numeric","numeric",
                     "integer", "numeric","character", "character")

df_full <- read.table("SternbergTaskSternbergTask_S1_Ses1exp.txt",
                      col.names = data_colnames,
                      colClasses = data_colclasses)

## Reading Filenames of RawData & read in all Experimental Data in one Dataframe ----


data_files <- list.files(pattern = "exp.txt")
data_files <- data.frame(data_files,stringsAsFactors = F)
df_full <- NULL
i <- 1
for(i in 1:nrow(data_files))
{
  datafile <- data_files[i,]
  
  
  dat <- read.table(datafile,col.names = data_colnames,
                    colClasses = data_colclasses)
  
  df_full  <- rbind(df_full,dat)
  rm(dat)
  i <- i + 1
}

##################################################
# Descriptives and transformations
##################################################

# Remove impossible values
fast <- df_full$RT < 0.1 # fast guesses
slow <- df_full$RT > 10 # slow

# sub with NA
df_full[fast & !is.na(fast)] <- NA
df_full[slow & !is.na(slow)] <- NA

# Transformation to reduce skewness
df_full$logRT <- log(df_full$RT) #log-transform RT

#Plot distributions
hist(df_full$logRT)
hist(df_full$RT)


######################################
# Write to csv
######################################

write.csv(df_full, "C:\\Users\\micha\\OneDrive\\Documents\\EEGvar\\Behaviour\\eegRT.csv", row.names = FALSE)



