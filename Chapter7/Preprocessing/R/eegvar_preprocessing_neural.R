######################################
# EEGVAR preprocessing neural data
######################################

#0. Load packages
library(data.table)
library(dplyr)
library(ggplot2)

###########################################
# Preprocessing global vs regional
###########################################

# select one
regions = "all"
regions = "posterior"

################
#1. Load data
################

# Make list to store data
all_subs <- list()

# Loop through subjects and load data
for (i in 1:151) {
  # Define the file path
  file_path <- sprintf("~\\EEGvar\\signal_variability\\sub-%03d_task-multiple_freq_var_single_trial.tsv", i)
  
  # Check if the file exists
  if (!file.exists(file_path)) {
    # Print a message if the file does not exist and skip to the next iteration
    message(sprintf("File for subject %03d does not exist. Skipping to next subject.", i))
    next
  }
  
  # Load the data and add subject identifier
  sub_data <- fread(file_path)
  sub_data$subject <- i
  
  # Store the data in the list
  all_subs[[i]] <- sub_data
}

# Stack dataframes
eeg_dat1 <- rbindlist(all_subs, use.names = TRUE, fill = TRUE)

##################
#2. Preprocess
##################

# if all regions
if (regions == "all"){
# get estimates averaged over all channels
eeg_dat1<- rename(eeg_dat1, pink_exponent = "1f_exponent", pink_offset = "1f_offset")
eeg_dat1<- eeg_dat1 %>% group_by(subject, epoch) %>% mutate(item_pink = mean(pink_exponent)) %>%
  mutate(item_entropy = mean(spectral_entropy)) %>%
  mutate(item_offset = mean(pink_offset))

eeg_dat1 <- eeg_dat1[which(eeg_dat1$cue_type == "0"),] #only keeps first row
hist(eeg_dat1$item_pink)

# remove participants with less than 500 items
eeg_dat1 <- eeg_dat1[!(eeg_dat1$subject %in% c("27","30")),]

# create trial counter
eeg_dat1 <- eeg_dat1 %>% group_by(subject) %>%
  mutate(trialnr = (rep(1:100, each = 5)))

# if posterior regions
# 20 & 21 parietal / 28 & 29 lateral occipital
} else if(regions == "posterior"){
  # get estimates averaged over 5 posterior channels
  eeg_dat1 <- rename(eeg_dat1, pink_exponent = "1f_exponent", pink_offset = "1f_offset")
  eeg_dat1 <- eeg_dat1 %>% filter(cue_type %in% c(22, 19, 11, 21, 20, 28, 29, 30)) %>% 
    group_by(subject, epoch) %>% 
    mutate(item_pink = mean(pink_exponent)) %>%
    mutate(item_entropy = mean(spectral_entropy)) %>%
    mutate(item_offset = mean(pink_offset))
  
  eeg_dat1 <- eeg_dat1[which(eeg_dat1$cue_type == "30"),] #only keeps first row
  hist(eeg_dat1$item_pink)
  
  # remove participants with less than 500 items
  eeg_dat1 <- eeg_dat1[!(eeg_dat1$subject %in% c("27","30")),]
  
  # create trial counter
  eeg_dat1 <- eeg_dat1 %>% group_by(subject) %>%
    mutate(trialnr = (rep(1:100, each = 5)))
}

####################
#3. Visualization
####################

sub1 = eeg_dat1 %>% 
        filter(subject %in% c(1))

# Distributions
hist(sub1$item_entropy)
hist(sub1$item_pink)
hist(sub1$item_offset)

# Traceplot
p1 <- ggplot(data = sub1, aes(y=item_offset, x=epoch, group = as.factor(trialnr))) +
  geom_point(color = sub1$trialnr)+
  geom_line(color = sub1$trialnr)
p1 <- p1 + geom_hline(yintercept=mean(sub1$item_offset, na.rm = TRUE), size=1, linetype = 'dotted', col = 'red')

#Add mean lines for each day, for 1 subject
p1 <- p1 + stat_smooth(data = sub1, aes(y= item_offset, x = epoch, group = trialnr), method="lm", formula=y~1, se=FALSE, color = "black")
p1 <- p1 + labs(x="Item number", y = "1/f offset")+
  theme_classic(base_size =  50)
p1

# save image of 1f noise per item
ggsave(p1, file = "~\\EEGvar\\Visualization\\trace_sub1_offset.png",
       dpi = 300, height = 10 , width = 40)

####################################################
#4. Write to csv
####################################################
write.csv(eeg_dat1, "~\\EEGvar\\Data\\eegnoise_posterior.csv", row.names = FALSE)
