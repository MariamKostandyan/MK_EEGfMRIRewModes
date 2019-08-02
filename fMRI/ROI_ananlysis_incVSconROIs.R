# Author: MK
###############################################################################
###############################################################################
###############################################################################

rm(list=ls()) # clear the working directory

#--------------importing libraries---------------
library(plyr)
library(dplyr)
library(base)
library(gtools)
library(ggplot2)
library(useful) # for shifting the columns 
library(stringr)
library(data.table)
library(grid)
library(Rmisc)
library(nlme)
library(lme4)
library(ez)
library(schoRsch)

# if the packages don't work 
# remove.packages(c("stringr", "data.table"))
# install.packages('stringr', dependencies = TRUE)
# install.packages('data.table', dependencies = TRUE)

#---------------setting up the path---------------

myPath <- "C:/Users/mariam/Documents/PhD/EXPERIMENTS/MK_MREEGRewModes_2016/Analysis/fMRIdata/ROI/inc_con/v6_SPMcluster/rIFG"
setwd(myPath)

#--------------create a list of files-------------
fs <- list.files(myPath, pattern = glob2rx("*.txt"))

CSRA_Rinc_mean = NULL
CSRA_Rcon_mean = NULL
CSRA_Uinc_mean = NULL
CSRA_Ucon_mean = NULL

SRA_Rinc_mean = NULL
SRA_Rcon_mean = NULL
SRA_Uinc_mean = NULL
SRA_Ucon_mean = NULL

MID_Rinc_mean = NULL
MID_Rcon_mean = NULL
MID_Uinc_mean = NULL
MID_Ucon_mean = NULL

CSRA_con_mean = NULL
CSRA_inc_mean = NULL

SRA_con_mean = NULL
SRA_inc_mean = NULL

MID_con_mean = NULL
MID_inc_mean = NULL

Neut_con_mean = NULL
Neut_inc_mean = NULL

CSRA_cue_mean = NULL
MID_cueR_mean = NULL
MID_cueU_mean = NULL

Cingulum_Ant_L = NULL


#-------------creating a loop for the cleaning of files------------
for (f in fs) {
  fname <- file.path(myPath, f)               ## current file name
  df <- read.table(fname, sep=",")     # read in file
  df <- str_split_fixed(df$V1, "\t", 2) # split into two colums 
  df <- as.data.frame(df) # bring back to the data frame
  
  CSRA_R_con_match <- c("^CSRA_Rcon") # choosing specific pattern to search through conditions
  CSRA_R_con <- df[grepl(paste(CSRA_R_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_R_con$V1 <- as.numeric(as.character(CSRA_R_con$V1)) 
  mean_CSRA_R_con <- mean(CSRA_R_con$V1) # calculating mean
  CSRA_Rcon_mean = rbind(CSRA_Rcon_mean, data.frame(mean_CSRA_R_con))
  
  CSRA_R_inc_match <- c("^CSRA_Rinc") # choosing specific pattern to search through conditions
  CSRA_R_inc <- df[grepl(paste(CSRA_R_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_R_inc$V1 <- as.numeric(as.character(CSRA_R_inc$V1)) 
  mean_CSRA_R_inc <- mean(CSRA_R_inc$V1) # calculating mean
  CSRA_Rinc_mean = rbind(CSRA_Rinc_mean, data.frame(mean_CSRA_R_inc))

  CSRA_U_con_match <- c("^CSRA_Ucon") # choosing specific pattern to search through conditions
  CSRA_U_con <- df[grepl(paste(CSRA_U_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_U_con$V1 <- as.numeric(as.character(CSRA_U_con$V1)) 
  mean_CSRA_U_con <- mean(CSRA_U_con$V1) # calculating mean
  CSRA_Ucon_mean = rbind(CSRA_Ucon_mean, data.frame(mean_CSRA_U_con))
  
  CSRA_U_inc_match <- c("^CSRA_Uinc") # choosing specific pattern to search through conditions
  CSRA_U_inc <- df[grepl(paste(CSRA_U_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_U_inc$V1 <- as.numeric(as.character(CSRA_U_inc$V1)) 
  mean_CSRA_U_inc <- mean(CSRA_U_inc$V1) # calculating mean
  CSRA_Uinc_mean = rbind(CSRA_Uinc_mean, data.frame(mean_CSRA_U_inc))

  
  SRA_R_con_match <- c("^SRA_Rcon") # choosing specific pattern to search through conditions
  SRA_R_con <- df[grepl(paste(SRA_R_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_R_con$V1 <- as.numeric(as.character(SRA_R_con$V1)) 
  mean_SRA_R_con <- mean(SRA_R_con$V1) # calculating mean
  SRA_Rcon_mean = rbind(SRA_Rcon_mean, data.frame(mean_SRA_R_con))
  
  SRA_R_inc_match <- c("^SRA_Rinc") # choosing specific pattern to search through conditions
  SRA_R_inc <- df[grepl(paste(SRA_R_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_R_inc$V1 <- as.numeric(as.character(SRA_R_inc$V1)) 
  mean_SRA_R_inc <- mean(SRA_R_inc$V1) # calculating mean
  SRA_Rinc_mean = rbind(SRA_Rinc_mean, data.frame(mean_SRA_R_inc))
  
  SRA_U_con_match <- c("^SRA_Ucon") # choosing specific pattern to search through conditions
  SRA_U_con <- df[grepl(paste(SRA_U_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_U_con$V1 <- as.numeric(as.character(SRA_U_con$V1)) 
  mean_SRA_U_con <- mean(SRA_U_con$V1) # calculating mean
  SRA_Ucon_mean = rbind(SRA_Ucon_mean, data.frame(mean_SRA_U_con))
  
  SRA_U_inc_match <- c("^SRA_Uinc") # choosing specific pattern to search through conditions
  SRA_U_inc <- df[grepl(paste(SRA_U_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_U_inc$V1 <- as.numeric(as.character(SRA_U_inc$V1)) 
  mean_SRA_U_inc <- mean(SRA_U_inc$V1) # calculating mean
  SRA_Uinc_mean = rbind(SRA_Uinc_mean, data.frame(mean_SRA_U_inc))
  
  
  MID_R_con_match <- c("^MID_R_rewcol_con", "^MID_R_nrcol_con") # choosing specific pattern to search through conditions
  MID_R_con <- df[grepl(paste(MID_R_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_R_con$V1 <- as.numeric(as.character(MID_R_con$V1)) 
  mean_MID_R_con <- mean(MID_R_con$V1) # calculating mean
  MID_Rcon_mean = rbind(MID_Rcon_mean, data.frame(mean_MID_R_con))
  
  MID_R_inc_match <- c("^MID_R_rewcol_inc", "^MID_R_nrcol_inc") # choosing specific pattern to search through conditions
  MID_R_inc <- df[grepl(paste(MID_R_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_R_inc$V1 <- as.numeric(as.character(MID_R_inc$V1)) 
  mean_MID_R_inc <- mean(MID_R_inc$V1) # calculating mean
  MID_Rinc_mean = rbind(MID_Rinc_mean, data.frame(mean_MID_R_inc))
  
  MID_U_con_match <- c("^MID_U_rewcol_con", "^MID_U_nrcol_con") # choosing specific pattern to search through conditions
  MID_U_con <- df[grepl(paste(MID_U_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_U_con$V1 <- as.numeric(as.character(MID_U_con$V1)) 
  mean_MID_U_con <- mean(MID_U_con$V1) # calculating mean
  MID_Ucon_mean = rbind(MID_Ucon_mean, data.frame(mean_MID_U_con))
  
  MID_U_inc_match <- c("^MID_U_rewcol_inc", "^MID_U_nrcol_inc") # choosing specific pattern to search through conditions
  MID_U_inc <- df[grepl(paste(MID_U_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_U_inc$V1 <- as.numeric(as.character(MID_U_inc$V1)) 
  mean_MID_U_inc <- mean(MID_U_inc$V1) # calculating mean
  MID_Uinc_mean = rbind(MID_Uinc_mean, data.frame(mean_MID_U_inc))
  
  
  Neut_con_match <- c("^NEUT_rewcol_con", "^NEUT_nrcol_con") # choosing specific pattern to search through conditions
  Neut_con <- df[grepl(paste(Neut_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  Neut_con$V1 <- as.numeric(as.character(Neut_con$V1)) 
  mean_Neut_con <- mean(Neut_con$V1) # calculating mean
  Neut_con_mean = rbind(Neut_con_mean, data.frame(mean_Neut_con))
  
  Neut_inc_match <- c("^NEUT_rewcol_inc", "^NEUT_nrcol_inc") # choosing specific pattern to search through conditions
  Neut_inc <- df[grepl(paste(Neut_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  Neut_inc$V1 <- as.numeric(as.character(Neut_inc$V1)) 
  mean_Neut_inc <- mean(Neut_inc$V1) # calculating mean
  Neut_inc_mean = rbind(Neut_inc_mean, data.frame(mean_Neut_inc))
  
  
  CSRA_cue_match <- c("^CSRA_cue") # choosing specific pattern to search through conditions
  CSRA_cue <- df[grepl(paste(CSRA_cue_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_cue$V1 <- as.numeric(as.character(CSRA_cue$V1)) 
  mean_CSRA_cue <- mean(CSRA_cue$V1) # calculating mean
  CSRA_cue_mean = rbind(CSRA_cue_mean, data.frame(mean_CSRA_cue))
  
  MID_cueR_match <- c("MID_cueR") # choosing specific pattern to search through conditions
  MID_cueR <- df[grepl(paste(MID_cueR_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_cueR$V1 <- as.numeric(as.character(MID_cueR$V1)) 
  mean_MID_cueR <- mean(MID_cueR$V1) # calculating mean
  MID_cueR_mean = rbind(MID_cueR_mean, data.frame(mean_MID_cueR))
  
  MID_cueU_match <- c("MID_cueU") # choosing specific pattern to search through conditions
  MID_cueU <- df[grepl(paste(MID_cueU_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_cueU$V1 <- as.numeric(as.character(MID_cueU$V1)) 
  mean_MID_cueU <- mean(MID_cueU$V1) # calculating mean
  MID_cueU_mean = rbind(MID_cueU_mean, data.frame(mean_MID_cueU))
  
  
  ROI <- data.frame(CSRA_Rcon_mean, CSRA_Uinc_mean, CSRA_Ucon_mean, CSRA_Rinc_mean, 
                    SRA_Rcon_mean, SRA_Uinc_mean, SRA_Ucon_mean, SRA_Rinc_mean, 
                    MID_Rcon_mean, MID_Uinc_mean, MID_Ucon_mean, MID_Rinc_mean,
                    Neut_con_mean, Neut_inc_mean, 
                    CSRA_cue_mean, MID_cueR_mean, MID_cueU_mean)
  
  ROI <- plyr::rename(ROI, c("mean_CSRA_R_con"="CSRA_R_con", "mean_CSRA_U_con"="CSRA_U_con","mean_CSRA_U_inc"="CSRA_U_inc", "mean_CSRA_R_inc"="CSRA_R_inc",
                             "mean_SRA_R_con"="SRA_R_con", "mean_SRA_U_inc"="SRA_U_inc","mean_SRA_U_con"="SRA_U_con", "mean_SRA_R_inc"="SRA_R_inc",
                             "mean_MID_R_con"="MID_R_con", "mean_MID_U_inc"="MID_U_inc","mean_MID_U_con"="MID_U_con", "mean_MID_R_inc"="MID_R_inc",
                             "mean_Neut_con"="Neut_con", "mean_Neut_inc"="Neut_inc",
                             "mean_CSRA_cue"="CSRA_cue","mean_MID_cueR"="MID_cueR", "mean_MID_cueU"="MID_cueU"))
  
}



CSRA_Rcon_mean$condition <- "CSRA"
CSRA_Rcon_mean <- plyr::rename(CSRA_Rcon_mean, c("mean_CSRA_R_con"="betas"))
CSRA_Rcon_mean$Subject <- 1:nrow(CSRA_Rcon_mean) 
CSRA_Rcon_mean$Reward <- "reward"
CSRA_Rcon_mean$Congruency <- "congruent"

CSRA_Uinc_mean$condition <- "CSRA"
CSRA_Uinc_mean <- plyr::rename(CSRA_Uinc_mean, c("mean_CSRA_U_inc"="betas"))
CSRA_Uinc_mean$Subject <- 1:nrow(CSRA_Uinc_mean) 
CSRA_Uinc_mean$Reward <- "noreward"
CSRA_Uinc_mean$Congruency <- "incongruent"

CSRA_Ucon_mean$condition <- "CSRA"
CSRA_Ucon_mean <- plyr::rename(CSRA_Ucon_mean, c("mean_CSRA_U_con"="betas"))
CSRA_Ucon_mean$Subject <- 1:nrow(CSRA_Ucon_mean) 
CSRA_Ucon_mean$Congruency <- "congruent"
CSRA_Ucon_mean$Reward <- "noreward"


CSRA_Rinc_mean$condition <- "CSRA"
CSRA_Rinc_mean <- plyr::rename(CSRA_Rinc_mean, c("mean_CSRA_R_inc"="betas"))
CSRA_Rinc_mean$Subject <- 1:nrow(CSRA_Rinc_mean) 
CSRA_Rinc_mean$Congruency <- "incongruent"
CSRA_Rinc_mean$Reward <- "reward"


SRA_Rcon_mean$condition <- "SRA"
SRA_Rcon_mean <- plyr::rename(SRA_Rcon_mean, c("mean_SRA_R_con"="betas"))
SRA_Rcon_mean$Subject <- 1:nrow(SRA_Rcon_mean) 
SRA_Rcon_mean$Reward <- "reward"
SRA_Rcon_mean$Congruency <- "congruent"

SRA_Uinc_mean$condition <- "SRA"
SRA_Uinc_mean <- plyr::rename(SRA_Uinc_mean, c("mean_SRA_U_inc"="betas"))
SRA_Uinc_mean$Subject <- 1:nrow(SRA_Uinc_mean) 
SRA_Uinc_mean$Reward <- "noreward"
SRA_Uinc_mean$Congruency <- "incongruent"

SRA_Ucon_mean$condition <- "SRA"
SRA_Ucon_mean <- plyr::rename(SRA_Ucon_mean, c("mean_SRA_U_con"="betas"))
SRA_Ucon_mean$Subject <- 1:nrow(SRA_Ucon_mean) 
SRA_Ucon_mean$Congruency <- "congruent"
SRA_Ucon_mean$Reward <- "noreward"


SRA_Rinc_mean$condition <- "SRA"
SRA_Rinc_mean <- plyr::rename(SRA_Rinc_mean, c("mean_SRA_R_inc"="betas"))
SRA_Rinc_mean$Subject <- 1:nrow(SRA_Rinc_mean) 
SRA_Rinc_mean$Congruency <- "incongruent"
SRA_Rinc_mean$Reward <- "reward"

MID_Rcon_mean$condition <- "MID"
MID_Rcon_mean <- plyr::rename(MID_Rcon_mean, c("mean_MID_R_con"="betas"))
MID_Rcon_mean$Subject <- 1:nrow(MID_Rcon_mean) 
MID_Rcon_mean$Reward <- "reward"
MID_Rcon_mean$Congruency <- "congruent"

MID_Uinc_mean$condition <- "MID"
MID_Uinc_mean <- plyr::rename(MID_Uinc_mean, c("mean_MID_U_inc"="betas"))
MID_Uinc_mean$Subject <- 1:nrow(MID_Uinc_mean) 
MID_Uinc_mean$Reward <- "noreward"
MID_Uinc_mean$Congruency <- "incongruent"

MID_Ucon_mean$condition <- "MID"
MID_Ucon_mean <- plyr::rename(MID_Ucon_mean, c("mean_MID_U_con"="betas"))
MID_Ucon_mean$Subject <- 1:nrow(MID_Ucon_mean) 
MID_Ucon_mean$Congruency <- "congruent"
MID_Ucon_mean$Reward <- "noreward"


MID_Rinc_mean$condition <- "MID"
MID_Rinc_mean <- plyr::rename(MID_Rinc_mean, c("mean_MID_R_inc"="betas"))
MID_Rinc_mean$Subject <- 1:nrow(MID_Rinc_mean) 
MID_Rinc_mean$Congruency <- "incongruent"
MID_Rinc_mean$Reward <- "reward"


Neut_con_mean$condition <- "Neut"
Neut_con_mean <- plyr::rename(Neut_con_mean, c("mean_Neut_con"="betas"))
Neut_con_mean$Subject <- 1:nrow(Neut_con_mean) 
Neut_con_mean$Congruency <- "congruent"

Neut_inc_mean$condition <- "Neut"
Neut_inc_mean <- plyr::rename(Neut_inc_mean, c("mean_Neut_inc"="betas"))
Neut_inc_mean$Subject <- 1:nrow(Neut_inc_mean) 
Neut_inc_mean$Congruency <- "incongruent"

CSRA_cue_mean$condition <- "CSRA_cue"
CSRA_cue_mean <- plyr::rename(CSRA_cue_mean, c("mean_CSRA_cue"="betas"))
CSRA_cue_mean$Subject <- 1:nrow(CSRA_cue_mean) 

MID_cueR_mean$condition <- "MID_cueR"
MID_cueR_mean <- plyr::rename(MID_cueR_mean, c("mean_MID_cueR"="betas"))
MID_cueR_mean$Subject <- 1:nrow(MID_cueR_mean) 

MID_cueU_mean$condition <- "MID_cueU"
MID_cueU_mean <- plyr::rename(MID_cueU_mean, c("mean_MID_cueU"="betas"))
MID_cueU_mean$Subject <- 1:nrow(MID_cueU_mean) 

#------------------- creating boxplots of the signal change means for all conditions and all ROIs --------------------------

boxplot(ROI, use.cols = TRUE, 
        # main = "Signal Change Values in rACCl\nROI from contrast R>U", 
        col="darkgreen", ylim=c(-3,5), las=2)

# ----------------------- STATS----------------
##### block  (Rew*Block)
MID_cueR_mean$condition <- "MID_cue"
MID_cueR_mean$Reward <- 'reward'
MID_cueU_mean$condition <- "MID_cue"
MID_cueU_mean$Reward <- 'noreward'

block <- rbind(CSRA_Rcon_mean, CSRA_Ucon_mean, CSRA_Rinc_mean, CSRA_Uinc_mean, SRA_Rcon_mean, SRA_Ucon_mean, SRA_Rinc_mean, SRA_Uinc_mean, MID_Rcon_mean, MID_Ucon_mean, MID_Rinc_mean, MID_Uinc_mean)
#block <- aggregate(data=block, betas~condition*Congruency*Subject, mean) 
#block <- rbind(block, Neut_inc_mean, Neut_con_mean)

demoAnova <- ezANOVA(block, # specify data frame
                     dv = betas, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(condition, Reward, Congruency), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

sum1 = summarySE(block, measurevar="betas", groupvars=c("Reward"))
sum2 = summarySE(block, measurevar="betas", groupvars=c("Reward", "condition"))


pairwise.t.test(block$betas, block$condition, p.adjust.method = "none"  )

boxplot(betas~condition, data=block)

# 2-way
df_R <- subset(block, Reward=="reward")
df_U <- subset(block, Reward=="noreward")
df_R$rew_effect <- df_R$betas-df_U$betas

dataCON <- subset(block, Congruency == "congruent")
dataINC <- subset(block, Congruency == "incongruent")
dataCON$con_effect <- dataCON$betas-dataINC$betas

dataCON <- subset(df_R, Congruency == "congruent")
dataINC <- subset(df_R, Congruency == "incongruent")
t.test(dataCON$rew_effect, dataINC$rew_effect, paired = TRUE) # variables are numeric 
sum3 = summarySE(df_R, measurevar="rew_effect", groupvars=c("Reward", "Congruency"))


dataSRA <- subset(df_R, condition == "SRA")
dataCSRA <- subset(df_R, condition == "CSRA")
dataMID <- subset(df_R, condition == "MID")
t.test(dataSRA$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataSRA$rew_effect, dataMID$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataMID$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 

#3-way
CON <- subset(block, Congruency == "congruent")
INC <- subset(block, Congruency == "incongruent")

df_R <- subset(INC, Reward=="reward")
df_U <- subset(INC, Reward=="noreward")
df_R$rew_effect <- df_R$betas-df_U$betas


dataSRA <- subset(df_R, condition == "SRA")
dataCSRA <- subset(df_R, condition == "CSRA")
dataMID <- subset(df_R, condition == "MID")
t.test(dataSRA$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataSRA$rew_effect, dataMID$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataMID$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 



# barplot
Rcon <- subset(block, Reward=="reward" & Congruency=="congruent")
Rinc <- subset(block, Reward=="reward" & Congruency=="incongruent")
Ucon <- subset(block, Reward=="noreward" & Congruency=="congruent")
Uinc <- subset(block, Reward=="noreward" & Congruency=="incongruent")


Rcon$Condition <- 'Rew_con'
Rinc$Condition <- 'Rew_inc'
Ucon$Condition <- "NoRew_con"
Uinc$Condition <- 'NoRew_inc'


df3 <- rbind(Rcon, Rinc, Ucon, Uinc)

# C-SRA + MID + SRA
sum = summarySEwithin(df3, measurevar="betas", withinvars=c("condition", "Condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)

sum$condition <- revalue(sum$condition , c("CSRA"="C-SRA"))
sum$condition <- factor(sum$condition, levels = c("C-SRA", "SRA", "MID"))
sum$Condition <- factor(sum$Condition, levels = c("Rew_con", "NoRew_con", "Rew_inc", "NoRew_inc"))


ggplot(sum, 
       aes(x=condition, y=betas, group=Condition, 
           ymax=betas+se, ymin=betas-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Condition)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black") +
  labs(x = "",
       y = "Parameter estimate ?? [AU]")  +
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5)) +
  coord_cartesian(ylim = c(0, 2))+
  # scale_fill_brewer()+
  ggtitle("right Inferior Frontal gyrus") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
    theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_fill_grey()
  scale_fill_manual(name="", # Legend label, use darker colors
                    values=c("#CC6666", "#9999CC", "#66CC99", "blanchedalmond"))
  
# barplot 2
sum = summarySE(dataCON, measurevar="con_effect", groupvars=c("condition"))

sum$condition <- revalue(sum$condition , c("CSRA"="C-SRA"))

ggplot(sum, 
       aes(x=condition, y=con_effect, 
           ymax=con_effect+se, ymin=con_effect-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=condition)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Parameter estimate ?? [AU]")  +
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5)) +
  coord_cartesian(ylim = c(-1,0))+
  #scale_fill_brewer()+
  scale_fill_grey() +
  ggtitle("left Precentral gyrus") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



################# 4x2 ANOVA #######################
# re-run everyting before 3x2x2 ANOVA

# Reward betas 
CSRA_R <- rbind(CSRA_Rcon_mean, CSRA_Rinc_mean)
SRA_R <- rbind(SRA_Rcon_mean, SRA_Rinc_mean)
MID_R <- rbind(MID_Rcon_mean, MID_Rinc_mean)
Neut_R <- rbind(Neut_inc_mean, Neut_con_mean)

beta_CSRA <- aggregate(data=CSRA_R, betas~condition*Congruency*Subject, mean) 
beta_SRA <- aggregate(data=SRA_R, betas~condition*Congruency*Subject, mean) 
beta_MID <- aggregate(data=MID_R, betas~condition*Congruency*Subject, mean) 
beta_Neut <- aggregate(data=Neut_R, betas~condition*Congruency*Subject, mean) 

df <- rbind(beta_CSRA, beta_SRA, beta_MID, beta_Neut)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = betas, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .( condition, Congruency), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

mean(df$betas[df$Congruency == "congruent"])
sd(df$betas[df$Congruency == "congruent"])

mean(df$betas[df$Congruency == "incongruent"])
sd(df$betas[df$Congruency == "incongruent"])

pairwise.t.test(df$betas, df$condition, p.adjust.method = "hochberg"  )

# barplot
sum = summarySEwithin(df, measurevar="betas", withinvars=c("condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)


ggplot(sum, 
       aes(x=condition, y=betas, 
           ymax=betas+se, ymin=betas-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=condition)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Block",
       y = "Parameter estimate [??]")  +
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5)) +
  coord_cartesian(ylim = c(-1,2.5))+
  #scale_fill_brewer()+
  scale_fill_grey() +
  ggtitle("left Inferior Frontal Gyrus triangular") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



# Reward betas 
CSRA_U <- rbind(CSRA_Ucon_mean, CSRA_Uinc_mean)
SRA_U <- rbind(SRA_Ucon_mean, SRA_Uinc_mean)
MID_U <- rbind(MID_Ucon_mean, MID_Uinc_mean)
Neut_U <- rbind(Neut_inc_mean, Neut_con_mean)

beta_CSRA <- aggregate(data=CSRA_U, betas~condition*Congruency*Subject, mean) 
beta_SRA <- aggregate(data=SRA_U, betas~condition*Congruency*Subject, mean) 
beta_MID <- aggregate(data=MID_U, betas~condition*Congruency*Subject, mean) 
beta_Neut <- aggregate(data=Neut_U, betas~condition*Congruency*Subject, mean) 

df <- rbind(beta_CSRA, beta_SRA, beta_MID, beta_Neut)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = betas, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(Congruency, condition), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

mean(df$betas[df$Congruency == "congruent"])
sd(df$betas[df$Congruency == "congruent"])

mean(df$betas[df$Congruency == "incongruent"])
sd(df$betas[df$Congruency == "incongruent"])

pairwise.t.test(df$betas, df$condition, p.adjust.method = "hochberg"  )


df_Con <- subset(df, Congruency=="congruent")
df_Inc <- subset(df, Congruency=="incongruent")
df_Con$inc_effect <- df_Con$betas-df_Inc$betas

dataSRA <- subset(df_Con, condition == "SRA")
dataCSRA <- subset(df_Con, condition == "CSRA")
dataMID <- subset(df_Con, condition == "MID")
dataNEUT <- subset(df_Con, condition == "Neut")
t.test(dataMID$inc_effect, dataSRA$inc_effect, paired = TRUE) # variables are numeric
t.test(dataMID$inc_effect, dataCSRA$inc_effect, paired = TRUE) # variables are numeric
t.test(dataMID$inc_effect, dataNEUT$inc_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$inc_effect, dataCSRA$inc_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$inc_effect, dataNEUT$inc_effect, paired = TRUE) # variables are numeric
t.test(dataCSRA$inc_effect, dataNEUT$inc_effect, paired = TRUE) # variables are numeric

# barplot
sum = summarySEwithin(df, measurevar="betas", withinvars=c("condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)


ggplot(sum, 
       aes(x=condition, y=betas, 
           ymax=betas+se, ymin=betas-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=condition)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Block",
       y = "Parameter estimate [??]")  +
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5)) +
  coord_cartesian(ylim = c(-1,2.5))+
  #scale_fill_brewer()+
  scale_fill_grey() +
  ggtitle("left Inferior Frontal Gyrus triangular") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



