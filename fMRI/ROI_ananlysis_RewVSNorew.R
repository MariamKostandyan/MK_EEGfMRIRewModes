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

myPath <- "C:/Users/mariam/Documents/PhD/EXPERIMENTS/MK_MREEGRewModes_2016/Analysis/fMRIdata/ROI/dopaminergic_midbrain/VTA_SN_bilat"
setwd(myPath)

#--------------create a list of files-------------
fs <- list.files(myPath, pattern = glob2rx("*.txt"))

CSRA_R_mean = NULL
CSRA_U_mean = NULL

SRA_R_mean = NULL
SRA_U_mean = NULL

MID_R_mean = NULL
MID_U_mean = NULL

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

  CSRA_R_match <- c("^CSRA_Rcon", "^CSRA_Rinc") # choosing specific pattern to search through conditions
  CSRA_R <- df[grepl(paste(CSRA_R_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_R$V1 <- as.numeric(as.character(CSRA_R$V1)) 
  mean_CSRA_R <- mean(CSRA_R$V1) # calculating mean
  CSRA_R_mean = rbind(CSRA_R_mean, data.frame(mean_CSRA_R))
  
  CSRA_U_match <- c("^CSRA_Ucon", "^CSRA_Uinc") # choosing specific pattern to search through conditions
  CSRA_U <- df[grepl(paste(CSRA_U_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_U$V1 <- as.numeric(as.character(CSRA_U$V1)) 
  mean_CSRA_U <- mean(CSRA_U$V1) # calculating mean
  CSRA_U_mean = rbind(CSRA_U_mean, data.frame(mean_CSRA_U))
  
  CSRA_con_match <- c("^CSRA_Rcon", "^CSRA_Ucon") # choosing specific pattern to search through conditions
  CSRA_con <- df[grepl(paste(CSRA_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_con$V1 <- as.numeric(as.character(CSRA_con$V1)) 
  mean_CSRA_con <- mean(CSRA_con$V1) # calculating mean
  CSRA_con_mean = rbind(CSRA_con_mean, data.frame(mean_CSRA_con))
  
  CSRA_inc_match <- c("^CSRA_Rinc", "^CSRA_Uinc") # choosing specific pattern to search through conditions
  CSRA_inc <- df[grepl(paste(CSRA_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  CSRA_inc$V1 <- as.numeric(as.character(CSRA_inc$V1)) 
  mean_CSRA_inc <- mean(CSRA_inc$V1) # calculating mean
  CSRA_inc_mean = rbind(CSRA_inc_mean, data.frame(mean_CSRA_inc))
  
  SRA_R_match <- c("^SRA_Rcon", "^SRA_Rinc") # choosing specific pattern to search through conditions
  SRA_R <- df[grepl(paste(SRA_R_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_R$V1 <- as.numeric(as.character(SRA_R$V1)) 
  mean_SRA_R <- mean(SRA_R$V1) # calculating mean
  SRA_R_mean = rbind(SRA_R_mean, data.frame(mean_SRA_R))
  
  SRA_U_match <- c("^SRA_Ucon", "^SRA_Uinc") # choosing specific pattern to search through conditions
  SRA_U <- df[grepl(paste(SRA_U_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_U$V1 <- as.numeric(as.character(SRA_U$V1)) 
  mean_SRA_U <- mean(SRA_U$V1) # calculating mean
  SRA_U_mean = rbind(SRA_U_mean, data.frame(mean_SRA_U))
  
  SRA_con_match <- c("^SRA_Rcon", "^SRA_Ucon") # choosing specific pattern to search through conditions
  SRA_con <- df[grepl(paste(SRA_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_con$V1 <- as.numeric(as.character(SRA_con$V1)) 
  mean_SRA_con <- mean(SRA_con$V1) # calculating mean
  SRA_con_mean = rbind(SRA_con_mean, data.frame(mean_SRA_con))
  
  SRA_inc_match <- c("^SRA_Rinc", "^SRA_Uinc") # choosing specific pattern to search through conditions
  SRA_inc <- df[grepl(paste(SRA_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  SRA_inc$V1 <- as.numeric(as.character(SRA_inc$V1)) 
  mean_SRA_inc <- mean(SRA_inc$V1) # calculating mean
  SRA_inc_mean = rbind(SRA_inc_mean, data.frame(mean_SRA_inc))
  
  
  MID_R_match <- c("^MID_R_rewcol_con", "^MID_R_nrcol_con", "^MID_R_rewcol_inc", "^MID_R_nrcol_inc") # choosing specific pattern to search through conditions
  MID_R <- df[grepl(paste(MID_R_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_R$V1 <- as.numeric(as.character(MID_R$V1)) 
  mean_MID_R <- mean(MID_R$V1) # calculating mean
  MID_R_mean = rbind(MID_R_mean, data.frame(mean_MID_R))
  
  MID_U_match <- c("^MID_U_rewcol_con", "^MID_U_nrcol_con", "^MID_U_rewcol_inc", "^MID_U_nrcol_inc") # choosing specific pattern to search through conditions
  MID_U <- df[grepl(paste(MID_U_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_U$V1 <- as.numeric(as.character(MID_U$V1)) 
  mean_MID_U <- mean(MID_U$V1) # calculating mean
  MID_U_mean = rbind(MID_U_mean, data.frame(mean_MID_U))
  
  MID_con_match <- c("^MID_R_rewcol_con", "^MID_R_nrcol_con", "^MID_U_rewcol_con", "^MID_U_nrcol_con") # choosing specific pattern to search through conditions
  MID_con <- df[grepl(paste(MID_con_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_con$V1 <- as.numeric(as.character(MID_con$V1)) 
  mean_MID_con <- mean(MID_con$V1) # calculating mean
  MID_con_mean = rbind(MID_con_mean, data.frame(mean_MID_con))
  
  MID_inc_match <- c("^MID_R_rewcol_inc", "^MID_R_nrcol_inc", "^MID_U_rewcol_inc", "^MID_U_nrcol_inc") # choosing specific pattern to search through conditions
  MID_inc <- df[grepl(paste(MID_inc_match, collapse="|"), df$V2), ] # grouping across this pattern
  MID_inc$V1 <- as.numeric(as.character(MID_inc$V1)) 
  mean_MID_inc <- mean(MID_inc$V1) # calculating mean
  MID_inc_mean = rbind(MID_inc_mean, data.frame(mean_MID_inc))
  
  
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
  
 
    ROI <- data.frame(CSRA_R_mean, CSRA_U_mean, CSRA_con_mean, CSRA_inc_mean, 
                         SRA_R_mean, SRA_U_mean, SRA_con_mean, SRA_inc_mean, 
                         MID_R_mean, MID_U_mean, MID_con_mean, MID_inc_mean,
                         Neut_con_mean, Neut_inc_mean, 
                         CSRA_cue_mean, MID_cueR_mean, MID_cueU_mean)
  
  ROI <- plyr::rename(ROI, c("mean_CSRA_R"="CSRA_R", "mean_CSRA_U"="CSRA_U","mean_CSRA_con"="CSRA_con", "mean_CSRA_inc"="CSRA_inc",
                                 "mean_SRA_R"="SRA_R", "mean_SRA_U"="SRA_U","mean_SRA_con"="SRA_con", "mean_SRA_inc"="SRA_inc",
                                 "mean_MID_R"="MID_R", "mean_MID_U"="MID_U","mean_MID_con"="MID_con", "mean_MID_inc"="MID_inc",
                                 "mean_Neut_con"="Neut_con", "mean_Neut_inc"="Neut_inc",
                                 "mean_CSRA_cue"="CSRA_cue","mean_MID_cueR"="MID_cueR", "mean_MID_cueU"="MID_cueU"))

  }
  


CSRA_R_mean$condition <- "CSRA"
CSRA_R_mean <- plyr::rename(CSRA_R_mean, c("mean_CSRA_R"="betas"))
CSRA_R_mean$Subject <- 1:nrow(CSRA_R_mean) 
CSRA_R_mean$Reward <- "reward"

CSRA_U_mean$condition <- "CSRA"
CSRA_U_mean <- plyr::rename(CSRA_U_mean, c("mean_CSRA_U"="betas"))
CSRA_U_mean$Subject <- 1:nrow(CSRA_U_mean) 
CSRA_U_mean$Reward <- "noreward"

CSRA_con_mean$condition <- "CSRA"
CSRA_con_mean <- plyr::rename(CSRA_con_mean, c("mean_CSRA_con"="betas"))
CSRA_con_mean$Subject <- 1:nrow(CSRA_con_mean) 
CSRA_con_mean$congruency <- "congruent"

CSRA_inc_mean$condition <- "CSRA"
CSRA_inc_mean <- plyr::rename(CSRA_inc_mean, c("mean_CSRA_inc"="betas"))
CSRA_inc_mean$Subject <- 1:nrow(CSRA_inc_mean) 
CSRA_inc_mean$congruency <- "incongruent"

SRA_R_mean$condition <- "SRA"
SRA_R_mean <- plyr::rename(SRA_R_mean, c("mean_SRA_R"="betas"))
SRA_R_mean$Subject <- 1:nrow(SRA_R_mean) 
SRA_R_mean$Reward <- "reward"

SRA_U_mean$condition <- "SRA"
SRA_U_mean <- plyr::rename(SRA_U_mean, c("mean_SRA_U"="betas"))
SRA_U_mean$Subject <- 1:nrow(SRA_U_mean) 
SRA_U_mean$Reward <- "noreward"

SRA_con_mean$condition <- "SRA"
SRA_con_mean <- plyr::rename(SRA_con_mean, c("mean_SRA_con"="betas"))
SRA_con_mean$Subject <- 1:nrow(SRA_con_mean) 
SRA_con_mean$congruency <- "congruent"

SRA_inc_mean$condition <- "SRA"
SRA_inc_mean <- plyr::rename(SRA_inc_mean, c("mean_SRA_inc"="betas"))
SRA_inc_mean$Subject <- 1:nrow(SRA_inc_mean) 
SRA_inc_mean$congruency <- "incongruent"

MID_R_mean$condition <- "MID"
MID_R_mean <- plyr::rename(MID_R_mean, c("mean_MID_R"="betas"))
MID_R_mean$Subject <- 1:nrow(MID_R_mean) 
MID_R_mean$Reward <- "reward"

MID_U_mean$condition <- "MID"
MID_U_mean <- plyr::rename(MID_U_mean, c("mean_MID_U"="betas"))
MID_U_mean$Subject <- 1:nrow(MID_U_mean) 
MID_U_mean$Reward <- "noreward"

MID_con_mean$condition <- "MID"
MID_con_mean <- plyr::rename(MID_con_mean, c("mean_MID_con"="betas"))
MID_con_mean$Subject <- 1:nrow(MID_con_mean) 
MID_con_mean$congruency <- "congruent"

MID_inc_mean$condition <- "MID"
MID_inc_mean <- plyr::rename(MID_inc_mean, c("mean_MID_inc"="betas"))
MID_inc_mean$Subject <- 1:nrow(MID_inc_mean) 
MID_inc_mean$congruency <- "incongruent"

Neut_con_mean$condition <- "Neut"
Neut_con_mean <- plyr::rename(Neut_con_mean, c("mean_Neut_con"="betas"))
Neut_con_mean$Subject <- 1:nrow(Neut_con_mean) 
Neut_con_mean$congruency <- "congruent"

Neut_inc_mean$condition <- "Neut"
Neut_inc_mean <- plyr::rename(Neut_inc_mean, c("mean_Neut_inc"="betas"))
Neut_inc_mean$Subject <- 1:nrow(Neut_inc_mean) 
Neut_inc_mean$congruency <- "incongruent"

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

block <- rbind(CSRA_R_mean, CSRA_U_mean, SRA_R_mean, SRA_U_mean, MID_R_mean, MID_U_mean, MID_cueR_mean, MID_cueU_mean)
block2 <- rbind(CSRA_R_mean, CSRA_U_mean, SRA_R_mean, SRA_U_mean, MID_R_mean, MID_U_mean)

demoAnova <- ezANOVA(block, # specify data frame
                     dv = betas, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(condition, Reward), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(block$betas, block$condition, p.adjust.method = "none"  )
sum = summarySEwithin(block, measurevar="betas", withinvars=c("condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)

df_R <- subset(block, Reward=="reward")
df_U <- subset(block, Reward=="noreward")
df_R$rew_effect <- df_R$betas-df_U$betas

sum = summarySEwithin(df_R, measurevar="rew_effect", withinvars=c("condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)


dataSRA <- subset(df_R, condition == "SRA")
dataCSRA <- subset(df_R, condition == "CSRA")
dataMID <- subset(df_R, condition == "MID")
dataMID_cue <- subset(df_R, condition == "MID_cue")
t.test(dataSRA$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataSRA$rew_effect, dataMID$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataSRA$rew_effect, dataMID_cue$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataMID$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataMID$rew_effect, dataMID_cue$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataCSRA$rew_effect, dataMID_cue$rew_effect, paired = TRUE) # variables are numeric 

# barplot
sum = summarySE(block, measurevar="betas", groupvars=c("condition", 'Reward'))
sum$condition <- revalue(sum$condition , c("CSRA"="C-SRA"))

ggplot(sum, 
       aes(x=condition, y=betas, group = Reward,
           ymax=betas+se, ymin=betas-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Reward)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Block",
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
  coord_cartesian(ylim = c(-1.5,1 ))+
  #scale_fill_brewer()+
  scale_fill_grey() +
  ggtitle("left Insula") +
  #theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# barplot 2
sum = summarySEwithin(df_R, measurevar="rew_effect", withinvars=c("condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)

sum$condition <- revalue(sum$condition , c("CSRA"="C-SRA"))

ggplot(sum, 
       aes(x=condition, y=rew_effect, 
           ymax=rew_effect+se, ymin=rew_effect-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=condition)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "R-minus-U difference\nParameter estimate ?? [AU]")  +
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
  coord_cartesian(ylim = c(-0.5,0.5))+
  #scale_fill_brewer()+
  scale_fill_manual(name="", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "lightblue1", "#009E73")) +
 # scale_fill_grey() +
  ggtitle("left Dopaminergic midbrain") +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# # CSRA cue barplot
# sum = summarySE(CSRA_cue_mean, measurevar="betas", groupvars=c("condition"))
# sum$condition <- revalue(sum$condition , c("CSRA_cue"="C-SRA_cue"))
# 
# ggplot(sum, 
#        aes(x=condition, y=betas, 
#            ymax=betas+se, ymin=betas-se))  +
#   geom_bar(stat="identity", position = "dodge", aes(fill=condition)) +
#   geom_errorbar(position=position_dodge(width=0.7), 
#                 width=0.0, size=0.5, color="black")  +
#   labs(x = "Block",
#        y = "Parameter estimate [??]")  +
#   theme_bw()  +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey50"),
#         plot.title = element_text(size = rel(1.5), 
#                                   face = "bold", vjust = 1.5),
#         axis.title = element_text(face = "bold"),
#         legend.key.size = unit(0.4, "cm"),
#         legend.key = element_rect(fill = "black"),
#         axis.title.y = element_text(vjust= 1.8),
#         axis.title.x = element_text(vjust= -0.5)) +
#   coord_cartesian(ylim = c(-3, 3))+
#   #scale_fill_brewer()+
#   scale_fill_grey() +
#   #theme(legend.position="none") +
#   theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 


# # cues: old analysis
# cues <- rbind(CSRA_cue_mean, MID_cueR_mean, MID_cueU_mean)
# 
# demoAnova <- ezANOVA(cues, # specify data frame
#                      dv = betas, # specify dependent variable 
#                      wid = Subject, # specify the subject variable
#                      within = .(condition), # specify within-subject variables
#                      detailed = TRUE, # get a detailed table that includes SS
#                      type = 3
# )
# demoAnova
# demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)
# 
# pairwise.t.test(cues$betas, cues$condition, p.adjust.method = "none"  )
# 
# 
# # barplot
# sum = summarySE(cues, measurevar="betas", groupvars=c("condition"))
# sum$condition <- revalue(sum$condition , c("CSRA_cue"="C-SRA_cue"))
# 
# ggplot(sum, 
#        aes(x=condition, y=betas, group = condition,
#            ymax=betas+se, ymin=betas-se))  +
#   geom_bar(stat="identity", position = "dodge", aes(fill=condition)) +
#   geom_errorbar(position=position_dodge(width=0.7), 
#                 width=0.0, size=0.5, color="black")  +
#   labs(x = "Block",
#        y = "Parameter estimate [???]")  +
#   theme_bw()  +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey50"),
#         plot.title = element_text(size = rel(1.5), 
#                                   face = "bold", vjust = 1.5),
#         axis.title = element_text(face = "bold"),
#         legend.key.size = unit(0.4, "cm"),
#         legend.key = element_rect(fill = "black"),
#         axis.title.y = element_text(vjust= 1.8),
#         axis.title.x = element_text(vjust= -0.5)) +
#   coord_cartesian(ylim = c(-0.8, -1.8))+
#   #scale_fill_brewer()+
#   scale_fill_grey() +
#     theme(legend.position="none") +
#   theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# 
# # targets - reward
# # MID_cueR_mean$condition <- "MID_cue"
# # MID_cueU_mean$condition <- "MID_cue"
# # MID_cueR_mean$Reward <- "reward"
# # MID_cueU_mean$Reward <- "noreward"
# 
# rew <- rbind(CSRA_R_mean, CSRA_U_mean, SRA_R_mean, SRA_U_mean, MID_R_mean, MID_U_mean)
# 
# rew$condition <- as.factor(rew$condition)
# rew$Subject <- as.factor(rew$Subject)
# rew$Reward <- as.factor(rew$Reward)
# rew$betas <- as.numeric(as.character(rew$betas))
# 
# demoAnova <- ezANOVA(rew, # specify data frame
#                      dv = betas, # specify dependent variable 
#                      wid = Subject, # specify the subject variable
#                      within = .(condition, Reward), # specify within-subject variables
#                      detailed = TRUE, # get a detailed table that includes SS
#                      type = 3
# )
# demoAnova
# demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)
# 
# pairwise.t.test(rew$betas, rew$condition, p.adjust.method = "none"  )
# pairwise.t.test(rew$betas, rew$Reward, p.adjust.method = "none"  )
# 
# data_rew <- subset(rew, Reward == "reward")
# demoAnova <- ezANOVA(data_rew, # specify data frame
#                      dv = betas, # specify dependent variable 
#                      wid = Subject, # specify the subject variable
#                      within = .(condition), # specify within-subject variables
#                      detailed = TRUE, # get a detailed table that includes SS
#                      type = 3
# )
# demoAnova
# demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)
# pairwise.t.test(data_rew$betas, data_rew$condition, p.adjust.method = "none"  )
# data_nr <- subset(rew, Reward == "noreward")
# demoAnova <- ezANOVA(data_nr, # specify data frame
#                      dv = betas, # specify dependent variable 
#                      wid = Subject, # specify the subject variable
#                      within = .(condition), # specify within-subject variables
#                      detailed = TRUE, # get a detailed table that includes SS
#                      type = 3
# )
# demoAnova
# demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)
# 
# rew <- subset(dataSRA, Reward =="0")
# norew <- subset(dataSRA, Reward=="1")
# t.test(rew$Freq, norew$Freq, paired = TRUE) # variables are numeric 
# 
# # barplot
# sum = summarySE(rew, measurevar="betas", groupvars=c("condition", "Reward"))
# sum$condition <- revalue(sum$condition , c("CSRA"="C-SRA"))
# 
# ggplot(sum, 
#        aes(x=condition, y=betas, group = Reward,
#            ymax=betas+se, ymin=betas-se))  +
#   geom_bar(stat="identity", position = "dodge", aes(fill=Reward)) +
#   geom_errorbar(position=position_dodge(width=0.7), 
#                 width=0.0, size=0.5, color="black")  +
#   labs(x = "Block",
#        y = "Parameter estimate [???]")  +
#   theme_bw()  +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey50"),
#         plot.title = element_text(size = rel(1.5), 
#                                   face = "bold", vjust = 1.5),
#         axis.title = element_text(face = "bold"),
#         legend.key.size = unit(0.4, "cm"),
#         legend.key = element_rect(fill = "black"),
#         axis.title.y = element_text(vjust= 1.8),
#         axis.title.x = element_text(vjust= -0.5)) +
#   #coord_cartesian(ylim = c(-0.5, -1.5))+
#   #scale_fill_brewer()+
#   scale_fill_grey() +
#   #theme(legend.position="none") +
#     theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))

