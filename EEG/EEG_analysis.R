# Author: MK
###############################################################################
###############################################################################
###############################################################################
################# EEG analysis ###############################################

rm(list=ls()) # clear the working directory

#---------------setting up the path---------------
myPath <- "C:/Users/mariam/Documents/PhD/EXPERIMENTS/MK_MREEGRewModes_2016/Analysis/EEGdata/stats"  
setwd(myPath)

#--------------import all libraries---------------
library(dplyr)
library(plyr)
library(useful) # for shifting the columns 
library(stringr)
library(data.table)
library(ggplot2) 
library(Rmisc)
library(ez)
library(schoRsch)
library(car)
library(compute.es)
library(effects)
library(multcomp)
library(pastecs)
library(WRS2)

Rew <- read.table("N1_REW_160-200.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                           dec = ".",
                           na.strings = "NaN", colClasses = NA, nrows = -1,
                           skip = 0, check.names = TRUE,
                           strip.white = FALSE, blank.lines.skip = TRUE,
                           comment.char = "#",
                           allowEscapes = FALSE, flush = FALSE)

Norew <- read.table("N1_NoRew_160-200.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                           dec = ".",
                           na.strings = "NaN", colClasses = NA, nrows = -1,
                           skip = 0, check.names = TRUE,
                           strip.white = FALSE, blank.lines.skip = TRUE,
                           comment.char = "#",
                           allowEscapes = FALSE, flush = FALSE)

MID_cueU <- read.table("CNV_700_2300_MID_cueU.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                           dec = ".",
                           na.strings = "NaN", colClasses = NA, nrows = -1,
                           skip = 0, check.names = TRUE,
                           strip.white = FALSE, blank.lines.skip = TRUE,
                           comment.char = "#",
                           allowEscapes = FALSE, flush = FALSE)
MID_cueR <- read.table("CNV_700_2300_MID_cueR.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                       dec = ".",
                       na.strings = "NaN", colClasses = NA, nrows = -1,
                       skip = 0, check.names = TRUE,
                       strip.white = FALSE, blank.lines.skip = TRUE,
                       comment.char = "#",
                       allowEscapes = FALSE, flush = FALSE)

CSRA_cue <- read.table("CNV_700_2300_CSRA_cue.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                       dec = ".",
                       na.strings = "NaN", colClasses = NA, nrows = -1,
                       skip = 0, check.names = TRUE,
                       strip.white = FALSE, blank.lines.skip = TRUE,
                       comment.char = "#",
                       allowEscapes = FALSE, flush = FALSE)

CSRA_cue_Even <- read.table("CSRA_cue_Even.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                       dec = ".",
                       na.strings = "NaN", colClasses = NA, nrows = -1,
                       skip = 0, check.names = TRUE,
                       strip.white = FALSE, blank.lines.skip = TRUE,
                       comment.char = "#",
                       allowEscapes = FALSE, flush = FALSE)

CSRA_cue_Odd <- read.table("CSRA_cue_Odd.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                            dec = ".",
                            na.strings = "NaN", colClasses = NA, nrows = -1,
                            skip = 0, check.names = TRUE,
                            strip.white = FALSE, blank.lines.skip = TRUE,
                            comment.char = "#",
                            allowEscapes = FALSE, flush = FALSE)

Con <- read.table("Negativity_Con_550_750.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                  dec = ".",
                  na.strings = "NaN", colClasses = NA, nrows = -1,
                  skip = 0, check.names = TRUE,
                  strip.white = FALSE, blank.lines.skip = TRUE,
                  comment.char = "#",
                  allowEscapes = FALSE, flush = FALSE)

Inc <- read.table("Negativity_Inc_550_750.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                    dec = ".",
                    na.strings = "NaN", colClasses = NA, nrows = -1,
                    skip = 0, check.names = TRUE,
                    strip.white = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#",
                    allowEscapes = FALSE, flush = FALSE)


THETA_Rcon <- read.table("THETA_Rcon_300_800.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                            dec = ".",
                            na.strings = "NaN", colClasses = NA, nrows = -1,
                            skip = 0, check.names = TRUE,
                            strip.white = FALSE, blank.lines.skip = TRUE,
                            comment.char = "#",
                            allowEscapes = FALSE, flush = FALSE)

THETA_Rinc <- read.table("THETA_Rinc_300_800.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                           dec = ".",
                           na.strings = "NaN", colClasses = NA, nrows = -1,
                           skip = 0, check.names = TRUE,
                           strip.white = FALSE, blank.lines.skip = TRUE,
                           comment.char = "#",
                           allowEscapes = FALSE, flush = FALSE)

THETA_Uinc <- read.table("THETA_Uinc_300_800.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                  dec = ".",
                  na.strings = "NaN", colClasses = NA, nrows = -1,
                  skip = 0, check.names = TRUE,
                  strip.white = FALSE, blank.lines.skip = TRUE,
                  comment.char = "#",
                  allowEscapes = FALSE, flush = FALSE)

THETA_Ucon <- read.table("THETA_Ucon_300_800.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                  dec = ".",
                  na.strings = "NaN", colClasses = NA, nrows = -1,
                  skip = 0, check.names = TRUE,
                  strip.white = FALSE, blank.lines.skip = TRUE,
                  comment.char = "#",
                  allowEscapes = FALSE, flush = FALSE)


# #-------------- FFTs ---------------------
# # outlier rejection
# Rew <- subset(Rew, select=c("File", "Pz.FFT_RBlockNorm", "P3.FFT_RBlockNorm", "P4.FFT_RBlockNorm", "POz.FFT_RBlockNorm",  "PO3.FFT_RBlockNorm",  "PO4.FFT_RBlockNorm" ))
# Norew <- subset(Norew, select=c("File", "Pz.FFT_NRBlockNorm", "P3.FFT_NRBlockNorm", "P4.FFT_NRBlockNorm", "POz.FFT_NRBlockNorm", "PO3.FFT_NRBlockNorm", "PO4.FFT_NRBlockNorm"))
# 
# Rew$MeanR <- rowMeans(Rew[,2:7])
# Norew$MeanNR <- rowMeans(Norew[,2:7])
# 
# Rew <- Rew[,-(2:7),drop=FALSE]
# Norew <- Norew[,-(2:7),drop=FALSE]
# 
# qnt <- quantile(Rew$MeanR, probs=c(.25, .75))
# H <- 1.5 * IQR(Rew$MeanR)
# outlierCheck <- (Rew$MeanR) > qnt[1]-H & (Rew$MeanR<qnt[2]+H)
# noOutliers <- Rew[outlierCheck,]
# 
# Rew <- Rew[-c(7, 9), ]
# 
# qnt <- quantile(Norew$MeanNR, probs=c(.25, .75))
# H <- 1.5 * IQR(Norew$MeanNR)
# outlierCheck <- (Norew$MeanNR) > qnt[1]-H & (Norew$MeanNR<qnt[2]+H)
# noOutliers <- Norew[outlierCheck,]
# 
# Norew <- Norew[-c(7, 9, 15, 17), ]
# 
# # with outliers
# exp.3 <- merge(Rew, Norew, all = TRUE)
# with(exp.3, t.test(MeanNR ,MeanR, equal.var=TRUE) )
#   
# # no outliers
# t.test(Rew$MeanR, Norew$MeanNR, paired = T)
# 
# # ANCOVA
# df$Bl <- "Rew" 
# df2$Bl <- "NoRew" 
# df <- rename(df, c("Pz.FFT_RBlockNorm"="FFT"))
# df2 <- rename(df2, c("Pz.FFT_NRBlockNorm"="FFT"))
# 
# expKl <- rbind(df, df2)
# 
# expKl$BlOrder <- 0 # no rew - rew - no rew - rew
# #expKl$BlOrder[expKl$File %in% c("PP01_expC", "PP02_expC", "PP03_expC", "PP07_expC", "PP08_expC", "PP09_expC", "PP13_expC", "PP14_expC", "PP15_expC", "PP18_expC", "PP20_expC", "PP21_expC")] <- 1 # rew - no rew - rew - no rew
# expKl$BlOrder[expKl$File %in% c("Mariam_Subject_03", "Mariam_Subject_05", "Mariam_Subject_08", "Mariam_Subject_10", "Mariam_Subject_12", "Mariam_Subject_13", "Mariam_Subject_15", "Mariam_Subject_17", "Mariam_Subject_19", "Mariam_Subject_20", "Mariam_Subject_22", "Mariam_Subject_24")] <- 1 # rew - no rew - rew - no rew
# 
# options(contrasts = c("contr.treatment", "contr.poly"))
# model.1 = lm (FFT ~ Bl + BlOrder + Bl:BlOrder, data = expKl)
# Anova(model.1, type="II")
#        
# # many electrodes
# Rew <- subset(Rew, select=c("File", "Pz.FFT_RBlockNorm", "P3.FFT_RBlockNorm", "P4.FFT_RBlockNorm", "POz.FFT_RBlockNorm",  "PO3.FFT_RBlockNorm",  "PO4.FFT_RBlockNorm" ))
# Norew <- subset(Norew, select=c("File", "Pz.FFT_NRBlockNorm", "P3.FFT_NRBlockNorm", "P4.FFT_NRBlockNorm", "POz.FFT_NRBlockNorm", "PO3.FFT_NRBlockNorm", "PO4.FFT_NRBlockNorm"))
# 
# Rew$Mean <- rowMeans(Rew[,2:7])
# Norew$Mean <- rowMeans(Norew[,2:7])
# 
# qnt <- quantile(Rew$Mean, probs=c(.25, .75))
# H <- 1.5 * IQR(Rew$Mean)
# outlierCheck <- (Rew$Mean) > qnt[1]-H & (Rew$Mean<qnt[2]+H)
# noOutliers <- Rew[outlierCheck,]
# 
# Rew <- Rew[-c(7, 9), ]
# 
# qnt <- quantile(Norew$Mean, probs=c(.25, .75))
# H <- 1.5 * IQR(Norew$Mean)
# outlierCheck <- (Norew$Mean) > qnt[1]-H & (Norew$Mean<qnt[2]+H)
# noOutliers <- Norew[outlierCheck,]
# 
# Norew <- Norew[-c(7, 9, 15), ]
# 
# Rew$Bl <- "Rew" 
# Norew$Bl <- "NoRew" 
# Rew <- rename(Rew, c("Mean"="FFT"))
# Norew <- rename(Norew, c("Mean"="FFT"))
# 
# Rew <- Rew[,-(2:7),drop=FALSE]
# Norew <- Norew[,-(2:7),drop=FALSE]
# 
# expKl <- rbind(Rew, Norew)
# 
# expKl$BlOrder <- 0 # no rew - rew - no rew - rew
# expKl$BlOrder[expKl$File %in% c("Mariam_Subject_03", "Mariam_Subject_05", "Mariam_Subject_08", "Mariam_Subject_10", "Mariam_Subject_12", "Mariam_Subject_13", "Mariam_Subject_15", "Mariam_Subject_17", "Mariam_Subject_19", "Mariam_Subject_20", "Mariam_Subject_22", "Mariam_Subject_24")] <- 1 # rew - no rew - rew - no rew
# 
# options(contrasts = c("contr.treatment", "contr.poly"))
# model.1 = lm (FFT ~ Bl + BlOrder + Bl:BlOrder, data = expKl)
# Anova(model.1, type="II")
# 




# --------------------- THETA -------------
CSRA_Rinc <- subset(THETA_Rinc, select=c("File", "Cz.AverageTHETA_CSRA_Rinc"))
CSRA_Rcon <- subset(THETA_Rcon, select=c("File", "Cz.AverageTHETA_CSRA_Rcon"))
CSRA_Uinc <- subset(THETA_Uinc, select=c("File", "Cz.AverageTHETA_CSRA_Uinc"))
CSRA_Ucon <- subset(THETA_Ucon, select=c("File", "Cz.AveragTHETA_CSRA_Ucon"))
SRA_Rinc <- subset(THETA_Rinc, select=c("File", "Cz.AverageTHETA_SRA_Rinc"))
SRA_Rcon <- subset(THETA_Rcon, select=c("File", "Cz.AverageTHETA_SRA_Rcon"))
SRA_Uinc <- subset(THETA_Uinc, select=c("File", "Cz.AverageTHETA_SRA_Uinc"))
SRA_Ucon <- subset(THETA_Ucon, select=c("File", "Cz.AverageTHETA_SRA_Ucon"))
MID_Rinc <- subset(THETA_Rinc, select=c("File", "Cz.AverageTHETA_MID_Rinc"))
MID_Rcon <- subset(THETA_Rcon, select=c("File", "Cz.AverageTHETA_MID_Rcon"))
MID_Uinc <- subset(THETA_Uinc, select=c("File", "Cz.AverageTHETA_MID_Uinc"))
MID_Ucon <- subset(THETA_Ucon, select=c("File", "Cz.AverageTHETA_MID_Ucon"))


CSRA_Rinc$Block <- "C-SRA"
CSRA_Rcon$Block <- "C-SRA"
CSRA_Uinc$Block <- "C-SRA"
CSRA_Ucon$Block <- "C-SRA"
SRA_Rinc$Block <- "SRA"
SRA_Rcon$Block <- "SRA"
SRA_Uinc$Block <- "SRA"
SRA_Ucon$Block <- "SRA"
MID_Rinc$Block <- "MID"
MID_Rcon$Block <- "MID"
MID_Uinc$Block <- "MID"
MID_Ucon$Block <- "MID"


CSRA_Rinc$Reward <- "rew"
CSRA_Rcon$Reward <- "rew"
CSRA_Uinc$Reward <- "no rew"
CSRA_Ucon$Reward <- "no rew"
SRA_Rinc$Reward <- "rew"
SRA_Rcon$Reward <- "rew"
SRA_Uinc$Reward <- "no rew"
SRA_Ucon$Reward <- "no rew"
MID_Rinc$Reward <- "rew"
MID_Rcon$Reward <- "rew"
MID_Uinc$Reward <- "no rew"
MID_Ucon$Reward <- "no rew"


CSRA_Rinc$Congruency <- "inc"
CSRA_Rcon$Congruency <- "con"
CSRA_Uinc$Congruency <- "inc"
CSRA_Ucon$Congruency <- "con"
SRA_Rinc$Congruency <- "inc"
SRA_Rcon$Congruency <- "con"
SRA_Uinc$Congruency <- "inc"
SRA_Ucon$Congruency <- "con"
MID_Rinc$Congruency <- "inc"
MID_Rcon$Congruency <- "con"
MID_Uinc$Congruency <- "inc"
MID_Ucon$Congruency <- "con"

#rename columns 
dfs <- c("CSRA_Rinc", "CSRA_Rcon", "CSRA_Uinc", "CSRA_Ucon", "SRA_Rinc", "SRA_Rcon", "SRA_Uinc", "SRA_Ucon", "MID_Rinc", "MID_Rcon", "MID_Uinc", "MID_Ucon")
for(df in dfs) {
  df.tmp <- get(df)
  names(df.tmp) <- c("File", "Amp", "Block", "Reward", "congruency") 
  assign(df, df.tmp)
}

df <- rbind(CSRA_Rinc, CSRA_Rcon, CSRA_Uinc, CSRA_Ucon, SRA_Rinc, SRA_Rcon, SRA_Uinc, SRA_Ucon, MID_Rinc, MID_Rcon, MID_Uinc, MID_Ucon)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Amp, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Block, Reward, congruency), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$Mean, df$Cond, p.adjust.method="none", paired=TRUE)


# barplot
sum = summarySEwithin(df, measurevar="Amp", withinvars=c("Block", "Reward", "congruency"), idvar="File", na.rm=FALSE, conf.interval=.95)

ggplot(sum, 
       aes(x = Reward, y = Amp, fill = Block, 
           ymax=Amp+se, ymin=Amp-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  facet_wrap(  ~ congruency) +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Cue type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(4, 8))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))



#-------------------- CUE ---------------
MID_cueU <- subset(MID_cueU, select=c("File", "Cz.Average_MID_cueU_NEW", 'FCz.Average_MID_cueU_NEW'))
MID_cueR <- subset(MID_cueR, select=c("File", "Cz.Average_MID_cueR_NEW", 'FCz.Average_MID_cueR_NEW'))
CSRA_cue <- subset(CSRA_cue, select=c("File", "Cz.Average_CSRA_cue_NEW", 'FCz.Average_CSRA_cue_NEW'))

MID_cueU$Mean <- rowMeans(MID_cueU[,2:3])
MID_cueR$Mean <- rowMeans(MID_cueR[,2:3])
CSRA_cue$Mean <- rowMeans(CSRA_cue[,2:3])

MID_cueU$Cond <- "MID_cueU"
MID_cueR$Cond <- "MID_cueR"
CSRA_cue$Cond <- "CSRA_cue"

#MID_cueU <- plyr::rename(MID_cueU, c("Cz.Average_MID_cueU"="Ampl"))
#MID_cueR <- plyr::rename(MID_cueR, c("Cz.Average_MID_cueR"="Ampl"))
#CSRA_cue <- plyr::rename(CSRA_cue, c("Cz.Average_CSRA_cue"="Ampl"))

MID_cueU$Cond <- "MID non-rew"
MID_cueR$Cond <- "MID rew"
CSRA_cue$Cond <- "C-SRA"

MID_cueU <- MID_cueU[,-(2:3),drop=FALSE]
MID_cueR <- MID_cueR[,-(2:3),drop=FALSE]
CSRA_cue <- CSRA_cue[,-(2:3),drop=FALSE]


df <- rbind(CSRA_cue, MID_cueR, MID_cueU)

df$File <- as.factor(df$File)
df$Cond <- as.factor(df$Cond)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Cond), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$Mean, df$Cond, p.adjust.method="none", paired=TRUE)

mean(CSRA_cue$Mean)
sd(CSRA_cue$Mean)
mean(MID_cueR$Mean)
sd(MID_cueR$Mean)
mean(MID_cueU$Mean)
sd(MID_cueU$Mean)

# barplot
sum = summarySEwithin(df, measurevar="Mean", withinvars=c("Cond"), idvar="File", na.rm=FALSE, conf.interval=.95)

ggplot(sum, 
       aes(x = Cond, y = Mean, fill = Cond, 
           ymax=Mean+se, ymin=Mean-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Cue type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(0, -1.2))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))


# P300
MID_cueU <- subset(MID_cueU, select=c("File", "P1.Average_MID_cueU_NEW", 'P2.Average_MID_cueU_NEW', "PO3.Average_MID_cueU_NEW", 'PO4.Average_MID_cueU_NEW', "Pz.Average_MID_cueU_NEW", 'POz.Average_MID_cueU_NEW'))
MID_cueR <- subset(MID_cueR, select=c("File", "P1.Average_MID_cueR_NEW", 'P2.Average_MID_cueR_NEW', "PO3.Average_MID_cueR_NEW", 'PO4.Average_MID_cueR_NEW', "Pz.Average_MID_cueR_NEW", 'POz.Average_MID_cueR_NEW'))
CSRA_cue <- subset(CSRA_cue, select=c("File", "P1.Average_CSRA_cue_NEW", 'P2.Average_CSRA_cue_NEW', "PO3.Average_CSRA_cue_NEW", 'PO4.Average_CSRA_cue_NEW', "Pz.Average_CSRA_cue_NEW", 'POz.Average_CSRA_cue_NEW'))

MID_cueU$Mean <- rowMeans(MID_cueU[,2:7])
MID_cueR$Mean <- rowMeans(MID_cueR[,2:7])
CSRA_cue$Mean <- rowMeans(CSRA_cue[,2:7])

MID_cueU$Cond <- "MID_cueU"
MID_cueR$Cond <- "MID_cueR"
CSRA_cue$Cond <- "CSRA_cue"


MID_cueU$Cond <- "MID non-rew"
MID_cueR$Cond <- "MID rew"
CSRA_cue$Cond <- "C-SRA"

MID_cueU <- MID_cueU[,-(2:7),drop=FALSE]
MID_cueR <- MID_cueR[,-(2:7),drop=FALSE]
CSRA_cue <- CSRA_cue[,-(2:7),drop=FALSE]


df <- rbind(CSRA_cue, MID_cueR, MID_cueU)

df$File <- as.factor(df$File)
df$Cond <- as.factor(df$Cond)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Cond), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$Mean, df$Cond, p.adjust.method="none", paired=TRUE)

mean(CSRA_cue$Mean)
sd(CSRA_cue$Mean)
mean(MID_cueR$Mean)
sd(MID_cueR$Mean)
mean(MID_cueU$Mean)
sd(MID_cueU$Mean)

# barplot
sum = summarySEwithin(df, measurevar="Mean", withinvars=c("Cond"), idvar="File", na.rm=FALSE, conf.interval=.95)

ggplot(sum, 
       aes(x = Cond, y = Mean, fill = Cond, 
           ymax=Mean+se, ymin=Mean-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Cue type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(1.5, 3))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))

#--------------------- TARGET _ REWARD --------------------------------
# P300
CSRA_R <- subset(Rew, select=c("File", "Fz.Average_CSRA_R", "FCz.Average_CSRA_R", "FC1.Average_CSRA_R", "FC2.Average_CSRA_R", "F1.Average_CSRA_R", "F2.Average_CSRA_R"))
SRA_R <- subset(Rew, select=c("File", "Fz.Average_SRA_R", "FCz.Average_SRA_R", "FC1.Average_SRA_R", "FC2.Average_SRA_R", "F1.Average_SRA_R", "F2.Average_SRA_R"))
MID_R <- subset(Rew, select=c("File", "Fz.Average_MID_R", "FCz.Average_MID_R",  "FC1.Average_MID_R",  "FC2.Average_MID_R", "F1.Average_MID_R", "F2.Average_MID_R"))

CSRA_U <- subset(Norew, select=c("File", "Fz.Average_CSRA_U", "FCz.Average_CSRA_U", "FC1.Average_CSRA_U", "FC2.Average_CSRA_U", "F1.Average_CSRA_U", "F2.Average_CSRA_U"))
SRA_U<- subset(Norew, select=c("File", "Fz.Average_SRA_U", "FCz.Average_SRA_U", "FC1.Average_SRA_U", "FC2.Average_SRA_U", "F1.Average_SRA_U", "F2.Average_SRA_U"))
MID_U <- subset(Norew, select=c("File", "Fz.Average_MID_U", "FCz.Average_MID_U",  "FC1.Average_MID_U",  "FC2.Average_MID_U", "F1.Average_MID_U", "F2.Average_MID_U"))

CSRA_R$Mean <- rowMeans(CSRA_R[,2:7])
SRA_R$Mean <- rowMeans(SRA_R[,2:7])
MID_R$Mean <- rowMeans(MID_R[,2:7])

CSRA_U$Mean <- rowMeans(CSRA_U[,2:7])
SRA_U$Mean <- rowMeans(SRA_U[,2:7])
MID_U$Mean <- rowMeans(MID_U[,2:7])


CSRA_R$Block <- "CSRA"
MID_R$Block <- "MID"
SRA_R$Block <- "SRA"
CSRA_R$Reward <- "Rew"
MID_R$Reward <- "Rew"
SRA_R$Reward <- "Rew"

CSRA_U$Block <- "CSRA"
MID_U$Block <- "MID"
SRA_U$Block <- "SRA"
CSRA_U$Reward <- "NoRew"
MID_U$Reward <- "NoRew"
SRA_U$Reward <- "NoRew"

CSRA_R <- CSRA_R[,-(2:7),drop=FALSE]
SRA_R <- SRA_R[,-(2:7),drop=FALSE]
MID_R <- MID_R[,-(2:7),drop=FALSE]

CSRA_U <- CSRA_U[,-(2:7),drop=FALSE]
SRA_U <- SRA_U[,-(2:7),drop=FALSE]
MID_U <- MID_U[,-(2:7),drop=FALSE]



df <- rbind(CSRA_R, MID_R, SRA_R, CSRA_U, MID_U, SRA_U)

# aov <- with(df, aov(Mean ~ Cond * Block + Error(File / (Cond * Block)))) # two way ANOVA repeated measures
# summary(aov)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Reward, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$Mean, df$Block, p.adjust.method="none", paired=TRUE)

mean(df$Mean[df$Reward == "Rew"])
sd(df$Mean[df$Reward == "Rew"])
mean(df$Mean[df$Reward == "NoRew"])
sd(df$Mean[df$Reward == "NoRew"])

mean(df$Mean[df$Block == "CSRA"])
sd(df$Mean[df$Block == "CSRA"])
mean(df$Mean[df$Block == "SRA"])
sd(df$Mean[df$Block == "SRA"])

# post hoc interactions
dataMID <- subset(df, Cond == "MID")
dataSRA <- subset(df, Cond == "SRA")
dataCSRA <- subset(df, Cond == "CSRA")

dataMID$Cond <- as.factor(dataMID$Cond)
dataMID$Block <- as.factor(dataMID$Block)
dataMID$File <- as.factor(dataMID$File)
dataMID$Mean <- as.numeric(as.character(dataMID$Mean))
dataSRA$Cond <- as.factor(dataSRA$Cond)
dataSRA$Block <- as.factor(dataSRA$Block)
dataSRA$File <- as.factor(dataSRA$File)
dataSRA$Mean <- as.numeric(as.character(dataSRA$Mean))
dataCSRA$Cond <- as.factor(dataCSRA$Cond)
dataCSRA$Block <- as.factor(dataCSRA$Block)
dataCSRA$File <- as.factor(dataCSRA$File)
dataCSRA$Mean <- as.numeric(as.character(dataCSRA$Mean))

rew <- subset(dataSRA, Block =="Rew")
norew <- subset(dataSRA, Block=="NoRew")
t.test(rew$Mean, norew$Mean, paired = TRUE) # variables are numeric 

rew <- subset(dataCSRA, Block =="Rew")
norew <- subset(dataCSRA, Block=="NoRew")
t.test(rew$Mean, norew$Mean, paired = TRUE) # variables are numeric 

rew <- subset(dataMID, Block =="Rew")
norew <- subset(dataMID, Block=="NoRew")
t.test(rew$Mean, norew$Mean, paired = TRUE) # variables are numeric 

means <- aggregate(Mean ~  Cond*Block, df, mean)
means$Mean  <- round(means$Mean, digits = 2)

ggplot(data=df, aes(x=Cond, y=Mean, fill=Block)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE, position=position_dodge(width=0.75)) + 
  geom_text(data = means, aes(label = Mean, y = Mean + 0.3), position=position_dodge(width=0.8))+
  xlab("Conditions") +
  ylab("Amplitude [µV]") 
#  ggtitle("Reward effect\nP300 component (280-400 ms)\nPO3, POz, PO4, Pz, P2, P4 electrodes") # always change if needed


# barplot
sum = summarySEwithin(df, measurevar="Mean", withinvars=c("Block","Reward"), idvar="File", na.rm=FALSE, conf.interval=.95)

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Reward <- revalue(sum$Reward , c("Rew"="Reward"))
sum$Reward <- revalue(sum$Reward , c("NoRew"="No Reward"))

ggplot(sum, 
       aes(x = Reward, y = Mean, fill = Block, 
           ymax=Mean+se, ymin=Mean-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Block type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(-1.5, 0.5))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))

# N1 component
CSRA_R <- subset(Rew, select=c("File", "PO7.Average_CSRA_R", "PO8.Average_CSRA_R"))
SRA_R <- subset(Rew, select=c("File", "PO7.Average_SRA_R", "PO8.Average_SRA_R"))
MID_R <- subset(Rew, select=c("File", "PO7.Average_MID_R", "PO8.Average_MID_R"))

CSRA_U <- subset(Norew, select=c("File", "PO7.Average_CSRA_U", "PO8.Average_CSRA_U"))
SRA_U<- subset(Norew, select=c("File", "PO7.Average_SRA_U", "PO8.Average_SRA_U"))
MID_U <- subset(Norew, select=c("File", "PO7.Average_MID_U", "PO8.Average_MID_U"))

CSRA_R$Mean <- rowMeans(CSRA_R[,2:3])
SRA_R$Mean <- rowMeans(SRA_R[,2:3])
MID_R$Mean <- rowMeans(MID_R[,2:3])

CSRA_U$Mean <- rowMeans(CSRA_U[,2:3])
SRA_U$Mean <- rowMeans(SRA_U[,2:3])
MID_U$Mean <- rowMeans(MID_U[,2:3])


CSRA_R$Block <- "CSRA"
MID_R$Block <- "MID"
SRA_R$Block <- "SRA"
CSRA_R$Reward <- "Rew"
MID_R$Reward <- "Rew"
SRA_R$Reward <- "Rew"

CSRA_U$Block <- "CSRA"
MID_U$Block <- "MID"
SRA_U$Block <- "SRA"
CSRA_U$Reward <- "NoRew"
MID_U$Reward <- "NoRew"
SRA_U$Reward <- "NoRew"

CSRA_R <- CSRA_R[,-(2:3),drop=FALSE]
SRA_R <- SRA_R[,-(2:3),drop=FALSE]
MID_R <- MID_R[,-(2:3),drop=FALSE]

CSRA_U <- CSRA_U[,-(2:3),drop=FALSE]
SRA_U <- SRA_U[,-(2:3),drop=FALSE]
MID_U <- MID_U[,-(2:3),drop=FALSE]



df <- rbind(CSRA_R, MID_R, SRA_R, CSRA_U, MID_U, SRA_U)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Reward, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

# barplot
sum = summarySEwithin(df, measurevar="Mean", withinvars=c("Block","Reward"), idvar="File", na.rm=FALSE, conf.interval=.95)

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Reward <- revalue(sum$Reward , c("Rew"="Reward"))
sum$Reward <- revalue(sum$Reward , c("NoRew"="No Reward"))

ggplot(sum, 
       aes(x = Reward, y = Mean, fill = Block, 
           ymax=Mean+se, ymin=Mean-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Block type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(-2, 0))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))

# Positivity 450-800 ms 
CSRA_R <- subset(Rew, select=c("File", "FCz.Average_CSRA_R", "Cz.Average_CSRA_R"))
SRA_R <- subset(Rew, select=c("File", "FCz.Average_SRA_R", "Cz.Average_SRA_R"))
MID_R <- subset(Rew, select=c("File", "FCz.Average_MID_R", "Cz.Average_MID_R"))

CSRA_U <- subset(Norew, select=c("File", "FCz.Average_CSRA_U", "Cz.Average_CSRA_U"))
SRA_U<- subset(Norew, select=c("File", "FCz.Average_SRA_U", "Cz.Average_SRA_U"))
MID_U <- subset(Norew, select=c("File", "FCz.Average_MID_U", "Cz.Average_MID_U"))

CSRA_R$Mean <- rowMeans(CSRA_R[,2:3])
SRA_R$Mean <- rowMeans(SRA_R[,2:3])
MID_R$Mean <- rowMeans(MID_R[,2:3])

CSRA_U$Mean <- rowMeans(CSRA_U[,2:3])
SRA_U$Mean <- rowMeans(SRA_U[,2:3])
MID_U$Mean <- rowMeans(MID_U[,2:3])


CSRA_R$Block <- "CSRA"
MID_R$Block <- "MID"
SRA_R$Block <- "SRA"
CSRA_R$Reward <- "Rew"
MID_R$Reward <- "Rew"
SRA_R$Reward <- "Rew"

CSRA_U$Block <- "CSRA"
MID_U$Block <- "MID"
SRA_U$Block <- "SRA"
CSRA_U$Reward <- "NoRew"
MID_U$Reward <- "NoRew"
SRA_U$Reward <- "NoRew"

CSRA_R <- CSRA_R[,-(2:3),drop=FALSE]
SRA_R <- SRA_R[,-(2:3),drop=FALSE]
MID_R <- MID_R[,-(2:3),drop=FALSE]

CSRA_U <- CSRA_U[,-(2:3),drop=FALSE]
SRA_U <- SRA_U[,-(2:3),drop=FALSE]
MID_U <- MID_U[,-(2:3),drop=FALSE]



df <- rbind(CSRA_R, MID_R, SRA_R, CSRA_U, MID_U, SRA_U)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Reward, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

# barplot
sum = summarySEwithin(df, measurevar="Mean", withinvars=c("Block","Reward"), idvar="File", na.rm=FALSE, conf.interval=.95)

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Reward <- revalue(sum$Reward , c("Rew"="Reward"))
sum$Reward <- revalue(sum$Reward , c("NoRew"="No Reward"))

ggplot(sum, 
       aes(x = Reward, y = Mean, fill = Block, 
           ymax=Mean+se, ymin=Mean-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Block type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(-0.5, 1.5))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))

#--------------------- TARGET _ CONGRUENCY --------------------------------
#N450
CSRA_con <- subset(Con, select=c("File", "Cz.Average_CSRA_con",  "CP1.Average_CSRA_con", "CP2.Average_CSRA_con", "CPz.Average_CSRA_con"))
SRA_con <- subset(Con, select=c("File", "Cz.Average_SRA_con",  "CP1.Average_SRA_con", "CP2.Average_SRA_con", "CPz.Average_SRA_con"))
MID_con <- subset(Con, select=c("File", "Cz.Average_MID_con", "CP1.Average_MID_con", "CP2.Average_MID_con", "CPz.Average_MID_con"))
NEUT_con <- subset(Con, select=c("File", "Cz.Average_NEUT_con", "CP1.Average_NEUT_con", "CP2.Average_NEUT_con", "CPz.Average_NEUT_con"))

CSRA_inc <- subset(Inc, select=c("File", "Cz.Average_CSRA_inc", "CP1.Average_CSRA_inc", "CP2.Average_CSRA_inc", "CPz.Average_CSRA_inc"))
SRA_inc <- subset(Inc, select=c("File", "Cz.Average_SRA_inc", "CP1.Average_SRA_inc", "CP2.Average_SRA_inc", "CPz.Average_SRA_inc"))
MID_inc <- subset(Inc, select=c("File", "Cz.Average_MID_inc",  "CP1.Average_MID_inc", "CP2.Average_MID_inc", "CPz.Average_MID_inc"))
NEUT_inc <- subset(Inc, select=c("File", "Cz.Average_Neut_inc_correct", "CP1.Average_Neut_inc_correct", "CP2.Average_Neut_inc_correct", "CPz.Average_Neut_inc_correct"))

CSRA_con$Mean <- rowMeans(CSRA_con[,2:5])
SRA_con$Mean <- rowMeans(SRA_con[,2:5])
MID_con$Mean <- rowMeans(MID_con[,2:5])
NEUT_con$Mean <- rowMeans(NEUT_con[,2:5])

CSRA_inc$Mean <- rowMeans(CSRA_inc[,2:5])
SRA_inc$Mean <- rowMeans(SRA_inc[,2:5])
MID_inc$Mean <- rowMeans(MID_inc[,2:5])
NEUT_inc$Mean <- rowMeans(NEUT_inc[,2:5])

CSRA_con$Block <- "CSRA"
SRA_con$Block <- "SRA"
MID_con$Block <- "MID"
NEUT_con$Block <- "NEUT"
CSRA_con$Congruency <- "Con"
SRA_con$Congruency <- "Con"
MID_con$Congruency <- "Con"
NEUT_con$Congruency <- "Con"

CSRA_inc$Block <- "CSRA"
SRA_inc$Block <- "SRA"
MID_inc$Block <- "MID"
NEUT_inc$Block <- "NEUT"
CSRA_inc$Congruency <- "Inc"
SRA_inc$Congruency <- "Inc"
MID_inc$Congruency <- "Inc"
NEUT_inc$Congruency <- "Inc"

CSRA_con <- CSRA_con[,-(2:5),drop=FALSE] 
SRA_con <- SRA_con[,-(2:5),drop=FALSE] 
MID_con <- MID_con[,-(2:5),drop=FALSE] 
NEUT_con <- NEUT_con[,-(2:5),drop=FALSE] 

CSRA_inc <- CSRA_inc[,-(2:5),drop=FALSE] 
SRA_inc <- SRA_inc[,-(2:5),drop=FALSE] 
MID_inc <- MID_inc[,-(2:5),drop=FALSE] 
NEUT_inc <- NEUT_inc[,-(2:5),drop=FALSE] 

df <- rbind(CSRA_con, SRA_con, MID_con, NEUT_con, CSRA_inc, SRA_inc, MID_inc, NEUT_inc)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Congruency, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)


means <- aggregate(Mean ~  Cond*Block, df, mean)
means$Mean  <- round(means$Mean, digits = 2)

ggplot(data=df, aes(x=Cond, y=Mean, fill=Block)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show_guide = FALSE, position=position_dodge(width=0.75)) + 
  geom_text(data = means, aes(label = Mean, y = Mean + 0.3), position=position_dodge(width=0.8))+
  xlab("Conditions") +
  ylab("Amplitude [µV]") 
# ggtitle("Congruency effect\nN450 component (350-500 ms)\nCentro-parietal electrodes\n(CP5, C3, C1, CP3, CP1, CPz, CP2, Cz, C2, Fz, FC1, FC2, CP4, C4)")

pairwise.t.test(df$Mean, df$Block, p.adjust.method="none", paired=TRUE)

mean(df$Mean[df$Congruency == "Con"])
sd(df$Mean[df$Congruency == "Con"])
mean(df$Mean[df$Congruency == "Inc"])
sd(df$Mean[df$Congruency == "Inc"])

mean(df$Mean[df$Cond == "CSRA"])
sd(df$Mean[df$Cond == "CSRA"])
mean(df$Mean[df$Cond == "SRA"])
sd(df$Mean[df$Cond == "SRA"])

# post hoc interactions
dataMID <- subset(df, Cond == "MID")
dataSRA <- subset(df, Cond == "SRA")
dataCSRA <- subset(df, Cond == "CSRA")
dataNeut <- subset(df, Cond == "NEUT")

dataMID$Cond <- as.factor(dataMID$Cond)
dataMID$Block <- as.factor(dataMID$Block)
dataMID$File <- as.factor(dataMID$File)
dataMID$Mean <- as.numeric(as.character(dataMID$Mean))
dataSRA$Cond <- as.factor(dataSRA$Cond)
dataSRA$Block <- as.factor(dataSRA$Block)
dataSRA$File <- as.factor(dataSRA$File)
dataSRA$Mean <- as.numeric(as.character(dataSRA$Mean))
dataCSRA$Cond <- as.factor(dataCSRA$Cond)
dataCSRA$Block <- as.factor(dataCSRA$Block)
dataCSRA$File <- as.factor(dataCSRA$File)
dataCSRA$Mean <- as.numeric(as.character(dataCSRA$Mean))
dataNeut$Cond <- as.factor(dataNeut$Cond)
dataNeut$Block <- as.factor(dataNeut$Block)
dataNeut$File <- as.factor(dataNeut$File)
dataNeut$Mean <- as.numeric(as.character(dataNeut$Mean))


con <- subset(dataSRA, Block =="Con")
inc <- subset(dataSRA, Block=="Inc")
t.test(con$Mean, inc$Mean, paired = TRUE) # variables are numeric 

con <- subset(dataCSRA, Block =="Con")
inc <- subset(dataCSRA, Block=="Inc")
t.test(con$Mean, inc$Mean, paired = TRUE) # variables are numeric 

con <- subset(dataMID, Block =="Con")
inc <- subset(dataMID, Block=="Inc")
t.test(con$Mean, inc$Mean, paired = TRUE) # variables are numeric 

con <- subset(dataNeut, Block =="Con")
inc <- subset(dataNeut, Block=="Inc")
t.test(con$Mean, inc$Mean, paired = TRUE) # variables are numeric 

mean(dataMID$Mean[dataMID$Block == "Inc"])
sd(dataMID$Mean[dataMID$Block == "Inc"])
mean(dataMID$Mean[dataMID$Block == "Con"])
sd(dataMID$Mean[dataMID$Block == "Con"])

mean(dataNeut$Mean[dataNeut$Block == "Inc"])
sd(dataNeut$Mean[dataNeut$Block == "Inc"])
mean(dataNeut$Mean[dataNeut$Block == "Con"])
sd(dataNeut$Mean[dataNeut$Block == "Con"])


# barplot
sum = summarySEwithin(df, measurevar="Mean", withinvars=c("Block","Congruency"), idvar="File", na.rm=FALSE, conf.interval=.95)

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("Con"="Congruent"))
sum$Congruency <- revalue(sum$Congruency , c("Inc"="Incongruent"))

ggplot(sum, 
       aes(x = Congruency, y = Mean, fill = Block, 
           ymax=Mean+se, ymin=Mean-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Block type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9","#F0E442", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(0.5,2))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))


# LPC
CSRA_con <- subset(Con, select=c("File", "P3.Average_CSRA_con", "P1.Average_CSRA_con",  "Pz.Average_CSRA_con", "P2.Average_CSRA_con", "P4.Average_CSRA_con"))
SRA_con <- subset(Con, select=c("File", "P3.Average_SRA_con", "P1.Average_SRA_con","Pz.Average_SRA_con", "P2.Average_SRA_con", "P4.Average_SRA_con"))
MID_con <- subset(Con, select=c("File", "P3.Average_MID_con", "P1.Average_MID_con",  "Pz.Average_MID_con", "P2.Average_MID_con", "P4.Average_MID_con"))
NEUT_con <- subset(Con, select=c("File", "P3.Average_NEUT_con", "P1.Average_NEUT_con", "Pz.Average_NEUT_con", "P2.Average_NEUT_con", "P4.Average_NEUT_con"))

CSRA_inc <- subset(Inc, select=c("File", "P3.Average_CSRA_inc", "P1.Average_CSRA_inc", "Pz.Average_CSRA_inc", "P2.Average_CSRA_inc", "P4.Average_CSRA_inc"))
SRA_inc <- subset(Inc, select=c("File", "P3.Average_SRA_inc", "P1.Average_SRA_inc", "Pz.Average_SRA_inc", "P2.Average_SRA_inc", "P4.Average_SRA_inc"))
MID_inc <- subset(Inc, select=c("File", "P3.Average_MID_inc", "P1.Average_MID_inc", "Pz.Average_MID_inc", "P2.Average_MID_inc", "P4.Average_MID_inc"))
NEUT_inc <- subset(Inc, select=c("File", "P3.Average_Neut_inc_correct", "P1.Average_Neut_inc_correct", "Pz.Average_Neut_inc_correct", "P2.Average_Neut_inc_correct", "P4.Average_Neut_inc_correct"))


CSRA_con$Mean <- rowMeans(CSRA_con[,2:6])
SRA_con$Mean <- rowMeans(SRA_con[,2:6])
MID_con$Mean <- rowMeans(MID_con[,2:6])
NEUT_con$Mean <- rowMeans(NEUT_con[,2:6])

CSRA_inc$Mean <- rowMeans(CSRA_inc[,2:6])
SRA_inc$Mean <- rowMeans(SRA_inc[,2:6])
MID_inc$Mean <- rowMeans(MID_inc[,2:6])
NEUT_inc$Mean <- rowMeans(NEUT_inc[,2:6])

CSRA_con$Block <- "CSRA"
SRA_con$Block <- "SRA"
MID_con$Block <- "MID"
NEUT_con$Block <- "NEUT"
CSRA_con$Congruency <- "Con"
SRA_con$Congruency <- "Con"
MID_con$Congruency <- "Con"
NEUT_con$Congruency <- "Con"

CSRA_inc$Block <- "CSRA"
SRA_inc$Block <- "SRA"
MID_inc$Block <- "MID"
NEUT_inc$Block <- "NEUT"
CSRA_inc$Congruency <- "Inc"
SRA_inc$Congruency <- "Inc"
MID_inc$Congruency <- "Inc"
NEUT_inc$Congruency <- "Inc"

CSRA_con <- CSRA_con[,-(2:6),drop=FALSE] 
SRA_con <- SRA_con[,-(2:6),drop=FALSE] 
MID_con <- MID_con[,-(2:6),drop=FALSE] 
NEUT_con <- NEUT_con[,-(2:6),drop=FALSE] 

CSRA_inc <- CSRA_inc[,-(2:6),drop=FALSE] 
SRA_inc <- SRA_inc[,-(2:6),drop=FALSE] 
MID_inc <- MID_inc[,-(2:6),drop=FALSE] 
NEUT_inc <- NEUT_inc[,-(2:6),drop=FALSE] 

df <- rbind(CSRA_con, SRA_con, MID_con, NEUT_con, CSRA_inc, SRA_inc, MID_inc, NEUT_inc)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Congruency, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$Mean, df$Block, p.adjust.method="none", paired=TRUE)

# post hoc interactions
df_CON <- subset(df, Congruency=="Con")
df_INC <- subset(df, Congruency=="Inc")
df_CON$con_effect <- df_CON$Mean-df_INC$Mean

dataSRA <- subset(df_CON, Block == "SRA")
dataCSRA <- subset(df_CON, Block == "CSRA")
dataMID <- subset(df_CON, Block == "MID")
dataNEUT <- subset(df_CON, Block == "NEUT")
t.test(dataMID$con_effect, dataSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataMID$con_effect, dataCSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataMID$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$con_effect, dataCSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric
t.test(dataCSRA$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric



mean(df$Mean[df$Congruency == "Con"])
sd(df$Mean[df$Congruency == "Con"])
mean(df$Mean[df$Congruency == "Inc"])
sd(df$Mean[df$Congruency == "Inc"])

mean(df$Mean[df$Cond == "CSRA"])
sd(df$Mean[df$Cond == "CSRA"])
mean(df$Mean[df$Cond == "SRA"])
sd(df$Mean[df$Cond == "SRA"])


mean(dataMID$Mean[dataMID$Block == "Inc"])
sd(dataMID$Mean[dataMID$Block == "Inc"])
mean(dataMID$Mean[dataMID$Block == "Con"])
sd(dataMID$Mean[dataMID$Block == "Con"])

mean(dataNeut$Mean[dataNeut$Block == "Inc"])
sd(dataNeut$Mean[dataNeut$Block == "Inc"])
mean(dataNeut$Mean[dataNeut$Block == "Con"])
sd(dataNeut$Mean[dataNeut$Block == "Con"])


# barplot
sum = summarySEwithin(df, measurevar="Mean", withinvars=c("Block","Congruency"), idvar="File", na.rm=FALSE, conf.interval=.95)

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("Con"="Congruent"))
sum$Congruency <- revalue(sum$Congruency , c("Inc"="Incongruent"))

ggplot(sum, 
       aes(x = Congruency, y = Mean, fill = Block, 
           ymax=Mean+se, ymin=Mean-se))  +
  geom_bar(stat="identity", position = "dodge", width = 0.7) +
  geom_bar(stat="identity", position = "dodge", 
           colour = "black", width = 0.7, 
           show_guide = FALSE)  +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "",
       y = "Amplitude [µV]")  +
  scale_fill_manual(name="Block type", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9","#F0E442", "#009E73")) +
  #ggtitle("The Effect of Reward and Congruency on\nResponses in MID block") + 
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_rect(fill = "black"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5))+
  coord_cartesian(ylim = c(2,4))+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16)) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))


# 550-750
CSRA_con <- subset(Con, select=c("File", "FCz.Average_CSRA_con", "Cz.Average_CSRA_con",  "Pz.Average_CSRA_con"))
SRA_con <- subset(Con, select=c("File", "FCz.Average_SRA_con", "Cz.Average_SRA_con","Pz.Average_SRA_con"))
MID_con <- subset(Con, select=c("File", "FCz.Average_MID_con", "Cz.Average_MID_con",  "Pz.Average_MID_con"))
NEUT_con <- subset(Con, select=c("File", "FCz.Average_NEUT_con", "Cz.Average_NEUT_con", "Pz.Average_NEUT_con"))

CSRA_inc <- subset(Inc, select=c("File", "FCz.Average_CSRA_inc", "Cz.Average_CSRA_inc", "Pz.Average_CSRA_inc"))
SRA_inc <- subset(Inc, select=c("File", "FCz.Average_SRA_inc", "Cz.Average_SRA_inc", "Pz.Average_SRA_inc"))
MID_inc <- subset(Inc, select=c("File", "FCz.Average_MID_inc", "Cz.Average_MID_inc", "Pz.Average_MID_inc"))
NEUT_inc <- subset(Inc, select=c("File", "FCz.Average_Neut_inc_correct", "Cz.Average_Neut_inc_correct", "Pz.Average_Neut_inc_correct"))


CSRA_con$Mean <- rowMeans(CSRA_con[,2:4])
SRA_con$Mean <- rowMeans(SRA_con[,2:4])
MID_con$Mean <- rowMeans(MID_con[,2:4])
NEUT_con$Mean <- rowMeans(NEUT_con[,2:4])

CSRA_inc$Mean <- rowMeans(CSRA_inc[,2:4])
SRA_inc$Mean <- rowMeans(SRA_inc[,2:4])
MID_inc$Mean <- rowMeans(MID_inc[,2:4])
NEUT_inc$Mean <- rowMeans(NEUT_inc[,2:4])

CSRA_con$Block <- "CSRA"
SRA_con$Block <- "SRA"
MID_con$Block <- "MID"
NEUT_con$Block <- "NEUT"
CSRA_con$Congruency <- "Con"
SRA_con$Congruency <- "Con"
MID_con$Congruency <- "Con"
NEUT_con$Congruency <- "Con"

CSRA_inc$Block <- "CSRA"
SRA_inc$Block <- "SRA"
MID_inc$Block <- "MID"
NEUT_inc$Block <- "NEUT"
CSRA_inc$Congruency <- "Inc"
SRA_inc$Congruency <- "Inc"
MID_inc$Congruency <- "Inc"
NEUT_inc$Congruency <- "Inc"

CSRA_con <- CSRA_con[,-(2:4),drop=FALSE] 
SRA_con <- SRA_con[,-(2:4),drop=FALSE] 
MID_con <- MID_con[,-(2:4),drop=FALSE] 
NEUT_con <- NEUT_con[,-(2:4),drop=FALSE] 

CSRA_inc <- CSRA_inc[,-(2:4),drop=FALSE] 
SRA_inc <- SRA_inc[,-(2:4),drop=FALSE] 
MID_inc <- MID_inc[,-(2:4),drop=FALSE] 
NEUT_inc <- NEUT_inc[,-(2:4),drop=FALSE] 

df <- rbind(CSRA_con, SRA_con, MID_con, NEUT_con, CSRA_inc, SRA_inc, MID_inc, NEUT_inc)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Mean, # specify dependent variable 
                     wid = File, # specify the subject variable
                     within = .(Congruency, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)

demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$Mean, df$Block, p.adjust.method="none", paired=TRUE)

mean(df$Mean[df$Block == "CSRA"])
sd(df$Mean[df$Cond == "CSRA"])
mean(df$Mean[df$Block == "MID"])
sd(df$Mean[df$Cond == "SRA"])
mean(df$Mean[df$Block == "NEUT"])
