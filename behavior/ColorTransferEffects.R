# Author: MK
###############################################################################
###############################################################################
###############################################################################
################################################################
rm(list=ls()) # clear the working directory

#---------------setting up the path---------------
myPath <- "C:/Users/mariam/Documents/PhD/EXPERIMENTS/MK_MREEGRewModes_2016/Analysis/BehavData/Data_1000cutoff_ColTransfEff"
setwd(myPath)

#--------------import all libraries---------------
library(dplyr)
library(plyr)
library(useful) # for shifting the columns 
library(stringr)
library(data.table)
library(grid)
library(ggplot2)
library(Rmisc)


#----------------- Transition of rewarded colors in CRSA and SRA to MID and Neut blocks -----------------------

########## MID ##############
#--------------reading in text files------------------------------

df <- choose.files(default = "", caption = "Select files",multi = TRUE, filters = Filters, index = nrow(Filters))


df <- ldply(df,    ## where are the files I should read?
            read.table,  ## for every entry in fileList, read a table
            header = T, ## header = T means the variable names are in the txt files
            fill = TRUE,
            sep="\t")     ## sep = "t" means different values in the text-file are seperated by "tab"


df <-df[!df$Subject %in% c(30, 31),] # exclude participants that are bad in the EEG (#31) and the fMRI (#30)

#------------- subset all conditions ----------------------
CSRA.df <- subset(df, Description %in% c('CSRA_Rinc', "CSRA_Rcon", 'CSRA_Ucon', "CSRA_Uinc"))
SRA.df <- subset(df, Description %in% c('SRA_Rinc', "SRA_Rcon", 'SRA_Ucon', "SRA_Uinc"))
MID.df <- subset(df, Description %in% c("MID_redR", "MID_yellowY", "MID_blueB", "MID_greenG", 
                                        "MID_yellowR", "MID_greenR", "MID_blueY", "MID_redY", "MID_yellowB", "MID_greenB", "MID_redG", "MID_blueG",
                                        "MID_redR", "MID_yellowY", "MID_blueB", "MID_greenG",
                                        "MID_yellowR", "MID_greenR", "MID_redY", "MID_blueY","MID_greenB", "MID_yellowB", "MID_redG", "MID_blueG"))
NEUT.df <- subset(df, Description %in% c("NEUT_redR", "NEUT_yellowY", "NEUT_blueB", "NEUT_greenG",
                                         "NEUT_greenR", "NEUT_yellowR","NEUT_redY", "NEUT_blueY", "NEUT_yellowB", "NEUT_greenB", "NEUT_blueG", "NEUT_redG"))   



#--------------dividing by rew col------------------------------

MID_RY <- subset(MID.df, Subject %in% c(9,11,12,13,14,21,22,23,24,25,26,35,36))

MID_RY$RewCol <- "RY"

MID_RY$Description <- revalue(MID_RY$Description, c("MID_redR"="RewCol", "MID_blueY"="RewCol", "MID_yellowY"="RewCol", "MID_yellowR"="RewCol", "MID_greenR"="RewCol", "MID_redY"="RewCol",
                                      "MID_redG"="NonRewCol", "MID_greenG"="NonRewCol", "MID_greenB"="NonRewCol", "MID_blueB"="NonRewCol", "MID_yellowB"="NonRewCol", "MID_blueG"="NonRewCol"))




MID_BG <- subset(MID.df, Subject %in% c(15,16,17,18,19,20,27,28,29,32,33,34))

MID_BG$RewCol <- "BG"

MID_BG$Description <- revalue(MID_BG$Description, c("MID_blueB"="RewCol", "MID_blueG"="RewCol", "MID_yellowB"="RewCol", "MID_greenG"="RewCol", "MID_greenB"="RewCol", "MID_redG"="RewCol",
                                      "MID_yellowY"="NonRewCol", "MID_redR"="NonRewCol", "MID_greenR"="NonRewCol", "MID_redY"="NonRewCol", "MID_blueY"="NonRewCol", "MID_yellowR"="NonRewCol"))

MID <- rbind(MID_BG, MID_RY)

#---------------prepare for the stats-------
MID.rt <-MID[!MID$Accuracy %in% c("miss"),] # exclude misses 
MID.rt <-MID.rt[!MID.rt$Accuracy2 %in% c("incorrect"),]  # exclude incorrect 

MID.rt <- subset(MID.rt, TTime.Shifted <=1000)

av.rt <- aggregate(data=MID.rt, TTime.Shifted~Description*Subject, mean) 

#--------------boxplot ---------

means <- aggregate(TTime.Shifted ~ Description, av.rt, mean)
means$TTime.Shifted <- round(means$TTime.Shifted,digits=2 )

ggplot(data = av.rt, aes(x=Description, y=TTime.Shifted, fill=Description)) +
  geom_boxplot() +
  xlab("Description") + ylab("RTs [ms]")+
  #  stat_summary(fun.y=mean,col='black',geom='point', shape=18, size=3,show_guide = FALSE) +
  #  geom_text(data = means, aes(label = means$TTTime.Shifted, y = TTTime.Shifted + 0.08), vjust = -0.7) +
  scale_x_discrete(breaks=c("con", "inc" ),
                   labels=c("Congruent", "Incongruent")) +
  scale_fill_manual(values=c("#E53304", "#164A7C"))+
  ggtitle("The Effect of Rewarded Color on\nResponses in the MID block") +
  theme(plot.title = element_text(face="bold", size=10)) +
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(text = element_text(size=20,face="bold"), axis.title=element_text(size=25,face="bold"))+
  coord_cartesian(ylim = c(400,800)) 

#------------- t test----------------
rew <- subset(av.rt, Description=="RewCol")
norew <- subset(av.rt, Description=="NonRewCol")

t.test(rew$TTime.Shifted, norew$TTime.Shifted, paired = TRUE) # variables are numeric 

mean(rew$TTime.Shifted)
sd(rew$TTime.Shifted)
mean(norew$TTime.Shifted)
sd(norew$TTime.Shifted)

#---------- Accuracies -------------------
MID$Accuracy[MID$Accuracy2 == "incorrect"] <- "incorrect"
MID$Accuracy[is.na(MID$Accuracy)] <- 1
MID$Accuracy <- revalue(MID$Accuracy, c("incorrect"=0, "miss"=0))
MID$Description <- revalue(MID$Description, c("RewCol"=1, "NonRewCol"=0))

MID$Accuracy <- as.numeric(as.character(MID$Accuracy))
MID$Description <- as.numeric(as.character(MID$Description))

MID.Acc <- as.data.frame(xtabs(Accuracy~Description+Subject,MID)/xtabs(~Description+Subject,MID))
MID.Acc$Freq <- MID.Acc$Freq * 100

#--------------boxplot ---------

means <- aggregate(Freq ~ Description, MID.Acc, mean)
means$Freq <- round(means$Freq,digits=2 )

ggplot(data = means, aes(x=Description, y=Freq, fill=Description)) +
  geom_boxplot() +
  xlab("Description") + ylab("Accuracy [%]")+
  #  stat_summary(fun.y=mean,col='black',geom='point', shape=18, size=3,show_guide = FALSE) +
  #  geom_text(data = means, aes(label = means$TTTime.Shifted, y = TTTime.Shifted + 0.08), vjust = -0.7) +
  scale_x_discrete(breaks=c("0", "1" ),
                   labels=c("NonRewCol", "RewCol")) +
  scale_fill_manual(values=c("#E53304", "#164A7C"))+
  ggtitle("The Effect of Rewarded Color on\nAccuracy in the MID block") +
  theme(plot.title = element_text(face="bold", size=10)) +
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(text = element_text(size=20,face="bold"), axis.title=element_text(size=25,face="bold"))+
  coord_cartesian(ylim = c(50,100)) 

#------------- t test----------------
rew <- subset(MID.Acc, Description=="1")
norew <- subset(MID.Acc, Description=="0")

t.test(rew$Freq, norew$Freq, paired = TRUE) # variables are numeric 


########## NEUT ##############
#--------------dividing by rew col------------------------------

NEUT_RY <- subset(NEUT.df, Subject %in% c(9,11,12,13,14,21,22,23,24,25,26,35,36))

NEUT_RY$RewCol <- "RY"

NEUT_RY$Description <- revalue(NEUT_RY$Description, c("NEUT_yellowR"="RewCol", "NEUT_redR"="RewCol", "NEUT_blueY"="RewCol", "NEUT_yellowY"="RewCol", "NEUT_greenR"="RewCol", "NEUT_redY"="RewCol",
                                        "NEUT_redG"="NonRewCol", "NEUT_greenG"="NonRewCol", "NEUT_greenB"="NonRewCol", "NEUT_blueB"="NonRewCol", "NEUT_yellowB"="NonRewCol", "NEUT_blueG"="NonRewCol"))



NEUT_BG <- subset(NEUT.df, Subject %in% c(15,16,17,18,19,20,27,28,29,30,31,32,33,34))

NEUT_BG$RewCol <- "BG"

NEUT_BG$Description <- revalue(NEUT_BG$Description, c("NEUT_greenB"="RewCol", "NEUT_blueB"="RewCol", "NEUT_greenG"="RewCol", "NEUT_redG"="RewCol", "NEUT_yellowB"="RewCol", "NEUT_blueG"="RewCol",
                                                      "NEUT_yellowY"="NonRewCol", "NEUT_blueY"="NonRewCol", "NEUT_yellowR"="NonRewCol", "NEUT_redR"="NonRewCol", "NEUT_redY"="NonRewCol", "NEUT_greenR"="NonRewCol"))

NEUT <- rbind(NEUT_BG, NEUT_RY)

#---------------prepare for the stats-------
NEUT.rt <-NEUT[!NEUT$Accuracy %in% c("miss"),] # exclude misses 
NEUT.rt <-NEUT.rt[!NEUT.rt$Accuracy2 %in% c("incorrect"),]  # exclude incorrect 

NEUT.rt <- subset(NEUT.rt, TTime.Shifted <=1000)

av.rt <- aggregate(data=NEUT.rt, TTime.Shifted~Description*Subject, mean) 

#--------------boxplot ---------

means <- aggregate(TTime.Shifted ~ Description, av.rt, mean)
means$TTime.Shifted <- round(means$TTime.Shifted,digits=2 )

ggplot(data = av.rt, aes(x=Description, y=TTime.Shifted, fill=Description)) +
  geom_boxplot() +
  xlab("Description") + ylab("RTs [ms]")+
  #  stat_summary(fun.y=mean,col='black',geom='point', shape=18, size=3,show_guide = FALSE) +
  #  geom_text(data = means, aes(label = means$TTTime.Shifted, y = TTTime.Shifted + 0.08), vjust = -0.7) +
  scale_x_discrete(breaks=c("con", "inc" ),
                   labels=c("Congruent", "Incongruent")) +
  scale_fill_manual(values=c("#E53304", "#164A7C"))+
  ggtitle("The Effect of Rewarded Color on\nResponses in the NEUT block") +
  theme(plot.title = element_text(face="bold", size=10)) +
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(text = element_text(size=20,face="bold"), axis.title=element_text(size=25,face="bold"))+
  coord_cartesian(ylim = c(400,800)) 

#------------- t test----------------
rew <- subset(av.rt, Description=="RewCol")
norew <- subset(av.rt, Description=="NonRewCol")

t.test(rew$TTime.Shifted, norew$TTime.Shifted, paired = TRUE) # variables are numeric 

mean(rew$TTime.Shifted)
sd(rew$TTime.Shifted)
mean(norew$TTime.Shifted)
sd(norew$TTime.Shifted)

#---------- Accuracies -------------------
NEUT$Accuracy[NEUT$Accuracy2 == "incorrect"] <- "incorrect"
NEUT$Accuracy[is.na(NEUT$Accuracy)] <- 1
NEUT$Accuracy <- revalue(NEUT$Accuracy, c("incorrect"=0, "miss"=0))
NEUT$Description <- revalue(NEUT$Description, c("RewCol"=1, "NonRewCol"=0))

NEUT$Accuracy <- as.numeric(as.character(NEUT$Accuracy))
NEUT$Description <- as.numeric(as.character(NEUT$Description))

NEUT.Acc <- as.data.frame(xtabs(Accuracy~Description+Subject,NEUT)/xtabs(~Description+Subject,NEUT))
NEUT.Acc$Freq <- NEUT.Acc$Freq * 100


av.rt <- aggregate(data=NEUT.Acc, Freq~Description*Subject, mean) 

means <- aggregate(Freq ~ Description, av.rt, mean)
means$Freq <- round(means$Freq,digits=2 )

#--------------boxplot ---------

ggplot(data = av.rt, aes(x=Description, y=Freq, fill=Description)) +
  geom_boxplot() +
  xlab("Description") + ylab("Accuracy [%]")+
  #  stat_summary(fun.y=mean,col='black',geom='point', shape=18, size=3,show_guide = FALSE) +
  #  geom_text(data = means, aes(label = means$TTTime.Shifted, y = TTTime.Shifted + 0.08), vjust = -0.7) +
  scale_x_discrete(breaks=c("0", "1" ),
                   labels=c("NonRewCol", "RewCol")) +
  scale_fill_manual(values=c("#E53304", "#164A7C"))+
  ggtitle("The Effect of Rewarded Color on\nAccuracy in the NEUT block") +
  theme(plot.title = element_text(face="bold", size=10)) +
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(text = element_text(size=20,face="bold"), axis.title=element_text(size=25,face="bold"))+
  coord_cartesian(ylim = c(50,100)) 

#------------- t test----------------
rew <- subset(NEUT.Acc, Description=="1")
norew <- subset(NEUT.Acc, Description=="0")

t.test(rew$Freq, norew$Freq, paired = TRUE) # variables are numeric 

mean(rew$Freq)
sd(rew$Freq)
mean(norew$Freq)
sd(norew$Freq)


