# Author: MK
###############################################################################
###############################################################################
###############################################################################
################################################################
rm(list=ls()) # clear the working directory

#---------------setting up the path---------------
myPath <- "C:/Users/mariam/Documents/PhD/EXPERIMENTS/MK_MREEGRewModes_2016/Analysis/BehavData/Data_1000cutoff"
setwd(myPath)

#--------------import all libraries---------------
library(plyr)
library(dplyr)
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
MID.df <- subset(df, Description %in% c('MID_Rinc', "MID_Rcon", 'MID_Ucon', "MID_Uinc"))
NEUT.df <- subset(df, Description %in% c('Neut_inc', "Neut_con"))



#------------renaming the events-----------------
CSRA.df$CodeNew <- revalue(CSRA.df$Description , c("CSRA_Uinc"="U_inc", 
                                                   "CSRA_Ucon"="U_con",
                                                   "CSRA_Rinc"="R_inc", 
                                                   "CSRA_Rcon"="R_con"))

# ------------------split the data ------------------------
setDT(CSRA.df)[, paste0("CodeNew", 1:2) := tstrsplit(CodeNew, "_")]


# ---------------rename colums-----------------
CSRA.df <- plyr::rename(CSRA.df, c("CodeNew1"="Reward", "CodeNew2"="Congruency")) 


#---------------prepare for the stats-------
df_CSRA <-CSRA.df[!CSRA.df$Accuracy %in% c("miss"),] # exclude misses 
df_CSRA <-df_CSRA[!df_CSRA$Accuracy2 %in% c("incorrect"),]  # exclude incorrect 

df_CSRA$Block <- "CSRA"


#---------- Accuracies -------------------
df_resp <- CSRA.df[!CSRA.df$Accuracy %in% c("miss"),] # exclude misses
df_resp$Accuracy2 <- as.character(df_resp$Accuracy2)
df_resp$Accuracy2[is.na(df_resp$Accuracy2)] <- 1
df_resp$Accuracy2 <- revalue(df_resp$Accuracy2, c("incorrect"=0))
df_resp$Reward <- revalue(df_resp$Reward, c("R"=0, "U"=1))
df_resp$Congruency <- revalue(df_resp$Congruency, c("inc"=1, "con"=0))

df_resp$Accuracy2 <- as.numeric(as.character(df_resp$Accuracy2))
df_resp$Reward <- as.numeric(as.character(df_resp$Reward))
df_resp$Congruency <- as.numeric(as.character(df_resp$Congruency))

CSRA.Acc <- as.data.frame(xtabs(Accuracy2~Reward+Congruency+Subject,df_resp)/xtabs(~Reward+Congruency+Subject,df_resp))
CSRA.Acc$Freq <- CSRA.Acc$Freq * 100

CSRA.Acc$Block <- "CSRA"

########## SRA ##############
#------------renaming the events-----------------
SRA.df$CodeNew <- revalue(SRA.df$Description , c("SRA_Uinc"="U_inc", 
                                                 "SRA_Ucon"="U_con",
                                                 "SRA_Rinc"="R_inc", 
                                                 "SRA_Rcon"="R_con"))

# ------------------split the data ------------------------
setDT(SRA.df)[, paste0("CodeNew", 1:2) := tstrsplit(CodeNew, "_")]


# ---------------rename colums-----------------
SRA.df <- plyr::rename(SRA.df, c("CodeNew1"="Reward", "CodeNew2"="Congruency")) 


#---------------prepare for the stats-------
df_SRA <-SRA.df[!SRA.df$Accuracy %in% c("miss"),] # exclude misses 
df_SRA <-df_SRA[!df_SRA$Accuracy2 %in% c("incorrect"),]  # exclude incorrect 

df_SRA$Block <- "SRA"

#---------- Accuracies -------------------
df_resp <-SRA.df[!SRA.df$Accuracy %in% c("miss"),] # exclude misses
df_resp$Accuracy2 <- as.character(df_resp$Accuracy2)
df_resp$Accuracy2[is.na(df_resp$Accuracy2)] <- 1
df_resp$Accuracy2 <- revalue(df_resp$Accuracy2, c("incorrect"=0))
df_resp$Reward <- revalue(df_resp$Reward, c("R"=0, "U"=1))
df_resp$Congruency <- revalue(df_resp$Congruency, c("inc"=1, "con"=0))

df_resp$Accuracy2 <- as.numeric(as.character(df_resp$Accuracy2))
df_resp$Reward <- as.numeric(as.character(df_resp$Reward))
df_resp$Congruency <- as.numeric(as.character(df_resp$Congruency))

SRA.Acc <- as.data.frame(xtabs(Accuracy2~Reward+Congruency+Subject,df_resp)/xtabs(~Reward+Congruency+Subject,df_resp))
SRA.Acc$Freq <- SRA.Acc$Freq * 100

SRA.Acc$Block <- "SRA"


########## MID ##############
#------------renaming the events-----------------
MID.df$CodeNew <- revalue(MID.df$Description , c("MID_Uinc"="U_inc", 
                                                 "MID_Ucon"="U_con",
                                                 "MID_Rinc"="R_inc", 
                                                 "MID_Rcon"="R_con"))

# ------------------split the data ------------------------
setDT(MID.df)[, paste0("CodeNew", 1:2) := tstrsplit(CodeNew, "_")]


# ---------------rename colums-----------------
MID.df <- plyr::rename(MID.df, c("CodeNew1"="Reward", "CodeNew2"="Congruency")) 



#---------------prepare for the stats-------
df_MID <-MID.df[!MID.df$Accuracy %in% c("miss"),] # exclude misses 
df_MID <-df_MID[!df_MID$Accuracy2 %in% c("incorrect"),]  # exclude incorrect 

df_MID$Block <- "MID"

#---------- Accuracies -------------------
df_resp <-MID.df[!MID.df$Accuracy %in% c("miss"),] # exclude misses
df_resp$Accuracy2 <- as.character(df_resp$Accuracy2)
df_resp$Accuracy2[is.na(df_resp$Accuracy2)] <- 1
df_resp$Accuracy2 <- revalue(df_resp$Accuracy2, c("incorrect"=0))
df_resp$Reward <- revalue(df_resp$Reward, c("R"=0, "U"=1))
df_resp$Congruency <- revalue(df_resp$Congruency, c("inc"=1, "con"=0))

df_resp$Accuracy2 <- as.numeric(as.character(df_resp$Accuracy2))
df_resp$Reward <- as.numeric(as.character(df_resp$Reward))
df_resp$Congruency <- as.numeric(as.character(df_resp$Congruency))

MID.Acc <- as.data.frame(xtabs(Accuracy2~Reward+Congruency+Subject,df_resp)/xtabs(~Reward+Congruency+Subject,df_resp))
MID.Acc$Freq <- MID.Acc$Freq * 100

MID.Acc$Block <- "MID"

########## NEUT ##############
#------------renaming the events-----------------
NEUT.df$CodeNew <- revalue(NEUT.df$Description , c("Neut_inc"="inc", 
                                                   "Neut_con"="con"))

NEUT.df <- plyr::rename(NEUT.df, c("CodeNew"="Congruency")) 


#---------------prepare for the stats-------
df_NEUT <-NEUT.df[!NEUT.df$Accuracy %in% c("miss"),] # exclude misses 
df_NEUT <-df_NEUT[!df_NEUT$Accuracy2 %in% c("incorrect"),]  # exclude incorrect 

df_NEUT$Block <- "NEUT"

#---------- Accuracies -------------------
df_resp <-NEUT.df[!NEUT.df$Accuracy %in% c("miss"),] # exclude misses
df_resp$Accuracy2 <- as.character(df_resp$Accuracy2)
df_resp$Accuracy2 <- revalue(df_resp$Accuracy2, c("incorrect"=0))
df_resp$Accuracy2[is.na(df_resp$Accuracy2)] <- 1
df_resp$Congruency <- revalue(df_resp$Congruency, c("inc"=1, "con"=0))

df_resp$Accuracy2 <- as.numeric(as.character(df_resp$Accuracy2))
df_resp$Congruency <- as.numeric(as.character(df_resp$Congruency))

NEUT.Acc <- as.data.frame(xtabs(Accuracy2~Congruency+Subject,df_resp)/xtabs(~Congruency+Subject,df_resp))
NEUT.Acc$Freq <- NEUT.Acc$Freq * 100

NEUT.Acc$Block <- "NEUT"

################# 3x2x2 ANOVA #######################
#-------------- REWARD RT --------------
df_MID <- df_MID[which(df_MID$TTime.Shifted < 1000),]
df_SRA <- df_SRA[which(df_SRA$TTime.Shifted < 1000),]
df_CSRA <- df_CSRA[which(df_CSRA$TTime.Shifted < 1000),]
df_NEUT <- df_NEUT[which(df_NEUT$TTime.Shifted < 1000),]

av.rt_CSRA <- aggregate(data=df_CSRA, TTime.Shifted~Reward*Block*Congruency*Subject, mean) 
av.rt_SRA <- aggregate(data=df_SRA, TTime.Shifted~Reward*Block*Congruency*Subject, mean) 
av.rt_MID <- aggregate(data=df_MID, TTime.Shifted~Reward*Block*Congruency*Subject, mean) 
av.rt_NEUT <- aggregate(data=df_NEUT, TTime.Shifted~Block*Congruency*Subject, mean) 


df <- rbind(av.rt_CSRA, av.rt_SRA, av.rt_MID)

df$Reward <- as.factor(df$Reward)
df$Block <- as.factor(df$Block)
df$Subject <- as.factor(df$Subject)
df$TTime.Shifted <- as.numeric(as.character(df$TTime.Shifted))
df$Congruency <- as.factor(df$Congruency)

df <- subset(df, Reward=="R")

demoAnova <- ezANOVA(df, # specify data frame
                     dv = TTime.Shifted, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(Congruency, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$TTime.Shifted, df$Block, p.adjust.method = "none"  )

t.test(av.rt_MID$TTime.Shifted, av.rt_SRA$TTime.Shifted, paired = TRUE)

# post hoc interactions
df_R <- subset(df, Reward=="R")
df_U <- subset(df, Reward=="U")
df_R$rew_effect <- df_R$TTime.Shifted-df_U$TTime.Shifted

dataSRA <- subset(df_R, Block == "SRA")
dataCSRA <- subset(df_R, Block == "CSRA")
dataMID <- subset(df_R, Block == "MID")
t.test(dataSRA$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataSRA$rew_effect, dataMID$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataCSRA$rew_effect, dataMID$rew_effect, paired = TRUE) # variables are numeric 

sum = summarySEwithin(df, measurevar="TTime.Shifted", withinvars=c("Block", "Congruency"), idvar="Subject", na.rm=FALSE, conf.interval=.95)


df_Con <- subset(df, Congruency=="con")
df_Inc <- subset(df, Congruency=="inc")
df_Con$inc_effect <- df_Con$TTime.Shifted-df_Inc$TTime.Shifted

dataR <- subset(df_Con, Reward == "R")
dataU <- subset(df_Con, Reward == "U")
t.test(dataR$inc_effect, dataU$inc_effect, paired = TRUE) # variables are numeric 
mean(dataR$inc_effect)
mean(dataU$inc_effect)

dataSRA <- subset(df_Con, Block == "SRA")
dataCSRA <- subset(df_Con, Block == "CSRA")
dataMID <- subset(df_Con, Block == "MID")
t.test(dataSRA$inc_effect, dataCSRA$inc_effect, paired = TRUE) # variables are numeric 
t.test(dataSRA$inc_effect, dataMID$inc_effect, paired = TRUE) # variables are numeric 
t.test(dataCSRA$inc_effect, dataMID$inc_effect, paired = TRUE) # variables are numeric 


mean(df$TTime.Shifted[df$Reward == "R"])
sd(df$TTime.Shifted[df$Reward == "R"])

mean(df$TTime.Shifted[df$Reward == "U"])
sd(df$TTime.Shifted[df$Reward == "U"])

mean(df$TTime.Shifted[df$Congruency == "con"])
sd(df$TTime.Shifted[df$Congruency == "con"])

mean(df$TTime.Shifted[df$Congruency == "inc"])
sd(df$TTime.Shifted[df$Congruency == "inc"])

mean(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "R" & av.rt_CSRA$Congruency == "con"])
sd(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "R" & av.rt_CSRA$Congruency == "con"])

mean(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "U" & av.rt_CSRA$Congruency == "con"])
sd(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "U" & av.rt_CSRA$Congruency == "con"])

mean(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "R" & av.rt_CSRA$Congruency == "inc"])
sd(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "R" & av.rt_CSRA$Congruency == "inc"])

mean(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "U" & av.rt_CSRA$Congruency == "inc"])
sd(av.rt_CSRA$TTime.Shifted[av.rt_CSRA$Reward == "U" & av.rt_CSRA$Congruency == "inc"])

mean(av.rt_NEUT$TTime.Shifted[av.rt_NEUT$Congruency == "con"])
sd(av.rt_NEUT$TTime.Shifted[av.rt_NEUT$Congruency == "con"])

mean(av.rt_NEUT$TTime.Shifted[av.rt_NEUT$Congruency == "inc"])
sd(av.rt_NEUT$TTime.Shifted[av.rt_NEUT$Congruency == "inc"])

av.rt_NEUT -> av.rt_NEUT_Rew
av.rt_NEUT -> av.rt_NEUT_NoRew
av.rt_NEUT_Rew$Reward <- "R"
av.rt_NEUT_NoRew$Reward <- "U"

df2 <- rbind(av.rt_NEUT_Rew, av.rt_NEUT_NoRew)
df <- rbind(av.rt_CSRA, av.rt_SRA, av.rt_MID, df2)

# barplot
sum = summarySEwithin(df, measurevar="TTime.Shifted", withinvars=c("Block","Reward", "Congruency"), idvar="Subject", na.rm=FALSE, conf.interval=.95)


sum$Reward <- revalue(sum$Reward , c("U"="No Reward", "R"="Reward"))
sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("con"="Congruent", "inc"="Incongruent"))

ggplot(sum, 
       aes(x=Congruency, y=TTime.Shifted, group=Block, 
           ymax=TTime.Shifted+se, ymin=TTime.Shifted-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Block)) +
  facet_wrap(  ~ Reward) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Condition",
       y = "RTs [ms]")  +
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
  coord_cartesian(ylim = c(500, 700))+
 # scale_fill_brewer()+
  scale_fill_manual(name="Block", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#F0E442", "#009E73")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# barplot #2 
Rcon <- subset(df, Reward=="R" & Congruency=="con")
Rinc <- subset(df, Reward=="R" & Congruency=="inc")
Ucon <- subset(df, Reward=="U" & Congruency=="con")
Uinc <- subset(df, Reward=="U" & Congruency=="inc")
Con <- subset(av.rt_NEUT, Congruency=="con")
Inc <- subset(av.rt_NEUT, Congruency=="inc")

Rcon$Condition <- 'R_con'
Rinc$Condition <- 'R_inc'
Ucon$Condition <- "NR_con"
Uinc$Condition <- 'NR_inc'
Con$Condition <- 'con'
Inc$Condition <- 'inc'

df3 <- rbind(Rcon, Rinc, Ucon, Uinc)

# C-SRA + MID + SRA
sum = summarySEwithin(df3, measurevar="TTime.Shifted", withinvars=c("Block","Condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Block <- revalue(sum$Block , c("NEUT"="Neutral"))
#sum <-sum[!sum$Block=="NEUT",] 
sum$Block <- factor(sum$Block, levels = c("C-SRA", "SRA", "MID", "Neutral"))
sum$Condition <- factor(sum$Condition, levels = c("R_con", "NR_con", "R_inc", "NR_inc"))
sum <- sum[-9, ]
sum <- sum[-9, ]

ggplot(sum, 
       aes(x=Block, y=TTime.Shifted, group=Condition, 
           ymax=TTime.Shifted+se, ymin=TTime.Shifted-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Condition)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black") +
  labs(x = "",
       y = "RTs [ms]")  +
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
  coord_cartesian(ylim = c(500, 700))+
  # scale_fill_brewer()+
  scale_fill_manual(name="Condition\n(except Neutral)", # Legend label, use darker colors
                    values=c("#CC6666", "#9999CC", "#66CC99", "blanchedalmond")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
   

#-------------- REWARD ACCURACIES  --------------
av.er_CSRA <- aggregate(data=CSRA.Acc, Freq~Reward*Block*Congruency*Subject, mean) 
av.er_SRA <- aggregate(data=SRA.Acc, Freq~Reward*Block*Congruency*Subject, mean) 
av.er_MID <- aggregate(data=MID.Acc, Freq~Reward*Block*Congruency*Subject, mean) 
av.er_NEUT <- aggregate(data=NEUT.Acc, Freq~Block*Congruency*Subject, mean) 

df <- rbind(av.er_CSRA, av.er_SRA, av.er_MID)

df$Reward <- as.factor(df$Reward)
df$Block <- as.factor(df$Block)
df$Subject <- as.factor(df$Subject)
df$Freq <- as.numeric(as.character(df$Freq))
df$Congruency <- as.factor(df$Congruency)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Freq, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(Reward, Block, Congruency), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

mean(NEUT.Acc$Freq[NEUT.Acc$Congruency == "0"])
sd(NEUT.Acc$Freq[NEUT.Acc$Congruency == "0"])

mean(av.er_CSRA$Freq[av.er_CSRA$Reward == "0" & av.er_CSRA$Congruency == "0"])
sd(av.er_CSRA$Freq[av.er_CSRA$Reward == "0" & av.er_CSRA$Congruency == "0"])

mean(NEUT.Acc$Freq[NEUT.Acc$Reward == "0" & NEUT.Acc$Congruency == "1"])
sd(NEUT.Acc$Freq[NEUT.Acc$Reward == "0" & NEUT.Acc$Congruency == "1"])

mean(NEUT.Acc$Freq[NEUT.Acc$Reward == "1" & NEUT.Acc$Congruency == "1"])
sd(NEUT.Acc$Freq[NEUT.Acc$Reward == "1" & NEUT.Acc$Congruency == "1"])


mean(df$Freq[df$Reward == "0"])
sd(df$Freq[df$Reward == "0"])

mean(df$Freq[df$Reward == "1"])
sd(df$Freq[df$Reward == "1"])

mean(df$Freq[df$Congruency == "0"])
sd(df$Freq[df$Congruency == "0"])

mean(df$Freq[df$Congruency == "1"])
sd(df$Freq[df$Congruency == "1"])

# post hoc interactions
df_R <- subset(df, Reward=="0")
df_U <- subset(df, Reward=="1")
df_R$rew_effect <- df_R$Freq-df_U$Freq

dataSRA <- subset(df_R, Block == "SRA")
dataCSRA <- subset(df_R, Block == "CSRA")
dataMID <- subset(df_R, Block == "MID")
t.test(dataSRA$rew_effect, dataCSRA$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataSRA$rew_effect, dataMID$rew_effect, paired = TRUE) # variables are numeric 
t.test(dataCSRA$rew_effect, dataMID$rew_effect, paired = TRUE) # variables are numeric 
mean(dataMID$rew_effect)
mean(dataSRA$rew_effect)

av.er_NEUT -> av.er_NEUT_Rew
av.er_NEUT -> av.er_NEUT_NoRew
av.er_NEUT_Rew$Reward <- "0"
av.er_NEUT_NoRew$Reward <- "1"

df2 <- rbind(av.er_NEUT_Rew, av.er_NEUT_NoRew)
df <- rbind(av.er_CSRA, av.er_SRA, av.er_MID, df2)


# barplot
sum = summarySE(df, measurevar="Freq", groupvars=c("Block","Reward", "Congruency"))

sum$Reward <- revalue(sum$Reward , c("1"="No Reward", "0"="Reward"))
sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("1"="Incongruent", "0"="Congruent"))

sum$Reward <- factor(sum$Reward, levels = c("No Reward", "Reward"))

ggplot(sum, 
       aes(x=Congruency, y=Freq, group=Block, 
           ymax=Freq+se, ymin=Freq-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Block)) +
  facet_wrap(  ~ Reward) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Condition",
       y = "Acuracy [%]")  +
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
  coord_cartesian(ylim = c(70, 100))+
#  scale_fill_brewer()+
  scale_fill_manual(name="Block", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#F0E442", "#009E73")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

sum = summarySE(df, measurevar="Freq", groupvars=c("Block", "Congruency"))
sum = summarySE(df, measurevar="Freq", groupvars=c("Block", "Reward"))

# barplot #2 
Rcon <- subset(df, Reward=="0" & Congruency=="0")
Rinc <- subset(df, Reward=="0" & Congruency=="1")
Ucon <- subset(df, Reward=="1" & Congruency=="0")
Uinc <- subset(df, Reward=="1" & Congruency=="1")


Rcon$Condition <- 'R_con'
Rinc$Condition <- 'R_inc'
Ucon$Condition <- "NR_con"
Uinc$Condition <- 'NR_inc'


df3 <- rbind(Rcon, Rinc, Ucon, Uinc)

sum = summarySEwithin(df3, measurevar="Freq", withinvars=c("Block","Condition"), idvar="Subject", na.rm=FALSE, conf.interval=.95)

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Block <- revalue(sum$Block , c("NEUT"="Neutral"))
#sum <-sum[!sum$Block=="NEUT",] 
sum$Block <- factor(sum$Block, levels = c("C-SRA", "SRA", "MID", "Neutral"))
sum$Condition <- factor(sum$Condition, levels = c("R_con", "NR_con", "R_inc", "NR_inc"))
sum <- sum[-9, ]
sum <- sum[-9, ]

ggplot(sum, 
       aes(x=Block, y=Freq, group=Condition, 
           ymax=Freq+se, ymin=Freq-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Condition)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black") +
  labs(x = "Block",
       y = "Accuracy [%]")  +
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
  coord_cartesian(ylim = c(70, 100))+
  # scale_fill_brewer()+
  scale_fill_manual(name="Condition\n(except Neutral)", # Legend label, use darker colors
                    values=c("#CC6666", "#9999CC", "#66CC99", "blanchedalmond")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=25), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

################# 4x2 ANOVA #######################
# re-run everyting before 3x2x2 ANOVA

# Reward RTs
CSRA_R <- subset(df_CSRA, Reward == 'R')
SRA_R <- subset(df_SRA, Reward == 'R')
MID_R <- subset(df_MID,Reward == 'R')

av.rt_CSRA <- aggregate(data=CSRA_R, TTime.Shifted~Block*Congruency*Subject, mean) 
av.rt_SRA <- aggregate(data=SRA_R, TTime.Shifted~Block*Congruency*Subject, mean) 
av.rt_MID <- aggregate(data=MID_R, TTime.Shifted~Block*Congruency*Subject, mean) 
av.rt_Neut <- aggregate(data=df_NEUT, TTime.Shifted~Block*Congruency*Subject, mean) 

df <- rbind(av.rt_CSRA, av.rt_SRA, av.rt_MID, av.rt_Neut)

df$Block <- as.factor(df$Block)
df$Subject <- as.factor(df$Subject)
df$TTime.Shifted <- as.numeric(as.character(df$TTime.Shifted))
df$Congruency <- as.factor(df$Congruency)


demoAnova <- ezANOVA(df, # specify data frame
                     dv = TTime.Shifted, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(Congruency, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

mean(df$TTime.Shifted[df$Congruency == "con"])
sd(df$TTime.Shifted[df$Congruency == "con"])

mean(df$TTime.Shifted[df$Congruency == "inc"])
sd(df$TTime.Shifted[df$Congruency == "inc"])

pairwise.t.test(df$TTime.Shifted, df$Block, p.adjust.method = "hochberg"  )

# post hoc interactions
df_inc <- subset(df, Congruency == 'inc')
df_con <- subset(df, Congruency == 'con')
df_inc$con_effect <- df_inc$TTime.Shifted-df_con$TTime.Shifted

dataSRA <- subset(df_inc, Block == "SRA")
dataCSRA <- subset(df_inc, Block == "CSRA")
dataMID <- subset(df_inc, Block == "MID")
dataNEUT <- subset(df_inc, Block == "NEUT")
t.test(dataMID$con_effect, dataSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataMID$con_effect, dataCSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataMID$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$con_effect, dataCSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric
t.test(dataCSRA$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric
mean(dataNEUT$con_effect)
mean(dataMID$con_effect)
mean(dataSRA$con_effect)
mean(dataCSRA$con_effect)

# barplot
sum = summarySE(df, measurevar="TTime.Shifted", groupvars=c("Block","Congruency"))

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("con"="Congruent", "inc"="Incongruent"))

ggplot(sum, 
       aes(x=Congruency, y=TTime.Shifted, group=Block, 
           ymax=TTime.Shifted+se, ymin=TTime.Shifted-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Block)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Condition in Rewarded Trials",
       y = "RTs [ms]")  +
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
  coord_cartesian(ylim = c(400, 800))+
  scale_fill_manual(name="Block", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#F0E442", "#009E73")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



# Reward Acc

CSRA_R <- subset(CSRA.Acc, Reward == 0)
SRA_R <- subset(SRA.Acc, Reward == 0)
MID_R <- subset(MID.Acc,Reward == 0)

av.er_CSRA <- aggregate(data=CSRA_R, Freq~Block*Congruency*Subject, mean) 
av.er_SRA <- aggregate(data=SRA_R, Freq~Block*Congruency*Subject, mean) 
av.er_MID <- aggregate(data=MID_R, Freq~Block*Congruency*Subject, mean) 
av.er_Neut <- aggregate(data=NEUT.Acc, Freq~Block*Congruency*Subject, mean) 

df <- rbind(av.er_CSRA, av.er_SRA, av.er_MID, av.er_Neut)

df$Block <- as.factor(df$Block)
df$Subject <- as.factor(df$Subject)
df$Freq <- as.numeric(as.character(df$Freq))
df$Congruency <- as.factor(df$Congruency)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Freq, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(Block, Congruency), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

mean(df$Freq[df$Congruency == "0"])
sd(df$Freq[df$Congruency == "0"])

mean(df$Freq[df$Congruency == "1"])
sd(df$Freq[df$Congruency == "1"])


pairwise.t.test(df$Freq, df$Block, p.adjust.method = "none"  )

# post hoc interactions
df_inc <- subset(df, Congruency == '0')
df_con <- subset(df, Congruency == '1')
df_inc$con_effect <- df_inc$Freq-df_con$Freq

dataSRA <- subset(df_inc, Block == "SRA")
dataCSRA <- subset(df_inc, Block == "CSRA")
dataMID <- subset(df_inc, Block == "MID")
dataNEUT <- subset(df_inc, Block == "NEUT")
t.test(dataMID$con_effect, dataSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataMID$con_effect, dataCSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataMID$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$con_effect, dataCSRA$con_effect, paired = TRUE) # variables are numeric
t.test(dataSRA$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric
t.test(dataCSRA$con_effect, dataNEUT$con_effect, paired = TRUE) # variables are numeric

# barplot
sum = summarySE(df, measurevar="Freq", groupvars=c("Block", "Congruency"))

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("1"="Incongruent", "0"="Congruent"))

ggplot(sum, 
       aes(x=Congruency, y=Freq, group=Block, 
           ymax=Freq+se, ymin=Freq-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Block)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Condition in Rewarded Trials",
       y = "Acuracy [%]")  +
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
  coord_cartesian(ylim = c(70, 100))+
  scale_fill_manual(name="Block", # Legend label, use darker colors
                    values=c("#E69F00", "#56B4E9", "#F0E442", "#009E73")) +
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



# No Reward RTs
CSRA_R <- subset(df_CSRA, Reward == 'U')
SRA_R <- subset(df_SRA, Reward == 'U')
MID_R <- subset(df_MID,Reward == 'U')

av.rt_CSRA <- aggregate(data=CSRA_R, TTime.Shifted~Block*Congruency*Subject, mean) 
av.rt_SRA <- aggregate(data=SRA_R, TTime.Shifted~Block*Congruency*Subject, mean) 
av.rt_MID <- aggregate(data=MID_R, TTime.Shifted~Block*Congruency*Subject, mean) 
av.rt_Neut <- aggregate(data=df_NEUT, TTime.Shifted~Block*Congruency*Subject, mean) 

df <- rbind(av.rt_CSRA, av.rt_SRA, av.rt_MID, av.rt_Neut)

df$Block <- as.factor(df$Block)
df$Subject <- as.factor(df$Subject)
df$TTime.Shifted <- as.numeric(as.character(df$TTime.Shifted))
df$Congruency <- as.factor(df$Congruency)


demoAnova <- ezANOVA(df, # specify data frame
                     dv = TTime.Shifted, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(Congruency, Block), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)

pairwise.t.test(df$TTime.Shifted, df$Block, p.adjust.method = "hochberg"  )


mean(df$TTime.Shifted[df$Congruency == "con"])
sd(df$TTime.Shifted[df$Congruency == "con"])

mean(df$TTime.Shifted[df$Congruency == "inc"])
sd(df$TTime.Shifted[df$Congruency == "inc"])


# barplot
sum = summarySE(df, measurevar="TTime.Shifted", groupvars=c("Block","Congruency"))

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("con"="Congruent", "inc"="Incongruent"))

ggplot(sum, 
       aes(x=Congruency, y=TTime.Shifted, group=Block, 
           ymax=TTime.Shifted+se, ymin=TTime.Shifted-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Block)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Block",
       y = "RTs [ms]")  +
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
  coord_cartesian(ylim = c(400, 800))+
  scale_fill_brewer()+
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



# Reward Acc

CSRA_R <- subset(CSRA.Acc, Reward == 1)
SRA_R <- subset(SRA.Acc, Reward == 1)
MID_R <- subset(MID.Acc,Reward == 1)

av.er_CSRA <- aggregate(data=CSRA_R, Freq~Block*Congruency*Subject, mean) 
av.er_SRA <- aggregate(data=SRA_R, Freq~Block*Congruency*Subject, mean) 
av.er_MID <- aggregate(data=MID_R, Freq~Block*Congruency*Subject, mean) 
av.er_Neut <- aggregate(data=NEUT.Acc, Freq~Block*Congruency*Subject, mean) 

df <- rbind(av.er_CSRA, av.er_SRA, av.er_MID, av.er_Neut)

df$Block <- as.factor(df$Block)
df$Subject <- as.factor(df$Subject)
df$Freq <- as.numeric(as.character(df$Freq))
df$Congruency <- as.factor(df$Congruency)

demoAnova <- ezANOVA(df, # specify data frame
                     dv = Freq, # specify dependent variable 
                     wid = Subject, # specify the subject variable
                     within = .(Block, Congruency), # specify within-subject variables
                     detailed = TRUE, # get a detailed table that includes SS
                     type = 3
)
demoAnova
demoAnova2 <- anova_out(demoAnova,etasq="partial",print=TRUE)


mean(df$Freq[df$Congruency == "0"])
sd(df$Freq[df$Congruency == "0"])

mean(df$Freq[df$Congruency == "1"])
sd(df$Freq[df$Congruency == "1"])


# barplot
sum = summarySE(df, measurevar="Freq", groupvars=c("Block", "Congruency"))

sum$Block <- revalue(sum$Block , c("CSRA"="C-SRA"))
sum$Congruency <- revalue(sum$Congruency , c("1"="Incongruent", "0"="Congruent"))

ggplot(sum, 
       aes(x=Congruency, y=Freq, group=Block, 
           ymax=Freq+se, ymin=Freq-se))  +
  geom_bar(stat="identity", position = "dodge", aes(fill=Block)) +
  geom_errorbar(position=position_dodge(width=0.7), 
                width=0.0, size=0.5, color="black")  +
  labs(x = "Block",
       y = "Acuracy [%]")  +
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
  coord_cartesian(ylim = c(70, 100))+
  scale_fill_brewer()+
  theme(panel.grid.major = element_blank(), text = element_text(size=15), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

