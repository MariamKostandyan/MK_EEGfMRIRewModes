#--------------------------------------------------------------
rm(list=ls()) # clear the working directory

myPath <- "C:/Users/mariam/Documents/PhD/EXPERIMENTS/MK_MREEGRewModes_2016/Analysis/EEGdata/plots"  
setwd(myPath)

#--------------import all libraries---------------
library(dplyr)
library(plyr)
library(useful) # for shifting the columns 
library(stringr)
library(data.table)
library(ggplot2)
library(Rmisc)


#-------------- P300 reward ---------------------
Rew <- read.table("P300_R_300_400_Filters.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                  dec = ".",
                  na.strings = "NaN", colClasses = NA, nrows = -1,
                  skip = 0, check.names = TRUE,
                  strip.white = FALSE, blank.lines.skip = TRUE,
                  comment.char = "#",
                  allowEscapes = FALSE, flush = FALSE)

Diff <- read.table("P300_RminusU_300_400_Diff. Waves.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                   dec = ".",
                   na.strings = "NaN", colClasses = NA, nrows = -1,
                   skip = 0, check.names = TRUE,
                   strip.white = FALSE, blank.lines.skip = TRUE,
                   comment.char = "#",
                   allowEscapes = FALSE, flush = FALSE)

NoRew <- read.table("P300_U_300_400_Filters.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                    dec = ".",
                    na.strings = "NaN", colClasses = NA, nrows = -1,
                    skip = 0, check.names = TRUE,
                    strip.white = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#",
                    allowEscapes = FALSE, flush = FALSE)

Rew <- plyr::rename(Rew, c("FCz"="R"))
Diff <- plyr::rename(Diff, c("FCz"="R-U"))
NoRew <- plyr::rename(NoRew, c("FCz"="U"))

df <- cbind(Rew, NoRew)
df <- data.frame(df, check.names=F)
df <- data.frame(a = c(-200,-198,-196,-194,-192,-190,-188,-186,-184,-182,-180,-178,-176,-174,-172,-170,-168,-166,-164,-162,-160,-158,-156,-154,-152,-150,-148,-146,-144,-142,-140,-138,-136,-134,-132,-130,-128,-126,-124,-122,-120,-118,-116,-114,-112,-110,-108,-106,-104,-102,-100,-98,-96,-94,-92,-90,-88,-86,-84,-82,-80,-78,-76,-74,-72,-70,-68,-66,-64,-62,-60,-58,-56,-54,-52,-50,-48,-46,-44,-42,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,
                       140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,
                       476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,552,554,556,558,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600,602,604,606,608,610,612,614,616,618,620,622,624,626,628,630,632,634,636,638,640,642,644,646,648,650,652,654,656,658,660,662,664,666,668,670,672,674,676,678,680,682,684,686,688,690,692,694,696,698,700,702,704,706,708,710,712,714,716,718,720,722,724,726,728,730,732,734,736,738,740,742,744,746,748,750,752,754,756,758,760,762,764,766,768,770,772,774,776,778,780,782,784,786,788,790,792,794,796,798,800,802,804,806,808,810,
                       812,814,816,818,820,822,824,826,828,830,832,834,836,838,840,842,844,846,848,850,852,854,856,858,860,862,864,866,868,870,872,874,876,878,880,882,884,886,888,890,892,894,896,898,900,902,904,906,908,910,912,914,916,918,920,922,924,926,928,930,932,934,936,938,940,942,944,946,948,950,952,954,956,958,960,962,964,966,968,970,972,974,976,978,980,982,984,986,988,990,992,994,996,998,1000,1002,1004,1006,1008,1010,1012,1014,1016,1018,1020,1022,1024,1026,1028,1030,1032,1034,1036,1038,1040,1042,1044,1046,1048,1050,1052,1054,1056,1058,1060,1062,1064,1066,1068,1070,1072,1074,1076,1078,1080,1082,1084,1086,1088,1090,1092,1094,1096,1098,1100,1102,1104,1106,1108,1110,1112,1114,1116,
                       1118,1120,1122,1124,1126,1128,1130,1132,1134,1136,1138,1140,1142,1144,1146,1148,1150,1152,1154,1156,1158,1160,1162,1164,1166,1168,1170,1172,1174,1176,1178,1180,1182,1184,1186,1188,1190,1192,1194,1196,1198), df,  check.names=F)
d <- melt(df,  id.vars = 'a', variable.name = 'series')

rect <- data.frame(xmin=450, xmax=800, ymin=-Inf, ymax=Inf)
ggplot(d, aes(a,value)) +
  geom_line(aes(color = series), size = 2)+
  scale_y_continuous(breaks=seq(-2.5, 2.5, 1.5)) +
  coord_cartesian(ylim=c(-2.5, 2.5))+
  #scale_x_continuous(breaks=seq(-200, 1200, 200)) +
  scale_color_manual(values=c("red", "black", "blue")) +
  # scale_linetype_manual(values=c("solid", "solid", "twodash")) +
  theme(legend.direction = 'vertical', 
        legend.position = 'right',
        legend.key = element_rect(size = 7),
        legend.key.size = unit(3, 'lines')) +
  geom_vline(xintercept=c(0), linetype="dotted", size=1.5) +
  theme(panel.grid.major = element_blank(), text = element_text(size=30), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.y = element_line(colour = "black"))+
  geom_vline(xintercept=c(0), linetype="dotted", size=1.5)+
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "gray88",
            color="gray88",
            alpha=0.4,
            inherit.aes = FALSE) +
  labs(title = "", x = "", y = "", color = "") +
  labs(x = "Time [ms]",
       y = "Amplitude [µV]") +
  theme(axis.text.x=element_text(size=20, vjust=80),
        axis.title.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  geom_hline(yintercept=0,  
             color = "black", size=0.5)


#-------------- N1 reward ---------------------
Rew <- read.table("N1_Rew_Filters.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                  dec = ".",
                  na.strings = "NaN", colClasses = NA, nrows = -1,
                  skip = 0, check.names = TRUE,
                  strip.white = FALSE, blank.lines.skip = TRUE,
                  comment.char = "#",
                  allowEscapes = FALSE, flush = FALSE)


NoRew <- read.table("N1_NoRew_Filters.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                    dec = ".",
                    na.strings = "NaN", colClasses = NA, nrows = -1,
                    skip = 0, check.names = TRUE,
                    strip.white = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#",
                    allowEscapes = FALSE, flush = FALSE)

Rew$Mean <- rowMeans(Rew[,1:2])
NoRew$Mean <- rowMeans(NoRew[,1:2])

Rew <- plyr::rename(Rew, c("Mean"="R"))
NoRew <- plyr::rename(NoRew, c("Mean"="U"))

Rew <- Rew[,-(1:2),drop=FALSE]
NoRew <- NoRew[,-(1:2),drop=FALSE]

df <- cbind(Rew, NoRew)
df <- data.frame(df, check.names=F)
df <- data.frame(a = c(-200,-198,-196,-194,-192,-190,-188,-186,-184,-182,-180,-178,-176,-174,-172,-170,-168,-166,-164,-162,-160,-158,-156,-154,-152,-150,-148,-146,-144,-142,-140,-138,-136,-134,-132,-130,-128,-126,-124,-122,-120,-118,-116,-114,-112,-110,-108,-106,-104,-102,-100,-98,-96,-94,-92,-90,-88,-86,-84,-82,-80,-78,-76,-74,-72,-70,-68,-66,-64,-62,-60,-58,-56,-54,-52,-50,-48,-46,-44,-42,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,
                       140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,
                       476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,552,554,556,558,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600,602,604,606,608,610,612,614,616,618,620,622,624,626,628,630,632,634,636,638,640,642,644,646,648,650,652,654,656,658,660,662,664,666,668,670,672,674,676,678,680,682,684,686,688,690,692,694,696,698,700,702,704,706,708,710,712,714,716,718,720,722,724,726,728,730,732,734,736,738,740,742,744,746,748,750,752,754,756,758,760,762,764,766,768,770,772,774,776,778,780,782,784,786,788,790,792,794,796,798,800,802,804,806,808,810,
                       812,814,816,818,820,822,824,826,828,830,832,834,836,838,840,842,844,846,848,850,852,854,856,858,860,862,864,866,868,870,872,874,876,878,880,882,884,886,888,890,892,894,896,898,900,902,904,906,908,910,912,914,916,918,920,922,924,926,928,930,932,934,936,938,940,942,944,946,948,950,952,954,956,958,960,962,964,966,968,970,972,974,976,978,980,982,984,986,988,990,992,994,996,998,1000,1002,1004,1006,1008,1010,1012,1014,1016,1018,1020,1022,1024,1026,1028,1030,1032,1034,1036,1038,1040,1042,1044,1046,1048,1050,1052,1054,1056,1058,1060,1062,1064,1066,1068,1070,1072,1074,1076,1078,1080,1082,1084,1086,1088,1090,1092,1094,1096,1098,1100,1102,1104,1106,1108,1110,1112,1114,1116,
                       1118,1120,1122,1124,1126,1128,1130,1132,1134,1136,1138,1140,1142,1144,1146,1148,1150,1152,1154,1156,1158,1160,1162,1164,1166,1168,1170,1172,1174,1176,1178,1180,1182,1184,1186,1188,1190,1192,1194,1196,1198), df,  check.names=F)
d <- melt(df,  id.vars = 'a', variable.name = 'series')

rect <- data.frame(xmin=160, xmax=200, ymin=-Inf, ymax=Inf)
ggplot(d, aes(a,value)) +
  geom_line(aes(color = series), size = 2)+
  scale_y_continuous(breaks=seq(-2.5, 2.5, 1.5)) +
  coord_cartesian(ylim=c(-2.5, 2.5))+
  #scale_x_continuous(breaks=seq(-200, 1200, 200)) +
  scale_color_manual(values=c("red", "black", "blue")) +
  # scale_linetype_manual(values=c("solid", "solid", "twodash")) +
  theme(legend.direction = 'vertical', 
        legend.position = 'right',
        legend.key = element_rect(size = 7),
        legend.key.size = unit(3, 'lines')) +
  geom_vline(xintercept=c(0), linetype="dotted", size=1.5) +
  theme(panel.grid.major = element_blank(), text = element_text(size=30), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.y = element_line(colour = "black"))+
  geom_vline(xintercept=c(0), linetype="dotted", size=1.5)+
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "gray88",
            color="gray88",
            alpha=0.4,
            inherit.aes = FALSE) +
  labs(title = "", x = "", y = "", color = "") +
  labs(x = "Time [ms]",
       y = "Amplitude [µV]") +
  theme(axis.text.x=element_text(size=20, vjust=80),
        axis.title.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  geom_hline(yintercept=0,  
             color = "black", size=0.5)


#-------------- congruency ---------------------
Con <- read.table("NEW_GA_con_Filters.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                  dec = ".",
                  na.strings = "NaN", colClasses = NA, nrows = -1,
                  skip = 0, check.names = TRUE,
                  strip.white = FALSE, blank.lines.skip = TRUE,
                  comment.char = "#",
                  allowEscapes = FALSE, flush = FALSE)

Diff <- read.table("NEW_GA_inc_Diff. Waves.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                   dec = ".",
                   na.strings = "NaN", colClasses = NA, nrows = -1,
                   skip = 0, check.names = TRUE,
                   strip.white = FALSE, blank.lines.skip = TRUE,
                   comment.char = "#",
                   allowEscapes = FALSE, flush = FALSE)

Inc <- read.table("NEW_GA_inc_Filters.txt", header = TRUE, fill=TRUE, sep = "", quote = "\"'",
                    dec = ".",
                    na.strings = "NaN", colClasses = NA, nrows = -1,
                    skip = 0, check.names = TRUE,
                    strip.white = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#",
                    allowEscapes = FALSE, flush = FALSE)

Con$Mean <- rowMeans(Con[,1:3])
Diff$Mean <- rowMeans(Diff[,1:3])
Inc$Mean <- rowMeans(Inc[,1:3])

Con <- plyr::rename(Con, c("Mean"="Con"))
Diff <- plyr::rename(Diff, c("Mean"="Inc-Con"))
Inc <- plyr::rename(Inc, c("Mean"="Inc"))

Con <- Con[,-(1:3),drop=FALSE]
Diff <- Diff[,-(1:3),drop=FALSE]
Inc <- Inc[,-(1:3),drop=FALSE]

df <- cbind(Con, Inc)
df <- data.frame(df, check.names=F)
df <- data.frame(a = c(-200,-198,-196,-194,-192,-190,-188,-186,-184,-182,-180,-178,-176,-174,-172,-170,-168,-166,-164,-162,-160,-158,-156,-154,-152,-150,-148,-146,-144,-142,-140,-138,-136,-134,-132,-130,-128,-126,-124,-122,-120,-118,-116,-114,-112,-110,-108,-106,-104,-102,-100,-98,-96,-94,-92,-90,-88,-86,-84,-82,-80,-78,-76,-74,-72,-70,-68,-66,-64,-62,-60,-58,-56,-54,-52,-50,-48,-46,-44,-42,-40,-38,-36,-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,
                       140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,246,248,250,252,254,256,258,260,262,264,266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302,304,306,308,310,312,314,316,318,320,322,324,326,328,330,332,334,336,338,340,342,344,346,348,350,352,354,356,358,360,362,364,366,368,370,372,374,376,378,380,382,384,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,454,456,458,460,462,464,466,468,470,472,474,
                       476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,552,554,556,558,560,562,564,566,568,570,572,574,576,578,580,582,584,586,588,590,592,594,596,598,600,602,604,606,608,610,612,614,616,618,620,622,624,626,628,630,632,634,636,638,640,642,644,646,648,650,652,654,656,658,660,662,664,666,668,670,672,674,676,678,680,682,684,686,688,690,692,694,696,698,700,702,704,706,708,710,712,714,716,718,720,722,724,726,728,730,732,734,736,738,740,742,744,746,748,750,752,754,756,758,760,762,764,766,768,770,772,774,776,778,780,782,784,786,788,790,792,794,796,798,800,802,804,806,808,810,
                       812,814,816,818,820,822,824,826,828,830,832,834,836,838,840,842,844,846,848,850,852,854,856,858,860,862,864,866,868,870,872,874,876,878,880,882,884,886,888,890,892,894,896,898,900,902,904,906,908,910,912,914,916,918,920,922,924,926,928,930,932,934,936,938,940,942,944,946,948,950,952,954,956,958,960,962,964,966,968,970,972,974,976,978,980,982,984,986,988,990,992,994,996,998,1000,1002,1004,1006,1008,1010,1012,1014,1016,1018,1020,1022,1024,1026,1028,1030,1032,1034,1036,1038,1040,1042,1044,1046,1048,1050,1052,1054,1056,1058,1060,1062,1064,1066,1068,1070,1072,1074,1076,1078,1080,1082,1084,1086,1088,1090,1092,1094,1096,1098,1100,1102,1104,1106,1108,1110,1112,1114,1116,
                       1118,1120,1122,1124,1126,1128,1130,1132,1134,1136,1138,1140,1142,1144,1146,1148,1150,1152,1154,1156,1158,1160,1162,1164,1166,1168,1170,1172,1174,1176,1178,1180,1182,1184,1186,1188,1190,1192,1194,1196,1198), df,  check.names=F)
d <- melt(df,  id.vars = 'a', variable.name = 'series')

d$a <- as.numeric(as.character(d$a))
d$series <- as.factor(d$series)
d$value <- as.numeric(as.character(d$value))

rect2 <- data.frame(xmin=550, xmax=750, ymin=-Inf, ymax=Inf)
ggplot(d, aes(a,value)) +
  geom_line(aes(color = series), size = 2)+
  scale_y_continuous(breaks=seq(-1, 2.5, 1)) +
  coord_cartesian(ylim=c(-1, 2.5))+
  #scale_x_continuous(breaks=seq(-200, 1200, 200)) +
  scale_color_manual(values=c("red", "black", "blue")) +
  # scale_linetype_manual(values=c("solid", "solid", "twodash")) +
  theme(legend.direction = 'vertical', 
        legend.position = 'right',
        legend.key = element_rect(size = 7),
        legend.key.size = unit(3, 'lines')) +
  theme(panel.grid.major = element_blank(), text = element_text(size=30), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.y = element_line(colour = "black"))+
  geom_vline(xintercept=c(0), linetype="dotted", size=1.5)+
  #geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "gray88",
  #          color="gray88",
  #          alpha=0.4,
  #          inherit.aes = FALSE)  +
  geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = "gray88",
            color="gray88",
            alpha=0.4,
            inherit.aes = FALSE)  +
  #  labs(title = "Reward effect in all blocks\nP300 Grand Average Pooled Fz, FCz, F1, F2, FC1, FC2", 
  #      x = "Time from stimulus onset [ms]", y = "Amplitude [µV]", color = "Conditions\n")
  
  labs(title = "", 
       x = "", y = "", color = "") +
  labs(x = "Time [ms]",
       y = "Amplitude [µV]"
  ) +
  theme(axis.text.x=element_text(size=20, vjust=50),
        axis.title.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  geom_hline(yintercept=0,  
             color = "black", size=0.5)



