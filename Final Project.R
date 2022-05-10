#Load libraries
library("sf")
library("ggspatial")
library("tidyverse")
library("dplyr")
library("lubridate")
library("ggmap")
library("ggstatsplot")
library("ggpubr")
library("corrgram")
library("ellipse")
library("mctest")
library("car")
library("ggplot2")
library("lme4")


#Set working directory and load data
setwd("~/MPS/Spring 2022/Data Manipulation Course/Final Project")
data <- read.csv("heat_stress_data.csv")
data
#Select out columns of interest
df1 <- data %>% dplyr::select(puck, AcclimationTankOriginal, TempRegimeTank, OriginalTreatment,
                       TreatmentTank, genotype, genotypeTag, region, experiment,
                       lat, long, daysUntilPulled, reasonPulled, bw1Date, bw1CalcMass..g., bw2Date,
                       bw2CalcMass..g., bw3Date, bw3CalcMass..g.)
df1$Difference <- (df1$bw3CalcMass..g. - df1$bw1CalcMass..g.)
df1$Difference_percentage <- ((df1$Difference/df1$bw1CalcMass..g.)*100)

ggbetweenstats(data=df1,
               x=OriginalTreatment,
               y=Difference_percentage,
               outlier.tagging=TRUE,
               outlier.label=puck)

#puck 4 and 235 seem like outliers so will not include them in the analysis

df1$Difference_percentage[which(df1$puck==4)]=NA
df1$Difference_percentage[which(df1$puck==235)]=NA
df1 <- na.omit(df1)

#Is the data normally distributed?
ggqqplot(df1$Difference_percentage) #looks normal
shapiro.test(df1$Difference_percentage) #Fail to reject null hypothesis so normal distribution

#Do we have homogeneity of variances?
resultBT_Treatment = bartlett.test(Difference_percentage~OriginalTreatment,data=df1)
resultBT_Treatment
#Reject null hypothesis that there is no difference in variances between groups.
summary(df1$OriginalTreatment)
#Although the treatment variable has unequal variances, the sample sizes are roughly equal so two-way anova will work (robust to unequal variances if sample sizes are equal)

resultBT_Genotype = bartlett.test(Difference_percentage~genotypeTag, data=df1)
resultBT_Genotype
#Fail to reject null hypothesis that all population variances are equal

df1$Reef <- sapply(df1$genotype, switch,
                   'Yungs-A' = "Yungs",
                   'Yungs-B' = "Yungs",
                   'BC-1'= "Broward County",
                   'BC-8b' = "Broward County",
                   'Kelseys-1'= "Kelseys",
                   'Kelseys-2'= "Kelseys")

#Exploring the data:
class(df1$genotypeTag)
class(df1$OriginalTreatment)
class(df1$Reef)
class(df1$Difference_percentage)

df1$OriginalTreatment <- as.factor(df1$OriginalTreatment)
df1$genotypeTag <- as.factor(df1$genotypeTag)
df1$Reef <- as.factor(df1$Reef)

#Site mapping
fl_map <-get_map(location=c(-80.5,25.25,-79.5,26.25), source = "google", maptype = "terrain")
fl_map2 <- get_map(location=c(-80.5,25.25,-79.5,26.25), source = "google", maptype = "satellite")
map2 <- ggmap(fl_map2) + 
  geom_point(data=df1, aes(x=long,y=lat, color=Reef, fill=Reef),
             shape =21, size = 2, color = "black", stroke = 0.75) +
  scale_fill_manual(values=c("darkseagreen2","thistle2","blanchedalmond"))+  
  ggtitle("Sample Collection Sites")+
  theme(plot.title=element_text(size=10))+
  xlab("")+
  ylab("")

#Days until fragments were pulled analysis:
df1_summarized <- df1 %>% group_by(genotypeTag, OriginalTreatment) %>% summarize(mean = mean(daysUntilPulled), sd = sd(daysUntilPulled))

ggplot(df1_summarized, aes(x=genotypeTag, y=mean, fill = OriginalTreatment)) + geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3, position=position_dodge(0.9)) +
  scale_fill_manual(values = c("steelblue4","powderblue")) +
  xlab("Coral Genotype") +
  ylab("Number of Days") + 
  ggtitle("Mean days before bleaching or tissue loss")

aov1 <-aov(daysUntilPulled ~ OriginalTreatment + genotypeTag + OriginalTreatment:genotypeTag, data=df1)
anova(aov1) #both significant but interaction is not significant

aov2 <- aov(daysUntilPulled ~ genotypeTag, data = df1)
anova(aov2) #genotype not significant

t.test(daysUntilPulled~OriginalTreatment, data=df1)

#Buoyant weight analysis
df1_differencepercentage_summarized <- df1 %>% group_by(genotypeTag, OriginalTreatment) %>% summarise(mean = mean(Difference_percentage), sd = sd(Difference_percentage))

par(mfrow = c(3, 2))
ggplot(df1, aes(x=OriginalTreatment, y=Difference_percentage)) + 
  geom_boxplot(fill="cornsilk") + 
  facet_wrap(~genotypeTag) +
  xlab("Treatment") +
  ylab("Percent Change in Growth") + 
  ggtitle("Comparing Differences in Buoyant Weight Changes for Treatment and Control")

#Anovas on Growth Difference
aov3 <- aov(Difference_percentage ~ OriginalTreatment + genotypeTag + OriginalTreatment:genotypeTag, data = df1)
anova(aov3) #interaction not significant

#Since the interaction was not significant, can separate out into genotype and treatment group
aov4 <- aov(Difference_percentage ~ genotypeTag, data = df1)
anova(aov4) #significant

post_hoc_aov4 <- TukeyHSD(aov4)
post_hoc_aov4
TK_data <- as.data.frame(post_hoc_aov4[1])
TK_data_significant <- filter(TK_data, genotypeTag.p.adj <= "0.05") 
TK_data_significant #two significant pairwise comparisons

t.test(Difference_percentage~OriginalTreatment, data=df1) #significant

#GLM:
library(visreg)
library(lme4)
library(MuMIn)
library(predictmeans)
library(mosaic)
library(car)
library(influence.ME)
library(DHARMa)
library(tidyverse)

par(mfrow = c(1,1))
qqnorm(df1$Difference)
qqline(df1$Difference)

#variables to include: genotypeTag, and OriginalTreatment, 
#Fixed = treatment assigned (OriginalTreatment, genotypeTag)
#Random = not assigning a treatment (Reef)
#Modeling
model_all <- lmer(Difference_percentage ~ OriginalTreatment + genotypeTag +
                     (1|Reef),
                      data=df1)

model_onlyfixed <- lm(Difference_percentage ~ OriginalTreatment + genotypeTag,
                      data=df1)

model_onlygeno <- lm(Difference_percentage ~ genotypeTag, data=df1)

model_onlytreatment <- lm(Difference_percentage ~ OriginalTreatment, data=df1)

model_originaltreatmentandreef <- lmer(Difference_percentage ~ OriginalTreatment +
                     (1|Reef),
                      data=df1)

model_genotypeTagandreef <- lmer(Difference_percentage ~ genotypeTag +
                     (1|Reef),
                      data=df1)

#Comparing AICs of models
anova(model_all, model_onlyfixed) #no difference
anova(model_all, model_onlygeno) #model all is better than only geno
anova(model_all, model_onlytreatment) #model all is better than only genotype
anova(model_all, model_originaltreatmentandreef) #model all is better
anova(model_all, model_genotypeTagandreef) #model all is better
anova(model_onlyfixed, model_onlygeno)
anova(model_onlyfixed, model_onlytreatment)

#Shows that there is no significant difference from including Reef as a random effect so best model is one that includes Original Treatment and Genotype tag

