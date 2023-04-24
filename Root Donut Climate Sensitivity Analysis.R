##This script is the final code for analysis of root donut data for 
## Vojdani et al. manuscript -- code authors: J. Rudgers, L. Baur, A. Vojdani
## Last updated 22 April 2023
rm(list=ls(all=TRUE)) #give R a blank slate

library(tidyverse)
library(nlme)
library(lme4)
library(car)
library(emmeans)
library(piecewiseSEM)
library(MuMIn)
library(visreg)
library(plotrix)
library(wiqid)

#set working directory
# examples ...
# JENN's COMPUTER
# setwd("C:/Users/jrcassvc/Desktop/Sev Data/Root donuts/")
# AZAD's COMPUTER
#setwd("C:/Users/vojda/OneDrive/Desktop/Root_Biomass_MS")
# LAUREN'S COMPUTER
#setwd("C:/Users/laure/Desktop/work temp/CLEAN DATA/Root donuts/Azads project/results_20230104")

## DATA STEPS ----------------------------
#import data from computer
roots_spei<-read.csv("roots_spei_20220415.csv", stringsAsFactors = T)
#1956 observations of 19 variables

#add treatment information
unique(roots_spei$Site)
roots_spei$treatment[roots_spei$Site=="F" & (roots_spei$Sample==2|roots_spei$Sample==4|roots_spei$Sample==6|roots_spei$Sample==8|roots_spei$Sample==10|roots_spei$Sample==11|roots_spei$Sample==14|roots_spei$Sample==16|roots_spei$Sample==17|roots_spei$Sample==19)]<-"fertilized"
roots_spei$treatment[roots_spei$Site=="M" & (roots_spei$Sample==1|roots_spei$Sample==3|roots_spei$Sample==5|roots_spei$Sample==6|roots_spei$Sample==9)]<-"small_addition"
roots_spei$treatment[roots_spei$Site=="M" & (roots_spei$Sample==2|roots_spei$Sample==4|roots_spei$Sample==7|roots_spei$Sample==10|roots_spei$Sample==12)]<-"large_addition"
roots_spei$treatment[roots_spei$Site=="DWB"]<-"burned"
roots_spei$treatment[is.na(roots_spei$treatment)]<-"control"
roots_spei$treatment<-as.factor(roots_spei$treatment)
summary(roots_spei$treatment)

#revise depth values and remove any data missing depth info
roots_spei$Depth2[roots_spei$Depth=="0-15"]<-1
roots_spei$Depth2[roots_spei$Depth=="15-30"]<-2
roots_spei<-subset(roots_spei, !is.na(Depth2)) # removes 1 obs
roots_spei$Depth2<-as.factor(roots_spei$Depth2)
roots_spei$Depth.f<-recode_factor(roots_spei$Depth,
                                     '0-15' = "Soil depth 0-15 cm",
                                     '15-30' = "Soil depth 15-30 cm")

#create a unique id for each root donut location for repeated measures analysis
roots_spei$sitesample<-as.factor(paste(roots_spei$Site, roots_spei$Sample, sep="_"))
summary(roots_spei$sitesample)
#create a unique id for each root donut location and depth for temporal autocorrelation
roots_spei$sitesampledpth<-as.factor(paste(roots_spei$Site, roots_spei$Sample,  roots_spei$Depth2, sep="_"))
summary(roots_spei$sitesampledepth)
#create year as a factor for analysis to include random effect of year
roots_spei$year.f<-as.factor(roots_spei$Date)

#exclude TWO very high outliers for root mass > 6, next highest value was ~4
qqnorm(roots_spei$wt_by_vol)
summary(roots_spei$wt_by_vol)
roots_spei<-subset(roots_spei,wt_by_vol<6)
# 1954 observations - removed 1 extreme outlier

##Removing additional outliers within the Fertilizer experiment
summary(roots_spei$Root_Weight)
subset(roots_spei, Root_Weight>4500 & Site=="F")
subset(roots_spei, Root_Weight>4500)
roots_spei<-subset(roots_spei, Root_Weight<=4500|Site!="F")
summary(roots_spei)
# removed an additional 4 extreme outliers, 1950 observations
 
# look at data distribution
hist(roots_spei$wt_by_vol)
max(roots_spei$wt_by_vol)
min(roots_spei$wt_by_vol)

# convert root mass to g per m3
# Upscale wt_by_vol by .002 number, converts cm^3 to m^3
# outer volume (cm3) 4854.82 (20.3 diameter, 15 cm height)
# inner volume (cm3) 2721.88 (15.2 diameter, 15 cm height)
# soil volume (cm3)  2132.94  (subtract inner from outer to get the ring volume)
# conversion cm3 to m2  2132.94/1000000 = 0.00213294 m2
# upscale our wt_by_vol to m2, you need to divide by 0.00213294
roots_spei$weight_volume<-(roots_spei$wt_by_vol/0.00213294)
summary(roots_spei$weight_volume)
# try log transformation
hist(roots_spei$weight_volume)
hist(log(roots_spei$weight_volume+10))
qqnorm(log(roots_spei$weight_volume+10))
#note: adding 10 as the constant gives the best normality plot

#look at levels of sites and treatments
summary(roots_spei$treatment)
summary(roots_spei$Site)

# METADATA Quick KEY ---------------
#C   Creosotebush Desert Shrubland Core Site
#DW  Deep Well Core site - Unburned
#DWB Deep Well Core site - Burned
#F   Nitrogen Fertilization site near DW - includes N addition plots
#M   Monsoon Rainfall Manipulation Desert Grassland (black grama grass) Site - includes rainfall treatments

####   Question (1) Compare ecosystems for root biomass vs. aridity index: Control plots ############################################################
## Question (1) Do dryland ecosystem types differ in the sensitivity of root biomass to the mean or interannual variance in aridity?
# Answer: Yes, shrubland root biomass was more sensitive to SPEI aridity than grasslands
#subset the data to control root donuts only
root_controls<-subset(roots_spei,treatment=="control")
summary(root_controls$Site)
summary(root_controls$Site:root_controls$year.f)

## model selection with log-transformed data
# L = linear, Q = quadratic, C = cubic
#data are insufficient to get convergence for a random slopes model (random =~ Date|sitesampledepth)
mL<-lme(log(weight_volume+10)~SPEI_12_Sept*Site*Depth2,data=root_controls,random=~1|sitesample/sitesampledepth,method="ML", na.action=na.omit)
mLAR1<-lme(log(weight_volume+10)~SPEI_12_Sept*Site*Depth2,data=root_controls,random=~1|sitesample/sitesampledepth,correlation=corAR1(form=~Date|sitesample/sitesampledepth),method="ML", na.action=na.omit)
mLAR2<-lme(log(weight_volume+10)~SPEI_12_Sept*Site*Depth2,data=root_controls,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2,q=0),method="ML", na.action=na.omit)
mQ<-update(mL, .~. +I(SPEI_12_Sept^2)*Site*Depth2)
mQAR1<-update(mLAR1, .~. +I(SPEI_12_Sept^2)*Site*Depth2)
mQAR2<-update(mLAR2, .~. +I(SPEI_12_Sept^2)*Site*Depth2)
mC<-update(mQ, .~. +I(SPEI_12_Sept^3)*Site*Depth2)
mCAR1<-update(mQAR1, .~. +I(SPEI_12_Sept^3)*Site*Depth2)
mCAR2<-update(mQAR2, .~. +I(SPEI_12_Sept^3)*Site*Depth2)

#compare models to pick the best one
AICcTable<-data.frame(MuMIn::AICc(mL,mLAR1,mLAR2,
                      mQ,mQAR1,mQAR2,
                      mC,mCAR1,mCAR2))

AICcTable[order(AICcTable$AICc),]

## Best model is mCAR2 -- better than mQAR2 by 5 AICc
## Test is cubic term significant?
## Results from best model
anova(mCAR2, type="marginal")
table<-data.frame(anova(mCAR2, type="marginal"))
table[order(table$p.value),]
## Significant: Depth, Site, Site*Depth, SPEI, SPEI^2.
## SPEI^3 is not significant nor are any interactions with it
## So cubic model is not really better, choose the more conservative Q model
## mQAR2 is next best model
anova(mQAR2, type="marginal")
table<-data.frame(anova(mQAR2, type="marginal"))
# write.csv(table,"statstable.site.roots.csv")
## Significant: SPEI, Depth, Site, SPEI*Site, Site*Depth, SPEI^2
## This means sites have different climate sensitivities
## CSFs - linear SPEI*Site, but similar CSF quadratic (SPEI^2*Site is n.s.)
## and root biomass changes with depth differently among sites, Site*Depth
## Check model assumptions
hist(resid(mQAR2))
qqnorm(resid(mQAR2)) #log fixed the normality of residuals assumption
plot(mQAR2) #homogeneity of variances assumption
#everything looks good for model assumptions

## Which Sites differed in parameters of the relationship with aridity?
## First, take the best model and standardize the SPEI_12_Sept variable using the poly function
m.poly.controls<- lme(weight_volume~poly(SPEI_12_Sept,degree=2)*Site*Depth.f,data=root_controls,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2,q=0),method="ML", na.action=na.omit)
#next get slopes and quadratic terms estimated from the poly model
#which has only one continuous predictor, a poly function of SPEI
em_site<-data.frame(emtrends(m.poly.controls, ~ Site|degree, var="SPEI_12_Sept", max.degree = 2, options = list()))
#note: 'degree' provides estimates of linear and nonlinear
em_site
#model R2
rsquared(m.poly.controls)
#write out parameters
# write.csv(em_site,"parameters.site.roots.csv")
# next contrast pairs of sites
c_site<-(emtrends(m.poly.controls, pairwise ~ Site|degree, var="SPEI_12_Sept", max.degree = 2, options = list()))
##"Note: Use 'contrast(regrid(object), ...)' to obtain contrasts of back-transformed estimates"
c_site
## RESULTS
## All sites have overall positive linear relationship 
## between root biomass and aridity - less arid, more roots. (Previously--data through 2018--all were negative except creosote)
## All sites have negative quadratic term
## CONTRASTS: Creosotebush significantly differed from DW and F (linear) with a stronger positive relationship with aridity
## CONTRASTS: No significant contrasts for quadratic term, 
## which is negative= concave downward and significantly so in all sites but M=MRME

## Controls only - Graph across sites without depth  ---------------------------------
## use recode_factor to make the values of Site more meaningful for the graph
root_controls$ecosystem<-recode_factor(root_controls$Site,
                                       C = "Desert shrubland",
                                       M = "Desert grassland",
                                       DW = "Mixed grassland",
                                       'F' = "Plains grassland")
#determine graph axis limits
#xaxis
min(root_controls$SPEI_12_Sept)
max(root_controls$SPEI_12_Sept)
#yaxis
min(log(root_controls$weight_volume+10))
max(log(root_controls$weight_volume+10))

# create the graph
ggplot(root_controls,aes(x=SPEI_12_Sept,y=log(weight_volume+10),colour=ecosystem))+
  scale_color_manual(values=c("olivedrab4","grey30","khaki3","steelblue3"))+
  facet_grid(rows = vars(ecosystem))+ #scales="free_y")+
  geom_point()+
  geom_smooth(formula=y~poly(x,2),method="lm")+
  xlab("SPEI aridity index")+
  ylab(bquote('ln (Belowground net primary production) (g '*~m^-3~y^-1*')'))+
  theme_minimal()+
  ylim(0,7.25)+
  xlim(-2,2)+
  theme(legend.position = "none")
#export graph
ggsave("root_controls.tiff",dpi=500,device=NULL,width=4,height=10)

## Effect size by depth --------------
# How do depths differ in root biomass?
em_depth<-data.frame(emmeans(m.poly.controls,~Depth.f))
em_depth
# There is more root biomass deeper, but how much more?
# Get treatment means in measurement scale (g/m3)
exp(summary(emmeans(m.poly.controls,~Depth.f))[c("emmean")])-10
#back transform the confidence limits so they are asymmetrical in original units
exp(summary(emmeans(m.poly.controls,~Depth.f))[c("lower.CL")])-10
exp(summary(emmeans(m.poly.controls,~Depth.f))[c("upper.CL")])-10
#*** These are the m_raw numbers
(130.6-82.3)/82.3*100
##59% more biomass on average at deeper depth across all sites combined

## Effect sizes for DEPTH : SITE ------------------
# How do depths differ in root biomass among ecosystems?
em_ds<-data.frame(emmeans(m.poly.controls,~Depth.f|Site))
em_ds
pairs(emmeans(m.poly.controls,~Depth.f|Site)) # Significantly different at C, F and M
# more root biomass deeper, but how much more?
#get treatment means in measurement scale (g/L)
emmeans(m.poly.controls,~Depth.f|Site)
m_raw_ds<-exp(summary(emmeans(m.poly.controls,~Depth.f|Site))[c("emmean")])-10
m_raw_ds$site_depth<-c("C0-15","C15-30","DW0-15","DW15-30","F0-15","F15-30","M0-15","M15-30")
# back transform the confidence limits so they are asymmetrical in original units
m_raw_ds$lower.CL<-exp(summary(emmeans(m.poly.controls,~Depth.f|Site))[c("lower.CL")])-10
m_raw_ds$upper.CL<-exp(summary(emmeans(m.poly.controls,~Depth.f|Site))[c("upper.CL")])-10
m_raw_ds
# Creosotebush Desert shrubland had 131% more root biomass at deeper depth
(103.2-44.6)/44.6*100
# Mixed grassland ecosystem had 12% LESS root biomass at deeper depth, n.s.
(100.2-113.6)/113.6*100
# Plains grassland ecosystem had 37% more roots at deeper depth
# interesting that the grasslands differ in root distribution!
(134.0-97.5)/97.5*100
# Desert black grama grassland ecosystem had 131% more roots at deeper depth
(207.5-89.9)/89.9*100
#write.csv(m_raw_ds,"depth_means_by_site.csv")

##  Supplement:   Time series analysis for all ecosystems ------
# Results are in supplement, not main text
mtime<-lme(log(weight_volume+10)~Date*Site*Depth2,data=root_controls,random=~1|sitesampledepth,correlation=corAR1(form=~Date|sitesampledepth),method="ML", na.action=na.omit)
#results
anova(mtime, type="marginal")
#temporal patterns differed among the sites: Date:Site P=0.0100
#check model assumptions
hist(resid(mtime))
qqnorm(resid(mtime)) # log transformation was required
plot(mtime)
#compare temporal trends among sites
emtrends(mtime, "Site", var = "Date")
# get average per site per year
root_controls_by<-root_controls%>% group_by(Date,ecosystem,Depth.f)
root_controls_mean<-data.frame(root_controls_by %>% summarise_at(vars(weight_volume),funs(mean,std.error)))
root_controls_mean
#make a time plot
ggplot(root_controls_mean,aes(x=Date,y=mean,colour=ecosystem))+
  scale_color_manual(values=c("olivedrab4","grey30","khaki3","steelblue3"))+
  facet_grid(rows = vars(ecosystem),cols = vars(Depth.f))+ #scales="free_y")+
  geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error,colour=ecosystem),width=0)+
  geom_line(aes(colour=ecosystem),size=0.5,na.rm = TRUE)+
  geom_point()+
  xlab("Time")+
  ylab(bquote('Belowground net primary production (g '*~m^-3~y^-1*')'))+
  theme_minimal()+
  theme(legend.position = "none")
ggsave("time_root_controls.tiff",dpi=300,device=NULL,width=10,height=10)

##### Question (2)  Monsoon Rainfall Manipulation Experiment -----------------
## Question (2): Does the intraannual rainfall regime alter the sensitivity 
##  of root biomass to the mean or variance in interannual aridity?
## Answer: No. (Only SPEI and depth are significant.)
#subset the data to Monsoon Rainfall Manipulation Experiment  site only
root_mrme<-subset(roots_spei,Site=="M")
##model selection with log-transformed data
#L = linear, Q = quadratic, C = cubic
#data are insufficient to get convergence for a random slopes model (random =~ Date|sitesampledepth)
mrmeL<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_mrme,random=~1|sitesample/sitesampledepth,method="ML", na.action=na.omit)
mrmeLAR1<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_mrme,random=~1|sitesample/sitesampledepth,correlation=corAR1(form=~Date|sitesample/sitesampledepth),method="ML", na.action=na.omit)
mrmeLAR2<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_mrme,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2),method="ML", na.action=na.omit)
mrmeQ<-update(mrmeL, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
mrmeQAR1<-update(mrmeLAR1, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
mrmeQAR2<-update(mrmeLAR2, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
mrmeC<-update(mrmeQ, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
mrmeCAR1<-update(mrmeQAR1, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
mrmeCAR2<-update(mrmeQAR2, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
#compare models to pick the best one
AICcTable_mrme<-data.frame(MuMIn::AICc(mrmeL,mrmeLAR1,mrmeLAR2,
                           mrmeQ,mrmeQAR1,mrmeQAR2,
                           mrmeC,mrmeCAR1,mrmeCAR2))
AICcTable_mrme[order(AICcTable_mrme$AICc),]
#results from best model which was 6.5 AICc better than linear
(623.0223-616.4840)
anova(mrmeQAR2, type="marginal")
## SEPI^2 is not significant, but next best model is mrmeLAR2, not very close
# and also does not match results for Controls from Q1.
anova(mrmeLAR2, type="marginal")
table<-data.frame(anova(mrmeQAR2, type="marginal"))
## significant: SPEI, Depth, no significant interactions
#write.csv(table, "mrme_anova_table.csv")
#check model assumptions
hist(resid(mrmeQAR2))
qqnorm(resid(mrmeQAR2)) 
plot(mrmeQAR2)
#everything looks good for model assumptions

#which treatments, if any, differ in parameters of the relationship with aridity?
#main stats results indicate no differences, look at parameter estimates for each treatment too
mrme.poly<- lme(log(weight_volume+10)~poly(SPEI_12_Sept,degree=2)*treatment*Depth2,data=root_mrme,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2),method="ML", na.action=na.omit)
em_mrme<-data.frame(emtrends(mrme.poly, ~degree|treatment, var="SPEI_12_Sept", max.degree = 2, options = list()))
em_mrme
#write.csv(em_mrme,"parameters.mrme.roots.csv")
#next contrast pairs of mrme treatment levels
c_mrme<-(emtrends(mrmeQAR2, pairwise ~ treatment|degree, var="SPEI_12_Sept", max.degree = 2, options = list()))
c_mrme #confirms no treatment differences in root ~ SPEI
# interpretation: MRME treatments do not affect the relationship of root biomass with SPEI (linear/quadratic estimates are all positive)
## Depth X Treatment trends (interaction was n.s.)
em_mrme.depth<-data.frame(emtrends(mrmeQAR2, ~ Depth2|treatment, var="SPEI_12_Sept", options = list()))
em_mrme.depth
## Treatment effects on root biomass
em_mrme.trt<-data.frame(emmeans(mrmeQAR2, ~ treatment, options = list()))
em_mrme.trt
emmeans(mrmeQAR2, pairwise ~ treatment, options = list())
# Graph results
#give treatments better names and reorder them
root_mrme$treatment<-recode_factor(root_mrme$treatment,
                                   control="Control",
                                   small_addition="Many small rains",
                                   large_addition="Few large rains")
# plot raw data, log transformed                                   
ggplot(root_mrme,aes(x=SPEI_12_Sept,y=log(weight_volume+10),colour=treatment))+
   scale_color_manual(values=c("black","slategray3","slateblue1"))+
   facet_grid(rows = vars(treatment))+
   geom_point()+
   geom_smooth(formula=y~poly(x,2),method="lm")+
   xlab("SPEI aridity index")+
   ylab(bquote('ln (Belowground net primary production) (g '*~m^-3~y^-1*')'))+
   theme_minimal()+
  ylim(0,7.5)+
  theme(legend.position = "none")
ggsave("root_mrme.tiff",dpi=300,device=NULL,width=3,height=10)

####  (3)  Nitrogen Fertilization ############################################################
##  Does nitrogen addition alter the sensitivity 
##  of root biomass to the mean or variance in interannual aridity?
##No (only SPEI and depth are significant)
#subset the data to NFert site only
root_fert<-subset(roots_spei,Site=="F")

##model selection with log-transformed data
#L = linear, Q = quadratic, C = cubic
fertL<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_fert,random=~1|sitesample/sitesampledepth,method="ML", na.action=na.omit)
fertLAR1<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_fert,random=~1|sitesample/sitesampledepth,correlation=corAR1(form=~Date|sitesample/sitesampledepth),method="ML", na.action=na.omit)
fertLAR2<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_fert,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2),method="ML", na.action=na.omit)
fertQ<-update(fertL, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
fertQAR1<-update(fertLAR1, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
fertQAR2<-update(fertLAR2, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
fertC<-update(fertQ, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
fertCAR1<-update(fertQAR1, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
fertCAR2<-update(fertQAR2, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)

#compare models to pick the best one
AICcTable_fert<-data.frame(MuMIn::AICc(fertL,fertLAR1,fertLAR2,
                                fertQ,fertQAR1,fertQAR2,
                                fertC,fertCAR1,fertCAR2))
AICcTable_fert[order(AICcTable_fert$AICc),]
#results from best model 
##fertQAR2 is best model, beating fertQAR2 by >3 AICc
anova(fertQAR2, type="marginal")
table<-data.frame(anova(fertQAR2, type="marginal"))
table[order(table$p.value),]
#write.csv(table, "fert_anova_table.csv")
#Significant:SPEI^2, Depth
#N fertilization did not affect root biomass
#check model assumptions
hist(resid(fertQAR2))
qqnorm(resid(fertQAR2)) #try log
plot(fertQAR2)
#first, take the best model and standardize the SPEI_12_Sept variable using poly function
fert.poly<-lme(log(weight_volume+10)~poly(SPEI_12_Sept,degree=2)*treatment*Depth2,data=root_fert,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2),method="ML", na.action=na.omit)
#next get slopes, quadratic, and/or cubic terms estimated from the poly model,
#which has only one continuous predictor, a poly function of SPEI
em_fert<-data.frame(emtrends(fert.poly, ~ treatment|degree, var="SPEI_12_Sept", max.degree = 2, options = list()))
em_fert
#write.csv(em_fert,"parameters.fert.roots.csv")
#next contrast N addition against control - do relationships with SPEI differ?
c_fert<-(emtrends(fert.poly, pairwise ~ treatment|degree, var="SPEI_12_Sept", max.degree = 2, options = list()))
c_fert 
#results: fert treatments do not affect the relationship with SPEI
##depth trends
emmeans(fert.poly, ~ Depth2, options = list())
exp(summary(emmeans(fert.poly, ~ Depth2))[c("emmean")])-10
exp(summary(emmeans(fert.poly, ~ Depth2))[c("lower.CL")])-10
exp(summary(emmeans(fert.poly, ~ Depth2))[c("upper.CL")])-10
#depth effect size
(128.4696-102.8908)/128.4696*100
#give treatments better names and reorder them
root_fert$treatment<-recode_factor(root_fert$treatment,
                                   control="Control",
                                   fertilized= "Fertilized")

ggplot(root_fert,aes(x=SPEI_12_Sept,y=log(weight_volume+10),colour=treatment))+
   scale_color_manual(values=c("black","gray50"))+
   facet_grid(rows = vars(treatment))+ #,cols = vars(Depth.f)
   geom_point()+
   geom_smooth(formula=y~poly(x,2),method="lm")+
   xlab("SPEI aridity index")+
  ylab(bquote('ln (Belowground net primary production) (g '*~m^-3~y^-1*')'))+
   theme_minimal()+
   ylim(0,7.5)+
  xlim(-2,2)+
   theme(legend.position = "none")
#export graph
 ggsave("root_fert.tiff",dpi=300,device=NULL,width=3,height=10)

#### Question  (4)  Fire Experiment --------------
##  Does disturbance by fire alter the sensitivity 
##  of root biomass to the mean or variance in interannual aridity?
##  Ye, SPEI*treatment, treatment*depth, treatment all significant
#subset the data to DWB and DW sites only
root_burn<-subset(roots_spei,Site=="DW"|Site=="DWB")
##model selection with log-transformed data
#L = linear, Q = quadratic, C = cubic
#data are insufficient to get convergence for a random slopes model (random =~ Date|sitesample/sitesampledepth)
burnL<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_burn,random=~1|sitesample/sitesampledepth,method="ML", na.action=na.omit)
burnLAR1<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_burn,random=~1|sitesample/sitesampledepth,correlation=corAR1(form=~Date|sitesample/sitesampledepth),method="ML", na.action=na.omit)
burnLAR2<-lme(log(weight_volume+10)~SPEI_12_Sept*treatment*Depth2,data=root_burn,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2),method="ML", na.action=na.omit)
burnQ<-update(burnL, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
burnQAR1<-update(burnLAR1, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
burnQAR2<-update(burnLAR2, .~. +I(SPEI_12_Sept^2)*treatment*Depth2)
burnC<-update(burnQ, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
burnCAR1<-update(burnQAR1, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
burnCAR2<-update(burnQAR2, .~. +I(SPEI_12_Sept^3)*treatment*Depth2)
#compare models to pick the best one
AICcTable_burn<-data.frame(MuMIn::AICc(burnL,burnLAR1,burnLAR2,
                                burnQ,burnQAR1,burnQAR2,
                                burnC,burnCAR1,burnCAR2))
AICcTable_burn[order(AICcTable_burn$AICc),]
#burnQAR2 is best model, beating burnCAR2 by 1 point, and the best linear model, burnLAR1, by 51 points
anova(burnQAR2, type="marginal")
table<-data.frame(anova(burnQAR2, type="marginal"))
table[order(table$p.value),]
#write.csv(table, "burn_anova_table.csv")
#check model assumptions
hist(resid(burnQAR2))
qqnorm(resid(burnQAR2)) #try log
plot(burnQAR2)
#everything looks ok for model assumptions
#the low residual variation at low predicteds isn't great - could investigate
#take the best model and standardize the SPEI_12_Sept variable using poly function
burn.poly<- lme(log(weight_volume+10)~poly(SPEI_12_Sept,degree=2)*treatment*Depth.f,data=root_burn,random=~1|sitesample/sitesampledepth,correlation=corARMA(form=~Date|sitesample/sitesampledepth,p=2),method="ML", na.action=na.omit)
#next get slopes, quadratic, and cubic terms estimated from the poly model,
#which has only one continuous predictor, a poly function of SPEI
em_burn<-data.frame(emtrends(burn.poly, ~ treatment|degree, var="SPEI_12_Sept", max.degree = 2, options = list()))
em_burn
#write.csv(em_burn,"parameters.burn.roots.csv")
#effect size - 116% stronger linear sensitivity to aridity in burned than control
(0.41-0.19)/0.19*100
#next contrast pairs of burns
c_burn<-(emtrends(burn.poly, pairwise ~ treatment|degree, var="SPEI_12_Sept", max.degree = 2, options = list()))
c_burn 
##how much did treatment effects vary with depth?
emmeans(burn.poly, ~ treatment|Depth.f)
emmeans(burn.poly, pairwise ~ treatment|Depth.f)
#graph raw data with a cubic function
#give treatments better names and reorder them
root_burn$treatment<-recode_factor(root_burn$treatment,
                                   control="Control",
                                   burned="Burned")
ggplot(root_burn,aes(x=SPEI_12_Sept,y=log(weight_volume+10),colour=treatment))+
   scale_color_manual(values=c("black","orangered2"))+
   facet_grid(rows = vars(treatment))+ #,cols = vars(Depth.f)
   geom_point()+
   geom_smooth(formula=y~poly(x,2),method="lm")+
   xlab("SPEI aridity index")+
  ylab(bquote('ln (Belowground net primary production) (g '*~m^-3~y^-1*')'))+
   theme_minimal()+
  ylim(0,7.5)+
  xlim(-2,2)+
   theme(legend.position = "none")
# #export graph
ggsave("root_burn.tiff",dpi=300,device=NULL,width=3,height=10)
