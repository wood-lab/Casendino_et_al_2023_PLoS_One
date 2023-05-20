###############
# Mixed models of sea star catch (within susceptibility categories) from 1997 to 2019. 
# Includes code for Figures 2, 3 and 4 in Casendino et al. (2023)
# Also includes code testing for serial depletion over successive trawls.  

#####======DEPENDENCIES=======
library(readr)
library(tidyverse)
library(dplyr)
library(lmerTest)
library(lme4)
library(janitor)
library(ggeffects)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(gtable)
library(grid)
library(DHARMa) 
# This is a package for assessing goodness of fit of generalized linear models
library(mgcv) # package for fitting generalized additive models
library(glmmTMB) #using glmmTMB because gmler.nb wouldn't converge
library(here)
library(stringi)
library(MuMIn)
library(magrittr)

######======Loading in data, getting it formatted to run mods=========

source("src/data_processing_function.R")
dat_all<-process_dat() %>% 
  mutate(est_depth=as.factor(est_depth))

  # make year a continuous variable (first year is 0, second year is 1, etc)
  dat_all$year <- as.numeric(dat_all$year)
  dat_all<- dat_all %>% mutate("newcol" = year) %>% 
    mutate(year = year - 1997) %>% 
    rename("year_cont" = year) %>%
    rename("year" = newcol) 
  
dat_all$year <- as.character(dat_all$year)     
dat_all %<>% mutate(total_catch_SM_UM = sum(total_catch_SM + total_catch_UM)) 
dat_all %<>% mutate(total_density_SM_UM = sum(total_density_SM + total_density_UM)) 
# creating extra column combining UM ("likely affected" species) 
# and SM (species with "some mortality"). This column represents moderate-susceptibility spp
# HM columns represent high susceptibility spp

# Function to plot observed and predicted density for each tested model
plot_fun<-function(mod,plot_rand=NULL){
  #generate model predictions
  pred<-predict(mod,se=TRUE, re.form=plot_rand)
  #merge model predictions (log scale) and standard errors with observation data
  new_dat<-dat_all %>% left_join(mod$frame %>% mutate(fit=pred$fit,se.fit=pred$se.fit) 
  %>% rename(test=1)) %>% 
  # calculate confidence intervals and fit on the normal scale
  mutate(pred_d=exp(fit),l95=exp(fit-1.96*se.fit),u95=exp(fit+1.96*se.fit)) %>% 
    group_by(est_depth,year) %>% summarize(Observed=mean(test), Predicted=mean(pred_d),u95=mean(u95),l95=mean(l95)) %>% 
    arrange(year) %>%   pivot_longer(c(Observed,Predicted)) %>% 
    mutate(est_depth=paste(est_depth,"m"))

  ggplot(new_dat, aes(x=as.numeric(year),y=value,color=name)) +
    facet_wrap(~est_depth)+geom_line()+geom_point()+ 
    scale_color_manual("",values=c('black','blue'))+ 
    geom_ribbon(data = new_dat, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.2,.2,.7,.4),color="NA")+
    guides(fill=guide_legend(title="New Legend Title"))+xlab("")+
    ylab(expression("Catch / 1000 m"^2))+ 
    theme(text = element_text(size = 14),panel.spacing = unit(1, "lines"))  
}

######======Models of total catch, moderate susceptibility group=========

# Depth
m1<-glmmTMB(total_catch_SM_UM~est_depth+
                ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
                family="nbinom2",data=dat_all)
summary(m1)
plot_fun(m1,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m1,plot=T) 

# Depth + PrePost (i.e., PrePost effect the same across depths)
m2<-glmmTMB(total_catch_SM_UM~est_depth +
              prepost+
              ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m2)
plot_fun(m2,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m2,plot=T) 

# Depth + PrePost * Depth
m3<-glmmTMB(total_catch_SM_UM~ est_depth+
              prepost:est_depth+
              ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)), 
              family="nbinom2",data=dat_all)
summary(m3)
plot_fun(m3,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m3,plot=T) 

# Depth + Year (i.e., year effect the same across depths)
m4<-glmmTMB(total_catch_SM_UM~est_depth + 
                              year_cont+
                              ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
                              family="nbinom2",data=dat_all)
summary(m4)
plot_fun(m4,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m4,plot=T) 

# Depth + Year * Depth 
m5<-glmmTMB(total_catch_SM_UM~est_depth + 
                              year_cont:est_depth+
                              ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
                                family="nbinom2",data=dat_all)
summary(m5)
plot_fun(m5,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m5,plot=T) 

# Depth + PrePost + Year (prepost and year effects common across depths)
m6<-glmmTMB(total_catch_SM_UM~est_depth+
                year_cont+
                prepost+
                ar1(year+0|est_depth)+
                offset(log(trawlarea/1e3)),family="nbinom2",data=dat_all)
summary(m6)
plot_fun(m6,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m6,plot=T) 

# Depth + PrePost * Year * Depth
m7<-glmmTMB(total_catch_SM_UM~ est_depth+ 
              year_cont:est_depth:prepost+
              ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m7)
plot_fun(m7,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m7,plot=T) 

# Depth + PrePost * Year (i.e., unique slope for pre and post, but same across depths)
m8<-glmmTMB(total_catch_SM_UM~est_depth + 
                  year_cont:prepost+
                ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)),
                  family="nbinom2",data=dat_all)
summary(m8)
plot_fun(m8,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m8,plot=T) 

# Depth + PrePost * Depth + Year * Depth 
m9<-glmmTMB(total_catch_SM_UM~ est_depth + 
              year_cont:est_depth+
              prepost:est_depth+
              ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m9)
plot_fun(m9,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m9,plot=T) 


# Model selection 
AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9) #m5 is ranked best

#add temp to best model fit
m.temp<-glmmTMB(total_catch_SM_UM~est_depth +
                            meanTemp + 
                            year_cont:est_depth+
                            ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
                            family="nbinom2",data=dat_all)
summary(m.temp) # temp is significant positive correlation w/catch
plot_fun(m.temp,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m.temp,plot=T) 

######======plotting chosen model of moderate susceptibility catch (Fig 4)=========

pred<-predict(m5,se=TRUE, re.form=~0)
#merge model predictions (log scale) and standard errors with observation data

new_dat<-dat_all %>% left_join(m5$frame %>% mutate(fit=pred$fit,se.fit=pred$se.fit)) %>% 
  # calculate confidence intervals and fit on the normal scale
  mutate(pred_d=exp(fit),l95=exp(fit-1.96*se.fit),u95=exp(fit+1.96*se.fit)) %>%
  group_by(est_depth,year) %>%
  summarize(Observed=mean(total_density_SM_UM)*1000, Predicted=mean(pred_d),u95=mean(u95),l95=mean(l95)) %>%
  arrange(year) %>%   pivot_longer(c(Observed,Predicted)) %>% 
  mutate(est_depth=paste(est_depth,"m"))

# Because this is best fit model, let's get estimates + CI for catch at each depth in 1997
new_dat %>% filter(year == 1997 & name == "Predicted") %>% group_by(est_depth) %>% summarise(mean(value))
new_dat %>% filter(year == 1997 & name == "Predicted") %>% group_by(est_depth) %>% summarise(mean(l95))
new_dat %>% filter(year == 1997 & name == "Predicted") %>% group_by(est_depth) %>% summarise(mean(u95))

#10 
new_dat_10 <- new_dat %>% filter(est_depth == "10 m")
modelPlot_10<- ggplot(new_dat_10, aes(x=as.numeric(year),y=value,color=name)) + 
  geom_line() + 
  geom_point(data = new_dat_10[new_dat_10$name == "Observed",c(2,5,6)] ) + 
  scale_color_manual("",values=c('black','darkslateblue'), guide = guide_legend(override.aes = list(
    linetype = c("solid", "solid"), shape=c(19, NA)))) + 
  geom_ribbon(data = new_dat_10, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.1,.3,.8,.4),color="NA")+ 
  xlab("Year")+ylab(expression("Catch / 1000 m"^2)) +  
  theme(text = element_text(size = 13))  +theme_minimal() + 
  coord_cartesian(ylim = c(0,9))+ 
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015)) + 
  annotate("text", x=2015, y=7, label= "10 m", size=5, col="black") + 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.text=element_text(size=13)) 
modelPlot_10

#25
new_dat_25 <- new_dat %>% filter(est_depth == "25 m")
modelPlot_25<- ggplot(new_dat_25, aes(x=as.numeric(year),y=value,color=name))+ 
  geom_line()+geom_point(data = new_dat_25[new_dat_25$name == "Observed",c(2,5,6)] )+
  scale_color_manual("",values=c('black','darkslateblue'))+ 
  geom_ribbon(data = new_dat_25, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.1,.3,.8,.4),color="NA")+ 
  xlab("Year")+ylab(expression("Catch / 1000 m"^2))+  
 theme_minimal()+ coord_cartesian(ylim = c(0,9))+
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015))+
  annotate("text", x=2015, y=7, label= "25 m", size=5, col="black") + 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

#50
new_dat_50 <- new_dat %>% filter(est_depth == "50 m")
modelPlot_50<- ggplot(new_dat_50, aes(x=as.numeric(year),y=value,color=name))+ 
  geom_line()+geom_point(data = new_dat_50[new_dat_50$name == "Observed",c(2,5,6)] )+
  scale_color_manual("",values=c('black','darkslateblue'))+ 
  geom_ribbon(data = new_dat_50, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.1,.3,.8,.4),color="NA")+
  xlab("Year")+ylab(expression("Catch / 1000 m"^2))+ 
 theme_minimal()+ coord_cartesian(ylim = c(0,9)) +
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015)) + 
  annotate("text", x=2015, y=7, label= "50 m", size=5, col="black") +
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

#70
new_dat_70 <- new_dat %>% filter(est_depth == "70 m")
modelPlot_70<- ggplot(new_dat_70, aes(x=as.numeric(year),y=value,color=name))+ 
  geom_line()+geom_point(data = new_dat_70[new_dat_70$name == "Observed",c(2,5,6)] )+
  scale_color_manual("",values=c('black','darkslateblue'))+
  geom_ribbon(data = new_dat_70, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.1,.3,.8,.4),color="NA")+ 
  xlab("Year")+ylab(expression("Catch / 1000 m"^2))+ 
theme_minimal() + coord_cartesian(ylim = c(0,9)) + 
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015))+ 
  annotate("text", x=2015, y=7, label= "70 m", size=5, col="black")+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

#full 
legend <- get_legend(modelPlot_10)
modelPlot_10 <- modelPlot_10 + theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
modelPlot_25 <- modelPlot_25 + theme(legend.position="none") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
modelPlot_50 <- modelPlot_50 + theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
modelPlot_70 <- modelPlot_70 + theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

top_row <- plot_grid(NULL, modelPlot_10, NULL, modelPlot_25, legend, 
                     labels = c('a','','b','', ''), label_size = 15, 
                     rel_widths = c(.1, 1, 0.15, 1, .6), nrow=1)
bottom_row <- plot_grid(NULL, modelPlot_50, NULL, modelPlot_70,  NULL, 
                        labels = c('c','','d' ,'', ''), label_size = 15, 
                        rel_widths = c(.1, 1, 0.15, 1, .6),nrow=1)
all_row<- plot_grid(top_row, bottom_row, ncol = 1)


line_1 <- "Moderate Susceptibility Species"

y.lab <- textGrob(expression(paste("Sea star density (number of sea stars / 1000 m"^2,")")),
                  gp=gpar(fontsize=14,font=8), rot = 90, vjust = 1)
x.lab <- textGrob("Year", gp=gpar(fontsize=14))
add_x_lab<- plot_grid(all_row, x.lab, ncol = 1,  rel_heights = c(1, 0.07))
add_y_lab <- plot_grid(y.lab, add_x_lab, nrow = 1, rel_widths = c(0.07, 1))                        
add_y_lab <- add_y_lab + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
add_y_lab<- add_y_lab + draw_label(line_1, x = 0.01, y = 0.5, angle = 90, fontface = "bold")
figure_4 <- add_y_lab 
figure_4

######======Models of total catch, high susceptibility group =========

# Depth
m1<-glmmTMB(total_catch_HM~est_depth+
              ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m1)
plot_fun(m1,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m1,plot=T) 

# Depth + PrePost (i.e., PrePost effect the same across depths)
m2<-glmmTMB(total_catch_HM~est_depth +
              prepost+
              ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m2)
plot_fun(m2,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m2,plot=T) 

# Depth + PrePost * Depth
m3<-glmmTMB(total_catch_HM~ est_depth+
              prepost:est_depth+
              ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)), 
            family="nbinom2",data=dat_all)
summary(m3)
plot_fun(m3,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m3,plot=T) 

# Depth + Year (i.e., year effect the same across depths)
m4<-glmmTMB(total_catch_HM~est_depth + 
              year_cont+
              ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m4)
plot_fun(m4,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m4,plot=T) 

# Depth + Year * Depth 
m5<-glmmTMB(total_catch_HM~est_depth + 
              year_cont:est_depth+
              ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m5)
plot_fun(m5,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m5,plot=T) 

# Depth + PrePost + Year (prepost and year effects common across depths)
m6<-glmmTMB(total_catch_HM~est_depth+
              year_cont+
              prepost+
              ar1(year+0|est_depth)+
              offset(log(trawlarea/1e3)),family="nbinom2",data=dat_all)
summary(m6)
plot_fun(m6,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m6,plot=T) 

# Depth + PrePost * Year * Depth
m7<-glmmTMB(total_catch_HM~ est_depth+ 
              year_cont:est_depth:prepost+
              ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m7)
plot_fun(m7,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m7,plot=T) 

# Depth + PrePost * Year (i.e., unique slope for pre and post, but same across depths)
m8<-glmmTMB(total_catch_HM~est_depth + 
              year_cont:prepost+
              ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m8)
plot_fun(m8,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m8,plot=T) 

# Depth + PrePost * Depth + Year * Depth 
m9<-glmmTMB(total_catch_HM~ est_depth + 
              year_cont:est_depth+
              prepost:est_depth+
              ar1(year+0|est_depth)+ offset(log(trawlarea/1e3)),
            family="nbinom2",data=dat_all)
summary(m9)
plot_fun(m9,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m9,plot=T) 


# Model selection 
AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9) #m2 is ranked best

#add temp to best model fit
m.temp<-glmmTMB(total_catch_HM~est_depth +
                  meanTemp + 
                  prepost +
                   ar1(year+0|est_depth)+offset(log(trawlarea/1e3)),
                family="nbinom2",data=dat_all)
summary(m.temp) # temp is n.s.
plot_fun(m.temp,plot_rand=~0)
sim<-DHARMa::simulateResiduals(m.temp,plot=T) 

######======plotting chosen model of high susceptibility catch (Fig 3)=========

pred<-predict(m2,se=TRUE, re.form=~0)

new_dat<-dat_all %>% 
  left_join(m2$frame %>% mutate(fit=pred$fit,se.fit=pred$se.fit)) %>% 
  # calculate confidence intervals and fit on the normal scale
  mutate(pred_d=exp(fit),l95=exp(fit-1.96*se.fit),u95=exp(fit+1.96*se.fit)) %>% 
  group_by(est_depth,year) %>% 
  summarize(Observed=mean(total_density_HM)*1000, Predicted=mean(pred_d),u95=mean(u95),l95=mean(l95)) %>%
  arrange(year) %>%   
  pivot_longer(c(Observed,Predicted)) %>% mutate(est_depth=paste(est_depth,"m")) 

# Because this is best fit model, let's get estimates + CI for catch at each depth in 1997
 new_dat %>% filter(year == 1997 & name == "Predicted") %>% group_by(est_depth) %>% summarise(mean(value))
 new_dat %>% filter(year == 1997 & name == "Predicted") %>% group_by(est_depth) %>% summarise(mean(l95))
 new_dat %>% filter(year == 1997 & name == "Predicted") %>% group_by(est_depth) %>% summarise(mean(u95))

#10 
new_dat_10 <- new_dat %>% filter(est_depth == "10 m")
modelPlot_10<- ggplot(new_dat_10, aes(x=as.numeric(year),y=value,color=name)) +
  geom_line()+geom_point(data = new_dat_10[new_dat_10$name == "Observed",c(2,5,6)] )+ 
  scale_color_manual("",values=c('black','darkorange4'),
                      guide = guide_legend(override.aes = list(
                       linetype = c("solid", "solid"), shape=c(19, NA)))) + 
  geom_ribbon(data = new_dat_10, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.9,.5,0,.4),color="NA")+
  xlab("Year")+ylab(expression("Catch / 1000 m"^2)) +
  theme_minimal() +
  coord_cartesian(ylim = c(0,4.5)) + 
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015))+
  annotate("text", x=2015, y=3.7, label= "10 m", size=5, col="black")+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text=element_text(size=13)) 

#25
new_dat_25 <- new_dat %>% filter(est_depth == "25 m")
modelPlot_25<- ggplot(new_dat_25, aes(x=as.numeric(year),y=value,color=name)) +
  geom_line()+geom_point(data = new_dat_25[new_dat_25$name == "Observed",c(2,5,6)] )+ 
  scale_color_manual("",values=c('black','darkorange4'))+ 
  geom_ribbon(data = new_dat_25, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.9,.5,0,.4),color="NA")+
  xlab("Year")+ylab(expression("Catch / 1000 m"^2))+ 
  theme_minimal() +
  coord_cartesian(ylim = c(0,4.5)) + 
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015))+
  annotate("text", x=2015, y=3.7, label= "25 m", size=5, col="black")+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

#50
new_dat_50 <- new_dat %>% filter(est_depth == "50 m")
modelPlot_50<- ggplot(new_dat_50, aes(x=as.numeric(year),y=value,color=name)) +
  geom_line()+geom_point(data = new_dat_50[new_dat_50$name == "Observed",c(2,5,6)] )+ 
  scale_color_manual("",values=c('black','darkorange4'))+ 
  geom_ribbon(data = new_dat_50, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.9,.5,0,.4),color="NA")+
  xlab("Year")+ylab(expression("Catch / 1000 m"^2))+
  theme_minimal() +
  coord_cartesian(ylim = c(0,4.5)) + 
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015))+
  annotate("text", x=2015, y=3.7, label= "50 m", size=5, col="black")+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

#70
new_dat_70 <- new_dat %>% filter(est_depth == "70 m")
modelPlot_70<- ggplot(new_dat_70, aes(x=as.numeric(year),y=value,color=name)) +
  geom_line()+geom_point(data = new_dat_70[new_dat_70$name == "Observed",c(2,5,6)] )+ 
  scale_color_manual("",values=c('black','darkorange4'))+ 
  geom_ribbon(data = new_dat_70, aes(x=as.numeric(year),ymax=u95,ymin=l95),fill=rgb(.9,.5,0,.4),color="NA")+
  xlab("Year")+ylab(expression("Catch / 1000 m"^2))+ 
  theme_minimal() +
  coord_cartesian(ylim = c(0,4.5)) + 
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015))+
  annotate("text", x=2015, y=3.7, label= "70 m", size=5, col="black")+ 
  theme(axis.text=element_text(size=10), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 


#full 
legend <- get_legend(modelPlot_10)
modelPlot_10 <- modelPlot_10 + theme(legend.position="none") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
modelPlot_25 <- modelPlot_25 + theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
modelPlot_50 <- modelPlot_50 + theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
modelPlot_70 <- modelPlot_70 + theme(legend.position="none") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

top_row <- plot_grid(NULL, modelPlot_10, NULL, modelPlot_25, legend, 
                     labels = c('a','','b','', ''), label_size = 15, 
                     rel_widths = c(.1, 1, 0.15, 1, .6), nrow=1)
bottom_row <- plot_grid(NULL, modelPlot_50, NULL, modelPlot_70,  NULL, 
                        labels = c('c','','d' ,'', ''), label_size = 15, 
                        rel_widths = c(.1, 1, 0.15, 1, .6),nrow=1)
all_row<- plot_grid(top_row, bottom_row, ncol = 1)


line_1 <- "High Susceptibility Species"

y.lab <- textGrob(expression(paste("Sea star density (number of sea stars / 1000 m"^2,")")), 
                  gp=gpar(fontsize=14,font=8), rot = 90, vjust = 1)
x.lab <- textGrob("Year", gp=gpar(fontsize=14))
add_x_lab<- plot_grid(all_row, x.lab, ncol = 1,  rel_heights = c(1, 0.07))
add_y_lab <- plot_grid(y.lab, add_x_lab, nrow = 1, rel_widths = c(0.07, 1))                        
add_y_lab <- add_y_lab + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
add_y_lab<- add_y_lab + draw_label(line_1, x = 0.01, y = 0.5, angle = 90, fontface = "bold")
figure_3<- add_y_lab
figure_3

######======plotting species-level trends in catch overtime (Fig 2)====

dat_all_pivot <- dat_all %>% 
  ungroup() %>%
  select(!c(year_cont, distance, trawlarea, prepost, actual_depth, 
            total_catch, total_density, total_density_HM, total_density_SM, 
            total_density_UM, meanTemp, total_catch_HM,
            total_catch_SM,total_catch_UM, total_catch_SM_UM , 
            total_density_SM_UM)) %>%
  pivot_longer(!c(est_depth, time, year),names_to = "species_name", 
               values_to = "catch") %>% 
  unique() %>%
  mutate(susceptibility = 
           case_when(species_name == "S.stimpsoni" | 
                       species_name == "P.helianthoides"|
                       species_name == "S.dawsoni"| 
                       species_name == "P.brevispinus"| 
                       species_name == "E.troschelli" ~ "High Susceptibility",
                       species_name == "H.leviuscula" | 
                       species_name == "L.foliolata"| 
                       species_name == "H.spinosa"| 
                       species_name == "D.imbricata"| 
                       species_name == "M.aequalis" | 
                       species_name == "C.papposus" ~ "Moderate Susceptibility",
                                    TRUE ~ "")) %>% 
  group_by(est_depth, year, species_name, susceptibility) %>%
  summarize(summed_catch_across_times = sum(catch)) 

dat_all_pivot$est_depth <- paste(dat_all_pivot$est_depth, "m")
dat_all_pivot$year <- as.numeric(dat_all_pivot$year)
dat_all_pivot$species_name <- gsub("[.]", ". ", dat_all_pivot$species_name)

#unique(ggplot_build(plot)$data[[1]]$fill) 
#extract colors of paired palette so I can manually change 1 color

#plot
species_plot <-dat_all_pivot[which(dat_all_pivot$summed_catch_across_times >0),] %>% 
  ggplot(aes(fill = species_name, x = year, y = summed_catch_across_times)) +
  geom_bar(position="stack", stat="identity", color="black", lwd=0.1) +
  facet_wrap( ~ est_depth, ncol = 2) + 
  scale_x_continuous(name="Year",  breaks=c(2000,2005,2010,2015)) + 
  theme_bw() +  
  labs(y = "Total Catch", x = "Year") + 
  scale_fill_manual(values = 
                      c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C","#FB9A99",
                        "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6","#6A3D9A",
                        "yellow1")) +
  guides(fill=guide_legend(title="Species")) + 
 theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
       legend.text = element_text(size = 8, face = "italic"),
       legend.title = element_text(size = 10),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank())

figure_2 <- species_plot
figure_2

######======Testing for serial depletion during trawls========

# Test whether or not serial depletion occurred throughout the day during 
# sampling.  

trawltime_rank <- dat_all %>% mutate(trawl_order = case_when(time == "afternoon" ~ 1,
                                                             time == "evening"~ 2,
                                                             time == "night" ~ 3,
                                                             time == "early morning" ~ 4,
                                                             time == "morning"~ 5))

mod<-glmmTMB(total_catch~as.numeric(trawl_order) + as.numeric(trawl_order):est_depth+
               est_depth,family="poisson",data=trawltime_rank)
summary(mod)

ggplot(trawltime_rank, aes(x=trawl_order,y=total_catch,color=est_depth)) + 
  geom_point() + 
  xlab("time")+ylab(expression("Catch ")) +  
  facet_wrap( ~ est_depth) + 
  theme_classic()  