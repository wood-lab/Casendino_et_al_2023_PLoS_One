###############
# Uses the cpm package to identify changepoints in moderate-susceptibility and 
# high susceptibility sea star catch from 1997 to 2019.
# Includes code for Figures 5 and 6 in Casendino et al. (2023)

#####======DEPENDENCIES=======

library(readr)
library(tidyverse)
library(dplyr)
library(janitor)
library(ggeffects)
library(cpm)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(gtable)
library(grid)
library(magrittr)

########===Cleaning up data=====

source("src/data_processing_function.R")
dat_all<-process_dat() %>% 
  mutate(est_depth=as.factor(est_depth))

# make year a continuous variable (first year is 0, second year is 1, etc)
dat_all$year <- as.numeric(dat_all$year)
dat_all[,ncol(dat_all) + 1] <- dat_all[,8]
dat_all[1:430,8] <-  dat_all[1:430,8] - 1997   
colnames(dat_all)[8] <- "year_cont"
dat_all$year <- as.character(dat_all$year) 

dat_all %<>% mutate(total_density_SM_UM = sum(total_density_SM + total_density_UM)) 
# combine species assigned "some mortality" with species assigned "unknown mortality" (likely affected) 
# to produce a group of moderate-susceptibility species

####=====CPM: High susceptibility subset===============

#average across depths for each year.
dat_all$year <- as.numeric(dat_all$year)
density_across_depth<- dat_all %>%
  group_by(year) %>%
  summarise(n = mean(total_density_HM))

densityvector <- select(density_across_depth,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
# NONE

# Plot 

cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000)
cpPlot <- ggplot(data=cpPlotData) + 
  geom_line(aes(year,density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", size = 0.9) +
  annotate("text", x=2001, y=5.5, label= "all", size=6, col="black") + 
  annotate("text", x=2001, y=4.9, label= "depths", size=6, col="black") + 
  theme_bw()   + ylim(0, 7) +  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_blank()) +
  scale_linetype_manual(name = element_blank() , values = c(2),
                        guide = guide_legend(override.aes = list(color = c( "ivory4")))) +
  theme(axis.text=element_text(size=13))

## 10m
density_10m<- dat_all %>%
  filter( est_depth == "10")

density_10m<- density_10m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_HM))

densityvector <- select(density_10m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)
plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
# NONE

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_10 <- ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", size = 0.9) +
  annotate("text", x=2001, y=5.5, label= "10 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7) +  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() , 
                        values = c(2),
                        guide = guide_legend(override.aes = list(color = c("ivory4")))) +
  theme(axis.text=element_text(size=13))

## 25m
density_25m<- dat_all %>%
  filter( est_depth == "25")

density_25m<- density_25m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_HM))

densityvector <- select(density_25m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
# NONE

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_25 <- ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", size = 0.9) +
  annotate("text", x=2001, y=5.5, label= "25 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7) +  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() , values = c( 2),
                        guide = guide_legend(override.aes = list(color = c("ivory4")))) +
  theme(axis.text=element_text(size=13))

## 50m
density_50m<- dat_all %>%
  filter( est_depth == "50")

density_50m<- density_50m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_HM))

densityvector <- select(density_50m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
#NONE

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_50 <- ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD outbreak"),  colour="ivory4", size = 0.9) +
  annotate("text", x=2001, y=5.5, label= "50 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) +
  theme(legend.text = element_text(size=17)) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() , 
                        values = c( 2),guide = guide_legend(override.aes = list(color = c("ivory4"))))+
  theme(axis.text=element_text(size=13))


#######  70m
density_70m<- dat_all %>%
  filter( est_depth == "70")

density_70m<- density_70m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_HM))

densityvector <- select(density_70m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
#none

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_70 <- ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", size = 0.9) +
  annotate("text", x=2001, y=5.5, label= "70 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7)+ theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015",""))  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() , 
                        values = c( 2),
                        guide = guide_legend(keyheight=2,
                              override.aes = list(color = c("ivory4")))) + 
   theme(axis.text=element_text(size=13),
         legend.text=element_text(size=13))

# Code for Figure 5 (changepoint analysis for high-susceptibility species)
legend <- get_legend(cpPlot_70)

cpPlot <- cpPlot + theme(legend.position="none")
cpPlot_10 <- cpPlot_10 + theme(legend.position="none")
cpPlot_25 <- cpPlot_25 + theme(legend.position="none")
cpPlot_50 <- cpPlot_50 + theme(legend.position="none")
cpPlot_70 <- cpPlot_70 + theme(legend.position="none")

top_row <- plot_grid(cpPlot_10, cpPlot_25, legend, labels = c('a', 'b', ''), 
                     label_size = 16, nrow=1, label_x = 0, label_y = 0,
                     hjust = -0.5, vjust = -0.5)
bottom_row <- plot_grid(cpPlot_50, cpPlot_70,  cpPlot, labels = c('c', 'd', 'e'), 
                        label_size = 16, nrow=1, label_x = 0, label_y = 0,
                        hjust = -0.5, vjust = -0.5)
all_row<- plot_grid(top_row, bottom_row, ncol = 1)
y.lab <- textGrob(expression(paste("Sea star density (number of sea stars / 1000 m"^2,")")),
                  gp=gpar(fontsize=15,font=8), rot = 90, vjust = 1)

x.lab <- textGrob("Year", 
                  gp=gpar(fontsize=15))

add_x<- plot_grid(all_row, x.lab, ncol = 1,  rel_heights = c(1, 0.07))
complete_plot <- plot_grid(y.lab, add_x, nrow = 1, rel_widths = c(0.07, 1))                        

line_1 <- "High Susceptibility Species"
complete_plot<- complete_plot + draw_label(line_1, x = 0.01, y = 0.5, angle = 90, 
                                           fontface = "bold")
figure_5 <- complete_plot
figure_5

####=====CPM: Moderate susceptibility subset===============

dat_all$year <- as.numeric(dat_all$year)
density_across_depth<- dat_all %>%
  group_by(year) %>%
  summarise(n = mean(total_density_SM_UM))

densityvector <- select(density_across_depth,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
density_across_depth$year[9]
# 2006

# Plot 

cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000)
cpPlot <- ggplot(data=cpPlotData) + 
  geom_line(aes(year,density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", size = 0.9) +
  geom_vline(aes(xintercept = 2006, linetype = "Changepoint"),  colour="royalblue2", size = 0.9)  +
  annotate("text", x=2000, y=5.5, label= "all", size=6, col="black") + 
  annotate("text", x=2000, y=4.9, label= "depths", size=6, col="black") + 
  theme_bw()   + ylim(0, 7) +  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_blank()) +
  scale_linetype_manual(name = element_blank() , values = c(1, 2),
                        guide = guide_legend(override.aes = list(color = c("royalblue2", "ivory4")))) +
  theme(axis.text=element_text(size=13))
cpPlot

##### 10m
density_10m<- dat_all %>%
  filter( est_depth == "10")

density_10m<- density_10m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_SM_UM))

densityvector <- select(density_10m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)
plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
# none

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_10 <- ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", size = 0.9) +
  annotate("text", x=2002, y=5.5, label= "10 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7) +  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() ,
                        values = c(2),guide = guide_legend(override.aes = list(color = c( "ivory4")))) +
  theme(axis.text=element_text(size=13))

##### 25m
density_25m<- dat_all %>%
  filter( est_depth == "25")

density_25m<- density_25m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_SM_UM))

densityvector <- select(density_25m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
# none

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_25 <-ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", linetype="dashed", size = 0.9) +
   annotate("text", x=2005, y=5.5, label= "25 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7) + theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() , values = c(2),
                        guide = guide_legend(override.aes = list(color = c( "ivory4")))) +
  theme(axis.text=element_text(size=13))

##### 50m
density_50m<- dat_all %>%
  filter( est_depth == "50")

density_50m<- density_50m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_SM_UM))

densityvector <- select(density_50m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
#none

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_50 <- ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD outbreak"),  colour="ivory4", size = 0.9) +
   annotate("text", x=2001, y=5.5, label= "50 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015","")) +
  theme(legend.text = element_text(size=17)) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())  +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() , 
                        values = c( 2),guide = guide_legend(override.aes = list(color = c( "ivory4"))))+
  theme(axis.text=element_text(size=13))

#######  70m
density_70m<- dat_all %>%
  filter( est_depth == "70")

density_70m<- density_70m %>%
  group_by(year) %>%
  summarise(n = mean(total_density_SM_UM))

densityvector <- select(density_70m,n)
densityvector<- as.vector(densityvector)
densityvector<- as.numeric(unlist(densityvector))

fit_cpm = processStream(densityvector, cpmType="Mann-Whitney", ARL0=370, startup=22)

plot(densityvector,type="l")
for (i in 1:length(fit_cpm$changePoints)){
  abline(v=fit_cpm$changePoints[i], lty=2)}
fit_cpm$changePoints
density_across_depth$year[9]
#2006

# Plot 
cpPlotData<- data.frame(year=c(1997,1999:2019),
                        density_vec=densityvector*1000 )
cpPlot_70 <- ggplot(data=cpPlotData) + 
  geom_line(aes(year, density_vec), lty=1, color="black") + 
  geom_vline(aes(xintercept = 2013, linetype = "SSWD Outbreak"),  colour="ivory4", size = 0.9) +
  geom_vline(aes(xintercept = 2006, linetype = "Changepoint"),  colour="royalblue2", size = 0.9)  +
  annotate("text", x=2001, y=5.5, label= "70 m", size=6, col="black") + 
  theme_bw() + ylim(0, 7)+ theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  scale_x_continuous(name="Year", limits=c(1997, 2019), labels=c("1997","2000","2005","2010","2015",""))  +
   theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.title = element_blank())  +
  scale_linetype_manual(name = element_blank() , values = c(1, 2),
                        guide = guide_legend(keyheight=2, 
                                             override.aes = list(color = c("royalblue2", "ivory4")))) + 
  theme(axis.text=element_text(size=13),
        legend.text=element_text(size=13)) 

# ALL PLOTS 
legend <- get_legend(cpPlot_70)

cpPlot <- cpPlot + theme(legend.position="none")
cpPlot_10 <- cpPlot_10 + theme(legend.position="none")
cpPlot_25 <- cpPlot_25 + theme(legend.position="none")
cpPlot_50 <- cpPlot_50 + theme(legend.position="none")
cpPlot_70 <- cpPlot_70 + theme(legend.position="none")

top_row <- plot_grid(cpPlot_10, cpPlot_25, legend, labels = c('a', 'b', ''), 
                     label_size = 16, nrow=1, label_x = 0, label_y = 0,
                     hjust = -0.5, vjust = -0.5)
bottom_row <- plot_grid(cpPlot_50, cpPlot_70,  cpPlot, labels = c('c', 'd', 'e'), 
                        label_size = 16, nrow=1, label_x = 0, label_y = 0,
                        hjust = -0.5, vjust = -0.5)
all_row<- plot_grid(top_row, bottom_row, ncol = 1)
y.lab <- textGrob(expression(paste("Sea star density (number of sea stars / 1000 m"^2,")")), 
                  gp=gpar(fontsize=15,font=8), rot = 90, vjust = 1)

x.lab <- textGrob("Year", 
                  gp=gpar(fontsize=15))

add_x<- plot_grid(all_row, x.lab, ncol = 1,  rel_heights = c(1, 0.07))
complete_plot <- plot_grid(y.lab, add_x, nrow = 1, rel_widths = c(0.07, 1))                        

line_1 <- "Moderate Susceptibility Species"
complete_plot<- complete_plot + draw_label(line_1, x = 0.01, y = 0.5, angle = 90, 
                                           fontface = "bold")
complete_plot

figure_6 <- complete_plot
figure_6
