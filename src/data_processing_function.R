###############
# Data processing function used in Casendino et al. (2023)
# Standardizes trawling data, invertebrate data, and environmental data from 
# 1997 to 2019 to be used for modeling and changepoint analyses (analyses in separate script).
####################

# DEPENDENCIES: 

library(readr)
library(tidyverse)
library(dplyr)
library(janitor)
library(stringi)

process_dat<-function(){
# Read in trawling data, rename relevant columns, remove empty rows and columns. 
rawdata <- read_csv("data/trawling_dataset.csv")
data<-select(rawdata, colnames(rawdata)[2], colnames(rawdata)[4],colnames(rawdata)[5],colnames(rawdata)[11])
data<- data %>% rename(date=colnames(rawdata)[2]) %>% rename(time=colnames(rawdata)[4]) %>% 
          rename(depth=colnames(rawdata)[5]) %>% rename(distance=colnames(rawdata)[11])
data <- janitor::remove_empty(data)
data<- data[-c(1:3),]
attach(data)

######=========2019 TRAWL DATA ONLY ==========
# Depth was organized for 2019 differently than other years, 
# so it was processed separately.

# Removing irrelevant (empty) data 
data_only2019<-select(rawdata, colnames(rawdata)[1],colnames(rawdata)[2],colnames(rawdata)[4],colnames(rawdata)[11])
data_only2019<- data_only2019 %>% rename(date=colnames(rawdata)[2]) %>% 
  rename(time=colnames(rawdata)[4]) %>% rename(depth=colnames(rawdata)[1]) %>% 
  rename(distance=colnames(rawdata)[11])
data_only2019 <- janitor::remove_empty(data_only2019)
data_only2019<- data_only2019[-c(1:4),]

# Isolate 2019 data
data_only2019<-data_only2019 %>% slice(1506:1585)

### Adding column for date (2019)
data_only2019 <- data_only2019 %>% fill(date) 
data_only2019$date <- as.Date(data_only2019$date,  "%d-%B-%y")

### Each trawl tow was 370 m.
data_only2019<- data_only2019 %>% fill(distance) 

#### Arrange rows by time to determine which belong to which time of day. 
data_only2019 <- data_only2019 %>% arrange(time) 
data_night<- data_only2019 %>% slice(1:8) %>% select(date,depth,distance) %>% 
  mutate(time="night") %>% arrange(date)
data_early_morning<- data_only2019 %>% slice(9:16) %>% select(date,depth,distance) %>% 
  mutate(time="early morning") %>% arrange(date)
data_morning<- data_only2019%>% slice(17:24) %>% select(date,depth,distance) %>% 
  mutate(time="morning") %>% arrange(date)
data_afternoon<- data_only2019 %>% slice(25:32) %>% select(date,depth,distance) %>% 
  mutate(time="afternoon") %>% arrange(date)
data_evening<- data_only2019 %>% slice(33:40) %>% select(date,depth,distance) %>%
  mutate(time="evening") %>% arrange(date)

#### Rejoin the above tables. 
data_only2019 <- do.call("rbind", list(data_night,data_early_morning,data_morning, data_afternoon,data_evening))

### Repeat to assign depth. 
data_only2019<-data_only2019 %>% arrange(depth)
data_10<- data_only2019 %>% slice(1:5) %>% select(date,time,distance) %>% 
  mutate(depth="10") %>% arrange(date)
data_25<- data_only2019 %>% slice(6:10) %>% select(date,time,distance) %>% 
  mutate(depth="25") %>% arrange(date)
data_50<-data_only2019 %>% slice(11:15) %>% select(date,time,distance) %>% 
  mutate(depth="50") %>% arrange(date)
data_70<- data_only2019 %>% slice(16:20) %>% select(date,time,distance) %>% 
  mutate(depth="70") %>% arrange(date)
data_only2019<- do.call("rbind", list(data_10,data_25,data_50, data_70))

#####===========ORGANIZING REMAINING TRAWL DATA (2000-2018)============

### Dropped the NAs in the depth category to get rid of irrelevant entries. 
# Filled in NAs with respective years.
data<- drop_na(data, depth)
data<- data %>% fill(date)
data<-data %>% arrange(time)

# Arrange rows by time to determine which belong to which time of day. 

data_night<- data %>% slice(1:207) %>% select(date,depth,distance) %>% mutate(time="night") %>% arrange(date)
data_early_morning<- data %>% slice(208:421) %>% select(date,depth,distance) %>% mutate(time="early morning") %>% arrange(date)
data_morning<- data %>% slice(422:637) %>% select(date,depth,distance) %>% mutate(time="morning") %>% arrange(date)
data_afternoon<- data %>% slice(638:839) %>% select(date,depth,distance) %>% mutate(time="afternoon") %>% arrange(date)
data_evening<- data %>% slice(840:1063) %>% select(date,depth,distance) %>% mutate(time="evening") %>% arrange(date)
data_night_2<- data %>% slice(1064:1072) %>% select(date,depth,distance) %>% mutate(time="night") %>% arrange(date)
data_unknown <- data %>% slice(1083:1084, 1086:1087) %>% select(date,depth,distance) %>% mutate(time="unknown") %>% arrange(date)

### Join the time-grouped dataframes. 
data<- do.call("rbind", list(data_night,data_night_2,data_early_morning,data_morning,data_afternoon, data_evening,data_unknown))

# Change date column's class to date class, and depth column's class to numeric class.
data$date <- as.Date(data$date, "%d-%B-%y")
data <- data %>% arrange(date)
data$depth <- as.numeric(data$depth)

### Including continuous depths (as opposed to 10, 25,50 and 70).

# Fill in "unknown" times. Make sure sets of 3 depths are together by deleting the row 
#in the wrong position and putting it in the correct place.  
which(data$date == "2005-05-13" & data$depth=="10" & data$time=="unknown")
data<- data[-309,]
data <- data %>% add_row(date =as.Date("2005-05-13", "%Y-%m-%d"), depth=as.numeric(10), 
                         distance="0.20 nm", time ="afternoon"    
                         ,.before = 290)
data[310,4] <- "evening"

data <- data[-442,]
data <- data %>% add_row(date =as.Date("2007-05-11", "%Y-%m-%d"), depth=as.numeric(24), 
                         distance="0.20 nm", time ="afternoon" ,.before = 423)

data <- data[-466,]
data <- data %>% add_row(date =as.Date("2007-05-12", "%Y-%m-%d"), depth=as.numeric(52), 
                         distance="0.20 nm", time ="early morning" ,.before = 461)

data <- data[-289,]
data <- data %>% add_row(date =as.Date("2004-05-15", "%Y-%m-%d"), depth=as.numeric(13), 
                         distance="166  degrees mag.", time ="early morning" ,.before = 268)

data <- data[-181,]
data <- data %>% add_row(date =as.Date("2003-05-16", "%Y-%m-%d"), depth=as.numeric(48),
                         distance="0.20 nm", time ="evening" ,.before = 187)

data <- data %>% add_row(date =as.Date("2003-05-16", "%Y-%m-%d"), depth=as.numeric(46), 
                         distance="Hang up, Rocks!", time ="evening" ,.before = 190)

data <- data %>% add_row(date =as.Date("2003-05-16", "%Y-%m-%d"), depth=as.numeric(46), 
                         distance="Hang up, Rocks!",  time ="evening" ,.before = 190)

data <- data %>% add_row(date =as.Date("2013-05-17", "%Y-%m-%d"), depth=as.numeric(12), 
                         distance=NA,  time ="evening" ,.before = 769)

data <- data %>% add_row(date =as.Date("2018-05-11", "%Y-%m-%d"), depth=as.numeric(12),
                         distance=NA, time ="afternoon" ,.before = 1021)

### For each set of 3 depth observations, take the average of the start tow and end two 
#(the last two of the three), put in  avg_depth column. 
data <- data %>% mutate(avg_depth = NA)
data$avg_depth <- as.numeric(data$avg_depth)
data<- data %>% arrange(date)

start_index <- 1 
for(i in 1:360){
  if(i ==1) {data[1,5] <- (data[[start_index + 1 ,2]]+ data[[start_index + 2,2]])/2
  }
  else {
    start_index <- start_index + 3
    data[start_index,5] <- ( data[[start_index + 1, 2]] + data[[start_index + 2, 2]])/2
  }
}


### Put in the standard depths, based on how they are ordered in the trawling dataset (10, 25, 50 and 70). 
data<-data %>% arrange(depth)
data_10<- data %>% slice(1:270) %>% select(date,time,distance,avg_depth) %>% 
  mutate(depth="10") %>% arrange(date)
data_25<- data %>% slice(271:540) %>% select(date,time,distance,avg_depth) %>%
  mutate(depth="25") %>% arrange(date)
data_50<- data %>% slice(543, 547, 549:788, 790:814, 816:818) %>% 
  select(date,time,distance,avg_depth) %>% mutate(depth="50") %>% arrange(date)
data_70<- data %>% slice(541:542, 544:546, 548, 789, 815, 819:1080) %>% 
  select(date,time,distance,avg_depth) %>% mutate(depth="70") %>% arrange(date)
data<- do.call("rbind", list(data_10,data_25, data_50,data_70))

# Make separate table with average continuous depths for each trawl tow. 
actualDepth<- data %>% select(date, time, avg_depth,depth) %>% unite(trawl_ID, date, time) 
actualDepth<- drop_na(actualDepth)

# Include 2019 depths (which are just the standard, discrete depths). 
avg_depth_2019 <- data_only2019 %>% unite(trawl_ID, date, time) %>% select(-distance) %>%
  rename(avg_depth= depth) 
avg_depth_2019$avg_depth <- as.numeric(avg_depth_2019$avg_depth)
actualDepth <- full_join(actualDepth, avg_depth_2019)
standard_depths<- c(rep(10, times=5), rep(25, times=5),rep(50, times=5),rep(70, times=5))
actualDepth[361:380,3] <- as.character(standard_depths) 
# actualDepth includes trawl tow ID, average depth for each trawl (continuous), 
#and discrete depth for each trawl. 

# Drop NAs from main table.
data <- data %>% select(-avg_depth)
data<- drop_na(data)
data<-data %>% arrange(date)

####===========# COMBINE 2019, 2000-2018, AND 1999 TRAWL DATA ================

### Join the 2019 dataframe with the rest of the data, and add 1999 trawl data. 
data<-rbind(data,data_only2019)
rawdata_1999 <- read_csv("data/trawling_dataset_1999.csv")

### Isolate 1999 from the raw data.  
data_1999<- rawdata_1999 %>% slice(176:195) %>% 
  select(`Puget Sound trawl survey times, locations and depths`,
         colnames(rawdata_1999)[3],colnames(rawdata_1999)[5], 
         colnames(rawdata_1999)[7],colnames(rawdata_1999)[13])
data_1999<- data_1999 %>% rename(date=`Puget Sound trawl survey times, locations and depths`) %>%
  rename(time=colnames(rawdata_1999)[5]) %>% rename(depth=colnames(rawdata_1999)[7]) %>%
  rename(distance=colnames(rawdata_1999)[13]) %>% 
  mutate(date= as.Date("1999-May-14", "%Y-%B-%d")) %>% select(-colnames(rawdata_1999)[3])

### Replace "dawn" with "early morning". 
data_early_morning_1999<- data_1999 %>% slice(13:16) %>% select(date,depth,distance) %>%
  mutate(time="early morning")
data_1999 <- data_1999 %>% slice(1:12,17:20)
data_1999<- full_join(data_1999, data_early_morning_1999)

### Add 1999 depths to actualDepth dataframe. 
avg_depth_1999 <- data_1999 %>% unite(trawl_ID, date, time) %>% select(-distance) %>% 
  rename(avg_depth= depth) 
avg_depth_1999 <- avg_depth_1999 %>% mutate(depth=avg_depth)
avg_depth_1999$avg_depth <- as.numeric(avg_depth_1999$avg_depth)
actualDepth <- full_join(avg_depth_1999, actualDepth)

### Join all years together. 
data<- full_join(data, data_1999)

### To sort out distance, sort the distance column in ascending order and 
#remove section with non-distance measurements. 
# Then, add in a column of distances based on the distances we see in specific 
#rows to standardize distance measurements.   

data<- data %>% arrange(distance)
data<- data %>% slice(1:201, 755:946, 947:953)

data_370<- data %>% slice(1:193, 214:307) %>% select(date,depth,time) %>% mutate(distance=370)
data_371<- data %>% slice(308:337) %>% select(date,depth,time) %>% mutate(distance=371)
data_372<- data %>% slice(194, 338:359) %>% select(date,depth,time) %>% mutate(distance=372)
data_373<- data %>% slice(195:196, 360:368) %>% select(date,depth,time) %>% mutate(distance=373)
data_374<- data %>% slice(369:376) %>% select(date,depth,time) %>% mutate(distance=374)
data_375<- data %>% slice(377:380) %>% select(date,depth,time) %>% mutate(distance=375)
data_389<- data %>% slice(395,199) %>% select(date,depth,time) %>% mutate(distance=389)
data_400<- data %>% slice(200,397) %>% select(date,depth,time) %>% mutate(distance=400)
data_445<- data %>% slice(201) %>% select(date,depth,time) %>% mutate(distance=445)
data_367<- data %>% slice(202) %>% select(date,depth,time) %>% mutate(distance=367)
data_368<- data %>% slice(203:205) %>% select(date,depth,time) %>% mutate(distance=368)
data_369<- data %>% slice(206:213) %>% select(date, depth, time) %>% mutate(distance=369)
data_421<- data %>% slice(397) %>% select(date,depth,time) %>% mutate(distance=400)
data_482<- data %>% slice(398) %>% select(date,depth,time) %>% mutate(distance=421)
data_538<- data %>% slice(399) %>% select(date,depth,time) %>% mutate(distance=482)
data_376<- data %>% slice(381:382) %>% select(date,depth,time) %>% mutate(distance=376)
data_378<- data %>% slice(383:384) %>% select(date,depth,time) %>% mutate(distance=378)
data_379<- data %>% slice(385) %>% select(date,depth,time) %>% mutate(distance=379)
data_380<- data %>% slice(386:389) %>% select(date,depth,time) %>% mutate(distance=380)
data_382<- data %>% slice(390) %>% select(date,depth,time) %>% mutate(distance=382)
data_385<- data %>% slice(391:392) %>% select(date,depth,time) %>% mutate(distance=385)
data_386<- data %>% slice(394) %>% select(date,depth,time) %>% mutate(distance=386)
data_395<- data %>% slice(396) %>% select(date,depth,time) %>% mutate(distance=395)
data_538<- data %>% slice(400) %>% select(date,depth,time) %>% mutate(distance=538)

#Then, join all distance dataframes. 
data <- do.call("rbind", list(data_370,data_371,data_372, data_373,
                              data_374,data_375,data_389,data_400,data_445,data_367,
                              data_368,data_369, data_421,data_482,data_538, data_376,
                              data_378,data_379, data_380,data_382, data_385, data_386, data_395, data_538))
data <- data %>% arrange(date)

## Add trawl area column. (3.5 = net width)
data<- data %>% mutate(trawlarea=3.5*distance)

#####========ADDING INVERTS SPREADSHEET==============

rawdata_inv<- read_csv("data/invertebrates_dataset.csv", col_types = cols( .default=col_character(), number = col_double()))
data_inv<-rawdata_inv %>% select(year,time,depth,`genus species`,number, group) %>% 
  filter(group=='sea star') %>% select(year,time,depth,`genus species`,number) %>% 
  filter(`genus species` !='NA')

### Error in data: the star labelled "Crossaster papposus" 
#caught in 2006 was actually the species "Hippasteria spinosa".
data1<- data_inv %>% slice(1:58,60:110)
data2<- data_inv %>% slice(59)
data2<- data2 %>% select(year,time,depth,number) %>% mutate("genus species"="Hippasteria spinosa") %>% 
  select(year,time,depth, "genus species",number)
data_inv<- full_join(data1,data2)
data_inv<- arrange(data_inv,year)

######===ORGANIZING INVERTS DATA, MERGE WITH TRAWL DATA ========

# To get species-specific counts:
# Create dataframes for each species of sea star, then merge them with the inverts dataframe. 

MA <-  data_inv %>% filter(`genus species`== 'Mediaster aequalis') %>% 
  select(year,time,depth, number) %>% rename(M.aequalis=number)
CP <-  data_inv  %>% filter(`genus species`== 'Crossaster papposus') %>% 
  select(year,time,depth, number) %>% rename(C.papposus=number)
DI <-  data_inv %>% filter(`genus species`== 'Dermasterias imbricata') %>% 
  select(year,time,depth, number) %>% rename(D.imbricata=number)
HS <-   data_inv  %>% filter(`genus species`== 'Hippasteria spinosa') %>% 
  select(year,time,depth, number) %>% rename(H.spinosa=number)
M <- data_inv  %>% filter(`genus species`== 'menricia') %>% 
  select(year,time,depth, number) %>% rename(H.leviuscula=number)
PB<-  data_inv%>% filter(`genus species`== 'Pisaster brevispinus') %>% 
  select(year,time,depth, number) %>% rename(P.brevispinus=number)
PH<-  data_inv%>% filter(`genus species`== 'Pycnopodia helianthoides') %>% 
  select(year,time,depth, number) %>% rename(P.helianthoides=number)
SD<-  data_inv %>% filter(`genus species`== 'Solaster dawsoni') %>% 
  select(year,time,depth, number) %>% rename(S.dawsoni=number)
SS<-  data_inv %>% filter(`genus species`== 'Solaster stimpsoni') %>% 
  select(year,time,depth, number) %>% rename(S.stimpsoni=number)
AP <-  data_inv  %>% slice(99:100,105:106) %>% select(year,time,depth, number) %>% 
  rename(A.pugetana=number) %>% mutate("A.pugetana" = 0) 
## removing A. pugetana from analysis becuase it's an ophiuroid. 
#Changing counts to zero to avoid any hiccups in rest of script. 
LF <-  data_inv %>% slice(30:31,35:39, 81) %>% select(year,time,depth, number) %>% 
  rename(L.foliolata=number) ## These didn't subset correctly, subset by hand

#### Combine duplicate rows of trawls (specific year/time/depths) and summarise their TOTAL number of sea star catches. 
data_inv<- data_inv %>% select(year, time,depth, number) %>% group_by(year,time,depth) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)
data_inv<- data_inv %>% select(year, time,depth, number) %>% group_by(year,time,depth) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)

### Combine species-specific counts with invert dataframe. 
data_inv<- left_join(data_inv, MA)
data_inv<- left_join(data_inv, AP)
data_inv<- left_join(data_inv, CP)
data_inv<- left_join(data_inv, DI)
data_inv<- left_join(data_inv, HS)
data_inv<- left_join(data_inv, LF)
data_inv<- left_join(data_inv, M)
data_inv<- left_join(data_inv, PB)
data_inv<- left_join(data_inv, PH)
data_inv<- left_join(data_inv, SD)
data_inv<- left_join(data_inv, SS)

### Replace NAs with 0s & include sum of sea star catch per trawl. 
data_inv<- data_inv %>% mutate(M.aequalis=replace_na(M.aequalis,0))%>% 
  mutate(A.pugetana=replace_na(A.pugetana,0))%>% mutate(C.papposus=replace_na(C.papposus,0))%>%
  mutate(D.imbricata=replace_na(D.imbricata,0))%>% mutate(H.spinosa=replace_na(H.spinosa,0))%>% 
  mutate(L.foliolata=replace_na(L.foliolata,0))%>% mutate(H.leviuscula=replace_na(H.leviuscula,0))%>%
  mutate(P.brevispinus=replace_na(P.brevispinus,0))%>% 
  mutate(P.helianthoides=replace_na(P.helianthoides,0))%>% mutate(S.dawsoni=replace_na(S.dawsoni,0))%>% 
  mutate(S.stimpsoni=replace_na(S.stimpsoni,0))
data_inv<- data_inv %>% 
  select(year,time,depth,M.aequalis,A.pugetana,C.papposus,D.imbricata,H.spinosa,L.foliolata,H.leviuscula,P.brevispinus,P.helianthoides,S.dawsoni,S.stimpsoni) %>% 
  group_by(year,time,depth) %>% 
  mutate(total_catch= sum(M.aequalis,A.pugetana ,C.papposus ,D.imbricata ,H.spinosa, L.foliolata ,H.leviuscula ,P.brevispinus, P.helianthoides ,S.dawsoni ,S.stimpsoni)) 

### Merge the tables (replace NAs in "number" with 0). Convert date column to year format. 
dataYears <-   format(data$date, format="%Y")
data <- data %>% select(-date) %>% mutate(year=dataYears)

data<- full_join(data,data_inv)
data <- data %>% arrange(year)  %>% relocate(year) # Bring year to the front. 
data<- data %>% mutate_all(~replace(., is.na(.), 0))

# Invertebrate data shows that 4 sea stars were caught in 2013, 
#but trawling data only recorded afternoon and evening trawls. 
# This is also the case for the 2005_70m_afternoon trawl. 
#Assume they were uniform trawls (370 m distance). 
data[140,4]<- 370 # 2005 trawl row 
data[140,5]<- 1295
data[279:282,4] <- 370 # 2013 trawl row 
data[279:282,5] <- 1295

### Add in pre (<2014) and post (>2014) designations.
data<- data %>% mutate_all(~replace(., is.na(.), 0))
data<- data %>% mutate (total_density = (total_catch/trawlarea)) %>% 
  mutate(prepost = if_else(year<2014, "pre", "post",missing=NULL))

# Merge actualDepth (continuous depth) with main dataframe for density analysis at depth.
actualDepth<- actualDepth %>% separate(trawl_ID, into=c("year", "month","day" ,"time"))
actualDepth <- actualDepth %>% select(-c(month,day))
data <- left_join(data, actualDepth, by = c("year","time","depth"))
data <- data %>% rename(est_depth=depth) %>% rename(actual_depth=avg_depth)

#Remove NAs from data$actual_depth 
data$est_depth <- as.double(data$est_depth)
for(i in 1:length(data$actual_depth)){
  if(is.na(data[i,20]) == TRUE){
    data[i,20] <- data[i,2]
  }
}

#######====ADD 1997 DATA: INVERTS SPREADSHEET====== 
# Including additional sea star catch data from 1997

# formatting 1997 data (column headings, convert NA to 0, numeric formatting, etc.)
data_1997 <- read.csv("data/invertebrates_dataset_1997.csv")
data_1997 <- data_1997 %>% rename(col1 = Counts.of.macroinvertebrates.from.spring.1997.Port.Madison.survey) 
data_1997 <- data_1997[c(2,6:10),1:21]
data_1997<- t(data_1997)
colnames(data_1997) <- data_1997[1,]
data_1997 <- data_1997[-1,]
data_1997 <- as.data.frame(data_1997)
data_1997[is.na(data_1997)] = 0
data_1997[,2:6] <- lapply(data_1997[,2:6], as.numeric)

# Include time designations, column names, etc..
data_1997 <- data_1997 %>% mutate(time = c( rep("afternoon", times=4), rep("evening", times=4), 
             rep("night", times=4),   rep("early morning", times=4), rep("morning", times=4))) %>% 
  mutate(year = rep(1997, times=20))%>% mutate(distance = rep(370, times=20)) %>% 
  mutate(trawlarea = rep(1295, times=20)) %>% mutate(prepost = rep("pre", times=20))
colnames(data_1997) <- c("est_depth", "S.stimpsoni", "M.aequalis", "C.papposus", "P.helianthoides", "E.troschelli", "time", "year", "distance",  "trawlarea",  "prepost")
data_1997$year <- as.character(data_1997$year )
data_1997$est_depth <- as.numeric(data_1997$est_depth)

# Join 1997 data with main dataframe.  
data<- full_join(data_1997, data)
data[is.na(data)] = 0
data[1:20,21] <- data[1:20,1] 
# putting standard (discrete) depths in the actual 
#depth column for 1997, no continuous depths. 

# Deleting duplicate rows in 2005 (only 4 sea stars were caught 2005_10m)
data <- data[-c(148,150,151),]

#remove and recalculate total catch and density for each trawl 
#(must repeat calculations because 1997 data was added after original calculations).
data <- data %>% select(-c(total_catch, total_density))
data <- data %>% group_by(year,time,est_depth) %>% 
  mutate(total_catch= sum(E.troschelli,M.aequalis,A.pugetana ,C.papposus ,D.imbricata ,
                          H.spinosa, L.foliolata ,H.leviuscula ,P.brevispinus, 
                          P.helianthoides ,S.dawsoni ,S.stimpsoni)) %>%
  mutate(total_catch_HM = sum(S.stimpsoni, P.helianthoides, P.brevispinus, S.dawsoni, E.troschelli)) %>%    #species with high mortality
  mutate(total_catch_SM = sum(D.imbricata, H.leviuscula)) %>%      #species with noticeable mortality
  mutate(total_catch_UM = sum(M.aequalis, L.foliolata, C.papposus, H.spinosa))     #species with likely/unknown mortality
data<- data %>% mutate (total_density = (total_catch/trawlarea)) %>%
  mutate (total_density_HM = (total_catch_HM/trawlarea)) %>%
  mutate (total_density_SM = (total_catch_SM/trawlarea)) %>%
  mutate (total_density_UM = (total_catch_UM/trawlarea))

data <- data %>% select(-A.pugetana) 
# lastly, lets cut this species (brittle star, already set its counts to 0)

#####===AVG SEA STAR ABUNDANCE TABLE==========

speciesAvg<- as.data.frame(Map(c, lapply(data[data$est_depth == 10, c(2:6,12:17)], mean),
      lapply(data[data$est_depth == 25, c(2:6,12:17)], mean),
    lapply(data[data$est_depth == 50, c(2:6,12:17)], mean),
    lapply(data[data$est_depth == 70, c(2:6,12:17)], mean)))
speciesAvg <- t(speciesAvg)

# Sort by which species has highest density 
speciesAvg  <- as.data.frame(speciesAvg)
speciesAvg <- speciesAvg %>% mutate(Species=row.names(speciesAvg)) %>% relocate(Species) 
minus_sp <- speciesAvg[,-1]
rowmeans <- rowMeans(minus_sp)
speciesAvg <- speciesAvg %>% mutate(means=rowmeans) 
speciesAvg <- speciesAvg %>%  arrange(desc(means)) %>% select(-means)

###=====Temp Data from WA Dept of Ecology CTD (currently using) ====########

# Link:  https://apps.ecology.wa.gov/eim/search/Eim/EIMSearchResults.aspx?ResultType=EIMTabs&StudySystemIds=99970619&StudySystemIds=99970618&StudyUserIdSearchType=Equals&StudyUserIds=MarineWater-P&StudyUserIds=MarineWater
# Location PSB003 (Lat: 47.66001 Long: -122.4417); 1999-2017

full.CTD.data.raw <- read.csv("data/WADeptEcology_CTD.csv")
full.CTD.data <- full.CTD.data.raw %>% select("Field_Collection_Start_Date_Time",
                                              "Field_Collection_Upper_Depth",
                                              "Field_Collection_Lower_Depth",
                                              "Reslt_Parameter_Name",
                                              "Result_Value",
                                              "Result_Unit",
                                              "Base_Parameter_Name",
                                              "Base_Result_Method",
                                              "Base_Result_Method_Description",
                                              "Statistical_Basis",
                                              "Sample_Size",
                                              "Calculated_Latitude_Decimal_Degrees_NAD83HARN",
                                              "Calculated_Longitude_Decimal_Degrees_NAD83HARN")

temp.CTD.1 <- full.CTD.data %>% 
  filter(Reslt_Parameter_Name == "Temperature, water (profile median)") %>% 
  select(-c(Reslt_Parameter_Name,Result_Unit, Base_Parameter_Name,Base_Result_Method,
            Base_Result_Method_Description, Statistical_Basis, Calculated_Latitude_Decimal_Degrees_NAD83HARN,
            Calculated_Longitude_Decimal_Degrees_NAD83HARN,Sample_Size))

# Range of Collection depths
c(min(temp.CTD.1$Field_Collection_Upper_Depth), max(temp.CTD.1$Field_Collection_Upper_Depth)) # upper: (0.5 - 2 m)
c(min(temp.CTD.1$Field_Collection_Lower_Depth), max(temp.CTD.1$Field_Collection_Lower_Depth)) # lower: (10 - 97 m)
# hist(temp.CTD$Field_Collection_Lower_Depth, xlab="Depth (m)", main="", col="darkgreen")

# Subset only CTD entries with maximum depths between 30 and 50 m
temp.CTD <- temp.CTD.1 %>% filter(Field_Collection_Lower_Depth > 20 &
                                  Field_Collection_Lower_Depth < 60 )

# Remove am/pm & convert to date class  
temp.CTD <- temp.CTD %>% mutate(Field_Collection_Start_Date_Time = stri_sub(Field_Collection_Start_Date_Time, from=1, to=-4)) 
temp.CTD$Field_Collection_Start_Date_Time <- as.POSIXct(temp.CTD$Field_Collection_Start_Date_Time, format = "%m/%d/%Y %H:%M:%S")

# Extract year and month, take the average temp of April-June (to capture spring variability)
# 1999-2017
temp.CTD.values <- temp.CTD %>% 
  mutate(year = format(temp.CTD$Field_Collection_Start_Date_Time, format = "%Y")) %>% 
  mutate(month = format(temp.CTD$Field_Collection_Start_Date_Time, format = "%m")) %>% 
  rename(temp = "Result_Value") %>% 
  select(-c(Field_Collection_Start_Date_Time, Field_Collection_Upper_Depth,Field_Collection_Lower_Depth)) %>% 
  filter( month > "02" & month < "06") %>%
  group_by(year) %>% summarise(meanTemp = mean(temp))
temp.CTD.values

###===== Finally, put environmental data into main dataset ========

# Final dataframe, with temp 
dat_all<-left_join(data,temp.CTD.values,by="year")
return(dat_all)
}
