## To select 30 lines to represent NA germpalsm diversity
# Decide to pick from the recent 3 years lines 
# line data is downloaded from Horizon https://breeding.ag/horizonclassic/#!/

library(tidyverse)

# load ancestral Haplotype data
pbh = readRDS("/Users/taozuo/Downloads/ancestralHaplo_20230303.rds")

setwd("/Users/taozuo/Library/CloudStorage/OneDrive-Bayer/Project/2023/misc/Graziana/diverseLines")

line_2020 = read.csv('export_2020.csv')
line_2020 = line_2020 %>% filter(Line.Name %in% colnames(pbh),Trait.Name == "CONV",str_detect(MARKET_ME,"RM")) %>% 
  select(Line.Name,RM,GENDER,YEAR_CODED,MARKET_ME,ORIGIN)
colnames(line_2020)[1] = "lineName"

line_2021 = read.csv('export_2021.csv')
line_2021 = line_2021 %>% filter(Line.Name %in% colnames(pbh),Trait.Name == "CONV",str_detect(MARKET_ME,"RM")) %>% 
  select(Line.Name,RM,GENDER,YEAR_CODED,MARKET_ME,ORIGIN)
colnames(line_2021)[1] = "lineName"

line_2022 = read.csv('export_2022.csv')
line_2022 = line_2022 %>% filter(Line.Name %in% colnames(pbh),Trait.Name == "CONV",str_detect(MARKET_ME,"RM")) %>% 
  select(Line.Name,RM,GENDER,YEAR_CODED,MARKET_ME,ORIGIN)
colnames(line_2022)[1] = "lineName"

# combine all three year's data together
lineData = unique(bind_rows(line_2020,line_2021,line_2022))
table(lineData$GENDER,lineData$MARKET_ME)
#RM075 RM080 RM085 RM090 RM095 RM100 RM105 RM110 RM115 RM120
#F    80   105   107   135   225   282   366   536   813   187
#M    76   103   106   134   240   299   338   522   845   185

# decide to select 1 from each RM and Gender, then 3 more for RM115, 1 more for RM110 and 1 more for RM105 because of the size
# Reasons: 1) Better representative for RMs; RN105, RM110 and RM115 have much more lines (about 60%), making the selection based toward these three RMs
# Reason 2: To select 15 liens from 3000 takes too much time to complete. 


# subset the ancHap
pbh = pbh[,lineData$lineName]

# select 1 per RM and gender for all RMs
allRMs = sort(unique(lineData$MARKET_ME))
selLineData = data.frame()
genders = c("F","M")
for(rm in allRMs){
  for(gender in genders){
    subLines = subset(lineData,MARKET_ME == rm & GENDER == gender )$lineName
    print(paste(rm,gender,length(subLines),Sys.time()))
    popPwSim = getPwHapSim(subLines,pbh)
    avgSim = apply(popPwSim,2,function(x){return(mean(x,na.rm=T))})
    selLineData = bind_rows(selLineData,data.frame(lineName = names(which.max(avgSim))[1],MARKET_ME = rm,GENDER = gender))
  }
}


# one more for RM105 and RM110
for(rm in c("RM105","RM110")){
  for(gender in genders){
    subLines = subset(lineData,MARKET_ME == rm & GENDER == gender )$lineName
    print(paste(rm,gender,length(subLines),Sys.time()))
    selLine = subset(selLineData,MARKET_ME == rm & GENDER == gender)$lineName
    tmpPbh = pbh[,subLines]
    tmpSelLinePbh = pbh[,selLine]
    for(i in 1:ncol(tmpPbh)){
      idx = which(tmpPbh[,i] == tmpSelLinePbh)
      tmpPbh[idx,i] = NA
    }
    popPwSim = getPwHapSim(subLines,tmpPbh)
    avgSim = apply(popPwSim,2,function(x){return(mean(x,na.rm=T))})
    selLineData = bind_rows(selLineData,data.frame(lineName = names(which.max(avgSim))[1],MARKET_ME = rm,GENDER = gender))
  }
}

# three more on RM115
rm = "RM115"

for(gender in genders){
  subLines = subset(lineData,MARKET_ME == rm & GENDER == gender )$lineName
  selLine = subset(selLineData,MARKET_ME == rm & GENDER == gender)$lineName
  tmpPbh = pbh[,subLines]
  for(rep in 1:3){
    print(paste(rep,rm,gender,length(subLines),Sys.time()))
    tmpSelLinePbh = pbh[,selLine]
    for(i in 1:ncol(tmpPbh)){
      idx = which(tmpPbh[,i] == tmpSelLinePbh)
      if(length(idx) > 0){
        tmpPbh[idx,i] = NA
      }
    }
    popPwSim = getPwHapSim(subLines,tmpPbh)
    avgSim = apply(popPwSim,2,function(x){return(mean(x,na.rm=T))})
    selLineData = bind_rows(selLineData,data.frame(lineName = names(which.max(avgSim))[1],MARKET_ME = rm,GENDER = gender))
    selLine = names(which.max(avgSim))[1]
  }

}

write.csv(selLineData,"selLineData.csv",row.names=F)


# pairwise similarities
getPwHapSim = function(subLines,pbh,numOfBins = 17116){
  if(!all(subLines %in% colnames(pbh))){stop("Not all lines have PBH data")}
  simOut = matrix(NA, nrow=length(subLines),ncol=length(subLines),dimnames = list(subLines,subLines))
  subPbh = pbh[,subLines]
  for(i in 1:(length(subLines)-1)){
    for(j in 2:length(subLines)){
      tmpSim = round(table(subPbh[,i] == subPbh[,j])["TRUE"] / numOfBins,2)
      simOut[i,j] = simOut[j,i] = tmpSim
    }
  }
  return(simOut)
}



