#Using Plum for modeling
install.packages("rplum")
library(rplum)
library(rbacon)

#### JPC05 Composite with 63 um slump threshold
#### slumps are determined in python 
Plum(ask=FALSE,
     core="JPC05_Composite",
     date.sample=2022,
     ra.case=0,
     n.supp=1,
     otherdates="JPC05_Composite_C14.csv",
     d.max=2055,
     acc.mean=c(2, 5, 10),
     boundary=c(580, 1600),
     mem.strength=10,
     ssize=5000,
     thick=10,
     slump=c(46,57, 252,253, 381,382, 399,403, 419,423, 425,427, 444,445, 447,451, 518,561, 565,580,
             585,588, 616,618, 626,632, 637,638, 739,742, 744,745, 814,815, 851,853,
             905,908, 936,938, 986,992, 1136,1137, 1142,1145, 1267,1279,
             1338,1339, 1391,1392, 1415,1479, 1655,1669, 1701,1703, 1796,1802, 1961,1976, 2030,2031),
     cc=1,
     slump.col=grey(0.9),
     mainpanel.tickfontsize=1.5)  


#### Read in Event Depths From Python output

#on a mac
event_depths = read.csv("~/Library/CloudStorage/OneDrive-WoodsHoleOceanographicInstitution/Documents/WHOI/Alaska/Codes/Coherent_Codes/JPC05/For_Publication/core_data/event_depth.csv",
                        header = FALSE)

#on a pc 
event_depths = read.csv("C:/Users/kelly/OneDrive - Woods Hole Oceanographic Institution/Documents/WHOI/Alaska/Codes/Coherent_Codes/JPC05/For_Publication/core_data/event_depth.csv", header=FALSE)


#### Get Age Distributions at Event Depths
x = as.numeric(event_depths[,1])
ens = matrix(data=NA, nrow=5000, ncol=length(x)) #make empty matrix, use 4K because that's the size bacon outputs for the distributions
count = 1
for (val in x){
  event = Bacon.Age.d(val) #bacon command for the age distributions
  m = matrix(event)
  ens[,count] = m
  count = count+1
}


#### Save Age Distributions of Event Beds
#on a mac
write.csv(ens, file = "/Users/kellymckeon/Library/CloudStorage/OneDrive-WoodsHoleOceanographicInstitution/Documents/WHOI/Alaska/Codes/Coherent_Codes/JPC05/For_Publication/core_data/JPC05C_Age_Distributions.csv")
#on a pc
write.csv(ens, file = "C:/Users/kelly/OneDrive - Woods Hole Oceanographic Institution/Documents/WHOI/Alaska/Codes/Coherent_Codes/JPC05/core_data/JPC05C_Age_Distributions.csv")

