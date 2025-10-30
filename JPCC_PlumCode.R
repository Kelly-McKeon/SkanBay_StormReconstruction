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
     mainpanel.tickfontsize=1.5,
     xaxt="n", yaxt="n")  # Suppress default axis to add custom ticks


#### Change Tick Marks so Easier to Edit the Figure Visuals in Illustrator Because I don't know how to in R

# Add custom x-axis with 1000-year tick marks
# Get the age range from the current plot
age_range <- par("usr")[1:2]  # Get x-axis limits from the plot
# Force specific 1000-year intervals - adjust range as needed for your data
# You may need to modify these values based on your actual age range
major_ticks <- seq(from = 0, to = 20000, by = 500)
# Filter ticks to only those within the actual plot range
major_ticks <- major_ticks[major_ticks >= age_range[1] & major_ticks <= age_range[2]]
# Create minor tick marks every 500 years (optional)
minor_ticks <- seq(from = 0, to = 20000, by = 100)
minor_ticks <- minor_ticks[minor_ticks >= age_range[1] & minor_ticks <= age_range[2]]
# Force the axis to use ALL our tick marks
axis(1, at = major_ticks, labels = major_ticks,
     cex.axis = 1.2,  # Slightly smaller text to fit more labels
     las = 2,         # Rotate labels vertically
     mgp = c(3, 0.5, 0))  # Adjust label positioning
# Add minor tick marks without labels
axis(1, at = minor_ticks, labels = FALSE, tcl = -0.25)
# If labels are still too crowded, try alternating labels:
# alt_labels <- ifelse(major_ticks %% 2000 == 0, major_ticks, "")
# axis(1, at = major_ticks, labels = alt_labels, cex.axis = 1.2, las = 2)
# Optional: Add vertical grid lines at major tick marks
abline(v = major_ticks, col = "gray90", lty = 1, lwd = 0.5)

# Add custom y-axis with tick marks every 1000 years
# Get the depth range from the current plot
depth_range <- par("usr")[3:4]  # Get y-axis limits from the plot
# Create major tick marks every 1000 years (assuming your depths are in years)
# If your depths are in cm, you'll want different intervals - adjust accordingly
major_ticks <- seq(from = ceiling(depth_range[1]/1000)*1000,
                   to = floor(depth_range[2]/1000)*1000,
                   by = 1000)
# Create minor tick marks every 500 years (optional)
minor_ticks <- seq(from = ceiling(depth_range[1]/500)*500,
                   to = floor(depth_range[2]/500)*500,
                   by = 500)
# Add major tick marks with labels
axis(2, at = major_ticks, labels = major_ticks,
     cex.axis = 1.5,  # Match your mainpanel.tickfontsize
     las = 1)         # Horizontal labels
# Add minor tick marks without labels
axis(2, at = minor_ticks, labels = FALSE, tcl = -0.25)
# Optional: Add horizontal grid lines at major tick marks
abline(h = major_ticks, col = "gray90", lty = 1, lwd = 0.5)
# Optional: Add vertical grid lines at major tick marks
abline(v = major_ticks, col = "gray90", lty = 1, lwd = 0.5)


#### Read in Event Depths From Python output
#### Using 25 pt mov median 85th percentile anomaly

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
write.csv(ens, file = "/Users/kellymckeon/Library/CloudStorage/OneDrive-WoodsHoleOceanographicInstitution/Documents/WHOI/Alaska/Codes/Coherent_Codes/JPC05/For_Publication/core_data/JPC05_Age_Distributions.csv")
#on a pc
write.csv(ens, file = "C:/Users/kelly/OneDrive - Woods Hole Oceanographic Institution/Documents/WHOI/Alaska/Codes/Coherent_Codes/JPC05/core_data/JPC05_Age_Distributions.csv")

