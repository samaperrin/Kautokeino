## Import distance data
setwd("~/Dropbox/WorkPersonal_transfers/PhD/Kautikeino/Compendium")
distances <- read.csv("./Data/location_matrix.csv",header=T,sep=";")

# Get rid of all columns that have too many observations
library(dplyr)
dcount <- distances %>% group_by(locationID) %>%
  count %>%
  filter(n==318)
distances.clean <- distances %>%
  filter(locationID %in% dcount$locationID,
         locationID.1 %in% dcount$locationID)

# Generate list of locations used
locations_distance <- distinct(distances.clean,locationID)

# Generate distance matrix
distance_matrix <- data.frame(locations_distance)
for (i in locations_distance$locationID) {
  c <- distances.clean %>%
    filter(locationID==i) %>%
    dplyr::select(2,3)
  colnames(c) <- c("locationID",i)
  distance_matrix <- merge(distance_matrix,c,by="locationID",all=TRUE) 
}

# Give row names, get rid of locationID column
rownames(distance_matrix) <- distance_matrix$locationID
distance_matrix <- distance_matrix[-1]


distance_matrix <- distance_matrix %>%
  dplyr::select(rownames(distance_matrix))

distance_matrix[is.na(distance_matrix)] <- 0
save(distance_matrix,file="distance_matrix.rda")

