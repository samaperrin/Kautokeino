library(dplyr)

connect.species <- function(species, connectivity_matrix, locationIDs) {
  step1 <- richness_all[,c(species,"locationID","eurolst_bio10")]                 # Isolate occurrence and location column
  colnames(step1) <- c("occurrence", "locationID","temp")
  step2 <- merge(step1,locationIDs,by="locationID",all.x=TRUE)              # Merge with lakeid
  step3 <- step2[,c(2:4)]
  
  step4 <- connectivity_matrix %>%                                          # Whittle down to downstream leakes with
    filter(lakeid %in% filter(step3,occurrence == "1")$waterBodyID,         # fish present
           upstream_lake %in% filter(step3,occurrence != "NA")$waterBodyID)
  step5  <- merge(step4,rename(step3,upstream_lake = waterBodyID),by="upstream_lake",all.x=TRUE)
  # Introduce upstream lake
  colnames(step5)  <- c("upstream_lake", "upstream_slope_max","downstream_lake","occurrenceStatus","temp")
  step5 <- distinct(step5,upstream_lake,upstream_slope_max,downstream_lake,occurrenceStatus,temp)
  # Get rid of any rows we have conflicting                                                                             # duplicates for
  step6 <- merge(step5,tally(group_by(step5,upstream_lake,downstream_lake)),
                 by=c("upstream_lake","downstream_lake"),all.x=TRUE)
  step6$occurs  <- ifelse(step6$occurrenceStatus == "1",1,ifelse(step6$n > 1,1,0))
  step6 <- distinct(step6,upstream_lake,upstream_slope_max,downstream_lake,occurs,temp)
  return(step6)
}

