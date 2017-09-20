setwd("~/Dropbox/WorkPersonal_transfers/PhD/Kautikeino/Compendium")
source("./R/db_connect.R")
load("./Data/richness_all.rda")

#############################
#### Step 1: Matching ID ----
#############################

### We need to create a table that can match locationID with lakeID

# Produce a list only for lakes we have occurrence data for 
locationIDs_con <- tbl(nofa_db, "location")
locationIDs <- locationIDs_con %>%
  filter(locationID %in% richness_all$locationID) %>%
  select(locationID,waterBodyID) %>%
  collect()
save(locationIDs,file="./Data/locationIDs.rda")

# And a second list for every lake
locationIDs_all <- locationIDs_con %>%
  select(locationID,waterBodyID) %>%
  filter(!is.na(waterBodyID)) %>%
  collect()


#######################################
#### Step 2: Get connectivity data ----
#######################################

# The following data is going to be used to isolate all connections between lakes we have fish data for,
# and to identify which of the lakes we have occurrence data for are connected.

# THe following returnns every lakeID and all of its upstream lakes
connectivity <- 'SELECT * FROM (SELECT lakeid, CAST(unnest(string_to_array(upstream_lakes, \', \')) AS integer) AS upstream_lake FROM temporary.connectivity_14695 WHERE first_upstream_lakes != \'\') AS x'
connectivity <- get_postgis_query(con,connectivity)

# Then we whittle this down so that all lakes have occurrence data
connectivity <- connectivity %>%
  filter(lakeid %in% locationIDs$waterBodyID
         & upstream_lake %in% locationIDs$waterBodyID) %>%
  mutate(ConID = paste(upstream_lake,"downto",lakeid,sep="-"))
# So 'connectivity' holds a list of connections between lakes for which we have fish data. From now on Iæll simply refer to any multi-stream connection between two lakes with occurrence data as a Connection.
# 'connectivity also contains a column ConID that we'll use later to categorise all Connections.
         

### Now we get links from all lakes to their lakes directly upstream
ind_connect_con <- 'SELECT * FROM (SELECT lakeid, CAST(unnest(string_to_array(first_upstream_lakes, \', \'))
AS integer) AS upstream_lake FROM temporary.connectivity_14695 WHERE first_upstream_lakes != \'\') AS x'
ind_connect <- get_postgis_query(con,ind_connect_con)
# So 'ind_connect' holds all the individual links between lakes

######################################################
#### Step 3: Link occurrence lakes to each other -----
######################################################

# We now create a series of data frames within one larger data frame which detail every link within each Connection.

# Create a data frame for all this to start off with
big_R <- data.frame("ConID"=integer(0),"Upstream_Lake"=integer(0),"Downstream_Lake"=integer(0))

# ConID is used to demarcate different Connections
for(i in 1:length(connectivity$ConID))              # first produce a loop that cycles through all the con IDs
{
  # So now that we know we're looking at the first Connection,
  # we start with the upstream lake, and identify the lake immediately downstream
  x <- connectivity[i,3]                                   # connection id
  y <- connectivity[i,2]                                   # upstream lake
  z <- filter(ind_connect,upstream_lake == y)[,1]             # lake directly downstream
  # This link is then put into a smaller data frame
  small_R <- data.frame("ConID"=integer(0),"Upstream_Lake"=integer(0),"Downstream_Lake"=integer(0))
  a <- 1
  repeat
  {
    r <- c(x,y,z)  
    small_R[a,] <- r
    # Now we keep cycling through, going further downstream. And if our downstream lake finally reachers the end of the connection, the cycle breaks and we move onto the next connection.
    if (z == connectivity[i,1]) break
    # above basically says if our downstream lake equals our initial downstream lake stop looping
    # if not, we need to go back to the start
    # x stays the same, it's still the connection ID
    y <- z    # our downstream lake now becomes our upstream lake
    z <- filter(ind_connect,upstream_lake == y)[,1]      # we now start looking for the next downstream lake
    a <- a + 1
  } 
  big_R <- rbind(big_R,small_R)
  # This adds the table from our recently finished Connection to the larger table
}

###############################################
#### Step 4: Detect and delete duplicates -----
###############################################

# There are larger Connections that are unnecessary, because there are smaller Connections within them. For instance if we go from lake A to lake B, which then connects to lake C, we only need the AB and BC connections, not the AC connection, as it's less informative.
no_conex <- tally(group_by(big_R,ConID))
connections2 <- arrange(merge(connectivity,no_conex,by="ConID",all.x=TRUE),upstream_lake,n)
# So we've created a tally from the number of links that make up every Connection. Seeing as each lake can only have on downstream lake, we can delete larger connections by using the upstream lakes. If a lake is listed more than once as an upstream lake, we can simply take the Connection with the fewest links.

unique_conex <- summarize(group_by(connections2,upstream_lake),
                          minC = min(n))    # takes the connections with minimum
connections3 <- merge(unique_conex,rename(connections2,minC = n),by=c("upstream_lake","minC"),all.x=TRUE)

big_Rclean <- filter(big_R,ConID %in% connections3$ConID)

################################
#### Step 5: Import slopes -----
################################

slopes_con <-'SELECT * FROM (SELECT lakeid, CAST(unnest(string_to_array(first_upstream_lakes, \', \')) AS integer) 
AS upstream_lake, CAST(unnest(string_to_array(first_upstream_lakes_slope_mean, \', \')) AS numeric) 
AS upstream_lake_slope_mean, CAST(unnest(string_to_array(first_upstream_lakes_slope_max, \', \')) AS numeric) 
AS upstream_lake_slope_max FROM temporary.connectivity_14695 WHERE first_upstream_lakes != \'\') AS x'
slopes <- get_postgis_query(con,slopes_con)

# So for now we have brought in just slope max and mean

# Filter this down to just the connections we have fish for
stats <- filter(slopes,
                lakeid %in% big_Rclean$Downstream_Lake,
                upstream_lake %in% big_Rclean$Upstream_Lake)
head(stats)

# And merge to match slope with links.
almost_there <- arrange(merge(stats,rename(big_Rclean,lakeid = Downstream_Lake, upstream_lake = Upstream_Lake),
                              by=c("lakeid","upstream_lake")),ConID)
head(almost_there)



########################################################################
#### Step 6: Apply functions to create max slope over Connections  -----
########################################################################

final  <- almost_there %>%
  group_by(ConID) %>%
  summarize(max(upstream_lake_slope_max))
head(final)
dim(final)

connectivity_matrix  <- select(merge(final,connectivity,by = "ConID", all.x=TRUE), 2:4)
head(connectivity_matrix)
save(connectivity_matrix,file="./Data/connectivity_matrix.rda")




##############################################
#### Step 7: Add in fish absence/presence ----
##############################################

for (i in c("Burbot","Arcticcharr","Whitefish","Pike","Perch","Browntrout","Grayling"))
{
  step1 <- richness_all[,c(i,"locationID","eurolst_bio10")]                 # Isolate occurrence and location column
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
  assign(paste(strtrim(i,4),"D",sep=""),step6)
}

