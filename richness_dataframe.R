source("./R/db_connect.R")

# Initialise connections
setwd("~/Dropbox/WorkPersonal_transfers/PhD/Kautikeino/Compendium")


###################################
#### Step 1: Compile dataframe ----
###################################

# Compile location, event, taxon and occurrence data
location <- tbl(nofa_db,"location")
event <- tbl(nofa_db,"event")
l_taxon <- tbl(nofa_db,"l_taxon")
occurrence <- tbl(nofa_db,"occurrence")


df_base <- left_join(occurrence,event,by="eventID")           # Joins occurrence data to individual events
df_base <- left_join(df_base,location,by="locationID")        # Joins events to locations
df_base <- left_join(df_base,l_taxon,by="taxonID")            # Joins taxon to occurrences
df <- df_base %>% filter(datasetID=="Kautokeino") %>% collect()   # Collects everything

eurolst <- tbl(env_db,"location_EuroLST_BioClim") %>%         # Collects temperature data
  dplyr::select(c(1:2,5,8)) %>%
  collect()
df <- left_join(df,eurolst,by="locationID")                   # Adds temp data to occurrences


ecco <- tbl(env_db,"ecco_biwa_lakes_v_0_1") %>%               # Adds geomorph data to dataset 
  dplyr::select(vatn_lnr,areas32n,perim32n,area,perimtot) %>% 
  filter(vatn_lnr %in% df$no_vatn_lnr) %>%
  collect()
df <- left_join(df,rename(ecco,no_vatn_lnr=vatn_lnr),by="no_vatn_lnr")
df$no_vatn_lnr <- as.factor(df$no_vatn_lnr)



#########################################
### Step 2: Turn data to wide form -----
#########################################

# Produce better occurrence row
species_rich<- data.frame(df %>% 
  dplyr::select(eventID,vernacularName,occurrenceStatus))
species_rich$occurYESNO<-ifelse(species_rich$occurrenceStatus=="present",1,0)
species_rich$occurrenceStatus<-NULL

# The following gives us a data frame with our different species in different columns
# We can do this by eventID, as all our data comes from one year
speciesRNN <- reshape(species_rich,direction="wide", 
                      idvar="eventID", timevar="vernacularName")

# Produce species richness data for each lake
speciesRNN$richness<-rowSums(speciesRNN[,2:8],na.rm=TRUE)

# Turn occurrence into a factor, rename columns for simplicity
for (i in 2:8)
{
  speciesRNN[,i]<-factor(speciesRNN[,i],levels=c("0","1"))
}
names(speciesRNN) <- gsub("occurYESNO.","",names(speciesRNN))
names(speciesRNN) <- gsub(" ","",names(speciesRNN))


##########################################################
### Step 3: Merge richness matrix with habitat data ------
##########################################################

# Consolidate habitat data
habitat_var <- df %>%
  dplyr::select("eventID","locationID","decimalLatitude","maximumElevationInMeters",
                     "eurolst_bio01","eurolst_bio05","eurolst_bio10","area","perimtot") %>%
  rename(elevation = maximumElevationInMeters)
habitat_var <- habitat_var[!duplicated(habitat_var),]


 # Merge with richness dataframe
richness_all<-merge(speciesRNN,habitat_var,
                    by="eventID",all.x=TRUE)

# Create shoreline complexity variable
richness_all <- richness_all %>%
  filter(!is.na(area) & !is.na(perimtot))
shoreline_mod <- lm(data=richness_all,log(area) ~ log(perimtot))
summary(shoreline_mod)
richness_all$shoreline_complex <- residuals(shoreline_mod)

# Save data
save(richness_all,file="./Data/richness_all.rda")
load("./Data/richness_all.rda")
richness_all[,2:8]
names(richness_all[,2:8])
summary(richness_all[,2:8])
