SUMMARY
This project is designed to look at the effects of changes in environmental variables and connectivity on fish species composition in the Kautokeino lakes in Northern Norway.

SCRIPTS

connectivity_matrix.R
	This script constructs a connectivity matrix based on geological and topgraphic data. The output is a
	list of connected lakes and some data concerning information on lake connections both up and downstream.

dispersal_barriers.R
	This script merges our fish absence/presence data with connectivity measures derived from the above matrix,
	with the final output being a table indicating whether or not a species is absnet/present upstream
	given a maximum upstream slope.

richness_dataframe.R
	This script produces a wide-form dataframe showing species presence/absence and overall species richness
	for 325 lakes based  on a 1981 survey against a variety of environmental variables.

distance_matrix.R
  This script takes the distances between each lake with fish data and produces a series of coordinates to allow 
  us to analyse connectivity between lakes.

Kautokeino_JSDM.R
	This script runs a joint species distribution model as based on Pollock et al.'s 2014 paper "Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM)"
	
Kautokeino_HMSC.R
	This script produces the HMSC model for the Kautokeino database. It includes simulation for a temperature 
	increase of approximately 2.5 degrees.


SOURCE SCRIPTS

db_connect.R  
	Contains database connectors and relevant libraries.

fit_JSDM.R  
	Contains script for running JSDM
	
connectspecies_function.R  
  Contains function which produces connectivity matrix for a single species. Input is full matrix, locationIDs 
  and species vernacular name.


DATA FILES

connectivity_matrix.rda  
  Contains maximum slope between any lake i and closest downstream lake for both of which we have fish occurrence data.
  
locationIDs.rda  
  Link table which matches lakeIDs used in connectivity matrix to location UUIDs used in occurrence table.
  
richness_all.rda  
  Full occurrence table with climate and hydrological data.
  
location_matrix.csv  
  Shows the distance from every lake for which we have fish data to every other lake for which we have fish data. Length calculated using pgRouting. Lakes not included in length measurement. For more information, ask Jakob Miksch.
  
kauto_hmsc.rda  
  Full HMSC model for the Kautokeino catchment.