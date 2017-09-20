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

kautokeino_JSDM.R
	This script runs a joint species distribution model as based on Pollock et al.'s 2014 paper "Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM)"


DATA FILES

db_connect.R
	Contains database connectors and relevant libraries.

fit_JSDM.R
	Contains script for running JSDM
