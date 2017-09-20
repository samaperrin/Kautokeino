
packages.needed <-
  setdiff(
    c('doMC', 'igraph', 'rgrass7', 'RPostgreSQL','dplyr', 'postGIStools', 'pool', 'dbplyr'),
    rownames(installed.packages())
  )
if(length(packages.needed)) install.packages(packages.needed)


library(postGIStools)
library(RPostgreSQL)
library(dplyr)
library(dbplyr)
library(pool)

#Set connection parameters
pg_drv<-dbDriver("PostgreSQL")
pg_db <- "nofa_sandbox_plugin_2017-08-08"
pg_schema <- "Hydrography"
pg_tmp_schema <- "temporary"
pg_user <- "_"
pg_password <- "_"
pg_host <- "vm-srv-finstad.vm.ntnu.no"


#Initialise connection
con<-dbConnect(pg_drv,dbname=pg_db,user=pg_user, password=pg_password,host=pg_host)

nofa_db <- dbPool(
  drv = RPostgreSQL::PostgreSQL(),
  user="_",
  password = "_",
  host = "vm-srv-finstad.vm.ntnu.no",
  dbname = "nofa_sandbox_plugin_2017-08-08",
  options="-c search_path=nofa" # set db schema from where to look
)

temp_db <- dbPool(
  drv = RPostgreSQL::PostgreSQL(),
  user="_",
  password = "_",
  host = "vm-srv-finstad.vm.ntnu.no",
  dbname = "nofa_sandbox_plugin_2017-08-08",
  options="-c search_path=temporary" # set db schema from where to look
)

env_db <- dbPool(
  drv = RPostgreSQL::PostgreSQL(),
  user="_",
  password = "_",
  host = "vm-srv-finstad.vm.ntnu.no",
  dbname = "nofa",
  options="-c search_path=environmental" # set db schema from where to look
)