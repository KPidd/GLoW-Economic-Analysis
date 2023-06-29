###kp,dp,pb,sb,ct,jm,AA,sg,ab
library(MASS)
library(VGAM)
library(doParallel)
library(parallel)
library(dplyr)
set.seed(429)

#Create the Global options matrix
source("Global Options.R")
#Load in population and data files
source("all_model_files.R")


####build the population
population <-read.csv("Populations/GLOW_Population.csv")
population_clean <- build_population(population, PopulationVariables, GlobalVars)
#Generate an array of common random numbers for each patient for each event for each year
random_numbers <- generate_random(length(population_clean[,"ID"]))
calibration<-read.table("data/rejection_0.05.txt", header=TRUE)
parameter[1:4000,"CALIBRATION_RR_MI"]<-calibration[,3]
parameter[1:4000,"CALIBRATION_RR_STR"]<-calibration[,4]


#rebuild test simulations
set.seed(1)
start.time <- Sys.time()
resultsBL <- run_model(population_clean, 
                            parameter, 
                            50, 
                            "bl", 
                            GlobalVars,
                            random_numbers,
                            LifeTables)
end.time <- Sys.time()
end.time - start.time
saveRDS(resultsBL, file="Results/Testing psa 1944 and 1945 for NA.RData")
write.csv(resultsBL, "Results/Testing Remission/BL - adding remission to met.csv")
rm(resultsBL)

resultsGLOW <- run_model(population_clean, 
                      parameter, 
                      50, 
                      "GLOW", 
                      GlobalVars,
                      random_numbers,
                      LifeTables)
end.time <- Sys.time()
end.time - start.time
#saveRDS(resultsGLOW, file="Results/2000PSA - GLOW - Summary.RData")
write.table(resultsGLOW, "Results/Testing Remission/GLOW- adding remission to met.csv", sep=",")
rm(resultsGLOW)


results3 <- run_model(population_clean, 
                      parameter, 
                      50, 
                      "GLOW_beta_diff", 
                      GlobalVars,
                      random_numbers,
                      LifeTables)
end.time <- Sys.time()
end.time - start.time
saveRDS(results3, file="Results/20000PSA - BETA - Summary.RData")
#write.table(results3, "X:/ScHARR/PR_Amy_Ahern_WLM/General/GLoW modelling/Model/Model KPidd Repository/Results/results -  determ beta.csv", sep=",")











GlobalVars["psa_count", "Value"] <- 29
#rebuild test simulations
set.seed(1)
start.time <- Sys.time()
results1 <- run_model(population_clean, 
                      parameter, 
                      50, 
                      "bl", 
                      GlobalVars,
                      random_numbers,
                      LifeTables)
end.time <- Sys.time()
end.time - start.time

GlobalVars["psa_count", "Value"] <- 30
#rebuild test simulations
set.seed(1)
start.time <- Sys.time()
results1 <- run_model(population_clean, 
                      parameter, 
                      50, 
                      "bl", 
                      GlobalVars,
                      random_numbers,
                      LifeTables)
end.time <- Sys.time()
end.time - start.time
#####crashes here######

G