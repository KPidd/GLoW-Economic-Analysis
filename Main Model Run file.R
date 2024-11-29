#    Sheffield Type 2 Diabetes Treatment Model version 4.0: with the GLoW Trial health economic analysis implemented.
#    Copyright (C) 2023 Pidd, Pollard, Breeze, Bates, Thomas, Mueller, Ahern, Griffin, Brennan

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#    Contact person: Katharine Pidd, Email: k.pidd@sheffield.ac.uk, 
#    Address: Regent Court, 30 Regent Court, Sheffield, United Kingdom, S1 4DA
library(MASS)
library(VGAM)
library(doParallel)
library(parallel)
library(dplyr)
set.seed(1)

#Create the Global options matrix
source("Global Options.R")
#Load in population and data files
source("all_model_files.R")


####read in populations
population <-read.csv("Populations/Model_Population.csv")

###subset populations by Global vars specified subset
subsetpopulation <-function(population){
  if(GlobalVars["Population Subset", "Value"]=="A"){population<-population[1:10000,]}
  if(GlobalVars["Population Subset", "Value"]=="B"){population<-population[10001:20000,]}
  if(GlobalVars["Population Subset", "Value"]=="C"){population<-population[20001:30000,]}
  if(GlobalVars["Population Subset", "Value"]=="D"){population<-population[30001:40000,]}
  if(GlobalVars["Population Subset", "Value"]=="E"){population<-population[40001:50000,]}
  if(GlobalVars["Population Subset", "Value"]=="F"){population<-population[50001:60000,]}
  if(GlobalVars["Population Subset", "Value"]=="G"){population<-population[60001:70000,]}
  if(GlobalVars["Population Subset", "Value"]=="H"){population<-population[70001:80000,]}
  if(GlobalVars["Population Subset", "Value"]=="I"){population<-population[80001:90000,]}
  if(GlobalVars["Population Subset", "Value"]=="J"){population<-population[90001:100000,]}
  return(population)
}
population <- subsetpopulation(population)


####build the population##colClasses=c(rep("NULL", 1),rep("numeric", 4)),
population_clean <- build_population(population, PopulationVariables, GlobalVars)
#Generate an array of common random numbers for each patient for each event for each year
random_numbers <- generate_random(length(population_clean[,"ID"]))
calibration<-read.table("data/rejection_0.05.txt", header=TRUE, row.names=NULL)[,2:5]
parameter[1:4000,"CALIBRATION_RR_MI"]<-calibration[,3]
parameter[1:4000,"CALIBRATION_RR_STR"]<-calibration[,4]
lambda<-20000 


###########################################################################################################################
###########################        Table 2:Main incremental cost-effectiveness analysis      ################################

###########           assumes mixed F2F and online delivery
start.time_Mixed_F2F  <- Sys.time()
resultsDE_Mixed_F2F <- run_model(population_clean, 
                            parameter, 
                            50, 
                            "DE", 
                            GlobalVars,
                            random_numbers,
                            LifeTables)
write.csv(resultsDE_Mixed_F2F, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_Mixed_F2F_populationE.csv")

resultsDEW_Mixed_F2F <- run_model(population_clean, 
                      parameter, 
                      50, 
                      "DEW", 
                      GlobalVars,
                      random_numbers,
                      LifeTables)
write.csv(resultsDEW_Mixed_F2F, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_Mixed_F2F_populationE.csv")

#summarize results
Table2_MIXED_F2F<-Summarize_results(resultsDE_Mixed_F2F,resultsDEW_Mixed_F2F)
write.csv(Table2_MIXED_F2F, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_MIXED_F2F_populationE.csv")
rm(resultsDE_Mixed_F2F)
rm(resultsDEW_Mixed_F2F)
rm(Table2_MIXED_F2F)

end.time_Mixed_F2F  <- Sys.time()
run.time_Mixed_F2F <-end.time_Mixed_F2F  - start.time_Mixed_F2F 
saveRDS(run.time_Mixed_F2F,"Results/Table 2 - Main Incremental Cost-effectiveness Analysis/run time - 20000 people, 2000 psa, both arms.csv")

##################       assumes F2F delivery ONLY
parameter_F2F<-parameter
#replace COST_DESMOND with F2F cost
parameter_F2F[,"COST_DESMOND"]<-264.88

start.time_F2F <- Sys.time()
resultsDE_F2F <- run_model(population_clean, 
                                 parameter_F2F, 
                                 50, 
                                 "DE", 
                                 GlobalVars,
                                 random_numbers,
                                 LifeTables)
write.csv(resultsDE_F2F, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_F2F_populationE.csv")


resultsDEW_F2F <- run_model(population_clean, 
                                   parameter_F2F, 
                                   50, 
                                   "DEW", 
                                   GlobalVars,
                                   random_numbers,
                                   LifeTables)
write.csv(resultsDEW_F2F, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_F2F_populationE.csv")


Table2_F2F<-Summarize_results(resultsDE_F2F,resultsDEW_F2F)
write.csv(Table2_F2F, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_F2F_populationE.csv")
rm(resultsDE_F2F)
rm(resultsDEW_F2F)
rm(Table2_F2F)
rm(parameter_F2F)

##################       assumes ONLINE delivery ONLY

#replace COST_DESMOND with ONLINE cost
parameter_ONLINE<-parameter
parameter_ONLINE[,"COST_DESMOND"]<-12.47

start.time_ONLINE <- Sys.time()
resultsDE_ONLINE <- run_model(population_clean, 
                           parameter_ONLINE, 
                           50, 
                           "DE", 
                           GlobalVars,
                           random_numbers,
                           LifeTables)
write.csv(resultsDE_ONLINE, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDE_ONLINE_populationE.csv")


resultsDEW_ONLINE <- run_model(population_clean, 
                             parameter_ONLINE, 
                             50, 
                             "DEW", 
                             GlobalVars,
                             random_numbers,
                             LifeTables)
write.csv(resultsDEW_ONLINE, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/resultsDEW_ONLINE_populationE.csv")

Table2_ONLINE<-Summarize_results(resultsDE_ONLINE,resultsDEW_ONLINE)
write.csv(Table2_ONLINE, "Results/Table 2 - Main Incremental Cost-effectiveness Analysis/Table2_ONLINE_populationE.csv")
rm(resultsDE_ONLINE)
rm(resultsDEW_ONLINE)
rm(Table2_ONLINE)
rm(parameter_ONLINE)

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
################################             Table 3: Sub-group analysis          ########################################





##################       Disease duration <=1 year

population_DM1 <-read.csv("Populations/Model_Population_DM1.csv")
population_DM1 <-subsetpopulation(population_DM1)
##rebuild new population
population_clean_DM1 <- build_population(population_DM1, PopulationVariables, GlobalVars)



start.time_popDM1 <- Sys.time()
resultsDE_popDM1 <- run_model(population_clean_DM1, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_popDM1, "Results/Table 3 - Sub-group analysis/resultsDE_popDM1_populationE.csv")


resultsDEW_popDM1 <- run_model(population_clean_DM1, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_popDM1, "Results/Table 3 - Sub-group analysis/resultsDEW_popDM1_populationE.csv")

Table3_popDM1<-Summarize_results(resultsDE_popDM1,resultsDEW_popDM1)
write.csv(Table3_popDM1, "Results/Table 3 - Sub-group analysis/Table3_popDM1_populationE.csv")
rm(resultsDE_popDM1)
rm(resultsDEW_popDM1)
rm(Table3_popDM1)
rm(population_DM1)
rm(population_clean_DM1)


##################       Disease duration 2-3 years

population_DM2 <-read.csv("Populations/Model_Population_DM2.csv")
population_DM2 <-subsetpopulation(population_DM2)
##rebuild new population
population_clean_DM2 <- build_population(population_DM2, PopulationVariables, GlobalVars)


start.time_popDM2 <- Sys.time()
resultsDE_popDM2 <- run_model(population_clean_DM2, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_popDM2, "Results/Table 3 - Sub-group analysis/resultsDE_popDM2_populationE.csv")


resultsDEW_popDM2 <- run_model(population_clean_DM2, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_popDM2, "Results/Table 3 - Sub-group analysis/resultsDEW_popDM2_populationE.csv")

Table3_popDM2<-Summarize_results(resultsDE_popDM2,resultsDEW_popDM2)
write.csv(Table3_popDM2, "Results/Table 3 - Sub-group analysis/Table3_popDM2_populationE.csv")
rm(resultsDE_popDM2)
rm(resultsDEW_popDM2)
rm(Table3_popDM2)
rm(population_DM2)
rm(population_clean_DM2)

end.time_popDM2  <- Sys.time()
run.time_popDM2 <-end.time_popDM2  - start.time_popDM2






##################       BMI groups: 28-30kg/m2

population_BMI2830 <-read.csv("Populations/Model_Population_BMI2830.csv")
population_BMI2830 <-subsetpopulation(population_BMI2830)
##rebuild new population
population_clean_BMI2830 <- build_population(population_BMI2830, PopulationVariables, GlobalVars)


start.time_popBMI2830 <- Sys.time()
resultsDE_popBMI2830 <- run_model(population_clean_BMI2830, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_popBMI2830, "Results/Table 3 - Sub-group analysis/resultsDE_popBMI2830_populationE.csv")


resultsDEW_popBMI2830 <- run_model(population_clean_BMI2830, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_popBMI2830, "Results/Table 3 - Sub-group analysis/resultsDEW_popBMI2830_populationE.csv")

Table3_popBMI2830<-Summarize_results(resultsDE_popBMI2830,resultsDEW_popBMI2830)
write.csv(Table3_popBMI2830, "Results/Table 3 - Sub-group analysis/Table3_popBMI2830_populationE.csv")
rm(resultsDE_popBMI2830)
rm(resultsDEW_popBMI2830)
rm(Table3_popBMI2830)
rm(population_BMI2830)
rm(population_clean_BMI2830)




##################       BMI groups: 30-35kg/m2

population_BMI3035 <-read.csv("Populations/Model_Population_BMI3035.csv")
population_BMI3035 <-subsetpopulation(population_BMI3035)
##rebuild new population
population_clean_BMI3035 <- build_population(population_BMI3035, PopulationVariables, GlobalVars)


start.time_popBMI3035 <- Sys.time()
resultsDE_popBMI3035 <- run_model(population_clean_BMI3035, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_popBMI3035, "Results/Table 3 - Sub-group analysis/resultsDE_popBMI3035_populationE.csv")


resultsDEW_popBMI3035 <- run_model(population_clean_BMI3035, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_popBMI3035, "Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3035_populationE.csv")

Table3_popBMI3035<-Summarize_results(resultsDE_popBMI3035,resultsDEW_popBMI3035)
write.csv(Table3_popBMI3035, "Results/Table 3 - Sub-group analysis/Table3_popBMI3035_populationE.csv")
rm(resultsDE_popBMI3035)
rm(resultsDEW_popBMI3035)
rm(Table3_popBMI3035)
rm(population_BMI3035)
rm(population_clean_BMI3035)




##################       BMI groups: 35-40kg/m2

population_BMI3540 <-read.csv("Populations/Model_Population_BMI3540.csv")
population_BMI3540 <-subsetpopulation(population_BMI3540)
##rebuild new population
population_clean_BMI3540 <- build_population(population_BMI3540, PopulationVariables, GlobalVars)


start.time_popBMI3540 <- Sys.time()
resultsDE_popBMI3540 <- run_model(population_clean_BMI3540, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_popBMI3540, "Results/Table 3 - Sub-group analysis/resultsDE_popBMI3540_populationE.csv")


resultsDEW_popBMI3540 <- run_model(population_clean_BMI3540, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_popBMI3540, "Results/Table 3 - Sub-group analysis/resultsDEW_popBMI3540_populationE.csv")

Table3_popBMI3540<-Summarize_results(resultsDE_popBMI3540,resultsDEW_popBMI3540)
write.csv(Table3_popBMI3540, "Results/Table 3 - Sub-group analysis/Table3_popBMI3540_populationE.csv")
rm(resultsDE_popBMI3540)
rm(resultsDEW_popBMI3540)
rm(Table3_popBMI3540)
rm(population_BMI3540)
rm(population_clean_BMI3540)





##################       BMI groups: 40+kg/m2

population_BMI40 <-read.csv("Populations/Model_Population_BMI40.csv")
population_BMI40 <-subsetpopulation(population_BMI40)
##rebuild new population
population_clean_BMI40 <- build_population(population_BMI40, PopulationVariables, GlobalVars)


start.time_popBMI40 <- Sys.time()
resultsDE_popBMI40 <- run_model(population_clean_BMI40, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_popBMI40, "Results/Table 3 - Sub-group analysis/resultsDE_popBMI40_populationE.csv")


resultsDEW_popBMI40 <- run_model(population_clean_BMI40, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_popBMI40, "Results/Table 3 - Sub-group analysis/resultsDEW_popBMI40_populationE.csv")

Table3_popBMI40<-Summarize_results(resultsDE_popBMI40,resultsDEW_popBMI40)
write.csv(Table3_popBMI40, "Results/Table 3 - Sub-group analysis/Table3_popBMI40_populationE.csv")
rm(resultsDE_popBMI40)
rm(resultsDEW_popBMI40)
rm(Table3_popBMI40)
rm(population_BMI40)
rm(population_clean_BMI40)





##################      IMD Socioeconomic quantiles (1: least deprived)

population_IMD1 <-read.csv("Populations/Model_Population_IMD1.csv")
population_IMD1 <-subsetpopulation(population_IMD1)
##rebuild new population
population_clean_IMD1 <- build_population(population_IMD1, PopulationVariables, GlobalVars)


start.time_popIMD1 <- Sys.time()
resultsDE_popIMD1 <- run_model(population_clean_IMD1, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_popIMD1, "Results/Table 3 - Sub-group analysis/resultsDE_popIMD1_populationE.csv")


resultsDEW_popIMD1 <- run_model(population_clean_IMD1, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_popIMD1, "Results/Table 3 - Sub-group analysis/resultsDEW_popIMD1_populationE.csv")

Table3_popIMD1<-Summarize_results(resultsDE_popIMD1,resultsDEW_popIMD1)
write.csv(Table3_popIMD1, "Results/Table 3 - Sub-group analysis/Table3_popIMD1_populationE.csv")
rm(resultsDE_popIMD1)
rm(resultsDEW_popIMD1)
rm(Table3_popIMD1)
rm(population_IMD1)
rm(population_clean_IMD1)






##################      IMD Socioeconomic quantiles 2

population_IMD2 <-read.csv("Populations/Model_Population_IMD2.csv")
population_IMD2 <-subsetpopulation(population_IMD2)
##rebuild new population
population_clean_IMD2 <- build_population(population_IMD2, PopulationVariables, GlobalVars)


start.time_popIMD2 <- Sys.time()
resultsDE_popIMD2 <- run_model(population_clean_IMD2, 
                               parameter, 
                               50, 
                               "DE", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDE_popIMD2, "Results/Table 3 - Sub-group analysis/resultsDE_popIMD2_populationE.csv")


resultsDEW_popIMD2 <- run_model(population_clean_IMD2, 
                                parameter, 
                                50, 
                                "DEW", 
                                GlobalVars,
                                random_numbers,
                                LifeTables)
write.csv(resultsDEW_popIMD2, "Results/Table 3 - Sub-group analysis/resultsDEW_popIMD2_populationE.csv")

Table3_popIMD2<-Summarize_results(resultsDE_popIMD2,resultsDEW_popIMD2)
write.csv(Table3_popIMD2, "Results/Table 3 - Sub-group analysis/Table3_popIMD2_populationE.csv")
rm(resultsDE_popIMD2)
rm(resultsDEW_popIMD2)
rm(Table3_popIMD2)
rm(population_IMD2)
rm(population_clean_IMD2)


##################      IMD Socioeconomic quantiles 3

population_IMD3 <-read.csv("Populations/Model_Population_IMD3.csv")
population_IMD3 <-subsetpopulation(population_IMD3)
##rebuild new population
population_clean_IMD3 <- build_population(population_IMD3, PopulationVariables, GlobalVars)


start.time_popIMD3 <- Sys.time()
resultsDE_popIMD3 <- run_model(population_clean_IMD3, 
                               parameter, 
                               50, 
                               "DE", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDE_popIMD3, "Results/Table 3 - Sub-group analysis/resultsDE_popIMD3_populationE.csv")


resultsDEW_popIMD3 <- run_model(population_clean_IMD3, 
                                parameter, 
                                50, 
                                "DEW", 
                                GlobalVars,
                                random_numbers,
                                LifeTables)
write.csv(resultsDEW_popIMD3, "Results/Table 3 - Sub-group analysis/resultsDEW_popIMD3_populationE.csv")

Table3_popIMD3<-Summarize_results(resultsDE_popIMD3,resultsDEW_popIMD3)
write.csv(Table3_popIMD3, "Results/Table 3 - Sub-group analysis/Table3_popIMD3_populationE.csv")
rm(resultsDE_popIMD3)
rm(resultsDEW_popIMD3)
rm(Table3_popIMD3)
rm(population_IMD3)
rm(population_clean_IMD3)



##################      IMD Socioeconomic quantiles 4

population_IMD4 <-read.csv("Populations/Model_Population_IMD4.csv")
population_IMD4 <-subsetpopulation(population_IMD4)
##rebuild new population
population_clean_IMD4 <- build_population(population_IMD4, PopulationVariables, GlobalVars)


start.time_popIMD4 <- Sys.time()
resultsDE_popIMD4 <- run_model(population_clean_IMD4, 
                               parameter, 
                               50, 
                               "DE", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDE_popIMD4, "Results/Table 3 - Sub-group analysis/resultsDE_popIMD4_populationE.csv")


resultsDEW_popIMD4 <- run_model(population_clean_IMD4, 
                                parameter, 
                                50, 
                                "DEW", 
                                GlobalVars,
                                random_numbers,
                                LifeTables)
write.csv(resultsDEW_popIMD4, "Results/Table 3 - Sub-group analysis/resultsDEW_popIMD4_populationE.csv")

Table3_popIMD4<-Summarize_results(resultsDE_popIMD4,resultsDEW_popIMD4)
write.csv(Table3_popIMD4, "Results/Table 3 - Sub-group analysis/Table3_popIMD4_populationE.csv")
rm(resultsDE_popIMD4)
rm(resultsDEW_popIMD4)
rm(Table3_popIMD4)
rm(population_IMD4)
rm(population_clean_IMD4)



##################      IMD Socioeconomic quantiles (5: most deprived)

population_IMD5 <-read.csv("Populations/Model_Population_IMD5.csv")
population_IMD5 <-subsetpopulation(population_IMD5)
##rebuild new population
population_clean_IMD5 <- build_population(population_IMD5, PopulationVariables, GlobalVars)


start.time_popIMD5 <- Sys.time()
resultsDE_popIMD5 <- run_model(population_clean_IMD5, 
                               parameter, 
                               50, 
                               "DE", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDE_popIMD5, "Results/Table 3 - Sub-group analysis/resultsDE_popIMD5_populationE.csv")


resultsDEW_popIMD5 <- run_model(population_clean_IMD5, 
                                parameter, 
                                50, 
                                "DEW", 
                                GlobalVars,
                                random_numbers,
                                LifeTables)
write.csv(resultsDEW_popIMD5, "Results/Table 3 - Sub-group analysis/resultsDEW_popIMD5_populationE.csv")

Table3_popIMD5<-Summarize_results(resultsDE_popIMD5,resultsDEW_popIMD5)
write.csv(Table3_popIMD5, "Results/Table 3 - Sub-group analysis/Table3_popIMD5_populationE.csv")
rm(resultsDE_popIMD5)
rm(resultsDEW_popIMD5)
rm(Table3_popIMD5)
rm(population_IMD5)
rm(population_clean_IMD5)










###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
################################             Table 4: Sensitivity analysis          ########################################



##################       Treatment Effect Maintenance 3 years

###replace treatment effect duration in global variables
#GlobalVars_TreatDur3<-GlobalVars
GlobalVars["Treatment_effect_duration_HBA", "Value"] <-3
GlobalVars["Treatment_effect_duration_BMI", "Value"] <-3

start.time_TreatDur3 <- Sys.time()
resultsDE_TreatDur3 <- run_model(population_clean, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_TreatDur3, "Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur3_populationE.csv")


resultsDEW_TreatDur3 <- run_model(population_clean, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_TreatDur3, "Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur3_populationE.csv")

Table4_TreatDur3<-Summarize_results(resultsDE_TreatDur3,resultsDEW_TreatDur3)
write.csv(Table4_TreatDur3, "Results/Table 4 - Sensitivity Analysis/Table4_TreatDur3_populationE.csv")
rm(resultsDE_TreatDur3)
rm(resultsDEW_TreatDur3)
rm(Table4_TreatDur3)




##################       Treatment Effect Maintenance 13 years

###replace treatment effect duration in global variables
#GlobalVars_TreatDur13<-GlobalVars
GlobalVars["Treatment_effect_duration_HBA", "Value"] <-13
GlobalVars["Treatment_effect_duration_BMI", "Value"] <-13


start.time_TreatDur13 <- Sys.time()
resultsDE_TreatDur13 <- run_model(population_clean, 
                              parameter, 
                              50, 
                              "DE", 
                              GlobalVars,
                              random_numbers,
                              LifeTables)
write.csv(resultsDE_TreatDur13, "Results/Table 4 - Sensitivity Analysis/resultsDE_TreatDur13_populationE.csv")


resultsDEW_TreatDur13 <- run_model(population_clean, 
                               parameter, 
                               50, 
                               "DEW", 
                               GlobalVars,
                               random_numbers,
                               LifeTables)
write.csv(resultsDEW_TreatDur13, "Results/Table 4 - Sensitivity Analysis/resultsDEW_TreatDur13_populationE.csv")

Table4_TreatDur13<-Summarize_results(resultsDE_TreatDur13,resultsDEW_TreatDur13)
write.csv(Table4_TreatDur13, "Results/Table 4 - Sensitivity Analysis/Table4_TreatDur13_populationE.csv")
rm(resultsDE_TreatDur13)
rm(resultsDEW_TreatDur13)
rm(Table4_TreatDur13)





##################       Cost of DESMOND face to face assuming optimistic utilisation (£95.05: £159 F2F and £7.47 online)

GlobalVars["Treatment_effect_duration_HBA", "Value"] <-5
GlobalVars["Treatment_effect_duration_BMI", "Value"] <-10
#replace COST_DESMOND with Optimistic cost
parameter_OptimisticUtil<-parameter
parameter_OptimisticUtil[,"COST_DESMOND"]<-96.06

start.time_OptimisticUtil <- Sys.time()
resultsDE_OptimisticUtil <- run_model(population_clean, 
                                  parameter_OptimisticUtil, 
                                  50, 
                                  "DE", 
                                  GlobalVars,
                                  random_numbers,
                                  LifeTables)
write.csv(resultsDE_OptimisticUtil, "Results/Table 4 - Sensitivity Analysis/resultsDE_OptimisticUtil_populationE.csv")


resultsDEW_OptimisticUtil <- run_model(population_clean, 
                                   parameter_OptimisticUtil, 
                                   50, 
                                   "DEW", 
                                   GlobalVars,
                                   random_numbers,
                                   LifeTables)
write.csv(resultsDEW_OptimisticUtil, "Results/Table 4 - Sensitivity Analysis/resultsDEW_OptimisticUtil_populationE.csv")

Table4_OptimisticUtil<-Summarize_results(resultsDE_OptimisticUtil,resultsDEW_OptimisticUtil)
write.csv(Table4_OptimisticUtil, "Results/Table 4 - Sensitivity Analysis/Table4_OptimisticUtil_populationE.csv")
rm(resultsDE_OptimisticUtil)
rm(resultsDEW_OptimisticUtil)
rm(Table4_OptimisticUtil)
rm(parameter_OptimisticUtil)


################################          Applying a zero utility decrement for BMI  ###########
parameter_zerobmidec<-parameter
parameter_zerobmidec[,"UTIL_BMI"]<-0

start.time_zerobmidec <- Sys.time()
resultsDE_zerobmidec <- run_model(population_clean, 
                                      parameter_zerobmidec, 
                                      50, 
                                      "DE", 
                                      GlobalVars,
                                      random_numbers,
                                      LifeTables)
write.csv(resultsDE_zerobmidec, "Results/Table 4 - Sensitivity Analysis/resultsDE_zerobmidec_populationE.csv")


resultsDEW_zerobmidec <- run_model(population_clean, 
                                       parameter_zerobmidec, 
                                       50, 
                                       "DEW", 
                                       GlobalVars,
                                       random_numbers,
                                       LifeTables)
write.csv(resultsDEW_zerobmidec, "Results/Table 4 - Sensitivity Analysis/resultsDEW_zerobmidec_populationE.csv")

Table4_zerobmidec<-Summarize_results(resultsDE_zerobmidec,resultsDEW_zerobmidec)
write.csv(Table4_zerobmidec, "Results/Table 4 - Sensitivity Analysis/Table4_zerobmidec_populationE.csv")
rm(resultsDE_zerobmidec)
rm(resultsDEW_zerobmidec)
rm(Table4_zerobmidec)
rm(parameter_zerobmidec)















