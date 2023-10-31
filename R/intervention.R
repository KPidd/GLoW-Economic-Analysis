##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@return INTE_A1c is a matrix that gives the reduction in A1c for each patient 
##'in each year
##'

  
##this function calls in the the correct intervention function depending on which treatment arm is being run
##  "DEW" is the intervention and "DE" will be ran if applying no treatment effect. 
##  "GLOW_beta_treat" or GLOW_beta_bl" apply the beta regression approach to find 12-month HbA1c just for intervention arm or baseline arm, while 
##  "GLOW_beta_diff" estimated both the treatment and baseline 12 month outcomes for HbA1c and finds the difference - this is the treatment effect that is applied.
initialise_intervention_dt_HbA1c <- function(n_,treatment_, parameters_,population_,endtime_,HBA1c_underlying_, random_numbs_) {
  INTE_A1c <- matrix(data=NA, nrow = n_, ncol =endtime_+2)
  INTE_A1c[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  if(treatment_ == "DEW"){
    INTE_A1c[,2] <- parameters_[,"intv_effect_GLOW_Hba1c"]
  }else if(treatment_ == "GLOW_beta_treat"){
    INTE_A1c[,2]<-INTE_HBA1c_12m_DEW(parameters_, population_, treatment_,random_numbs_)-HBA1c_underlying_[,"1"]
  }else if(treatment_ == "GLOW_beta_bl"){
    INTE_A1c[,2]<-INTE_HBA1c_12m_DE(parameters_, population_, treatment_,random_numbs_)-HBA1c_underlying_[,"1"]
  }else if(treatment_ == "GLOW_beta_diff"){
    INTE_A1c[,2]<-(INTE_HBA1c_12m_DEW(parameters_, population_, treatment_,random_numbs_)-HBA1c_underlying_[,"1"])-(INTE_HBA1c_12m_DE(parameters_, population_, treatment_,random_numbs_)-HBA1c_underlying_[,"1"])
  }else{ #if no treatment option is selected, leave them at baseline values
    INTE_A1c[,2:(endtime_+2)] <- 0
  }
  
  return(INTE_A1c)
}

##this function allows the initial 12 month treatment to decay linearly for the assumed duration of treatment effect.
INTE_HBA1c_decay<-function(GlobalVars_, treatment_,HBA1c_INTV,year_, endtime_){
  if(year_<(as.numeric(GlobalVars["Treatment_effect_duration_HBA","Value" ]))){
    HBA1c_INTV[,(year_+3):(as.numeric(GlobalVars["Treatment_effect_duration_HBA","Value" ])+2)] <-HBA1c_INTV[,2]-(HBA1c_INTV[,2]*((as.numeric(year_+1))/(as.numeric(GlobalVars["Treatment_effect_duration_HBA","Value" ]))))
  }else{
    HBA1c_INTV[,(as.numeric(GlobalVars["Treatment_effect_duration_HBA","Value"])+2):(endtime_+2)] <- 0
  }
  return(HBA1c_INTV)
  
} 

##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@return INTE_BMI is a matrix that gives the reduction in BMI for each patient 
##'in each year

##this function initiates the treatment effect on BMI
initialise_intervention_dt_BMI <- function(n_,population_,treatment_, parameters_,endtime_) {
  INTE_BMI <- matrix(data=NA, nrow = n_, ncol =endtime_+2)
  INTE_BMI[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  if(treatment_ == "DEW"|treatment_ == "GLOW_beta_diff"){
    INTE_BMI[,2] <- parameters_[,"intv_effect_GLOW_BMI"]
  }else{ #if no treatment option is selected, leave them at baseline values
    INTE_BMI[,2:(endtime_+2)] <- 0
  }
   return(INTE_BMI)
}
##allows for treatment effect of BMI to decay linearly for assumed treatment duration
INTE_BMI_decay<-function(GlobalVars_, treatment_,BMI_INTV,year_, endtime_){
  if(year_<(as.numeric(GlobalVars["Treatment_effect_duration_BMI","Value" ]))){
     BMI_INTV[,(year_+3):(as.numeric(GlobalVars["Treatment_effect_duration_BMI","Value" ])+2)] <-BMI_INTV[,2]-(BMI_INTV[,2]*((as.numeric(year_+1))/(as.numeric(GlobalVars["Treatment_effect_duration_BMI","Value" ]))))
  }else{
     BMI_INTV[,(as.numeric(GlobalVars["Treatment_effect_duration_BMI","Value"])+2):(endtime_+2)] <- 0
   }
  return(BMI_INTV)
    
}   
  

########## this is a alternative treatment effect function for  bmi which uses the treatment effect on weight and applied the bmi change based on the individuals height.---------------issue! i don't think BMI is equal to sampled weight/height^2
initialise_intervention_dt_BMI_2 <- function(n_,population_,treatment_, parameter_,endtime_) {
  INTE_BMI <- matrix(data=NA, nrow = n_, ncol =endtime_+2)
  INTE_BMI[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  if(treatment_ == "test"|treatment_ == "GLOW_BMI"){
    INTE_BMI[,2] <- (population_[,"WEIGHT"]-0.81)/((population_[,"HEIGHT"]/100)^2)-BMI_underlying[,2]
  }else{ #if no treatment option is selected, leave them at baseline values
    INTE_BMI[,2:(endtime_+2)] <- 0
  }
  return(INTE_BMI)
}










####beta regression functions
##function to find the 12 month treatment effect difference for hba1c from 0 from beta regression
##population_[,"RAND_beta_reg"]

####finds 12 month HbA1c had the individual received treatment DEW
INTE_HBA1c_12m_DEW<-function(parameters_, population_, treatment_, random_numbs_){
  
  beta_conmean_DEW<-GLOW_12m_intervention_effect_mean_DEW(parameters_, population_, treatment_)
  beta_convar_DEW<-GLOW_12m_intervention_effect_var_DEW(parameters_, population_, treatment_)
  beta_alpha_DEW<-(beta_conmean_DEW*(((beta_conmean_DEW*(1-beta_conmean_DEW))/beta_convar_DEW)-1))
  beta_beta_DEW<-((1-beta_conmean_DEW)/beta_conmean_DEW)*beta_alpha_DEW
  INTE_HBA1c_12m_DEW<-qbeta(0.5,beta_alpha_DEW,beta_beta_DEW)
  ##converting it back from 0,1 interval (using min and max of 12 month hba1c in GLOW)
  INTE_HBA1c_12m_DEW<-(INTE_HBA1c_12m_DEW*(17.15595-4.162993))+4.162993
  return(INTE_HBA1c_12m_DEW)
  
}
#random_numbs_[,"OST",75]
####finds 12 month HbA1c had the individual received usual care DE
INTE_HBA1c_12m_DE<-function(parameters_, population_, treatment_,random_numbs_){
  
  beta_conmean_DE<-GLOW_12m_intervention_effect_mean_DE(parameters_, population_, treatment_)
  beta_convar_DE<-GLOW_12m_intervention_effect_var_DE(parameters_, population_, treatment_)
  beta_alpha_DE<-beta_conmean_DE*(((beta_conmean_DE*(1-beta_conmean_DE))/beta_convar_DE)-1)
  beta_beta_DE<-((1-beta_conmean_DE)/beta_conmean_DE)*beta_alpha_DE
  INTE_HBA1c_12m_DE<-qbeta(0.5,beta_alpha_DE,beta_beta_DE)
  ##converting it back from 0,1 interval (using min and max of 12 month hba1c in GLOW)
  INTE_HBA1c_12m_DE<-(INTE_HBA1c_12m_DE*(17.15595-4.162993))+4.162993
  return(INTE_HBA1c_12m_DE)
  
}
##adding treatment effect
###y~B(U,Var)

###conditional mean of treatment group DEW
GLOW_12m_intervention_effect_mean_DEW<-function(parameters_, population_, treatment_){
  fitted_value<-matrix(parameters_[,"beta_intercept_model1"]
                       +parameters_[,"beta_hba1c_bl_model1"]*population_[,"HBA"]
                       +parameters_[,"beta_treatment_model1"]*1 
                       +parameters_[,"beta_male_model1"]*population_[,"MALE"]
                       +parameters_[,"beta_dm_duration_1_model1"]*ifelse(population_[,"DIAB_DUR"]==1,1,0))
  
  beta_conmean_DEW<-(exp(fitted_value)/(1+exp(fitted_value)))
  
  return(beta_conmean_DEW)
  
  
}
###conditional mean of baseline group DE 
GLOW_12m_intervention_effect_mean_DE<-function(parameters_, population_, treatment_){
  fitted_value<-matrix(parameters_[,"beta_intercept_model1"]
                       +parameters_[,"beta_hba1c_bl_model1"]*population_[,"HBA"]
                       +parameters_[,"beta_treatment_model1"]*0 
                       +parameters_[,"beta_male_model1"]*population_[,"MALE"]
                       +parameters_[,"beta_dm_duration_1_model1"]*ifelse(population_[,"DIAB_DUR"]==1,1,0))
  
  beta_conmean_DE<-(exp(fitted_value)/(1+exp(fitted_value)))
  
  return(beta_conmean_DE)
  
}
###conditional variance of treatment group DEW
GLOW_12m_intervention_effect_var_DEW<-function(parameters_, population_, treatment_){ 
  fitted_value<-matrix(parameters_[,"beta_intercept_model1"]
                       +parameters_[,"beta_hba1c_bl_model1"]*population_[,"HBA"]
                       +parameters_[,"beta_treatment_model1"]*1 
                       +parameters_[,"beta_male_model1"]*population_[,"MALE"]
                       +parameters_[,"beta_dm_duration_1_model1"]*ifelse(population_[,"DIAB_DUR"]==1,1,0))
  beta_conmean_DEW<-(exp(fitted_value)/(1+exp(fitted_value)))
  scale_factor<-matrix(parameters_[,"beta_phi_intercept_model1"]
                       +(parameters_[,"beta_phi_hba1c_bl_model1"]*population_[,"HBA"])
                      +(parameters_[,"beta_phi_treatment_model1"]*1)
                      +(parameters_[,"beta_phi_hba1c_bl_treatment_model1"]*(1*population_[,"HBA"])))
  beta_convar_DEW<-(beta_conmean_DEW*(1-beta_conmean_DEW))/(1+exp(scale_factor))
  
  return(beta_convar_DEW)
}
###conditional variance of treatment group DE
GLOW_12m_intervention_effect_var_DE<-function(parameters_, population_, treatment_){ 
  fitted_value<-matrix(parameters_[,"beta_intercept_model1"]
                       +parameters_[,"beta_hba1c_bl_model1"]*population_[,"HBA"]
                       +parameters_[,"beta_treatment_model1"]*0
                       +parameters_[,"beta_male_model1"]*population_[,"MALE"]
                       +parameters_[,"beta_dm_duration_1_model1"]*ifelse(population_[,"DIAB_DUR"]==1,1,0))
  beta_conmean_DE<-(exp(fitted_value)/(1+exp(fitted_value)))
  scale_factor<-matrix(parameters_[,"beta_phi_intercept_model1"]
                       +(parameters_[,"beta_phi_hba1c_bl_model1"]*population_[,"HBA"])
                       +(parameters_[,"beta_phi_treatment_model1"]*0)
                       +(parameters_[,"beta_phi_hba1c_bl_treatment_model1"]*(0*population_[,"HBA"]))) 

  beta_convar_DE<-(beta_conmean_DE*(1-beta_conmean_DE))/(1+exp(scale_factor))
  
  return(beta_convar_DE)
}



  













######interventions for other risk factors not considered in GLoW analysis - NOT IN USE
##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@return INTE_SBP is a matrix that gives the reduction in SBP for each patient 
##'in each year

initialise_intervention_dt_SBP <- function(n_,treatment_, parameter_,endtime_) {
  INTE_SBP <- matrix(data=NA, nrow = n_, ncol =endtime_+2)
  INTE_SBP[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  if(treatment_ == "test" | treatment_ == "Mt_HOOD_2018_QOL_SBP"){
    INTE_SBP[,2:(endtime_+2)] <- -10
  }else{ #if no treatment option is selected, leave them at baseline values
    INTE_SBP[,2:(endtime_+2)] <- 0
  }
  
  return(INTE_SBP)
}

##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@return INTE_HDL is a matrix that gives the reduction in HDL cholesterol 
##'for each patient in each year

initialise_intervention_dt_HDL <- function(n_,treatment_, parameter_,endtime_) {
  INTE_HDL <- matrix(data=NA, nrow = n_, ncol =endtime_+2)
  INTE_HDL[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  #store the trajectory by assumed treatment
  if(treatment_ == "test" | treatment_ == "Mt_HOOD_2018_QOL_SBP"){
    INTE_HDL[,2:(endtime_+2)] <- 0.5
  }else{ #if no treatment option is selected, leave them at baseline values
    INTE_HDL[,2:(endtime_+2)] <- 0
  }
  
  return(INTE_HDL)
}

##'@param n_ is the number of patients in the model
##'@param treatment_ is a text term indicating the current treatment option
##'@param parameter_ is the row of the parameter matrix
##'@param endtime_ is a number indicating how many years the simulation is being run for
##'@return INTE_LDL is a matrix that gives the reduction in HDL cholesterol 
##'for each patient in each year

initialise_intervention_dt_LDL <- function(n_,treatment_, parameter_,endtime_) {
  INTE_LDL <- matrix(data=NA, nrow = n_, ncol =endtime_+2)
  INTE_LDL[,1] <- 1:n_ #make the first column equivalent to the patient ID for matching later on
  #If it's a MT Hood scenario set a permanent benefit
  if(treatment_ == "test" | treatment_ == "Mt_HOOD_2018_QOL_LDL"){
    INTE_LDL[,2:(endtime_+2)] <- -0.5
  }else{ #if no treatment option is selected, leave them at baseline values
    INTE_LDL[,2:(endtime_+2)] <- 0
  }
  return(INTE_LDL)
}

##'@param Input_Prob_ is a numeric number between 0 and 1
##'@param HR_ is a hazard ratio 
##'@return output_prob is the new probability after applying the hazard ratio

Hazard_Ratio_Intervention <- function(Input_Prob_,
                                      HR_){
  
  input_rate <- -log(1-Input_Prob_)/1
  output_rate <- input_rate*HR_
  output_prob <- 1-exp(-output_rate*1)
  
  return(output_prob)
  
}


