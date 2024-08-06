library(drc)
library(tidyverse)

all.data <- read.csv("C:/Users/GMCDP/OneDrive - Bayer/Downloads/BLV_PA_RESULTS_4SPOTFI.csv")

time1 <- Sys.time()

###### make RES_VALUE numeric, maybe by default in Spotfire?

# all.data$RES_VALUE <- as.numeric(all.data$RES_VALUE) 

######

###### clean table to leave only Efficacies (not ED_50, not IC_50...), maybe not needed in Spotfire?

# all.data <- all.data[which(all.data$RES_VAR_NAME=="EFFICACY"),]

######


# create output table, select *efficacies* only on *testing samples (not controls)*
# then make a table with all compounds (PA_NAME) for all pathogens
# and add four columns: Model, ED50, UpperED50, lowerED50

testing.sample <- all.data[which(all.data$PA_TYPE=="T"), ]

unrepeated.testings <- unique(testing.sample[,c("SUBSTANCE", "PA_NAME", "PROC_NAME")])

ED.output <- data.frame(unrepeated.testings, Model = NA, ED_50 = NA, LowerED50 = NA, UpperED50 = NA, Top.Efficacy = NA)

model.params <- data.frame(SUBSTANCE = NA, PA_NAME = NA, PROC_NAME = NA, coef.b = NA, coef.c = NA, coef.d = NA, coef.e = NA)

# run your fit and estimate your ED50 line by line

for(i in 1:nrow(ED.output)) {
  
  produit <- ED.output$PA_NAME[i]
  
  case <- all.data[which(all.data$PA_NAME == ED.output$PA_NAME[i] & all.data$PROC_NAME == ED.output$PROC_NAME[i]),]
  
  ED.output$Top.Efficacy[i] <- case$RES_VALUE[which(case$DOSAGE==max(case$DOSAGE))]
  
  if(ED.output$Top.Efficacy[i] - case$RES_VALUE[which(case$DOSAGE==min(case$DOSAGE))] <10) {
    
    ED.output$Model[i] <- "Linear"
    
    next}
  
  ED.output$Model[i] <- "Log-logistic"
  
  if(case$RES_VALUE[which(case$DOSAGE==max(case$DOSAGE))] < 50) {next}
  
  model <- drm(RES_VALUE ~ DOSAGE, data = case, fct = LL.4())
  
  # faire cooks.distance()
  # if cooks.distance() > 3*mean(cooks.distance())
  # repeter modele sans outlier
  
  model.params[i,] <- data.frame(unique(case$SUBSTANCE), unique(case$PA_NAME), unique(case$PROC_NAME), t(model$coefficients))
  
  ED50 <- ED(model, 50, "delta", type = "absolute", display = F)
  
  ED.output$ED_50[i] <- ED50[1]
  ED.output$LowerED50[i] <- ED50[3]
  ED.output$UpperED50[i] <- ED50[4]
  
}

model.params <- na.omit(model.params)

time2 <- Sys.time()

###########################

unexplained <- ED.output %>%
  filter(Model == "Linear" & Top.Efficacy > 50)
