
#Perform data analysis

#Load data and remove CPR-number, birth date and hospital identification
filename <- file.path(data.out, "Colon.RData")
if(file.exists(filename)){
  load(filename)
}else{
  #Set working directory and load the datasets from the drives
  setwd("../StatCure/ExternalData/")
  Colon <- read.csv2("F_2017_05_05_DCCG_ekstra.csv", header = T, stringsAsFactors = F)
  
  #Derive birth date from CPR-number and format date of diagnosis. 
  #Patients born within the first 9 days of a month are missing a 0 in the beginning of the CPR-number.
  #Create an age variable and remove birth day. No hospital specification is a available
  wh <- nchar(Colon$Patient_CPR) == 9
  Colon$Patient_CPR[wh] <- paste0("0", Colon$Patient_CPR[wh])
  Colon$Patient_CPR <- substr(Colon$Patient_CPR, start = 1, stop = 6)
  Colon$birth_date <- format(as.Date(Colon$Patient_CPR, format = "%d%m%y"), "19%y-%m-%d")
  Colon$birth_date <- as.Date(Colon$birth_date)
  Colon$diag_date <- as.Date(Colon$Dato_diagnose)
  Colon$age <- as.numeric(Colon$diag_date - Colon$birth_date)
  Colon <- subset(Colon, select = -c(Patient_CPR))
  
  #Remove decrypted files from folder and save temporary file in project folder
  file.remove("F_2017_05_05_amlformat.csv", "F_2017_05_05_DCCGdel.csv", 
              "F_2017_05_05_DCCG_ekstra.csv", "F_2017_05_05_lyfodel.csv")
  setwd(project)
  save(Colon, file = filename)
}


#Clean relevant variables in the Colon cancer data
Colon$sex <- factor(Colon$Koen, levels = c("Mand", "Kvinde"), labels = c("male", "female"))
Colon <- Colon[!Colon$CPR_status %in% c("Inaktiv, men tildelt personnummer af skattehensyn", 
                                        "Speciel vejkode i dansk folkeregister", 
                                        ""),]
Colon$status <- as.numeric(Colon$CPR_status == ifelse(linux, "D\xf8d","Død"))
Colon$death_date <- Colon$CPR_status_dato
Colon$death_date[Colon$death_date == ""] <- Colon$opdat_dato[Colon$death_date == ""]
Colon$death_date <- as.Date(Colon$death_date)

wh <- which(Colon$death_date > Colon$opdat_dato)
Colon$death_date[wh] <- Colon$opdat_dato
Colon$status[wh] <- 0

Colon$FU <- as.numeric(Colon$death_date - Colon$diag_date)
Colon$age_years <- Colon$age / ayear
Colon$age_above_50 <- factor(as.numeric(Colon$age_years >= 50), levels = 0:1, labels = c("Young", "Old"))
Colon$FU_years <- Colon$FU / ayear
Colon <- Colon[Colon$FU_years > 0,]
Colon <- Colon[Colon$age_years >= 18, ] #& Colon$age_years <= 80,]


Colon <- Colon[Colon$Metastaser %in% c("Ja", "Nej"),]
Colon$Metastaser <- as.numeric(Colon$Metastaser == "Ja")
Colon$Charlson <- as.numeric(Colon$LPR_Charlson_score >= 2)

Colon$gender <- as.numeric(Colon$sex == "male")
Colon <- Colon[Colon$UICC %in% c("UICC stadium I", "UICC stadium II", 
                                 "UICC stadium III", "UICC stadium IV"),]
Colon$Stage <- as.numeric(Colon$UICC %in% c("UICC stadium III", "UICC stadium IV"))

#Extrac general population survival
Colon$exp_haz <- general.haz(time = "FU", age = "age", sex = "sex", year = "diag_date", 
                             data = Colon, ratetable = survexp.dk)

#Subset all datasets by the following few variables
variables <- c("age", "age_years", "age_above_50", "sex", "gender",  "diag_date", 
               "FU", "FU_years", "status", "exp_haz", "Charlson", "Metastaser", "Stage")
Colon <- Colon[, variables]
Colon$Metastases <- Colon$Metastaser
Colon$Age <- Colon$age_years
Colon$Gender <- Colon$gender



