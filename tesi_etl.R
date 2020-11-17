library("dplyr")
library("readxl")
library("stringr")
library("readxl")

#upload data
MLG <- read.table("MLG.profile", encoding="UTF-8", stringsAsFactors=FALSE)
#KO <- read.table("KO.profile", encoding="UTF-8", stringsAsFactors=FALSE)
sample <- read.delim("1301samples_config")
infoMLG <- read.delim("MLG.infor")

#MLG aka Metagenomics linkage group

#traspongo per avere soggetti in riga e cluster in colonna
#assegno alle colonne di MLG la rispettiva classe tassonomica
MLG <- as.data.frame(t(MLG))
colnames(MLG) <- infoMLG[,3]
#KO <- as.data.frame(t(KO))

# MLG_ACVD <- MLG %>%
#   filter(Project.Corhort=="ACVD")
# 
# spars_pre <- apply(X = MLG_ACVD[,-c(1,2,3)], 2, FUN = function(c) (sum(c > 0)/405*100))
# hist(spars_pre,
#      xlim = c(0,100),
#      nclass = 100,
#      xlab = "Percent of samples", ylim = c(0,40), ylab = "Number of MLGs",
#      col = "red",
#      main = "Numero di MLG presenti in una data percentuale dei campioni")


#aggiungo sample id, disease e status: 0 controllo 1 malato
MLG <- cbind(sample,MLG)
#KO <- cbind(sample,KO)

#estraggo il file intero indistinto per classificati e non classificati
#tot campioni: 1301 / tot MLG: 1043
MLG[,1] <- sapply(MLG[,1], as.character)
MLG[,3] <- sapply(MLG[,3], as.factor)
#write.csv(MLG, file = "MyData_all_1k.csv",row.names=FALSE)

#estraggo solo MLG classificati
# tot MLG 683
MLG_class <- MLG %>%
  select(-starts_with("unclass"))

MLG_class[,1] <- sapply(MLG_class[,1], as.character)
MLG_class[,3] <- sapply(MLG_class[,3], as.factor)

#write.csv(MLG_class, file = "MyData_class.csv",row.names=FALSE)

########################################
# #sommo tutti gli MLG non classificati
# MLG_unclass <- MLG %>%
#   select(unclass=starts_with("Unclassified-"))
# 
# MLG_unclass <- as.data.frame(t(MLG_unclass))
# 
# MLG_unclass <- MLG_unclass %>%
#   summarise_all(funs(sum))
# 
# MLG_unclass <- as.data.frame(t(MLG_unclass))
# 
# 
# MLG_ok<-cbind(MLG_class,MLG_unclass)
# names(MLG_ok)[names(MLG_ok) == "V1"] <- "Unclassified-MLG"
# MLG_ok[,1] <- sapply(MLG_ok[,1], as.character)
# MLG_ok[,3] <- sapply(MLG_ok[,3], as.factor)
# 
# write.csv(MLG_ok, file = "MyData_all.csv",row.names=FALSE)

#####################################################################
# #estraggo i casi di tutte le patologie
# 
# MLG_cases <- MLG_ok %>%
#   filter(Enrichment..0.control..1.case.==1)
# 
# write.csv(MLG_cases, file = "MyData_cases.csv",row.names=FALSE)
# 
# #as.h2o(MLG_cases, destination_frame = "MLG_cases")
# #h2o.flow()



#estraggo solo ACVD
MLG_ACVD <- MLG_class %>%
filter(Project.Corhort=="ACVD")

write.csv(MLG_ACVD, file = "MyData_ACVD.csv",row.names=FALSE)

# #as.h2o(MLG_ACVD, destination_frame = "MLG_ACVD")
# #h2o.flow()

###################################################################
#ACVD e clinical features
#aggiungo variabili cliniche
# MLG_ACVD <- MLG_ACVD[,-2]
# 
# fenotipo_ACVD <- read_excel("41467_2017_900_MOESM3_ESM.xlsx", 
#                             col_types = c("text", "numeric", "text", 
#                                           "numeric", "text", "numeric", "numeric", 
#                                           "numeric", "numeric", "numeric", 
#                                           "numeric", "text", "numeric", "numeric", 
#                                           "numeric", "numeric", "numeric", 
#                                           "numeric", "numeric", "numeric", "numeric", 
#                                           "numeric", "numeric", "numeric", "numeric", "numeric", 
#                                           "numeric", "numeric", "numeric", "numeric", "numeric", 
#                                           "numeric", "numeric", "numeric", "numeric", "numeric", 
#                                           "numeric", "text", "text", "text", "text", 
#                                           "text", "text", "text", "text", "text", 
#                                           "text", "text", "text", "text", "text", 
#                                           "text", "text", "text", "text"), 
#                             skip = 1)
# 
# fenotipo_ACVD <- as.data.frame(fenotipo_ACVD)
# 
# 
# #rinomino i nomi dei campioni in maniera uguale nelle due tabelle
# MLG_ACVD[,1] <- gsub('([A-Z])([0-9]+)A', '\\1\\2',MLG_ACVD[,1])
# MLG_ACVD[,1] <- gsub('(ZSL)([^0-9]*)([0-9]+)A', '\\1\\2\\3',MLG_ACVD[,1])
# MLG_ACVD[,1] <- gsub('(ZSL)([^0-9]*)([0-9]+)', '\\1\\3',MLG_ACVD[,1])
# MLG_ACVD[,1] <- paste0(gsub('([A-Z]+)([0-9]+)', '\\1', MLG_ACVD[,1]), str_pad(gsub('([A-Z]+)([0-9]+)', '\\2', MLG_ACVD[,1]),width = 4, pad='0'))
# 
# fenotipo_ACVD[,1] <- gsub('([A-Z])([0-9]+)A', '\\1\\2',fenotipo_ACVD[,1])
# fenotipo_ACVD[,1] <- gsub('(ZSL)([^0-9]*)([0-9]+)A', '\\1\\2\\3',fenotipo_ACVD[,1])
# fenotipo_ACVD[,1] <- gsub('(ZSL)([^0-9]*)([0-9]+)', '\\1\\3',fenotipo_ACVD[,1])
# fenotipo_ACVD[,1] <- paste0(gsub('([A-Z]+)([0-9]+)', '\\1', fenotipo_ACVD[,1]), str_pad(gsub('([A-Z]+)([0-9]+)', '\\2', fenotipo_ACVD[,1]),width = 4, pad='0'))
# 
# colnames(fenotipo_ACVD)[which(names(fenotipo_ACVD) == "Sample ID")] <- "SampleID"
# colnames(MLG_ACVD)[which(names(MLG_ACVD) == "Sample.ID")] <- "SampleID"
# 
# MLG_ACVD <- MLG_ACVD %>% arrange(SampleID)
# fenotipo_ACVD <- fenotipo_ACVD %>% arrange(SampleID)
# 
# all_ACVD <- cbind(fenotipo_ACVD, MLG_ACVD[,-c(1,2)])
# 
# write.csv(all_ACVD, file = "MyData_ACVDcl.csv", row.names=FALSE)


###################################
#confronto MLG classificati del file con MLG classificati del paper 

idMLG <- infoMLG %>% filter(!str_detect(Taxonomy, 'Unclassified')) %>%
  select(MLG.ID, Taxonomy)
idMLG$MLG.ID <- as.character(idMLG$MLG.ID)

sup_data <- read_excel("41467_2017_900_MOESM5_ESM - suppl.3.xlsx", 
                       sheet = "Sheet1", col_types = c("text", "numeric", "numeric", "text", "numeric", "numeric", "numeric",
                                                       "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                       "numeric", "numeric", "text", "numeric", "numeric"))

idMLG_paper <- sup_data %>% filter(!str_detect(Taxonomy, 'Unclassified')) %>%
  select(`MLG ID`, Taxonomy)

idMLG_paper_q <- sup_data %>%
  filter(`Q value` <= 0.05) %>%  select(`MLG ID`, Taxonomy)

#solo 360 combaciano, invece di 536 !
qq <- inner_join(idMLG, idMLG_paper_q, by = c("MLG.ID" = "MLG ID"))




