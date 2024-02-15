library(UniprotR)
library(protti)
library(GenomicAlignments)
library(r3dmol)

# reading imported data
accession1<-read.csv("/Users/arlinemartinez/Downloads/Accessions.csv")
accession1
ls(accession1)

#Removes unwanted variables in the original accession data frame. 
df = subset(accession1, select = -c(X,X.1))
df
class(df)

#Converts dataframe to list
List<- as.list(df)

List[["Accessions"]]
class(List[["Accessions"]])

PGI_obj<-GetProteinGOInfo(List[["Accessions"]],directorypath = NULL)

PlotGoInfo(PGI_obj,directorypath = NULL)
PlotGOAll(GOObj = PGI_obj, Top = 10, directorypath = getwd(), width = 8, height = 5)

Patho <-GetPathology_Biotech(List[["Accessions"]],directorypath=NULL)
Get.diseases(Patho,directory=NULL)

uni_protInfo<-fetch_uniprot("P0A799") 
pdb<- fetch_pdb (pdb_ids=c("1ZMR"))
head(pdb)                 

alphafold<-fetch_alphafold_prediction(uniprot_ids = c("P0A799"), return_data_frame = TRUE)
head(alphafold,n=10)

#Second accession
library(UniprotR)
library(protti)
library(GenomicAlignments)
library(r3dmol)

# reading imported data
accession2<-read.csv("/Users/arlinemartinez/Documents/GitHub/accessions 2.csv")
accession2
ls(accession2)

#Converts dataframe to list
List_2<- as.list(df)

List_2[["Accessions"]]
class(List[["Accessions"]])

PGI_obj2<-GetProteinGOInfo(List_2[["Accessions"]],directorypath = NULL)

PlotGoInfo(PGI_obj2,directorypath = NULL)
PlotGOAll(GOObj = PGI_obj2, Top = 10, directorypath = getwd(), width = 8, height = 5)

Patho2 <-GetPathology_Biotech(List[["Accessions"]],directorypath=NULL)
Get.diseases(Patho2,directory=NULL)

uni_protInfo2<-fetch_uniprot("P08839") 
pdb_2<- fetch_pdb (pdb_ids=c("2HWG"))
head(pdb_2)                 

alphafold2<-fetch_alphafold_prediction(uniprot_ids = c("P08839"), return_data_frame = TRUE)
head(alphafold,n=10)
