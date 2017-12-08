##Polar_Tract_Identification
##14th Nov 2017
#David Brown

#Load sequences
hi6<-read.table("Y:/DavidB/Data/Systematic_Search_for_Phase_Separating_Proteins/IDR_set_aa_properties.txt", sep="\t", header = TRUE) 

library(Biostrings)

        s<-hi6$IDR_seq
        
        ##Define Polar Tract Function for sequence vector 's'
        pt_fun<-function(s){
                      #may need to clean '*' from sequence
                    IDR_aa<-AAStringSet(as.character(s))
                    #names(IDR_aa)<-hi$ensembl_peptide_id
        
                    IDR_composition<-alphabetFrequency(IDR_aa)  #IDR_composition is a matrix and must be accessed via 'colnames' not 'names'
        
                    
        hi$proline_content<-IDR_composition[,"P"]/hi$width
        hi$proline_rich<-hi$proline_content>0.15
        
s<-IDR_aa

##Define Polar Tract Function for AAStringSet 's'
pt_fun<-function(s){
  #may need to clean '*' from sequence
  IDR_composition<-alphabetFrequency(s)  #IDR_composition is a matrix and must be accessed via 'colnames' not 'names'
  
  PT_frac<-rowSums(IDR_composition[,c("Q","N","S","G","P")]/width(s))                 #Neutral Polar AAs
  HP_frac<-rowSums(IDR_composition[,c("F","H","W","Y", "A", "I", "L", "V")]/width(s)) #Aromatic and Aliphatic AAs
  
  #Polar Tract
  PT_frac>0.5 & HP_frac<=0.5
}

hi6$PT<-pt_fun(s)
summary(hi6$PT)



##Define Polyelectolyte Function for AAStringSet 's'

        ##Avoid redundant calls to alphabetFreq!!!!!!!!!!!!!
pe_fun<-function(s){
  #may need to clean '*' from sequence
  IDR_composition<-alphabetFrequency(s)  #IDR_composition is a matrix and must be accessed via 'colnames' not 'names'
  
  PT_frac<-rowSums(IDR_composition[,c("Q","N","S","G","P")]/width(s))                 #Neutral Polar AAs
  HP_frac<-rowSums(IDR_composition[,c("F","H","W","Y", "A", "I", "L", "V")]/width(s)) #Aromatic and Aliphatic AAs
  
  #Polar Tract
  PT_frac>0.5 & HP_frac<=0.5
}

hi6$PT<-pt_fun(s)
summary(hi6$PT)

##Calculate the Charge Assymetry 

#Load dependancies
library(Peptides)

#Clean sequence 's'
cs_fun<-function(s) gsub(pattern = "[*]", replacement = "", x = s)

#Define fp and fn functios
fp_fun<-function(s) sapply(aaComp(s), function(x) x["Basic","Mole%"])/100
fn_fun<-function(s) sapply(aaComp(s), function(x) x["Acidic","Mole%"])/100

#Define charge asymmetry function
ca_fun<-function(s) ((fp_fun(s)-fn_fun(s))^2)/(fp_fun(s)+fn_fun(s))

hi6$CA<-ca_fun(s)

##Compare charge assymetry of PA and PE
plot(hi6[hi6$PE,]$NCPR, hi6[hi6$PE,]$CA, ylim=c(0,1))
points(hi6[hi6$PA,]$NCPR, hi6[hi6$PA,]$CA, col="blue")
points(hi6[hi6$PT,]$NCPR, hi6[hi6$PT,]$CA, col="green")

##Compare charge assymetry of PA and PE
plot(hi6[hi6$PE,]$Charge, hi6[hi6$PE,]$CA, xlim=c(0,1), ylim=c(0,1))
points(hi6[hi6$PA,]$Charge, hi6[hi6$PA,]$CA, col="blue")
points(hi6[hi6$PT,]$Charge, hi6[hi6$PT,]$CA, col="green")


##Compare charge assymetry of PA and PE
plot(hi6$Basic, hi6$Acidic, xlim=c(0,1), ylim=c(0,1))
points(hi6[hi6$PE,]$Basic, hi6[hi6$PE,]$Acidic, col="red")
points(hi6[hi6$PA,]$Basic, hi6[hi6$PA,]$Acidic, col="blue")
points(hi6[hi6$PT,]$Basic, hi6[hi6$PT,]$Acidic, col="green")

##Elastin Like Polypeptide (ELP)
elp_fun<-function(s) grepl(pattern="VPG.G", x = s)
library(stringr)
joint_IDR_dat2[elp_fun(joint_IDR_dat2$IDR_seq),]$IDR_seq

str_count(joint_IDR_dat2[elp_fun(joint_IDR_dat2$IDR_seq),]$IDR_seq, "VPG.G")
#No poly-ELPs

##ESX1 Like (ESXL)
esxl_fun<-function(s) grepl(pattern="PP..P.PP.", x = s)
joint_IDR_dat2[esxl_fun(joint_IDR_dat2$IDR_seq),]$IDR_seq

joint_IDR_dat2$ESXL_count<-str_count(joint_IDR_dat2$uniprot_sequence, "PP..P.PP.")
#2518 proteins with ESXL sequence

unique(joint_IDR_dat2[joint_IDR_dat2$ESXL_count>2,]$ensembl_gene_id)
#18 unique genes. Many important and many nuclear speckle associated. All HPA characterised.

##New candidates!!
nc<-unique(joint_IDR_dat2[joint_IDR_dat2$ESXL_count>2,]$ensembl_gene_id)

nc_dat<-joint_IDR_dat2[joint_IDR_dat2$ensembl_gene_id %in% nc,]


#Approximate matching
a_esxl_fun<-function(s) agrepl(pattern="PP..P.PP.", x = s)
joint_IDR_dat2[a_esxl_fun(joint_IDR_dat2$IDR_seq),]$IDR_seq

#Detect RS domains


##Detect scrambled ESXL
summary(str_count(joint_IDR_dat2$uniprot_sequence, "P..PP..PP")>2) #628
summary(str_count(joint_IDR_dat2$uniprot_sequence, ".P.PP..PP")>2) #882
summary(str_count(joint_IDR_dat2$uniprot_sequence, ".P.P.P.PP")>2) #502
summary(str_count(joint_IDR_dat2$uniprot_sequence, "P.P.P.P.P")>2) #417
summary(str_count(joint_IDR_dat2$uniprot_sequence, ".PP.P..PP")>2) #639

unique(joint_IDR_dat2[str_count(joint_IDR_dat2$uniprot_sequence, ".P.PP..PP")>2,]$ensembl_gene_id)
#28 genes

#scrambled esx (significant overlap with ESXL set)
sex<-unique(joint_IDR_dat2[str_count(joint_IDR_dat2$uniprot_sequence, ".P.PP..PP")>2,]$ensembl_gene_id)

unique(joint_IDR_dat2[str_count(joint_IDR_dat2$uniprot_sequence, ".P.PP..PP")>2,]$ensembl_gene_id)
#28 genes

#scrambled esx (significant overlap with ESXL set)
unique(joint_IDR_dat2[str_count(joint_IDR_dat2$uniprot_sequence, "P.P.P.P.P")>2,]$ensembl_gene_id)
#12 genes
