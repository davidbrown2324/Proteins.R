##Intrinsically_Disordered_Region_Bioinformatics_v4
#9th of November 2017
#David Brown

#v4 starts with the full list of D2P2 IDRs 'm'.
  
#Set up 1 - load IDR data
#m<-read.table("Y:/DavidB/Data/Systematic_Search_for_Phase_Separating_Proteins/p2p_conversion/D2P2_IDR_predictions_and_ensembl_peptide_id.txt", sep="\t", header=TRUE)

#Set up 2 - Fetch sequence information
  ##Load the BiomaRt package
  library(biomaRt)
  
  #Define mart and data set
  mart <- useDataset(dataset = "hsapiens_gene_ensembl", 
                     mart    = useMart("ENSEMBL_MART_ENSEMBL",
                                       host    = "www.ensembl.org"))
  
  #Fetch peptide sequences
  hpep <- getBM(attributes = c("peptide", "ensembl_peptide_id"),
                   filters = "ensembl_peptide_id", values = hd$ensembl_peptide_id,
                   mart = mart)
  
  hpep$Gene<-hd$hox.Gene.Name
  
  #Convert sequence data into AASet - Not sure if this is useful at this stage
  library(Biostrings)
  
  h_aa<-AAStringSet(as.character(hpep$peptide))
  names(h_aa)<-hpep$ensembl_peptide_id
  #Keep AAStringSet as a separate class for region based analysis!!
  
#Step 1 - Load IDR data (start with the N/C-terminals of HPA curated homeodomain proteins as an example dataset)

  #Use D2P2 IDR data, labelled with ensemble peptide id
  m<-read.table("Y:/DavidB/Data/Systematic_Search_for_Phase_Separating_Proteins/p2p_conversion/D2P2_IDR_predictions_and_ensembl_peptide_id.txt", sep="\t", header=TRUE)
  
  #Annotate sequences of interest with D2P2 data
  hi5<-merge(dis_pep, m, by="ensembl_peptide_id", all.x=TRUE)
  hi5$protein_length<-nchar(hi5$peptide)
  
# Step 2 - Generate IDR dataframe 'hi5' this is a list of all IDRs in our disordered protein set. 
  
  #a) extract IDR sequence
  hi5$IDR_seq<-substr(hi5$peptide, hi5$start, hi5$end)
  substr(hi5[1010,]$peptide, hi5[1010,]$start, hi5[1010,]$end)
  
  nchar(hi5[1010,]$peptide) # only 758aa
  #This calls into question the coherence of the data set!!
  hi5<-hi5[!hi5$IDR_seq=="",]

  #b) determine IDR position. I am counting IDRs within 3 aa of the termini as terminal!
  hi5$IDR_position<-"Unknown"
  hi5[hi5$start<3,]$IDR_position<-"N-terminal"
  hi5[hi5$end>(hi5$protein_length-3),]$IDR_position<-"C-terminal"
  hi5[hi5$IDR_position=="Unknown",]$IDR_position<-"Internal"
  #hi[Nter & Cter,]$IDR_position<-"IDP"  ##Sanity check that IDRs are not IDPs covering the whole sequence
  table(hi$IDR_position)
  
  #c) calculate the proline content/ count each aa
  IDR_set_aa<-AAStringSet(as.character(hi5$IDR_seq))
  names(IDR_set_aa)<-hi5$ensembl_peptide_id
  
  IDR_set_composition<-alphabetFrequency(IDR_set_aa)  #IDR_composition is a matrix and must be accessed via 'colnames' not 'names'
  hi5$proline_content<-IDR_set_composition[,"P"]/hi5$width
  hi5$proline_rich<-hi5$proline_content>0.15
  
  #d) calculate charged fraction
  library(Peptides)
  
  IDR_set_aa_properties<-aaComp(IDR_set_aa)            #IDR_aa_properties is a list

  ##This took a long time to figure out
  #sapply function works on each component of the list, to look up values yielding a named numeric vector, 
  #as.vector simplifies to an unamed vector
  
  #tiny<-as.vector(sapply(IDR_aa_properties, function(x) x["Tiny","Mole%"]))
  
  ##With a few adjustments (generalised to all properties, transposed for compatibility with dataframe, and converted from percent to fraction)
  IDR_set_prop<-t(sapply(IDR_set_aa_properties, function(x) x[,"Mole%"]))/100
  #This gives the properties of all 31,741 IDRs in our disordered protein set.

  ##Reorder matrix as desired
  IDR_set_props_srt <- IDR_set_prop[, order(colnames(IDR_set_prop))]
  
  hi6<-cbind(hi5,IDR_set_props_srt) 
  write.table(hi6, file="Y:/DavidB/Data/Systematic_Search_for_Phase_Separating_Proteins/IDR_set_aa_properties.txt", sep="\t", row.names = FALSE)
  #hi6<-read.table("Y:/DavidB/Data/Systematic_Search_for_Phase_Separating_Proteins/IDR_set_aa_properties.txt", sep="\t", header = TRUE) 
  