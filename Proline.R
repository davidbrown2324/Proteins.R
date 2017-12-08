##Proline
library("Biostrings", lib.loc="~/R/win-library/3.2")

s<-as.character(unique(foci_pep$uniprot_sequence))
s<-AAStringSet(s)

#Define the proline spacing function, which returns 'IDR_width', 'Proline_Count', 'Proline_Fraction', 'Mean_Spacing' and 'Spacing_StdDev'
proline_spacing_fun<-function(s){
  IDR_composition<-alphabetFrequency(s)
  P_cont<-IDR_composition[,"P"]
  P_frac<-P_cont/width(s)
  search<-vmatchPattern("P", s)
  P_space<-diff(start(search))
  list(width(s), P_cont, P_frac, mean(P_space), sd(P_space))
}

#Call 'proline_spacing_fun' on a sequence 's'
pro<-data.frame(proline_spacing_fun(s))
names(pro)<-c("IDR_width", "Proline_Count", "Proline_Fraction", "Mean_Spacing", "Spacing_StdDev")

#Subset proteins between 30 and 100 aa with >30% prolines
clean<-pro[pro$IDR_width>30 & pro$IDR_width<100 & pro$Proline_Fraction>0.3,]
plot(clean$Proline_Fraction, clean$Mean_Spacing)
plot(clean$Proline_Count, clean$Spacing_StdDev)
