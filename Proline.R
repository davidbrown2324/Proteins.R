##Proline
library("Biostrings", lib.loc="~/R/win-library/3.2")

s<-as.character(unique(foci_pep$uniprot_sequence))
s<-AAStringSet(s)


#Proline content

#Proline spacing

proline_spacing_fun<-function(s){
  IDR_composition<-alphabetFrequency(s)
  P_cont<-IDR_composition[,"P"]
  P_frac<-P_cont/width(s)
  search<-vmatchPattern("P", s)
  P_space<-diff(start(search))
  list(width(s), P_cont, P_frac, mean(P_space), sd(P_space))
}

pro<-data.frame(proline_spacing_fun(s))
names(pro)<-c("IDR_width", "Proline_Count", "Proline_Fraction", "Mean_Spacing", "Spacing_StdDev")

clean<-pro[pro$IDR_width>30 & pro$IDR_width<100 & pro$Proline_Fraction>0.3,]
plot(clean$Proline_Fraction, clean$Mean_Spacing)
plot(clean$Proline_Count, clean$Spacing_StdDev)

head(search, 2)
head(start(search), 2)
head(diff(start(search)), 2)
head(mode(diff(start(search))), 2)

head(s, 10)

mean(P_space)
