##Identify_tandem_repeats
#9th of November 2017
#David Brown


#f) Identify tandem repeats using rle()
s<-str2 #example sequence

#Set minimum tandem repeat length, I have used '3'
ml<-3

#Convert aa sequence 's' to numeric vector 'n'
# rle(s) does not work as s is a single string
rn<-rle(unlist(strsplit(s,"")))

#Calculate the position of each run
rn$pos<-cumsum(rn$lengths)

#Convert to data frame 'rep'
rep<-data.frame(rn$values, rn$pos, rn$pos+rn$lengths,rn$lengths)
names(rep)<-c("aa", "rep_start","rep_end","rep_width")

#Filter to tandem repeats above a chosen length
rep<-rep[rep$rep_width>ml,]

#Wrap the whole process in a function
trs<- function(x, ml) {
  rn<-rle(unlist(strsplit(s,"")))
  rn$pos<-cumsum(rn$lengths)
  rep<-data.frame(rn$values, rn$pos, rn$pos+rn$lengths,rn$lengths)
  names(rep)<-c("aa", "rep_start","rep_end","rep_width")
  rep<-rep[rep$rep_width>ml,]
}

#Usage
rep <- trs(s, ml)

#Is the a function to report the presence of a tandem repeat? Good place to start
tr<-function(x, ml){
  rn<-rle(unlist(strsplit(x,"")))
  trs<-max(rn$lengths)>ml
  trs_number<-sum(rn$lengths>ml)
  max_rep_width<-max(rn$lengths)[1]   #There is bias here, as the first of equal maxima will be taken
  max_rep_aa<-rn$values[rn$lengths==max_rep_width][1]
  return(list("tandem_repeat"=trs, "tandem_repeat_number"= trs_number, "max_repeat_width" = max_rep_width, "max_repeat_aa" =max_rep_aa))
}


#Usage
tr(x, ml)
test<-tr(s, 3)

#Batch screen for tandem repeat sequences
#transpose for compatibility with hi2 data.frame
t(sapply(hi2$IDR_seq, function(x) tr(x,3), USE.NAMES = FALSE))


