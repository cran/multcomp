# R version-1.7.0 on MacOSX ,and R-1.7.1 on WinNT
#multcomp version 0.3-11 ,mvtnorm version 0.6-3
library(multcomp)
Responses<-c(0,0,0,-0.3333333,0,0,0,0,0,0,0,0.3333333,-0.6666667,- 
0.3333333,0,0,0,0.3333333,0,0,-0.3333333,0,0,0,0,0,0,0,0,-0.3333333,0,- 
0.3333333,0,0,-0.3333333,0,0,0,0,0,0,0,0,0,0.6666667)
Events<- c("b","b","b","b","b","b","c","c","c","c","c","d","d","d","d","d","d","e","e","e","e","e","e","f","f","f","f","f","f","g","g","g","g","g","h","h","h","h","h","a","a","a","a","a","a")
testdata<-cbind(as.data.frame(Responses),Events)
simtest(Responses~Events,data=testdata,type="Dunnett")
