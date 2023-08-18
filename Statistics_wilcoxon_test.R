library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)


sig_test <- function(data,factors,values){
  
  compare_list <- unique(data[[factors]]) %>% combn(2) 
  compare_list <- t(compare_list )
  
  wil_test <- data.frame("comp1" = c(), "comp2" = c(), "p.value" = c(), "position" = c() )
  for(i in 1:dim(compare_list)[1] ){
    a = data[data["meta_tissue"] == compare_list[i,][1],values]
    b = data[data["meta_tissue"] == compare_list[i,][2],values]
    
    wilcox_test = wilcox.test(x = a, y = b)
    p.value = round(wilcox_test$p.value, 5)
    
    position = max(a,b, na.rm = T)
    
    tmp = data.frame("comp1" = compare_list[i,][1], "comp2" = compare_list[i,][2], "p.value" = p.value ,  "position" = position)
    
    wil_test = rbind(wil_test, tmp)
  }
  
  return(wil_test)
  
}
