
library(tidyverse)
library(readxl)
library(tinytex)

pdf(file = "Noise_distributions.pdf")
X129_IGFDose_T0Norm_Noise <- read_excel("129_IGFDose_T0Norm_Noise.xlsx") %>%
  rename(IGF_I_0_pM = "1",
        IGF_I_5_pM = "2",
        IGF_I_10_pM = "3",
        IGF_I_15_pM = "4",
        IGF_I_25_pM = "6",
        IGF_I_40_pM = "8",
        IGF_I_125_pM = "11")


qqplots <- function(x, col_name){
  title_text <- str_replace(col_name, "_", "-") %>%
    str_replace("_"," ") %>%
    str_remove("_") 
  qqnorm(x, pch = 20, main=paste("Normal QQ plot of", title_text))
  qqline(x, col = "blue")
  hist(x, main=paste("Histogram of",title_text), breaks = 100)

}

res <- map2(X129_IGFDose_T0Norm_Noise, names(X129_IGFDose_T0Norm_Noise), qqplots)

dev.off()


