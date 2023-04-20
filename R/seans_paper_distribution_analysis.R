
library(tidyverse)
library(readxl)
library(tinytex)

# pdf(file = "Noise_distributions_combined.pdf")
# X129_IGFDose_T0Norm_Noise <- read_excel("129_IGFDose_T0Norm_Noise.xlsx") %>%
#   rename(IGF_I_0_pM = "1",
#         IGF_I_5_pM = "2",
#         IGF_I_10_pM = "3",
#         IGF_I_15_pM = "4",
#         IGF_I_25_pM = "6",
#         IGF_I_40_pM = "8",
#         IGF_I_125_pM = "11")
# df <- gather(X129_IGFDose_T0Norm_Noise, key = "Dosage", value = "Noise") %>%
#   mutate(Dosage = factor(Dosage,ordered = TRUE,levels = c("IGF_I_0_pM", "IGF_I_5_pM",   "IGF_I_10_pM", "IGF_I_15_pM",  "IGF_I_25_pM", "IGF_I_40_pM", "IGF_I_125_pM")))
# unique(df$Dosage)
# cols <- c("IGF_I_0_pM" = "#1a1919",    "IGF_I_5_pM" = "#253781",   "IGF_I_10_pM"  = "#283378", "IGF_I_15_pM" = "#009647",  "IGF_I_25_pM"  = "#e25e28", "IGF_I_40_pM"  = "#ee2f36","IGF_I_125_pM" = "#3b3b38")
# p <- ggplot(df, aes(sample = Noise, colour = Dosage)) +
#   stat_qq(alpha = .6) +
#   stat_qq_line() +  
#   scale_colour_manual(
#     values = cols,
#     aesthetics = c("colour")
#   ) +
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         panel.grid = element_blank(),
#         axis.line = element_line(colour = "black"))
# 
# plot(p)


pdf(file = "Signal_distributions_combined.pdf")
X129_IGFDose_T0Norm_Signal <- read_excel("129_IGFDose_T0Norm_Signal90.xlsx") %>%
  rename(IGF_I_0_pM = "1",
         IGF_I_5_pM = "2",
         IGF_I_10_pM = "3",
         IGF_I_15_pM = "4",
         IGF_I_25_pM = "6",
         IGF_I_40_pM = "8",
         IGF_I_125_pM = "11")
df <- gather(X129_IGFDose_T0Norm_Signal, key = "Dosage", value = "Signal") %>%
  mutate(Dosage = factor(Dosage,ordered = TRUE,levels = c("IGF_I_0_pM", "IGF_I_5_pM",   "IGF_I_10_pM", "IGF_I_15_pM",  "IGF_I_25_pM", "IGF_I_40_pM", "IGF_I_125_pM")))
unique(df$Dosage)
cols <- c("IGF_I_0_pM" = "#1a1919",    "IGF_I_5_pM" = "#253781",   "IGF_I_10_pM"  = "#283378", "IGF_I_15_pM" = "#009647",  "IGF_I_25_pM"  = "#e25e28", "IGF_I_40_pM"  = "#ee2f36","IGF_I_125_pM" = "#3b3b38")
p <- ggplot(df, aes(sample = Signal, colour = Dosage)) +
  stat_qq(alpha = .6) +
  stat_qq_line() +  
  scale_colour_manual(
    values = cols,
    aesthetics = c("colour")
  ) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"))

plot(p)
# p <- ggplot(df, aes(sample = Noise)) +
#   stat_qq(alpha = .5) +
#   stat_qq_line() +
#   theme_bw() +
#   facet_wrap(~Dosage)
# plot(p)
# 
# p <- ggplot(df, aes(x = Noise)) +
# geom_histogram(bins = 100) +
#   theme_bw() +
#   facet_wrap(~Dosage)
# plot(p)


# qqplots <- function(x, col_name){
#   title_text <- str_replace(col_name, "_", "-") %>%
#     str_replace("_"," ") %>%
#     str_remove("_") 
#   qqnorm(x, pch = 20, main=paste("Normal QQ plot of", title_text))
#   qqline(x, col = "blue")
#   hist(x, main=paste("Histogram of",title_text), breaks = 100)
# 
# }
# 
# res <- map2(X129_IGFDose_T0Norm_Noise, names(X129_IGFDose_T0Norm_Noise), qqplots)

dev.off()


