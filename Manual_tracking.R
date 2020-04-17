library(tidyverse)
df_input <- read_csv("MTrackJ/MTrackJ- Points.csv") %>%
  as_tibble(.name_repair = "universal")

T1 <- df_input %>%
  filter(PID == 1)
  mutate(x_rel = )
p <- ggplot(df, aes(x..pixel., y..pixel., colour = factor(TID))) +
  geom_point(size = .8, alpha = .8) +
  theme_bw()
p

foo <- df %>%
  select(TID, PID, I..val., )

