

dat <- data.frame(a=1:75, b=rep(c("apple", "banana", "carrot"), each=25), c=runif(75))
library(ggplot2)
library(dplyr)
dat %>%
  mutate(fact_b=factor(b, levels=c("banana", "carrot", "apple"))) %>%
  ggplot() +
  geom_boxplot(aes(x=fact_b, y=c))
