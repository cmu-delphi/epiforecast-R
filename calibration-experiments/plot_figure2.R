library(dplyr)
library(ggplot2)

df = tibble(x=seq(-3,3,by=0.01)) %>%
  mutate(h = dnorm(x), f = dnorm(x,mean=0,sd=2), g = h/f) 

ggplot(df) +
  geom_line(aes(x=x,y=f,color="f(y)")) +
  geom_line(aes(x=x,y=g,color="g(F(y))")) +
  geom_line(aes(x=x,y=h,color="g(F(y))*f(y)"),size=2) +
  geom_line(aes(x=x,y=h,color="h(y)"),size=0.75) +
  labs(x="y",y="Density",color="Function") +
  theme(aspect.ratio=0.75)

ggsave("fig2.png",width=8,height=6)
