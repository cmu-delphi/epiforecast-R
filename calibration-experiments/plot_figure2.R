library(dplyr)
library(ggplot2)

df = tibble(x=seq(-3,3,by=0.01)) %>%
  mutate(h = dnorm(x), f = dnorm(x,mean=0,sd=2), g = h/f) 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(4)

ggplot(df) +
  geom_line(aes(x=x,y=f,color="f(y)")) +
  geom_line(aes(x=x,y=g,color="g(F(y))")) +
  geom_line(aes(x=x,y=h,color="g(F(y))*f(y)")) +
  geom_line(aes(x=x,y=h,color="h(y)"),linetype="dashed") +
  labs(x="y",y="Density",color="Function") +
  scale_color_manual(values=c(cols[4],cols[3],cols[2],cols[1])) +
  theme(aspect.ratio=0.75)

ggsave("fig2.png",width=8,height=6)
