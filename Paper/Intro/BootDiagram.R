library(tidyverse)

phase0 <- read.table("~/Documents/GitHub/logisticMLE/paper/section-5/code/phase0.txt", quote="\"", comment.char="")
colnames(phase0) <- c("gamma", "kappa")
# Plot the phase transition diagram

phase0Plot <- rbind(phase0, c(20, 0))
g <- ggplot(phase0Plot) + 
  geom_line(aes(x = kappa, y = gamma)) + 
  geom_ribbon(fill = "cyan", aes(x = kappa, ymax = gamma), ymin = -0.5) +
  coord_cartesian(xlim=c(0.05, 0.5),
                  ylim=c(0.8, 19)) +
  scale_y_continuous(breaks=seq(0, 15, 5)) + 
  xlab(expression(kappa)) + 
  ylab(expression(gamma)) + 
  theme_bw() + 
  theme(text = element_text(size = 18))
filename <- "/Users/zq/Documents/Simulation_Data/glm/glm_boot/fig/phase_diagram.png"
ggsave(filename = filename, plot = g, 
       units = "in",
       width = 6, height = 5,
       dpi = 150)

set.seed(143)
series <- data.frame(
  time = c(rep(2017, 4),rep(2018, 4), rep(2019, 4), rep(2020, 4)),
  type = rep(c('a', 'b', 'c', 'd'), 4),
  value = rpois(16, 10)
)


ggplot(series, aes(time, value)) +
  geom_area(aes(fill = type)) + 
  # coord_cartesian(xlim=c(2017.8, 2020)) +
  scale_x_continuous(breaks=seq(2018, 2021, 1))


