# Figure 7 Empirical inflation and inflation of the bootstrap MLE
# 1 -- Empirical inflation

# Compute the average MLE (this part is computed on the cluster)
dir_output <- ""

mle_names <- list.files(path = dir_output, pattern = "mle_")
indices <- sapply(mle_names, function(t) strsplit(t, "_")[[1]][2])

N <- length(indices)

mle <- matrix(NA, nrow = N, ncol = p)
for(i in 1:N){
  new_mle <- read.table(paste0(dir_output,"mle_",indices[i]))
  mle[i, ] <- new_mle[,1]
}

# Average MLE (Note I multiply both by sqrt(p), which corresponds to standardizing Xj tp have std 1/sqrt(p))
avgmle <- colMeans(mle) * sqrt(p)
# True coef
beta <- scan("beta.txt") * sqrt(p)


# store the results
result <- data.frame(avgmle = avgmle,  beta = beta)
write.table(result, "inflation-empirical.txt", row.names = F, col.names = F)

# Plot empirical inflation
tb <- read.table(paste0(loc, "inflation-empirical.txt"), header = F)
colnames(tb) <- c("avgmle", "beta")

a <- lm(avgmle ~ beta + 0, data = tb)$coef
ind <- which(tb$beta!=0)
fig_empirical_inflation <- ggplot(tb[ind, ]) +
  geom_point(aes(x = beta, y = avgmle)) +
  geom_abline(slope = a, intercept = 0, color = "red")+
  xlab("Coefficients") +
  ylab("Average MLE coefficients") +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave(plot = fig_empirical_inflation, 
       filename = paste0(figloc, "inflation-incorrect-empirical.png"),
       width = 8, height = 6, units = "cm", dpi = 300)

# 2 --  Inflation of the resized bootstrap MLE
# randomly pick one simulation, I pick the number 791
# Get beta_star
betas <- scan(paste0(loc, "betas_unknown.txt" )) * sqrt(p)
boot_sample <- read.table(paste0(loc, "boot_unknown.txt")) * sqrt(p)
avgboot <- rowMeans(boot_sample)
sd_boot <- apply(boot_sample, 1, sd)
aboot <- lm(avgboot~betas + 0, weights = 1/sd_boot^2)$coef

fig_boot_inflation <- ggplot() +
  geom_point(aes(x = betas, y = avgboot)) +
  geom_abline(slope = aboot, intercept = 0, color = "red")+
  xlab("Resized coefficients") +
  ylab("Average resized bootstrap MLE") +
  theme_bw() +
  theme(text = element_text(size = 10))
ggsave(plot = fig_boot_inflation, 
       filename = paste0(figloc, "inflation-incorrect-boot.png"),
       width = 8, height = 6, units = "cm", dpi = 300)




