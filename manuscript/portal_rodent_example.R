library(coda)
library(runjags)
library(portalr)
library(scoringRules)
library(e1071)
source("functions.R")

# only run if raw download is needed
# download_observations()

abunds <- plot_sp_abunds(plot = 19, species = "PP")
plot(abunds, xlim = c(200, 520), type = "l")

moon_dates <- load_trapping_data(clean = FALSE)$newmoons_table$newmoondate

summarize_abunds(abunds[200:length(abunds)])

in_ts <- training_ts(starts = 200, ends = 300:500)

# only run if models have not yet been fit
 mod1 <- mods(1, abunds, moon_dates, in_ts)
 mod2 <- mods(2, abunds, moon_dates, in_ts)
 mod3 <- mods(3, abunds, moon_dates, in_ts)

load_models(1:3)
eval_tab <- make_eval_tab(mod1, mod2, mod3)


# working space
# building out Fig 4 and the ensemble


set.seed(321)


tiff("fig4.tiff", width = 6, height = 7, units = "in", res = 200)



topoffset <- 0.94
par(fig= c(0, 0.5, topoffset, 1), mar = c(0, 0, 0, 0))
blank()
text(x = 1, y = 1.15, "Predictive distributions and true counts of", 
     cex = 0.8)
text(x = 1, y = 0.825, "C. penicillatus", cex = 0.8, font = 3)
par(fig= c(0.5, 0.75, topoffset, 1), new = TRUE)
blank()
text(x = 1, y = 1.15, "Observed (y) vs.", cex = 0.8, font = 1)
text(x = 1, y = 0.825, "Predicted (x)", cex = 0.8, font = 1)
par(fig= c(0.75, 1, topoffset, 1), new = TRUE)
blank()
text(x = 1, y = 1.15, "Probability Integral", cex = 0.8, font = 1)
text(x = 1, y = 0.825, "Transform (PIT)", cex = 0.8, font = 1)

m1 <- as.mcmc(combine.mcmc(as.mcmc.list(mod1[[201]]$model), 
                          collapse.chains = TRUE))
yy1 <- matrix(NA, nrow = nrow(m1), ncol = 12)
for(i in 302:313){
  yy1[,i-301] <- rpois(nrow(m1), m1[,paste0("predY[", i, "]")])
}

rs <- which(eval_tab$model == 1)
cs <- grep("PIT", colnames(eval_tab))
PIT_m1 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy1, PIT_m1, 1, title = "RW")


m2 <- as.mcmc(combine.mcmc(as.mcmc.list(mod2[[201]]$model), 
                          collapse.chains = TRUE))
yy2 <- matrix(NA, nrow = nrow(m2), ncol = 12)
for(i in 302:313){
  yy2[,i-301] <- rpois(nrow(m2), m2[,paste0("predY[", i, "]")])
}

rs <- which(eval_tab$model == 2)
cs <- grep("PIT", colnames(eval_tab))
PIT_m2 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy2, PIT_m2, 2, title = "AR(1)")


m3 <- as.mcmc(combine.mcmc(as.mcmc.list(mod3[[201]]$model), 
                          collapse.chains = TRUE))
yy3 <- matrix(NA, nrow = nrow(m3), ncol = 12)
for(i in 302:313){
  yy3[,i-301] <- rpois(nrow(m3), m3[,paste0("predY[", i, "]")])
}

rs <- which(eval_tab$model == 3)
cs <- grep("PIT", colnames(eval_tab))
PIT_m3 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy3, PIT_m3, 3, title = "cAR(1)")



yy4 <- rbind(yy1, yy2, yy3)

cs <- grep("PIT", colnames(eval_tab))
PIT_m4 <- apply(eval_tab[, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy4, PIT_m4, 4, title = "Naïve Ensemble")


dev.off()

