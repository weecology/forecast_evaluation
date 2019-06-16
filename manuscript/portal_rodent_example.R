library(coda)
library(runjags)
library(portalr)
library(scoringRules)
library(e1071)
source("functions.R")

# only run if raw download is needed
# download_observations()

abunds <- plot_sp_abunds(plot = 17, species = "DM")
moon_dates <- load_trapping_data(clean = FALSE)$newmoons_table$newmoondate

summarize_abunds(abunds)

in_ts <- training_ts(ends = 10:500)

# only run if models have not yet been fit
# mod1 <- mods(1, abunds, moon_dates, in_ts)
# mod2 <- mods(2, abunds, moon_dates, in_ts)
# mod3 <- mods(3, abunds, moon_dates, in_ts)

load_models(1:3)
eval_tab <- make_eval_tab(mod1, mod2, mod3)


# working space

# fits of the models from the test origin 500 indicate that theta has very
# marginal effects and that phi is very close to 1
# which thus gives us three models that are quite similar to each other in fit
# and it appears predictions as well



rs <- which(eval_tab$model == 1 & eval_tab$lead == 1 & eval_tab$origin < 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m1 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

rs <- which(eval_tab$model == 2 & eval_tab$lead == 1 & eval_tab$origin < 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m2 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

rs <- which(eval_tab$model == 3 & eval_tab$lead == 1 & eval_tab$origin < 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m3 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)


par(mfrow = c(3, 1))
plot(PIT_m1, type = "h", ylim = c(0, max(PIT_m1) * 1.25))
plot(PIT_m2, type = "h", ylim = c(0, max(PIT_m1) * 1.25))
plot(PIT_m3, type = "h", ylim = c(0, max(PIT_m1) * 1.25))

# PITs are similar across the three and tends to be imprecise even at step 1
# same out at step 12 basically
# although gets more negatively biased out at step 12

mod1[[491]]$summary
mod2[[491]]$summary
mod3[[491]]$summary


par(mfrow=c(3,1), mar = c(2, 2, 1, 1))
blank(xlim = c(490, 515), ylim = c(0, 20), bty = "L")

m1 <- as.mcmc(combine.mcmc(as.mcmc.list(mod1[[491]]$model), 
                          collapse.chains = TRUE))
for(i in 501:512){
  yy1 <- rpois(nrow(m1), m1[,paste0("predY[", i, "]")])
  violin(yy1, location = i)
}
points(abunds)


blank(xlim = c(490, 515), ylim = c(0, 20), bty = "L")

m2 <- as.mcmc(combine.mcmc(as.mcmc.list(mod2[[491]]$model), 
                          collapse.chains = TRUE))
for(i in 501:512){
  yy2 <- rpois(nrow(m2), m2[,paste0("predY[", i, "]")])
  violin(yy2, location = i)
}
points(abunds)




blank(xlim = c(490, 515), ylim = c(0, 20), bty = "L")

m3 <- as.mcmc(combine.mcmc(as.mcmc.list(mod3[[491]]$model), 
                          collapse.chains = TRUE))
for(i in 501:512){
  yy3 <- rpois(nrow(m3), m3[,paste0("predY[", i, "]")])
  violin(yy3, location = i)
}
points(abunds)

i <- 512
yy1 <- rpois(nrow(m1), m1[,paste0("predY[", i, "]")])
yy2 <- rpois(nrow(m2), m2[,paste0("predY[", i, "]")])
yy3 <- rpois(nrow(m3), m3[,paste0("predY[", i, "]")])
plot(table(yy1), type = "h", xlim = c(0, 25))
abline(v = mean(yy1))
plot(table(yy2), type = "h", xlim = c(0, 25))
abline(v = mean(yy2))
plot(table(yy3), type = "h", xlim = c(0, 25))
abline(v = mean(yy3))
