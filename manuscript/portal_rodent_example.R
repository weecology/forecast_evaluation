library(coda)
library(runjags)
library(portalr)
library(scoringRules)
library(e1071)
source("functions.R")

# only run if raw download is needed
# download_observations()

abunds <- plot_sp_abunds(plot = 19, species = "PP")
moon_dates <- load_trapping_data(clean = FALSE)$newmoons_table$newmoondate

summarize_abunds(abunds[200:length(abunds)])
fig4(abunds, moon_dates)

# only run if base models have not yet been fit
# in_ts <- training_ts(starts = 200, ends = 300:500)
# mod1 <- mods(1, abunds, moon_dates, in_ts)
# mod2 <- mods(2, abunds, moon_dates, in_ts)
# mod3 <- mods(3, abunds, moon_dates, in_ts)

# only run if stacking models have not yet been fit
# stack_df <- stack_mods(mods, 300:499, abunds, save = FALSE) 

load_models(1:3)
eval_tab <- make_eval_tab(mod1, mod2, mod3)

load_stack_df()

# working space



# building out Fig 5 and the ensembles


set.seed(321)


tiff("fig5.tiff", width = 6, height = 7, units = "in", res = 200)

maxv <- 5e3


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
  maxout <- which(yy1[,i-301] > maxv)
  yy1[maxout,i-301] <- maxv
}


mm1 <- as.mcmc(combine.mcmc(as.mcmc.list(mod1[[189]]$model), 
                          collapse.chains = TRUE))
yyy1 <- matrix(NA, nrow = nrow(mm1), ncol = 12)
for(i in 290:301){
  yyy1[,i-289] <- rpois(nrow(mm1), mm1[,paste0("predY[", i, "]")])
}


rs <- which(eval_tab$model == 1 & eval_tab$origin == 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m1 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy1, yyy1, PIT_m1, 1, title = "RW")


m2 <- as.mcmc(combine.mcmc(as.mcmc.list(mod2[[201]]$model), 
                          collapse.chains = TRUE))
yy2 <- matrix(NA, nrow = nrow(m2), ncol = 12)
for(i in 302:313){
  yy2[,i-301] <- rpois(nrow(m2), m2[,paste0("predY[", i, "]")])
}

mm2 <- as.mcmc(combine.mcmc(as.mcmc.list(mod2[[189]]$model), 
                          collapse.chains = TRUE))
yyy2 <- matrix(NA, nrow = nrow(mm2), ncol = 12)
for(i in 290:301){
  yyy2[,i-289] <- rpois(nrow(mm2), mm2[,paste0("predY[", i, "]")])
}


rs <- which(eval_tab$model == 2 & eval_tab$origin == 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m2 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy2, yyy2, PIT_m2, 2, title = "AR(1)")


m3 <- as.mcmc(combine.mcmc(as.mcmc.list(mod3[[201]]$model), 
                          collapse.chains = TRUE))
yy3 <- matrix(NA, nrow = nrow(m3), ncol = 12)
for(i in 302:313){
  yy3[,i-301] <- rpois(nrow(m3), m3[,paste0("predY[", i, "]")])
}

mm3 <- as.mcmc(combine.mcmc(as.mcmc.list(mod3[[189]]$model), 
                          collapse.chains = TRUE))
yyy3 <- matrix(NA, nrow = nrow(mm3), ncol = 12)
for(i in 290:301){
  yyy3[,i-289] <- rpois(nrow(mm3), mm3[,paste0("predY[", i, "]")])
}


rs <- which(eval_tab$model == 3 & eval_tab$origin == 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m3 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy3, yyy3, PIT_m3, 3, title = "cAR(1)")



yy4 <- rbind(yy1, yy2, yy3)
yyy4 <- rbind(yyy1, yyy2, yyy3)

rs <- which(eval_tab$origin == 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m4 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)

fig4_row(abunds, yy4, yyy4, PIT_m4, 4, title = "Naïve Ensemble")


dev.off()



  lead <- eval_tab$lead
  crps <- eval_tab$crps
  logs <- eval_tab$logs
  model <- eval_tab$model
  origin <- eval_tab$origin
  destin <- eval_tab$destin
  


m1 <-mean(logs[model==1&lead==12&destin<501], na.rm = TRUE)
m2 <-mean(logs[model==2&lead==12&destin<501], na.rm = TRUE)
m3 <-mean(logs[model==3&lead==12&destin<501], na.rm = TRUE)

(m2-m1)/(0-m1)
(m3-m1)/(0-m1)

m1 <-mean(crps[model==1&lead==12&destin<501], na.rm = TRUE)
m2 <-mean(crps[model==2&lead==12&destin<501], na.rm = TRUE)
m3 <-mean(crps[model==3&lead==12&destin<501], na.rm = TRUE)

(m2-m1)/(0-m1)
(m3-m1)/(0-m1)

mean(crps[model==1 & lead==1], na.rm = TRUE)
mean(crps[model==2 & lead==1], na.rm = TRUE)
mean(crps[model==3 & lead==1], na.rm = TRUE)

mean(logs[model==1 & lead==1], na.rm = TRUE)
mean(logs[model==2 & lead==1], na.rm = TRUE)
mean(logs[model==3 & lead==1], na.rm = TRUE)


mean(crps[model==1 & lead==12], na.rm = TRUE)
mean(crps[model==2 & lead==12], na.rm = TRUE)
mean(crps[model==3 & lead==12], na.rm = TRUE)

mean(logs[model==1 & lead==12], na.rm = TRUE)
mean(logs[model==2 & lead==12], na.rm = TRUE)
mean(logs[model==3 & lead==12], na.rm = TRUE)

plot(jitter(lead[origin>480], 1.5), crps[origin>480], col = model[origin>480])

plot(lead[origin==500], crps[origin==500], col = model[origin==500])


