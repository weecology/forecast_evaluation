library(coda)
library(runjags)
library(portalr)
library(scoringRules)
library(e1071)
library(nnet)
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

load_models(1:3)
modl <- list(mod1, mod2, mod3)

# only run if stacking models have not yet been fit
# stack1 <- stack_mods(modl, abunds, rule = "crps", linkname = "identity")
# stack2 <- stack_mods(modl, abunds, rule = "crps", linkname = "beta1")
# stack3 <- stack_mods(modl, abunds, rule = "crps", linkname = "beta2")
# stack4 <- stack_mods(modl, abunds, rule = "logscore", linkname = "identity")
# stack5 <- stack_mods(modl, abunds, rule = "logscore", linkname = "beta1")
# stack6 <- stack_mods(modl, abunds, rule = "logscore", linkname = "beta2")


load("model_output/stacks.RData")

eval_tab <- make_eval_tab(mod1, mod2, mod3)


# working space

wtses <- matrix(NA, nrow = 7, ncol = 3)
wtses[1,] <- round(rep(1/3, 3), 4)
wtses[2,] <- round(ialr(stack1$par[1:2]), 4)
wtses[3,] <- round(ialr(stack2$par[1:2]), 4)
wtses[4,] <- round(ialr(stack3$par[1:2]), 4)
wtses[5,] <- round(ialr(stack4$par[1:2]), 4)
wtses[6,] <- round(ialr(stack5$par[1:2]), 4)
wtses[7,] <- round(ialr(stack6$par[1:2]), 4)

write.csv(wtses, file = "model_output/model_weights.csv", row.names = FALSE)




  set.seed(123)
  nsamps <- 10000
  nsteps <- 12
  nmods <- 3
  origins <- 300:500
  ntests <- length(origins)
  pdists <- array(NA, dim = c(nsamps, nsteps, nmods, ntests))
  for(i in 1:ntests){
    inp <- draw_predictive_dists(modl, origin = origins[i], 
                                          n = nsamps)
    pdists[,1:dim(inp)[2],,i] <- inp
  }

  ensdist_wn <- combine(pdists, rep(1/3, 3), "cdf")
  tensdist_wn <- link(ensdist_wn, "identity")
  ensdraw_wn <- array(NA, dim = c(nsamps, nsteps, ntests)) 
  for(i in 1:ntests){
    ensdraw_wn[,,i] <- draw_from_cdftl(nsamps, tensdist_wn[[i]])
  }

  ensdist_wci <- combine(pdists, ialr(stack1$par[1:2]), "cdf")
  tensdist_wci <- link(ensdist_wci, "identity")
  ensdraw_wci <- array(NA, dim = c(nsamps, nsteps, ntests)) 
  for(i in 1:ntests){
    ensdraw_wci[,,i] <- draw_from_cdftl(nsamps, tensdist_wci[[i]])
  }

  ensdist_wcb1 <- combine(pdists, ialr(stack2$par[1:2]), "cdf")
  tensdist_wcb1 <- link(ensdist_wcb1, "beta1", stack2$par)
  ensdraw_wcb1 <- array(NA, dim = c(nsamps, nsteps, ntests)) 
  for(i in 1:ntests){
    ensdraw_wcb1[,,i] <- draw_from_cdftl(nsamps, tensdist_wcb1[[i]])
  }

  ensdist_wcb2 <- combine(pdists, ialr(stack3$par[1:2]), "cdf")
  tensdist_wcb2 <- link(ensdist_wcb2, "beta2", stack3$par)
  ensdraw_wcb2 <- array(NA, dim = c(nsamps, nsteps, ntests)) 
  for(i in 1:ntests){
    ensdraw_wcb2[,,i] <- draw_from_cdftl(nsamps, tensdist_wcb2[[i]])
  }

  ensdist_wli <- combine(pdists, ialr(stack4$par[1:2]), "pdf")
  tensdist_wli <- link(ensdist_wli, "identity")
  ensdraw_wli <- array(NA, dim = c(nsamps, nsteps, ntests)) 
  for(i in 1:ntests){
    ensdraw_wli[,,i] <- draw_from_pdftl(nsamps, tensdist_wli[[i]])
  }

  ensdist_wlb1 <- combine(pdists, ialr(stack5$par[1:2]), "pdf")
  tensdist_wlb1 <- link(ensdist_wlb1, "beta1", stack5$par, pdists)
  ensdraw_wlb1 <- array(NA, dim = c(nsamps, nsteps, ntests)) 
  for(i in 1:ntests){
    ensdraw_wlb1[,,i] <- draw_from_pdftl(nsamps, tensdist_wlb1[[i]])
  }

  ensdist_wlb2 <- combine(pdists, ialr(stack6$par[1:2]), "pdf")
  tensdist_wlb2 <- link(ensdist_wlb2, "beta2", stack6$par, pdists)
  ensdraw_wlb2 <- array(NA, dim = c(nsamps, nsteps, ntests)) 
  for(i in 1:ntests){
    ensdraw_wlb2[,,i] <- draw_from_pdftl(nsamps, tensdist_wlb2[[i]])
  }

dim(ensdraw_wn)


head(eval_tab)
ensdrawl <- list(ensdraw_wn = ensdraw_wn, ensdraw_wci = ensdraw_wci,
                 ensdraw_wcb1 = ensdraw_wcb1, ensdraw_wcb2 = ensdraw_wcb2, 
                 ensdraw_wli = ensdraw_wli, ensdraw_wlb1 = ensdraw_wlb1,
                 ensdraw_wlb2 = ensdraw_wlb2)


eeval_tab <- make_ens_eval_tab(abunds, origins, lead_time, ensdrawl)

etabs <- rbind(eval_tab, eeval_tab)

n_bins <- 10
last_in <- 500
umods <- unique(etabs$model)
numods <- length(umods)
PITs <- matrix(NA, nrow = numods, ncol = n_bins)
for(i in 1:numods){
  incl <- which(etabs$model == umods[i] & etabs$destin <= last_in)
  etabsi <- etabs[incl,]
  PITcols <- grepl("PIT", colnames(etabsi))
  etabsiPIT <- etabsi[,PITcols]
  for(j in 1:ncol(etabsiPIT)){
    etabsiPIT[,j] <- as.numeric(etabsiPIT[,j])
  }
  PITsum <- apply(etabsiPIT, 2, sum, na.rm = TRUE)
  PITmean <- PITsum / nrow(etabsiPIT)
  PITs[i,] <- PITmean
}

crpss <- rep(NA, numods)
logss <- rep(NA, numods)
for(i in 1:numods){
  incl <- which(etabs$model == umods[i] & etabs$destin <= last_in)
  etabsi <- etabs[incl,]
  crpscol <- 5
  logscol <- 6
  crpss[i] <- mean(as.numeric(etabsi[,5]), na.rm = TRUE)
  logss[i] <- mean(as.numeric(etabsi[,6]), na.rm = TRUE)
}

# ok now get the main model draws, and re-evaluate everything together

tiff("fig5.tiff", width = 6, height = 7, units = "in", res = 200)


mod_ids <- umods
mod_names <- c("RW", "AR(1)", "cAR(1)", "Naïve Ensemble",
               paste0("B", 0:2, "RS Ensemble"),
               paste0("B", 0:2, "LS Ensemble"))

par(mfrow=c(4,3), mar = c(1, 1.5, 2, 1))
for(i in 1:numods){
  blank(bty = "L", ylim = c(0, 1.75), xlim = c(0.5, 10.5))
  points(PITs[i,], type = "h", lwd = 3)
  abline(h = 1, lty = 3, lwd = 2)
  text(1, 1.7, mod_names[i], cex = 1, adj = 0)
}

max_crpss <- which.max(crpss)
max_logss <- which.max(logss)

par(mar=c(3, 3, 2, 1))

mod_names2 <- gsub(" Ensemble", "", mod_names)
blank(bty = "n", ylim = c(-2.5, -1.5), xlim = c(0.5, 10.5))
points(max_crpss, crpss[max_crpss], cex = 2.5, pch = 16, lwd = 0, 
       col = rgb(0.6, 0.6, 0.6, 0.6))
points(max_crpss, crpss[max_crpss], cex = 2.5, pch = 1, lwd = 1, 
       col = rgb(0.6, 0.6, 0.6))
points(crpss, pch = 16)
points(c(0, 0), c(-2.65, -1.4), lwd = 1, type = "l", xpd = TRUE)
points(c(0, 10.5), c(-2.65, -2.65), lwd = 1, type = "l", xpd = TRUE)
text(1:10, -2.7, mod_names2, xpd = TRUE, srt = 45, cex = 0.55, adj = 1)
mtext(side = 2, "Ranked Probability Score", line = 1, cex = 0.6)



blank(bty = "n", ylim = c(-2.5, -1.9), xlim = c(0.5, 10.5))
points(max_logss, logss[max_logss], cex = 2.5, pch = 16, lwd = 0, 
       col = rgb(0.6, 0.6, 0.6, 0.6))
points(max_logss, logss[max_logss], cex = 2.5, pch = 1, lwd = 1, 
       col = rgb(0.6, 0.6, 0.6))
points(logss, pch = 16)
points(c(0, 0), c(-2.575, -1.85), lwd = 1, type = "l", xpd = TRUE)
points(c(0, 10.5), c(-2.575, -2.575), lwd = 1, type = "l", xpd = TRUE)
text(1:10, -2.6, mod_names2, xpd = TRUE, srt = 45, cex = 0.55, adj = 1)
mtext(side = 2, "Log Score", line = 1, cex = 0.6)
dev.off()


etabsft <- etabs[etabs$origin == 500,]
etabsft_PIT <- matrix(NA, nrow = numods, ncol = n_bins)
etabsft_c <- rep(NA, numods)
etabsft_l <- rep(NA, numods)

for(i in 1:numods){
  incl <- which(etabsft$model == umods[i])
  PITcol <- grepl("PIT", colnames(etabsft))
  x <- etabsft[incl, PITcol]
  for(j in 1:n_bins){
    x[, j] <- as.numeric(x[,j])
  }

  etabsft_PIT[i, ] <- apply(x, 2, sum) / nrow(x)
  etabsft_c[i] <- mean(as.numeric(etabsft$crps[incl])) 
  etabsft_l[i] <- mean(as.numeric(etabsft$logs[incl])) 
}

which.max(etabsft_c)
which.max(etabsft_l)


head(eval_tab)





tiff("fig6.tiff", width = 6, height = 7, units = "in", res = 200)


mod_ids <- umods
mod_names <- c("RW", "AR(1)", "cAR(1)", "Naïve Ensemble",
               paste0("B", 0:2, "RS Ensemble"),
               paste0("B", 0:2, "LS Ensemble"))

par(mfrow=c(4,3), mar = c(1, 1.5, 2, 1))
for(i in 1:numods){
  blank(bty = "L", ylim = c(0, 4), xlim = c(0.5, 10.5))
  points(etabsft_PIT[i,], type = "h", lwd = 3)
  abline(h = 1, lty = 3, lwd = 2)
  text(1, 3.95, mod_names[i], cex = 1, adj = 0, xpd = TRUE)
}


max_etabsft_c <- which.max(etabsft_c)
max_etabsft_l <- which.max(etabsft_l)

par(mar=c(3, 3, 2, 1))

mod_names2 <- gsub(" Ensemble", "", mod_names)
blank(bty = "n", ylim = c(-2.5, -1.5), xlim = c(0.5, 10.5))
points(max_etabsft_c, etabsft_c[max_etabsft_c], cex = 2.5, pch = 16, 
       lwd = 0, col = rgb(0.6, 0.6, 0.6, 0.6))
points(max_etabsft_c, etabsft_c[max_etabsft_c], cex = 2.5, pch = 1, 
       lwd = 1, col = rgb(0.6, 0.6, 0.6))
points(etabsft_c, pch = 16)
points(c(0, 0), c(-2.65, -1.4), lwd = 1, type = "l", xpd = TRUE)
points(c(0, 10.5), c(-2.65, -2.65), lwd = 1, type = "l", xpd = TRUE)
text(1:10, -2.7, mod_names2, xpd = TRUE, srt = 45, cex = 0.55, adj = 1)
mtext(side = 2, "Ranked Probability Score", line = 1, cex = 0.6)

blank(bty = "n", ylim = c(-2.7, -2.15), xlim = c(0.5, 10.5))
points(max_etabsft_l, etabsft_l[max_etabsft_l], cex = 2.5, pch = 16, lwd = 0, 
       col = rgb(0.6, 0.6, 0.6, 0.6))
points(max_etabsft_l, etabsft_l[max_etabsft_l], cex = 2.5, pch = 1, lwd = 1, 
       col = rgb(0.6, 0.6, 0.6))
points(etabsft_l, pch = 16)
points(c(0, 0), c(-2.775, -2.05), lwd = 1, type = "l", xpd = TRUE)
points(c(0, 10.5), c(-2.775, -2.775), lwd = 1, type = "l", xpd = TRUE)
text(1:10, -2.8, mod_names2, xpd = TRUE, srt = 45, cex = 0.55, adj = 1)
mtext(side = 2, "Log Score", line = 1, cex = 0.6)
dev.off()















# to do
#   figure out what of them to graph
#   put summary stats together for all of the models for table 4

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

rs <- which(eval_tab$model == 1 & eval_tab$origin == 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m1 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)
fig5_row(abunds, pdists[,,1], PIT_m1, 1, title = "RW")

rs <- which(eval_tab$model == 2 & eval_tab$origin == 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m2 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)
fig5_row(abunds, pdists[,,2], PIT_m2, 2, title = "AR(1)")

rs <- which(eval_tab$model == 3 & eval_tab$origin == 500)
cs <- grep("PIT", colnames(eval_tab))
PIT_m3 <- apply(eval_tab[rs, cs], 2, mean, na.rm = TRUE)
fig5_row(abunds, pdists[,,3], PIT_m3, 3, title = "cAR(1)")



PIT <- matrix(NA, nrow = length(ensdist_wnt), ncol = 10)
for(i in 1:12){
  spot <- (501:512)[i]
  PIT[i,] <- nrPIT(abunds[spot], ecdf(ensdist_wnt[[i]]), n_bins = 10)
}
PIT_en1 <- apply(PIT, 2, mean, na.rm = TRUE)

fig5_row(abunds, listtomat(ensdist_wnt), PIT_en1, 4, title = "Naïve Ensemble")



PIT <- matrix(NA, nrow = ncol(ensdraw_wci), ncol = 10)
for(i in 1:12){
  spot <- (501:512)[i]
  PIT[i,] <- nrPIT(abunds[spot], ecdf(ensdraw_wci[,i]), n_bins = 10)
}
PIT_en2 <- apply(PIT, 2, mean, na.rm = TRUE)

fig5_row(abunds, ensdraw_wci, PIT_en2, 5, title = "Beta-0 RPS Stack")


PIT <- matrix(NA, nrow = ncol(ensdraw_wcb1), ncol = 10)
for(i in 1:12){
  spot <- (501:512)[i]
  PIT[i,] <- nrPIT(abunds[spot], ecdf(ensdraw_wcb1[,i]), n_bins = 10)
}
PIT_en3 <- apply(PIT, 2, mean, na.rm = TRUE)

fig5_row(abunds, ensdraw_wcb1, PIT_en3, 6, title = "Beta-1 RPS Stack")


PIT <- matrix(NA, nrow = ncol(ensdraw_wcb2), ncol = 10)
for(i in 1:12){
  spot <- (501:512)[i]
  PIT[i,] <- nrPIT(abunds[spot], ecdf(ensdraw_wcb2[,i]), n_bins = 10)
}
PIT_en4 <- apply(PIT, 2, mean, na.rm = TRUE)

fig5_row(abunds, ensdraw_wcb2, PIT_en4, 7, title = "Beta-2 RPS Stack")








ensdist_wlb1 <- listtomat(ensdist_wlb1)
PIT <- matrix(NA, nrow = ncol(ensdist_wlb1), ncol = 10)
for(i in 1:12){
  spot <- (501:512)[i]
  PIT[i,] <- nrPIT(abunds[spot], ecdf(ensdist_wlb1[,i]), n_bins = 10)
}
PIT_en <- apply(PIT, 2, mean, na.rm = TRUE)

fig5_row(abunds, ensdist_wlb1, PIT_en, 6, title = "B1LS Ensemble")



ensdist_wlb2 <- listtomat(ensdist_wlb2)
PIT <- matrix(NA, nrow = ncol(ensdist_wlb2), ncol = 10)
for(i in 1:12){
  spot <- (501:512)[i]
  PIT[i,] <- nrPIT(abunds[spot], ecdf(ensdist_wlb2[,i]), n_bins = 10)
}
PIT_en <- apply(PIT, 2, mean, na.rm = TRUE)

fig5_row(abunds, ensdist_wlb2, PIT_en, 7, title = "B2LS Ensemble")



ensdist_wli <- listtomat(ensdist_wli)
PIT <- matrix(NA, nrow = ncol(ensdist_wli), ncol = 10)
for(i in 1:12){
  spot <- (501:512)[i]
  PIT[i,] <- nrPIT(abunds[spot], ecdf(ensdist_wli[,i]), n_bins = 10)
}
PIT_en <- apply(PIT, 2, mean, na.rm = TRUE)

fig5_row(abunds, ensdist_wli, PIT_en, 5, title = "B0LS Ensemble")




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


