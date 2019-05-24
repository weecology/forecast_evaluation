library(coda)
library(runjags)
library(portalr)
library(scoringRules)
source("functions.R")

download_observations()

abunds <- plot_sp_abunds(plot = 17, species = "DM")
moon_dates <- load_trapping_data(clean = FALSE)$newmoons_table$newmoondate
in_ts <- training_ts(ends = 100:105)

#mod1 <- mods(model = 1, abunds = abunds, in_timeseries = in_ts)

#mod2 <- mods(model = 2, abunds = abunds, in_timeseries = in_ts)
#mod3 <- mods(model = 3, abunds = abunds, moon_dates = moon_dates, 
#             in_timeseries = in_ts)
#mod4 <- mods(model = 4, abunds = abunds, in_timeseries = in_ts)


mod1 <- vector("list", length = length(in_ts))
for(i in 1:length(in_ts)){
  mod1[[i]] <- tryCatch(
                modi(model = 1, abunds = abunds, 
                     in_timeseries = in_ts[[i]]),
                error = function(x){NA})
  names(mod1)[i] <- run_namer(1, in_ts[[i]], lead_time = 12)
}

save(mod1, file = "mod1_withmodels.RData")
save_without_models(mod1, file = "mod1_withoutmodels.RData")

mod2 <- vector("list", length = length(in_ts))
for(i in 1:length(in_ts)){
  mod2[[i]] <- tryCatch(
                modi(model = 2, abunds = abunds, 
                     in_timeseries = in_ts[[i]]),
                error = function(x){NA})
  names(mod2)[i] <- run_namer(2, in_ts[[i]], lead_time = 12)
}

save(mod2, file = "mod2_withmodels.RData")
save_without_models(mod2, file = "mod2_withoutmodels.RData")

mod3 <- vector("list", length = length(in_ts))
for(i in 1:length(in_ts)){
  mod3[[i]] <- tryCatch(
                modi(model = 3, abunds = abunds, moon_dates = moon_dates,
                     in_timeseries = in_ts[[i]]),
                error = function(x){NA})
  names(mod3)[i] <- run_namer(3, in_ts[[i]], lead_time = 12)
}

save(mod3, file = "mod3_withmodels.RData")
save_without_models(mod3, file = "mod3_withoutmodels.RData")

mod4 <- vector("list", length = length(in_ts))
for(i in 1:length(in_ts)){
  mod4[[i]] <- tryCatch(
                modi(model = 4, abunds = abunds, 
                     in_timeseries = in_ts[[i]]),
                error = function(x){NA})
  names(mod4)[i] <- run_namer(4, in_ts[[i]], lead_time = 12)
}

save(mod4, file = "mod4_withmodels.RData")
save_without_models(mod4, file = "mod4_withoutmodels.RData")
