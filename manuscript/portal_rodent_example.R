library(coda)
library(runjags)
library(portalr)
library(scoringRules)
source("functions.R")


###
# run this code if the portal observations need to be downloaded
download_observations()
###


###
# run this code if the large batch of rolling origin forecasts need to be run
# note that the models take a while to run, most especially model 4
# also, the returned model objects are quite large and so the output objects
#   that include the model components are not stored on github (they're in
#   the gitignore file)
abunds <- plot_sp_abunds(plot = 17, species = "DM")
moon_dates <- load_trapping_data(clean = FALSE)$newmoons_table$newmoondate
in_ts <- training_ts(ends = 10:500)

mod1 <- mods(model = 1, abunds = abunds, in_timeseries = in_ts)
mod2 <- mods(model = 2, abunds = abunds, in_timeseries = in_ts)
mod3 <- mods(model = 3, abunds = abunds, moon_dates = moon_dates, 
             in_timeseries = in_ts)
mod4 <- mods(model = 4, abunds = abunds, in_timeseries = in_ts)

save(mod1, file = "model_output/mod1_withmodels.RData")
save_without_models(mod1, file = "model_output/mod1_withoutmodels.RData")
save(mod2, file = "mode2_output/mod2_withmodels.RData")
save_without_models(mod2, file = "model_output/mod2_withoutmodels.RData")
save(mod3, file = "mode3_output/mod3_withmodels.RData")
save_without_models(mod3, file = "model_output/mod3_withoutmodels.RData")
save(mod4, file = "mode4_output/mod4_withmodels.RData")
save_without_models(mod4, file = "model_output/mod4_withoutmodels.RData")
###


###
# run this code to process the large batch of rolling origin forecasts
load("model_output/mod1_withoutmodels.RData")
load("model_output/mod2_withoutmodels.RData")
load("model_output/mod3_withoutmodels.RData")
load("model_output/mod4_withoutmodels.RData")
eval_tab <- make_eval_tab(mod1, mod2, mod3, mod4)


plot(eval_tab$origin, eval_tab$logs, col = eval_tab$model)