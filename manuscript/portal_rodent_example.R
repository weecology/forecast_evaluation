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




