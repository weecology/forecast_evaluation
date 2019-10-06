library(portalr)
library(coda)
library(runjags)
library(scoringRules)
source("functions.R")

abunds <- plot_sp_abunds(plot = 19, species = "PP")
moon_dates <- load_trapping_data(clean = FALSE)$newmoons_table$newmoondate

# if models have not yet been fit,
    in_ts <- training_ts(starts = 200, ends = 300:500)
    mod1 <- mods(1, abunds, moon_dates, in_ts)
    mod2 <- mods(2, abunds, moon_dates, in_ts)
    mod3 <- mods(3, abunds, moon_dates, in_ts)
# else, 
    load_models(1:3)

fig1()
fig2()
fig3()
fig4(abunds, moon_dates)
fig5(mod1, mod2, mod3)




