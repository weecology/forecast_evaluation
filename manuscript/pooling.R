# i'm not convinced that we're doing this right just yet
# i think what we actually want is a single set of beta parameters
# otherwise how do we get them for the true test?
# gneiting and ranjan use tons of one step aheads, and then give a single
# set of weights and a single set of beta params
#


mods <- list(mod1, mod2, mod3)

smods <- stack_mods(mods, abunds, test_origins = 300:309)

