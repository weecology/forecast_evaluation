source("functions.R")

load_dependencies()
data_set <- prep_data_set()
models <- run_models(data_set)
make_figures(data_set, models)

    






