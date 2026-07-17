## -----------------------------------------------------------------------------
#| label: setup
#| include: false

# Disable ANSI colours for vignette rendering
options(cli.num_colors = 1)
Sys.setenv("RSTUDIO" = "")
Sys.setenv("POSITRON" = "")
old_options <- options(digits = 3)


## ----loadlib, results = 'hide', include = FALSE-------------------------------
if (!requireNamespace("proximetricsR", quietly = TRUE)) {
  pkgload::load_all()
}
library("proximetricsR")


## ----loadlib2, eval = FALSE---------------------------------------------------
# library("proximetricsR")


## ----data---------------------------------------------------------------------
data("NIRsoil", package = "prospectr")

# NIRsoil is in nanometers and absorbance
mdata <- data.frame(
  sampleName = rownames(NIRsoil), 
  Ciso = NIRsoil$Ciso, 
  Nt = NIRsoil$Nt,
  CEC = NIRsoil$CEC
)

# ProxiScout data comes in reflectance in percentages, so we convert to 
# reflectance and scale to 0-100 range
mdata$spc <- 100 * 1 / 10^(NIRsoil$spc)

# Get wavelengths and convert from nm to wavenumbers (cm^-1)
wav_nm <- as.numeric(colnames(mdata$spc))
wav_cm <- 10000000 / wav_nm

# Update column names to wavenumbers for ProxiScout compatibility
colnames(mdata$spc) <- wav_cm

head(colnames(mdata$spc))


## -----------------------------------------------------------------------------
class(mdata) <- c("proxiscout_data", class(mdata))


## ----recipe-------------------------------------------------------------------
recipe_01 <- preprocess_recipe(
  prep_resample(grid = "proxiscout"), # necessary for almost all ProxiScout recipe
  prep_derivative(m = 2, w = 9, p = 2, algorithm = "savitzky-golay"),
  prep_snv(),
  device = "proxiscout"
)

recipe_01


## ----calibrate, results = 'hide'----------------------------------------------
model_c <- calibrate(
  Ciso ~ spc,
  data = mdata,
  preprocess = recipe_01,
  method = fit_plsr(ncomp = 12, type = "modified"),
  control = calibration_control(
    validation_type = "kfold",
    number = 5,
    seed = 42
  ),
  verbose = FALSE
)


## ----evaluate, results = 'hide'-----------------------------------------------
model_c


## -----------------------------------------------------------------------------
#| eval: false
# plot(model_c)


## -----------------------------------------------------------------------------
#| label: serialized
#| results: 'hide'
#| eval: false
# my_serialised_model_c <- proxiscout_write_model(model_c, file = NULL)
# my_serialised_model_c


## -----------------------------------------------------------------------------
#| eval: false
# proxiscout_write_model(model_c, file = "my_model_c.json")


## -----------------------------------------------------------------------------
my_formulas <- list(Ciso ~ spc, CEC ~ spc)


## -----------------------------------------------------------------------------
#| label: multi_calibrate1
#| results: 'hide'
recipe_02 <- preprocess_recipe(
  prep_resample(grid = "proxiscout"),
  prep_derivative(m = 1, w = 7, p = 2, algorithm = "savitzky-golay"),
  device = "proxiscout"
)

recipe_03 <- preprocess_recipe(
  prep_resample(grid = "proxiscout"),
  prep_wav_trim(
    band = c(4000, 7000),
    trim_constant_edges = TRUE
  ),
  prep_derivative(m = 1, w = 7, p = 2, algorithm = "savitzky-golay"),
  device = "proxiscout"
)

my_recipes <- list(recipe_01, recipe_02, recipe_03)


## -----------------------------------------------------------------------------
my_fitting <- list(
  fit_plsr(ncomp = 12, type = "modified"), 
  fit_plsr(ncomp = 10, type = "standard")
)


## -----------------------------------------------------------------------------
my_control <- calibration_control(
    validation_type = "kfold",
    number = 5,
    remove_outliers = 1,
    seed = 42
)


## -----------------------------------------------------------------------------
#| label: multi_calibrate2
#| results: 'hide'
multiple_models <- calibrate_models(
  formulas = my_formulas,
  data = mdata, 
  preprocess_recipes = my_recipes,
  methods = my_fitting,
  control = my_control,
  save_all = TRUE
)


## -----------------------------------------------------------------------------
#| eval: false
# multiple_models


## -----------------------------------------------------------------------------
#| eval: true
#| echo: false
df <- multiple_models$results_grid
num_cols <- sapply(df, is.numeric)
df[num_cols] <- lapply(df[num_cols], round, 2)
multiple_models$results_grid <- df
print(multiple_models)


## ----serialized2--------------------------------------------------------------
#| eval: false
# proxiscout_write_model(
#   multiple_models$final_models$`Ciso ~ spc`, file = "my_model_c2.json"
# )
# 
# proxiscout_write_model(
#   multiple_models$final_models$`CEC ~ spc`, file = "my_model_cec.json"
# )


## ----cleanup, include = FALSE-------------------------------------------------
options(old_options)

