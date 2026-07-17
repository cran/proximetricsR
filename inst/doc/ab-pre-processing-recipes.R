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


## ----data, results = 'hide'---------------------------------------------------
data("NIRcannabis")
X <- NIRcannabis$spc


## ----prep_resample_proximate--------------------------------------------------
prep_resample(grid = c(1001, 1700, 2))


## ----prep_resample_proxiscout-------------------------------------------------
prep_resample(grid = "proxiscout")


## ----prep_smooth_sg-----------------------------------------------------------
prep_smooth(w = 11, p = 3, algorithm = "savitzky-golay")


## ----prep_smooth_ma-----------------------------------------------------------
prep_smooth(w = 7, algorithm = "moving-average")


## ----prep_snv-----------------------------------------------------------------
prep_snv()


## ----prep_derivative_sg-------------------------------------------------------
prep_derivative(m = 1, w = 11, p = 3, algorithm = "savitzky-golay")


## ----prep_derivative_gap------------------------------------------------------
prep_derivative(m = 2, w = 9, p = 3, algorithm = "gap-segment")


## ----prep_derivative_nwp------------------------------------------------------
prep_derivative(m = 1, w = 5, p = 11, algorithm = "nwp")


## ----prep_detrend-------------------------------------------------------------
prep_detrend(p = 2)


## ----prep_transform-----------------------------------------------------------
prep_transform(to = "absorbance")


## ----prep_wav_trim------------------------------------------------------------
prep_wav_trim(band = c(1000, 1800), trim_constant_edges = TRUE)


## ----recipe_snv---------------------------------------------------------------
recipe_snv <- preprocess_recipe(prep_snv())
recipe_snv


## ----recipe_proxiscout--------------------------------------------------------
recipe_ps <- preprocess_recipe(
  prep_smooth(w = 7, p = 1, algorithm = "savitzky-golay"),
  prep_snv(),
  prep_derivative(m = 1, w = 5, p = 2, algorithm = "savitzky-golay"),
  device = "proxiscout"
)
recipe_ps


## ----recipe_proximate---------------------------------------------------------
recipe_pm <- preprocess_recipe(
  prep_smooth(w = 7, algorithm = "moving-average"),
  prep_snv(),
  prep_derivative(m = 1, w = 5, p = 11, algorithm = "nwp"),
  device = "proximate"
)
recipe_pm


## ----process_snv--------------------------------------------------------------
X_snv <- process(X, recipe_snv)
dim(X_snv)


## ----process_proxiscout-------------------------------------------------------
X_ps <- process(X, recipe_ps)
dim(X_ps)


## ----process_inspect----------------------------------------------------------
applied_recipe <- attr(X_ps, "preprocess_recipe")
applied_recipe


## ----example_proximate--------------------------------------------------------
recipe_pm_fat <- preprocess_recipe(
  prep_smooth(w = 7, algorithm = "moving-average"),
  prep_snv(),
  prep_derivative(m = 1, w = 5, p = 11, algorithm = "nwp"),
  device = "proximate"
)

X_fat_prep <- process(X, recipe_pm_fat)
head(X_fat_prep[, 1:5])


## ----example_proxiscout_full--------------------------------------------------
recipe_ps_full <- preprocess_recipe(
  prep_resample(grid = "proxiscout"),
  prep_smooth(w = 7, p = 1, algorithm = "savitzky-golay"),
  prep_snv(),
  prep_detrend(p = 2),
  prep_derivative(m = 1, w = 5, p = 2, algorithm = "savitzky-golay"),
  device = "proxiscout"
)

X_ps_full <- process(X, recipe_ps_full)
dim(X_ps_full)


## ----example_minimal----------------------------------------------------------
recipe_minimal <- preprocess_recipe(prep_snv())
X_minimal <- process(X, recipe_minimal)


## ----example_band-------------------------------------------------------------
recipe_band <- preprocess_recipe(
  prep_wav_trim(band = c(1100, 1600)),
  prep_smooth(w = 5, p = 1, algorithm = "savitzky-golay"),
  prep_snv(),
  device = "proxiscout"
)

X_band <- process(X, recipe_band)
colnames(X_band)[c(1, ncol(X_band))]


## ----cleanup, include = FALSE-------------------------------------------------
options(old_options)

