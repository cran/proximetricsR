## -----------------------------------------------------------------------------
#| label: setup
#| include: false
old_options <- options(digits = 3)


## ----loadlib, results = 'hide', include = FALSE-------------------------------
# Load package, or use devtools::load_all() if in development
if (!requireNamespace("proximetricsR", quietly = TRUE)) {
  devtools::load_all()
}
library("proximetricsR")


## ----loadlib2, eval = FALSE---------------------------------------------------
# library("prospectr")
# library("proximetricsR")


## ----read, results = 'show'---------------------------------------------------
# data location
data_loc <- "https://raw.githubusercontent.com/buchi-labortechnik-ag/demo_data/main/data/"

# location of the tsv file containing the spectral data (of soybean meal)
my_file_1 <- "SoybeanMeal_file1.tsv"
my_file_2 <- "SoybeanMeal_file2.tsv"

my_file_1 <- paste0(data_loc, my_file_1)
my_file_2 <- paste0(data_loc, my_file_2)

# read the files
mdata_1 <- proximate_read_data(my_file_1)
mdata_2 <- proximate_read_data(my_file_2)


# explore the structure of the objects of the imported data sets
str(mdata_1, vec.len = 0, give.attr = FALSE)

str(mdata_2, vec.len = 0, give.attr = FALSE)


## ----spc, results = 'hide'----------------------------------------------------
mdata_1$spc
mdata_2$spc


## ----wavelengths, results = 'hide'--------------------------------------------
original_wavs_1 <- as.numeric(colnames(mdata_1$spc))

original_wavs_2 <- as.numeric(colnames(mdata_2$spc))


## ----plotspectra, fig.cap = "Spectra loaded. Left: datastet 1; Right: Dataset 2",  fig.cap.style = "Image Caption", fig.align = "center", fig.width = 8, fig.height = 8, echo = TRUE, fig.retina = 0.85----
old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), mar = c(4, 4, 1, 4))

y_max <- max(c(max(mdata_1$spc), max(mdata_2$spc)))

matplot(
  x = original_wavs_1, 
  y = t(mdata_1$spc), 
  xlab = "Wavelengths, nm",
  ylab = "Absorbance", 
  xlim = c(400, 1700), 
  ylim = c(0, y_max),
  type = "l", 
  lty = 1, 
  col = rgb(1, 0.5, 0.5, 0.3)
)
grid()

matplot(
  x = original_wavs_2, 
  y = t(mdata_2$spc), 
  xlab = "Wavelengths, nm",
  ylab = "Absorbance", 
  xlim = c(400, 1700), 
  ylim = c(0, y_max),
  type = "l", 
  lty = 1, 
  col = rgb(0.9, 0.8, 0.5, 0.3)
)
grid()


## ----nwpmerge, results = 'hide'-----------------------------------------------
mdata <- proximate_merge(list(mdata_1, mdata_2))


## ----resample-----------------------------------------------------------------
# define the vector of new wavelengths (with constant resolution)
# (replace the psc matrix with the new resampled maxtrix)
# for working with proximetricsR we need to have our spectra in a spc matrix inside 
# the data matrix
mdata$spc <- process(mdata$spc, prep_resample(c(900, 1700, 2)), device = "proximate")


## ----finalwavs,  results = 'hide'---------------------------------------------
final_wavs <- as.numeric(colnames(mdata$spc))


## ----plotspectra2, fig.cap = "Merged spectra.",  fig.cap.style = "Image Caption", fig.align = "center", fig.width = 8, fig.height = 5, echo = TRUE, fig.retina = 0.85----
par(old_par)
matplot(
  x = final_wavs, 
  y = t(mdata$spc), 
  xlab = "Wavelengths, nm",
  ylab = "Absorbance", 
  type = "l", 
  lty = 1, 
  col = rgb(0.5, 0.5, 0.5, 0.3)
)
grid()


## ----settoflase,  results = 'hide'--------------------------------------------
# a vector with the row indices of the samples to be deactivated
deact <- c(5, 10, 50)
# deactivate
mdata$Check[deact] <- "False"


## ----settotrue,  results = 'hide'---------------------------------------------
mdata$Check <- "True"


## ----properties, results = 'hide'---------------------------------------------
# this gets the names of all the variables in the data
names(mdata)

# this returns the column names of all responses
y_names <- extract_property_names(mdata)

# the indices of all the response variables
ys_indices <- which(colnames(mdata) %in% y_names)

# samples with reference values
colSums(!is.na(mdata[, y_names]))


## ----properties3, results = 'hide'--------------------------------------------
#get the names of the response variables
y_names <- y_names[1:2]


## ----formula1, results = 'hide'-----------------------------------------------
f_simple <- Protein ~ spc

# or 

f_long <- as.formula(paste0("Protein ~ ", paste0(final_wavs, collapse = " + ")))


## ----formulas1, results = 'hide'----------------------------------------------
app_formulas <- list(Moisture ~ spc, Protein ~ spc)


## ----formulas2, results = 'hide'----------------------------------------------
app_formulas <- lapply(paste0(y_names, " ~ spc"), FUN = formula)
app_formulas


## ----metadata, results = 'hide'-----------------------------------------------
# Define the metadata of each model (in the same order they are in app_formulas)
models_metadata <- list(
  # for the first property 
  add_model_metadata(decimal_places = 2, unit = "%"),
  # for the second property
  add_model_metadata(decimal_places = 2, unit = "%") 
)


## ----recipe, results = 'show'-------------------------------------------------
# here a recipe with the following pre-processing steps is created:
# - spline/resampling: with spectral range between 900 and 1700 in steps of 4
# - StandardNormalVariate
# - First derivative: with a gap parameter of 5 and a smoothing factor of 9
my_precipe_1 <- preprocess_recipe(
  prep_resample(c(900, 1700, 4)),
  prep_snv(),
  prep_derivative(m = 1, w = 5, p = 9, algorithm = "nwp"), 
  device = "proximate"
)

my_precipe_1$preprocessing_order


## ----recipe2, results = 'show'------------------------------------------------
# here a recipe with the following pre-processing steps is created:
# - spline/resampling: with spectral range between 900 and 1700 in steps of 5
# - StandardNormalVariate
# - First derivative: with a gap parameter of 7 and a smoothing factor of 11
my_precipe_2 <- preprocess_recipe(
  prep_resample(c(900, 1700, 4)), 
  prep_snv(),
  prep_derivative(m = 2, w = 7, p = 11, algorithm = "nwp"),
  device = "proximate"
)

my_precipe_2$preprocessing_order

# here a recipe with the following pre-processing steps is created:
# - spline/resampling: with spectral range between 900 and 1700 in steps of 4
# - StandardNormalVariate
# - Second derivative: with a gap parameter of 9 and a smoothing factor of 13
my_precipe_3 <- preprocess_recipe(
  prep_resample(c(900, 1700, 4)), 
  prep_snv(),
  prep_derivative(m = 2, w = 9, p = 13, algorithm = "nwp"),
  device = "proximate"
)

my_precipe_3$preprocessing_order

# here a recipe with the following pre-processing steps is created:
# - spline/resampling: with spectral range between 900 and 1700 in steps of 8
# - StandardNormalVariate
# - Second derivative: with a gap parameter of 5 and a smoothing factor of 5
my_precipe_4 <- preprocess_recipe(
  prep_resample(c(900, 1700, 5)), 
  prep_snv(),
  prep_derivative(m = 2, w = 5, p = 5, algorithm = "nwp"),
  device = "proximate"
)

my_precipe_4$preprocessing_order


## ----recipelist, results = 'show'---------------------------------------------
my_precipes <- list(my_precipe_1, my_precipe_2, my_precipe_3, my_precipe_4)


## ----rmethod, results = 'show'------------------------------------------------
# here we use the modified pls regression method (equivalent to the one 
# implemented in NIRWise PLUS)
# We use a maximum 11 components
my_pls_method <- fit_plsr(ncomp = 11, type = "nwp")
my_pls_method


## ----rmethodxls, results = 'show'---------------------------------------------
# here we use the modified pls regression method (equivalent to the one 
# implemented in NIRWise PLUS)
# We use a maximum 11 components
my_xls_method <- fit_xlsr(ncomp = 11)
my_xls_method


## ----calcontrol, results = 'show'---------------------------------------------
# Control some aspects of how the calibration is built and optimized
# We use k-fold cross-validation with selection of folds based on the order of 
# the samples in the data table ("sequential")
# We also specify the number of times a model must be re-fitted after outlier 
# removal, e.g. 0 means no re-fiting i.e. no outlier removal; 1 means a model is 
# built, then it is used to identify and remove outliers and finally a the final 
# model is refitted; a value of 5 would mean that the model is refitted 4 times 
# for 4 outlier removal iterations. 
my_control <- calibration_control(
  validation_type = "kfold", 
  number = 4,
  folds = "sequential", 
  remove_outliers = 1 # the number of iterations of outlier removal
)


## ----calibration, results =  'hide'-------------------------------------------
optimized_app <- calibrate_models(
  formulas = app_formulas,
  data = mdata,
  preprocess_recipes = my_precipes,
  methods = list(my_pls_method, my_xls_method),
  control = my_control, 
  metadata_list = models_metadata, 
  save_all = TRUE
)


## ----showoptimize2, eval = FALSE----------------------------------------------
# optimized_app


## ----showoptimize1, eval = TRUE, echo = FALSE---------------------------------
print(optimized_app, separator = " >\n")


## ----finalmodels1, eval =  FALSE----------------------------------------------
# optimized_app$final_models


## ----finalmodels2, eval =  TRUE-----------------------------------------------
names(optimized_app$final_models)


## ----prepro1, echo = FALSE, eval =  TRUE, comment = NA------------------------
print(optimized_app$final_models[[1]]$preprocess, separator = " >\n")


## ----prepro2, echo = FALSE, eval =  TRUE, comment = NA------------------------
print(optimized_app$final_models[[2]]$preprocess, separator = " >\n")


## ----prepro3, eval =  FALSE---------------------------------------------------
# optimized_app$final_models$`Protein ~ spc`
# # and
# optimized_app$final_models$`Moisture ~ spc`


## ----modelplots, eval =  FALSE------------------------------------------------
# plot(optimized_app$final_models$`Protein ~ spc`, selection = "all")
# # and
# plot(optimized_app$final_models$`Moisture ~ spc`, selection = "all")


## ----summaryopt, eval =  FALSE------------------------------------------------
# optimized_app$results_grid


## ----summaryopt2, eval =  TRUE------------------------------------------------
optimized_app$results_grid[optimized_app$results_grid$selection, ]


## ----othermodels, eval =  FALSE-----------------------------------------------
# optimized_app$all_models


## ----othermodels2, eval =  FALSE----------------------------------------------
# my_other_models <- optimized_app$all_models[c(1, 6)]


## ----preds, eval =  TRUE, results = 'hide'------------------------------------
my_pred_file <- "SoybeanMeal_file3.tsv"
my_pred_file <- paste0(data_loc, my_pred_file)

pdata <- proximate_read_data(my_pred_file)
original_wavs_pred <- as.numeric(colnames(pdata$spc))


## ----plotspectra3, fig.cap = "Merged spectra.",  fig.cap.style = "Image Caption", fig.align = "center", fig.width = 8, fig.height = 5, echo = TRUE, fig.retina = 0.85----
matplot(
  x = original_wavs_pred, 
  y = t(pdata$spc), 
  xlab = "Wavelengths, nm",
  ylab = "Absorbance", 
  type = "l", 
  lty = 1, 
  col = rgb(0.3, 0.9, 0.5, 0.5)
)
grid()


## ----preds2, eval =  TRUE, results = 'hide'-----------------------------------
preds <- predict(optimized_app, newdata = pdata)

# matrix of predictions
preds$predictions

# info on the models used
preds$model_information


## ----preds3, fig.cap = "Validation plots for the preditions of the optimal models found",  fig.cap.style = "Image Caption", fig.align = "center", fig.width = 7, fig.height = 3.5, echo = TRUE, fig.retina = 0.85----
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

# Moisture plot
plot(
  x = preds$predictions[, "Moisture"],
  y = pdata$Moisture,
  xlab = "Predicted moisture, %",
  ylab = "Measured moisture, %",
  xlim = c(10.5, 12.5),
  ylim = c(10.5, 12.5),
  pch = 16,
  cex = 1.5,
  col = "grey"
  )
grid()
abline(0, 1, col = "red")

# Protein plot
plot(
  x = preds$predictions[, "Protein"],
  y = pdata$Protein,
  xlab = "Predicted protein, %",
  ylab = "Measured protein, %",
  xlim = c(44, 48),
  ylim = c(44, 48),
  pch = 16,
  cex = 1.5,
  col = "grey"
  )
grid()
abline(0, 1, col = "red")
par(old_par)


## ----preds4, eval =  FALSE, results = 'hide'----------------------------------
# # squared correlations
# cor(x = preds$predictions[, "Moisture"],  y = pdata$Moisture)^2
# cor(x = preds$predictions[, "Protein"], y = pdata$Protein)^2
# 
# # rmses
# mean((preds$predictions[, "Moisture"] - pdata$Moisture)^2)^0.5
# mean((x = preds$predictions[, "Protein"] - pdata$Protein)^2)^0.5


## ----mymodelslist, eval = TRUE, results = 'hide'------------------------------
my_models <- optimized_app$final_models

# add some important metadata to the application/model list
my_models <- add_application_metadata(
  my_models, 
  view = "Up", 
  name = "my_application", 
  description = "created with proximetricsR"
)


## ----appfile, eval = FALSE, results = 'hide'----------------------------------
# proximate_write_nax(
#   path = getwd(),
#   object = my_models
# )


## ----protein, eval = FALSE, results = 'hide'----------------------------------
# my_p_recipe <- preprocess_recipe(
#   prep_resample(c(900, 1700, 2)),
#   prep_snv(),
#   prep_derivative(m = 1, w = 5, p = 3, algorithm = "nwp"),
#   device = "proximate"
# )
# 
# protein_model <- calibrate(
#   Protein ~ spc,
#   data = mdata,
#   preprocess = my_p_recipe, # here we can use only one recipe
#   method = fit_xlsr(5, type = "nwp"),
#   control = my_control
# )
# protein_model
# 
# # add the units and the number of decimal places to the model
# protein_model <- add_model_metadata(
#   protein_model,
#   unit = "%",
#   decimal_places = 2
# )
# 
# plot(protein_model, selection = "all")


## ----predictions, results = 'hide', eval = FALSE------------------------------
# # add the units and the number of decimal places to the model
# my_protein_predictions <- predict(protein_model, newdata = pdata$spc)
# my_protein_predictions$model_information
# my_protein_predictions$predictions


## ----writenax, results = 'hide', eval = FALSE---------------------------------
# # create a list with all the models separately calibrated
# # (in this case we only have one)
# a_model <- list(protein_model)
# 
# a_model <- add_application_metadata(
#   a_model,
#   view = "Up",
#   name = "a_protein_application",
#   description = "created with proximetricsR"
# )
# 
# proximate_write_nax(
#   path = getwd(),
#   object = a_model
# )


## ----replacemodel, results = 'hide', eval = FALSE-----------------------------
# my_models <- optimized_app$final_models
# my_models$`Protein ~ spc` <- protein_model


## ----writetsv, eval = FALSE, results = 'hide'---------------------------------
# # let's write down the calibration data into a tsv file
# # define a name for the tsv to be written
# tsv_name <- "merged_cal_data.tsv"
# 
# # write the tsv
# proximate_write_data(
#   x = mdata,
#   file = tsv_name,
#   id = mdata$ID,
#   spc = "spc",
#   spc_round = 8,
#   barcode = "",
#   properties = y_names,
#   note = "merged_file_1_and_file_2",
#   recipe = "SoybeanMeal up-view",
#   created = mdata$Date,
#   snr = mdata$SNR
# )


## ----prjcal, eval = FALSE, results = 'hide'-----------------------------------
# proximate_write_model(
#   object = list(protein_model),
#   path = getwd(),
#   tsv_paths = tsv_name,
#   application_name = "Untitled",
#   cal = TRUE,
#   prj = TRUE,
#   rtf = TRUE
# )


## ----cleanup, include = FALSE-------------------------------------------------
options(old_options)

