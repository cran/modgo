## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8,
                      fig.height = 8,
                      out.height = "80%",
                      out.width = "80%",
                      dpi = 300)
counter <- 0

## ----setup, warning=FALSE, include = FALSE------------------------------------
library(survival)

## ----data---------------------------------------------------------------------
library(modgo)
data("Cleveland", package = "modgo")

## ----basic_arguments----------------------------------------------------------
# Specifying dichotomous and ordinal categorical variables
binary_variables <- c("Sex","HighFastBloodSugar","CAD","ExInducedAngina")
categorical_variables <- c("Chestpaintype","RestingECG")
nrep <- 500
plot_variables <- c("Age", "STDepression", binary_variables[c(1,3)], categorical_variables)

## ----default_test-------------------------------------------------------------
test <- modgo(data = Cleveland,
              bin_variables = binary_variables,
              categ_variables = categorical_variables,
              nrep = nrep)

## ----correlation_default,echo=FALSE, fig.cap=paste0("Figure ",counter,": Correlation plots for a default *modgo* run.")----
corr_plots(test, variables = plot_variables)
counter <- counter  + 1

## ----distr_default,echo=FALSE, fig.cap=paste0("Figure ",counter,": Distribution plots for a default *modgo* run.")----
distr_plots(test, variables = plot_variables)
counter <- counter  + 1

## -----------------------------------------------------------------------------
Variables <- c("Age")
thresh_left <- c(65)
thresh_right <- c(NA)
thresholds <- data.frame(Variables, thresh_left, thresh_right)

print(as.matrix(thresholds))

test_thresh <-  modgo(data = Cleveland,
                      bin_variables = binary_variables,
                      categ_variables = categorical_variables,
                      thresh_var = thresholds,
                      nrep = nrep,
                      thresh_force = TRUE)

## ----echo=FALSE, fig.cap=paste0("Figure ",counter,": Correlation plot for Age > 65 threshold *modgo* run")----
corr_plots(test_thresh, variables = plot_variables)
counter <- counter  + 1

## ----echo=FALSE, fig.cap=paste0("Figure ",counter,": Distribution plot for Age > 65 threshold *modgo* run")----
distr_plots(test_thresh, variables = plot_variables)
counter <- counter  + 1

## -----------------------------------------------------------------------------
#Create named vector
perturb_vector <- c(0.9,0.7)
names(perturb_vector) <- c("RestingBP","Cholsterol")

test_pertru <-  modgo(data = Cleveland,
                      bin_variables = binary_variables,
                      categ_variables = categorical_variables,
                      pertr_vec = perturb_vector,
                      nrep = nrep)

## ----echo=FALSE, fig.cap=paste0("Figure ",counter,": Correlation plot for Pertrubation Expansion *modgo* run")----
corr_plots(test_pertru, variables = c(plot_variables, names(perturb_vector)))
counter <- counter  + 1

## ----echo=FALSE, fig.cap=paste0("Figure ",counter,": Distribution plot for Pertrubation Expansion *modgo* run")----
distr_plots(test_pertru, variables = c(plot_variables, names(perturb_vector)))
counter <- counter  + 1

## ----GLD_run------------------------------------------------------------------
test_GLD <- modgo(data = Cleveland,
                  bin_variables = binary_variables,
                  categ_variables = categorical_variables,
                  generalized_mode = TRUE,
                  nrep = nrep)

## ----correlation_GLD, echo=FALSE, fig.cap=paste0("Figure ",counter,": Correlation plots for Generalized Lambda Distribtion *modgo* run")----
corr_plots(test_GLD, variables = plot_variables)
counter <- counter  + 1

## ----distr_GLD, echo=FALSE, fig.cap=paste0("Figure ",counter,": Distribution plot for Generalized Lambda Distribtion *modgo* run")----
distr_plots(test_GLD, variables = plot_variables)
counter <- counter  + 1

## ----arguments_GLD_def_model--------------------------------------------------
Variables <- c("Age","STDepression")
Model <- c("rprs", "star-rmfmkl")
model_matrix <- cbind(Variables,
                      Model)

## ----GLD_run_def_model--------------------------------------------------------
test_GLD_define_model <- modgo(data = Cleveland,
                  bin_variables = binary_variables,
                  categ_variables = categorical_variables,
                  generalized_mode = TRUE,
                  generalized_mode_model = model_matrix,
                  nrep = nrep)

## ----correlation_GLD_def_model, echo=FALSE, fig.cap=paste0("Figure ",counter,": Correlation plots for Generalized Lambda Distribtion *modgo* run with specified GLD models")----
corr_plots(test_GLD_define_model, variables = plot_variables)
counter <- counter  + 1

## ----distr_GLD_def_model,echo=FALSE, fig.cap=paste0("Figure ",counter,": Distribution plots for Generalized Lambda Distribtion *modgo* run with specified GLD models")----
distr_plots(test_GLD_define_model, variables = plot_variables)
counter <- counter  + 1

## ----GLD_run_def_model_set_lambdas--------------------------------------------
gener_lambdas_matrix <- generalizedMatrix(data = Cleveland,
                                          generalized_mode_model = model_matrix,
                                          bin_variables = binary_variables) 
test_GLD_define_model_set_lambdas <- modgo(data = Cleveland,
                  bin_variables = binary_variables,
                  categ_variables = categorical_variables,
                  generalized_mode = TRUE,
                  generalized_mode_lmbds = gener_lambdas_matrix,
                  nrep = nrep)

## ----GLD_run_no_data_set------------------------------------------------------
# Necessary arguments
gener_lambdas_matrix <- generalizedMatrix(data = Cleveland,
                                          generalized_mode_model = model_matrix,
                                          bin_variables = binary_variables)
sigma <- cor(Cleveland)
variables_names <- colnames(sigma)
sample_size <- 100

test_GLD_no_data_set <- modgo(data = NULL,
                              variables = variables_names,
                              bin_variables = binary_variables,
                              categ_variables = categorical_variables,
                              sigma = sigma,
                              generalized_mode = TRUE,
                              generalized_mode_lmbds = gener_lambdas_matrix,
                              n_samples = sample_size,
                              nrep = nrep)



## ----correlation_GLD_run_no_data_set, echo=FALSE, fig.cap=paste0("Figure ",counter,": Correlation plots for Generalized Lambda Distribtion *modgo* run without providing a data set")----
corr_plots(test_GLD_no_data_set, variables = plot_variables)
counter <- counter  + 1

## ----survival_data_set--------------------------------------------------------
# cancer prepare
data("cancer", package = "survival")

cancer <- na.omit(cancer)
cancer$sex <- cancer$sex - 1
cancer$status <- cancer$status - 1

time_var_cancer <- "time"
status_var_cancer <- "status"
bin_var_cancer <- c("status", "sex")
cat_var_list_cancer <- c("ph.ecog")

plot_variables_surv <- colnames(cancer)[1:6]


## ----survival_data_set_run----------------------------------------------------
# Survival run
test_surv <- modgo_survival(data = cancer,
               surv_method = 1,
               bin_variables = bin_var_cancer,
               categ_variables = cat_var_list_cancer,
               event_variable = status_var_cancer,
               time_variable = time_var_cancer,
               generalized_mode_model_no_event = "rmfmkl",
               generalized_mode_model_event = "rprs")


## ----survival_data_set_corr, echo=FALSE, fig.cap=paste0("Figure ",counter,": Correlation plots for modgo_survival run")----
corr_plots(test_surv, variables = plot_variables_surv)
counter <- counter  + 1

## ----survival_data_set_distr, echo=FALSE, fig.cap=paste0("Figure ",counter,": Distribution plots for modgo_survival run")----
distr_plots(test_surv, variables = plot_variables_surv)
counter <- counter  + 1

## ----survival_data_set_coxplots, fig.cap=paste0("Figure ",counter,": Survival fit curves plot for modgo_survival run")----
data_set_info <- c(rep("Original", dim(test_surv$original_data)[1]),
                   rep("Simulated", dim(test_surv$simulated_data[[1]])[1]))
combine_data_set <- rbind(test_surv$original_data,
                          test_surv$simulated_data[[1]])
combine_data_set <- cbind(combine_data_set,
                          data_set_info)
fit <- survfit(Surv(time, status) ~ data_set_info,
               data=combine_data_set)
plot(fit,
     fun = "F",
     col=1:2)

legend(700, 1,
       c("Original data set", "Simulated data set"),
       lty=c(1,1),
       col=c(1,2),
       bty='n',
       lwd=2)

