library(furrr)
library(tidymodels)

## Parallel Version
plan(multisession)
# plan(list(
#   tweak(multisession, workers = availableCores() %/% 4),
#   tweak(multisession, workers = 4)
# ))


### Function Definitions

# `object` will be an `rsplit` object from our `data_folds` tibble
# `penalty` is the tuning parameter
# `desired_outcome` is the string version of the desired value/class label (e.g., "good", "diseased", etc)
lasso_roc <- function(object, penalty = 0.1, desired_outcome) {
  
  #y_col <- ncol(object$data)
  
  lasso_spec <- logistic_reg(penalty = penalty, mixture = 1) %>% set_engine("glmnet")
  
  data_recipe <- recipe(y ~ ., data = analysis(object)) %>% 
    step_zv(all_numeric(), -all_outcomes()) %>% 
    step_normalize(all_numeric(), -all_outcomes()) 
  
  lasso_wflow <- workflow() %>% 
    add_model(lasso_spec) %>% 
    add_recipe(data_recipe)
  
  model <- lasso_wflow %>% fit(data = analysis(object))
  
  holdout_pred <- 
    predict(model, assessment(object) %>% dplyr::select(-y), type = "prob") %>% 
    bind_cols(assessment(object) %>% dplyr::select(y))
  
  roc_auc(holdout_pred, truth = y, any_of(paste0(".pred_", desired_outcome)), event_level = "second")$.estimate
  #roc_auc(holdout_pred, truth = y, .pred_20230824, event_level = "second")$.estimate
}



# In some cases, we want to parameterize the function over the tuning parameter:
roc_wrapper <- function(penalty, object, desired_outcome) lasso_roc(object, penalty, desired_outcome)


# `object` will be an `rsplit` object for the bootstrap samples
tune_over_penalty <- function(object, desired_outcome) {
  tibble(penalty = 10^seq(0, -2, length.out = 10)) %>% 
    mutate(roc = map_dbl(penalty, roc_wrapper, object = object, desired_outcome = desired_outcome))
}

# `object` is an `rsplit` object in `results$inner_resamples` 
summarize_tune_results <- function(object, desired_outcome) {
  
  # Return row-bound tibble that has the inner CV results
  
  # For each value of the tuning parameter, compute the 
  # average ROC which is the inner bootstrap estimate. 
  
  map_df(object$splits, tune_over_penalty, desired_outcome) %>%
    group_by(penalty) %>%
    summarize(mean_roc = mean(roc, na.rm = TRUE), 
              n = length(roc), 
              .groups = "drop")
}

run_lasso_nestedcv <- function(nestedcv_folds, desired_outcome) {
  
  # Step 1) Execute all the inner resampling loops:
  
  # The object tuning_results is a list of data frames for each of the 50 outer resamples.
  tuning_results <- future_map(nestedcv_folds$inner_resamples, 
                               summarize_tune_results, 
                               desired_outcome = desired_outcome, 
                               .options = furrr_options(seed = TRUE))
  
  ## Non-parallel Version
  #tuning_results <- map(nestedcv_folds$inner_resamples, summarize_tune_results)
  
  
  # Step 2) Collect the best penalties identified from the inner folds/resamples
  best_penalty <- function(results) results[which.max(results$mean_roc),]
  
  penalty_vals <- tuning_results %>% map_df(best_penalty) %>% select(penalty)
  
  folds_with_penalty <- bind_cols(nestedcv_folds, penalty_vals)
  
  # Step 3) Now that we have these estimates, we can compute the outer resampling results 
  # for each of the outer splits using the corresponding tuning parameter value:
  
  outer_results <- folds_with_penalty %>% 
    mutate(roc_auc = map2_dbl(splits, penalty, lasso_roc, desired_outcome = desired_outcome))
  
  return(outer_results)
}

## Helpful website: https://dionysus.psych.wisc.edu/iaml/unit_09.html#model-comparisons-feature-ablation
run_single_permutation <- function(data, num_folds = 5, num_repeats = 5, desired_outcome) {
  
  permuted_data <- data %>% mutate(y = sample(y))
  permuted_folds <- nested_cv(permuted_data, 
                              outside = vfold_cv(v = num_folds, repeats = num_repeats, strata = y), 
                              inside = vfold_cv(v = num_folds, strata = y))
  
  permuted_results <- run_lasso_nestedcv(permuted_folds, desired_outcome)
  
  return(mean(permuted_results$roc_auc))
}
