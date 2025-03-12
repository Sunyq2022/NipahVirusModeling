# Load required packages
library(sf)
library(tidyverse)
library(pROC)
library(ggspatial)
library(patchwork)

# Load preprocessed data ----------------------------------------------------------------
sasea_map <- readRDS('Data/sasea_map.rds')
ten_segment_line <- readRDS('Data/ten_segment_line.rds')
env_vars <- readRDS('Data/env_vars_processed.rds')

# Define evaluation and visualization function -----------------------------------------------
evaluate_model <- function(nipah_type) {
  # Set output directory based on model type
  if (nipah_type == "nipah_human") {
    output_dir <- 'Outputs/human_spillover_model'
  } else if (nipah_type == "nipah_zoonotic") {
    output_dir <- 'Outputs/zoonotic_model'
  }
  
  # Load model results
  training_predictions <- readRDS(paste0(output_dir, '/', nipah_type, "_training_predictions.rds"))
  test_predictions <- readRDS(paste0(output_dir, '/', nipah_type, "_test_predictions.rds"))
  all_predictions <- readRDS(paste0(output_dir, '/', nipah_type, "_all_predictions.rds"))
  variable_contributions <- readRDS(paste0(output_dir, '/', nipah_type, "_variable_contributions.rds"))
  response_curves <- readRDS(paste0(output_dir, '/', nipah_type, "_response_curves.rds"))
  combined_data <- readRDS(paste0(output_dir, '/', nipah_type, "_combined_data.rds"))
  
  # Calculate averages
  avg_variable_contributions <- data.frame(var_name = variable_contributions$var, 
                                           avg_rel_con = apply(variable_contributions[, -1], 1, mean))
  avg_training_predictions <- aggregate(cbind(case, pred) ~ index, data = training_predictions, mean)
  avg_test_predictions <- aggregate(cbind(case, pred) ~ index, data = test_predictions, mean)
  avg_all_predictions <- aggregate(cbind(pred) ~ OID_, data = all_predictions, mean)
  print(avg_variable_contributions)
  
  # Model evaluation with ROC curve
  auc_values_list <- numeric()
  for (i in 1:100) {
    dataset <- test_predictions[test_predictions$round == i, ]
    roc_obj <- roc(dataset$case, dataset$pred, auc = TRUE)
    auc_values_list <- c(auc_values_list, roc_obj$auc)
    if (i == 1) {
      plot(roc_obj, legacy.axes = TRUE, col = "#b2b2b2", lwd = 1, print.auc = FALSE, add = FALSE)
    } else {
      plot(roc_obj, legacy.axes = TRUE, col = "#b2b2b2", lwd = 1, print.auc = FALSE, add = TRUE)
    }
  }
  auc_summary <- data.frame(auc = auc_values_list, mean_auc = mean(auc_values_list), 
                            lower_ci = quantile(auc_values_list, 0.025), 
                            upper_ci = quantile(auc_values_list, 0.975))
  write_rds(auc_summary, paste0(output_dir, '/', nipah_type, "_auc_summary.rds"))
  
  avg_dataset <- avg_test_predictions
  roc_avg <- roc(avg_dataset$case, avg_dataset$pred)
  pdf(paste0(output_dir, '/', nipah_type, '_roc_curve.pdf'))
  plot(roc_avg, legacy.axes = TRUE, col = "red", thresholds = "best", print.thres = "best", 
       print.auc.col = "black", print.auc = TRUE, add = FALSE)
  dev.off()
  
  threshold_metrics <- coords(roc_avg, "best", best.method = "youden", 
                              ret = c("threshold", "specificity", "sensitivity", "accuracy", 
                                      "precision", "recall", "youden", "tn", "tp", "fn", "fp"), 
                              transpose = FALSE)
  write_rds(threshold_metrics, paste0(output_dir, '/', nipah_type, '_threshold_metrics.rds'))
  cutoff_value <- threshold_metrics[1, "threshold"]
  
  # Update prediction map with binary prediction
  prediction_map <- readRDS(paste0(output_dir, '/', nipah_type, '_prediction_map.rds')) %>%
    mutate(binary_pred = ifelse(probability > cutoff_value, 1, 0))
  write_rds(prediction_map, paste0(output_dir, '/', nipah_type, '_prediction_map.rds'))
  
  # Response curves for variables with contribution > 5%
  significant_vars <- avg_variable_contributions %>% filter(avg_rel_con > 5)
  significant_var_indices <- significant_vars$var_name
  n_significant_vars <- length(significant_var_indices)
  
  response_data <- list()
  n_predictors <- (ncol(response_curves) - 2) / 2  # Number of predictors
  for (i in 1:n_predictors) {
    var_names <- colnames(response_curves)[c(i * 2 + 1, i * 2 + 2)]
    response_data[[i]] <- response_curves %>% dplyr::group_by(get(var_names[1])) %>%
      summarise(mean = mean(get(var_names[2])), 
                min = quantile(get(var_names[2]), probs = 0.025), 
                max = quantile(get(var_names[2]), probs = 0.975))
    names(response_data[[i]]) <- c("x", "mean", "min", "max")
    names(response_data)[i] <- var_names[1]
  }
  
  ordered_var_names <- as.character(arrange(avg_variable_contributions, desc(avg_rel_con))[[1]])
  ordered_response_data <- list()
  for (i in 1:n_predictors) {
    var_name <- ordered_var_names[i]
    ordered_response_data[[i]] <- response_data[[var_name]]
  }
  names(ordered_response_data) <- ordered_var_names
  
  significant_response_data <- ordered_response_data[significant_var_indices]
  combined_data_df <- as.data.frame(combined_data)
  secondary_axes <- list()
  for (i in 1:n_significant_vars) {
    hist_max <- max(table(cut(combined_data_df[, names(ordered_response_data)[i]], 
                              breaks = seq(min(combined_data_df[, names(ordered_response_data)[i]]), 
                                           max(combined_data_df[, names(ordered_response_data)[i]]), length.out = 40), 
                              include.lowest = TRUE)))
    min_val <- abs(min(as.matrix(ordered_response_data[[i]][, 2:4])))
    max_val <- abs(max(as.matrix(ordered_response_data[[i]][, 2:4])))
    secondary_axes[[i]] <- c((hist_max / (min_val + max_val) / 100) * 100, min_val)
  }
  
  plot_response_curve <- function(i) {
    x_label <- paste0(names(ordered_response_data)[i], " (", 
                      round(avg_variable_contributions$avg_rel_con[avg_variable_contributions$var_name == names(ordered_response_data)[i]], 2), "%)")
    if (grepl("value", names(ordered_response_data)[i])) {
      curve <- ggplot(ordered_response_data[[i]]) +
        geom_ribbon(aes(x = x, ymin = (min + secondary_axes[[i]][2]) * secondary_axes[[i]][1], 
                        ymax = (max + secondary_axes[[i]][2]) * secondary_axes[[i]][1]), fill = "#3d3939", alpha = 0.5) +
        geom_histogram(data = combined_data_df, aes(x = get(names(ordered_response_data)[i])), 
                       binwidth = (max(combined_data_df[, names(ordered_response_data)[i]]) - min(combined_data_df[, names(ordered_response_data)[i]])) / 40, 
                       fill = "#423b3b", color = "#423b3b", size = 0.1, alpha = 0.5) +
        geom_line(aes(x = x, y = (mean + secondary_axes[[i]][2]) * secondary_axes[[i]][1]), color = "red", size = 0.35) +
        scale_x_continuous(name = x_label, labels = scales::percent) +
        scale_y_continuous(name = "Frequency", 
                           sec.axis = sec_axis(~ . / secondary_axes[[i]][1] - secondary_axes[[i]][2], name = "Fitted function")) +
        theme(axis.title = element_text(size = 7), axis.text.x = element_text(size = 5), 
              axis.text.y = element_text(size = 5))
      return(curve)
    } else {
      curve <- ggplot(ordered_response_data[[i]]) +
        geom_ribbon(aes(x = x, ymin = (min + secondary_axes[[i]][2]) * secondary_axes[[i]][1], 
                        ymax = (max + secondary_axes[[i]][2]) * secondary_axes[[i]][1]), fill = "#3d3939", alpha = 0.5) +
        geom_histogram(data = combined_data_df, aes(x = get(names(ordered_response_data)[i])), 
                       binwidth = (max(combined_data_df[, names(ordered_response_data)[i]]) - min(combined_data_df[, names(ordered_response_data)[i]])) / 40, 
                       fill = "#423b3b", color = "#423b3b", size = 0.1, alpha = 0.5) +
        geom_line(aes(x = x, y = (mean + secondary_axes[[i]][2]) * secondary_axes[[i]][1]), color = "red", size = 0.35) +
        scale_x_continuous(name = x_label, labels = scales::comma) +
        scale_y_continuous(name = "Frequency", 
                           sec.axis = sec_axis(~ . / secondary_axes[[i]][1] - secondary_axes[[i]][2], name = "Fitted function")) +
        theme(axis.title = element_text(size = 7), axis.text.x = element_text(size = 5), 
              axis.text.y = element_text(size = 5))
      return(curve)
    }
  }
  
  eval(parse(text = paste(paste0('plot_response_curve(', 1:n_significant_vars, ')'), collapse = '+'), '+ plot_layout(ncol=3)'))
  ggsave(paste0(output_dir, '/', nipah_type, "_response_curves.pdf"), device = "pdf", width = 10, height = 10)
  write_rds(ordered_response_data, paste0(output_dir, '/', nipah_type, '_response_data.rds'))
  
  # Uncertainty map
  prediction_uncertainty <- aggregate(cbind(pred) ~ OID_, data = all_predictions, sd)
  names(prediction_uncertainty) <- c('OID_', 'sd')
  env_vars %>% left_join(prediction_uncertainty, by = "OID_") %>% 
    ggplot() +
    geom_sf(color = "transparent", aes(fill = sd)) + 
    scale_fill_gradient2(low = "white", mid = "#fc8d59", high = "#b30000", 
                         midpoint = 0.5, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         na.value = "#368BB1", name = str_wrap("Uncertainty", width = 24)) +
    geom_sf(data = sasea_map, fill = 'transparent', colour = 'black', size = 0.5) +
    geom_sf(data = ten_segment_line) +
    theme_bw() + 
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white"), 
          axis.ticks = element_blank(), axis.title = element_blank(), 
          legend.title = element_blank(), legend.text = element_text(size = 10), 
          legend.background = element_blank(), legend.position = c(0.2, 0.1), 
          legend.direction = "horizontal", legend.key.height = unit(1, 'cm'), 
          legend.key.width = unit(1.5, 'cm'), legend.key = element_blank(), 
          legend.box.margin = margin(0, 0, 0, 0), legend.margin = margin(0, 0, 0, 0)) +
    coord_sf(xlim = c(60, 155), ylim = c(-45, 55)) +
    annotation_scale(location = "bl", width_hint = 0.4, text_cex = 1) + 
    annotation_north_arrow(location = "tr", which_north = "true", 
                           pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"), 
                           style = north_arrow_fancy_orienteering)
  ggsave(paste0(output_dir, '/', nipah_type, '_uncertainty_map.pdf'), device = "pdf", width = 10, height = 10, dpi = 300)
}

# Evaluate models
evaluate_model('nipah_human')
evaluate_model('nipah_zoonotic')