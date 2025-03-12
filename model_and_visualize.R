setwd('D:\\WPS\\699051479\\WPS云盘\\01历史数据备份\\南京疾控笔记本备份\\00011  技术路线\\01 尼帕病毒\\18 尼帕回复\\数据代码共享')

# Ensure required packages are installed -------------------------------------------------------
required_packages <- c("sf", "tidyverse", "dismo", "terra", "raster", "corrplot", 
                       "gbm", "pROC", "usdm", "doParallel", "ggspatial", "patchwork")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)
lapply(required_packages, library, character.only = TRUE)

# Set up parallel computing ------------------------------------------------------------------
no_cores <- detectCores() - 1
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)

# Load preprocessed data ----------------------------------------------------------------
sasea_map <- readRDS('Data/sasea_map.rds')
ten_segment_line <- readRDS('Data/ten_segment_line.rds')
nipah_occurrences <- readRDS('Data/nipah_occurrences.rds')
env_vars <- readRDS('Data/env_vars_processed.rds')
bat_pathogens_sasea <- readRDS('Data/bat_pathogens_sasea.rds')
chikv_dengue_sasea <- readRDS('Data/chikv_dengue_sasea.rds')
bat_deng_chikv <- readRDS('Data/bat_deng_chikv.rds')

# Define modeling function --------------------------------------------------------------------
run_model <- function(nipah_type) {
  # Set output directory based on model type
  if (nipah_type == "nipah_human") {
    output_dir <- 'Outputs/human_spillover_model'
  } else if (nipah_type == "nipah_zoonotic") {
    output_dir <- 'Outputs/zoonotic_model'
  }
  dir.create(output_dir, showWarnings = FALSE)
  
  # Load positive locations for this type
  positive_points <- st_join(nipah_occurrences[[nipah_type]], env_vars, left = TRUE) %>% na.omit()
  
  # Pseudo-negative point selection with 20000 buffer
  buffer_zones <- circles(as(st_geometry(nipah_occurrences[[nipah_type]]), 'Spatial'), d = 20000, lonlat = TRUE) %>% 
    polygons() %>% 
    st_as_sf() %>% 
    st_join(env_vars, left = TRUE)
  points_outside_buffer <- env_vars %>% filter(!OID_ %in% buffer_zones$OID_)
  
  # Background point selection based on nipah_type
  if (nipah_type == "nipah_human") {
    background_ids <- dplyr::intersect(chikv_dengue_sasea$OID_, points_outside_buffer$OID_) %>% unique() %>% as.character()
  } else if (nipah_type == "nipah_zoonotic") {
    background_ids <- dplyr::intersect(bat_deng_chikv$OID_, points_outside_buffer$OID_) %>% unique() %>% as.character()
  }
  
  background_points <- points_outside_buffer %>% 
    filter(OID_ %in% background_ids) %>% 
    filter(!OID_ %in% positive_points$OID_) %>% 
    unique()
  colnames(background_points)[52] <- 'geometry'
  st_geometry(background_points) <- 'geometry'
  
  # Sample background points (10 times positive points)
  set.seed(66948442)
  sampled_background <- background_points[sample(1:nrow(background_points), 10 * nrow(positive_points)), ]
  sampled_background$exist <- 0
  positive_points$exist <- 1
  combined_data <- rbind(sampled_background, positive_points)
  combined_centroids <- st_centroid(combined_data)
  
  # Visualize positive and background points
  ggplot() +
    geom_sf(data = sasea_map, fill = "white", size = 0.2) + 
    geom_sf(data = combined_centroids, aes(color = factor(exist)), size = 0.05) + 
    scale_color_manual(values = c("blue", "#ba1414"), name = "Occurrence", 
                       labels = c("Background point", "Occurrence point")) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white"), 
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
          legend.title = element_text(size = 8), legend.text = element_text(size = 7), 
          legend.background = element_blank(), legend.position = c(0.1, 0.18), 
          legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.5, 'cm'), 
          legend.key = element_blank(), legend.box.margin = margin(0, 0, 0, 0), 
          legend.margin = margin(0, 0, 0, 0))
  ggsave(paste0(output_dir, '/', nipah_type, "_occurrence_background.pdf"), device = "pdf", width = 10, height = 10)
  
  # Prepare data for modeling
  combined_data <- combined_data %>% mutate(index = 1:length(OID_))
  write_rds(combined_data, paste0(output_dir, '/', nipah_type, "_combined_data.rds"))
  prediction_data <- env_vars %>% mutate(index = 1:length(OID_), exist = 0)
  
  # Model preparation
  variable_contributions <- training_predictions <- test_predictions <- all_predictions <- NULL
  response_curves <- NULL
  predictor_vars <- colnames(env_vars)[2:51]  # Predictors as per original code
  predictor_indices <- which(colnames(combined_data) %in% predictor_vars)
  response_index <- which(colnames(combined_data) %in% "exist")
  n_predictors <- length(predictor_indices)
  gbm_models <- list()
  auc_values <- matrix(nrow = 100, ncol = 1, dimnames = list(NULL, "CV_AUC"))
  
  # Model training with 100 iterations
  for (i in 1:100) {
    print(i)
    set.seed(123456 + i)
    test_indices <- sample(1:nrow(combined_data), nrow(combined_data) / 5)
    test_data <- as.data.frame(combined_data[test_indices, ])
    training_data <- as.data.frame(combined_data[-test_indices, ])
    
    gbm_models[[i]] <- gbm.step(data = training_data, gbm.x = predictor_indices, gbm.y = response_index, 
                                family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, 
                                bag.fraction = 0.75, n.trees = 100, step.size = 50, n.folds = 10, 
                                max.trees = 5000, silent = FALSE, plot.main = TRUE)
    auc_values[i, "CV_AUC"] <- gbm_models[[i]]$cv.statistics$discrimination.mean
    
    tmp <- cbind(rep(i, nrow(training_data)), training_data[, c("index", 'OID_', "exist")], 
                 predict.gbm(gbm_models[[i]], training_data, n.trees = gbm_models[[i]]$n.trees, type = "response")) %>% 
      as.data.frame()
    colnames(tmp) <- c('round', "index", 'OID_', 'case', 'pred')
    if (is.null(training_predictions)) training_predictions <- tmp else training_predictions <- rbind(training_predictions, tmp)
    
    tmp <- cbind(rep(i, nrow(test_data)), test_data[, c("index", 'OID_', "exist")], 
                 predict.gbm(gbm_models[[i]], test_data, n.trees = gbm_models[[i]]$n.trees, type = "response")) %>% 
      as.data.frame()
    colnames(tmp) <- c('round', "index", 'OID_', 'case', 'pred')
    if (is.null(test_predictions)) test_predictions <- tmp else test_predictions <- rbind(test_predictions, tmp)
    
    tmp <- cbind(rep(i, nrow(prediction_data)), prediction_data[, c("index", 'OID_', "exist")], 
                 predict.gbm(gbm_models[[i]], prediction_data, n.trees = gbm_models[[i]]$n.trees, type = "response")) %>% 
      as.data.frame()
    colnames(tmp) <- c('round', "index", 'OID_', 'case', 'pred')
    if (is.null(all_predictions)) all_predictions <- tmp else all_predictions <- rbind(all_predictions, tmp)
    
    tmp <- summary(gbm_models[[i]])
    colnames(tmp) <- c('var', paste0('rel.inf.', i))
    if (is.null(variable_contributions)) variable_contributions <- tmp else variable_contributions <- merge(variable_contributions, tmp, by = 'var', sort = FALSE)
    
    mean_train <- NULL
    response_curve_tmp <- NULL
    for (vi in predictor_indices) {
      if (is.null(mean_train)) mean_train <- as.data.frame(combined_data)
      mean_train[, vi] <- mean(mean_train[, vi])
    }
    for (rsi in predictor_indices) {
      tmp <- mean_train[1:250, ]
      tmp[rsi] <- seq(from = min(as.data.frame(combined_data)[rsi]), to = max(as.data.frame(combined_data)[rsi]), length.out = 250)
      x_i <- tmp[rsi]
      y_i <- predict.gbm(gbm_models[[i]], tmp, n.trees = gbm_models[[i]]$n.trees, type = "response")
      y_i <- scale(y_i, center = TRUE, scale = FALSE)
      tmp_pred <- cbind(x_i, y_i)
      names(tmp_pred) <- c(names(as.data.frame(combined_data))[rsi], paste0("y.", rsi))
      if (rsi == predictor_indices[1]) response_curve_tmp <- cbind(cbind(rep(i, nrow(tmp_pred)), 1:nrow(tmp_pred), tmp_pred)) else {
        response_curve_tmp <- cbind(response_curve_tmp, tmp_pred)
      }
    }
    colnames(response_curve_tmp)[1:2] <- c('round', 'point')
    if (is.null(response_curves)) response_curves <- response_curve_tmp else response_curves <- rbind(response_curves, response_curve_tmp)
  }
  
  # Save modeling results
  all_predictions <- all_predictions %>% st_drop_geometry()
  training_predictions <- training_predictions %>% st_drop_geometry()
  test_predictions <- test_predictions %>% st_drop_geometry()
  
  write_rds(response_curves, paste0(output_dir, '/', nipah_type, "_response_curves.rds"))
  write_rds(variable_contributions, paste0(output_dir, '/', nipah_type, "_variable_contributions.rds"))
  write_rds(training_predictions, paste0(output_dir, '/', nipah_type, "_training_predictions.rds"))
  write_rds(test_predictions, paste0(output_dir, '/', nipah_type, "_test_predictions.rds"))
  write_rds(all_predictions, paste0(output_dir, '/', nipah_type, "_all_predictions.rds"))
  
  # Calculate average predictions for map
  avg_all_predictions <- aggregate(cbind(pred) ~ OID_, data = all_predictions, mean)
  
  # Prediction map
  prediction_map <- data.frame(OID_ = avg_all_predictions$OID_, probability = avg_all_predictions$pred)
  env_vars %>% left_join(prediction_map, by = "OID_") %>% 
    ggplot() +
    geom_sf(color = "transparent", aes(fill = probability)) + 
    scale_fill_gradient2(low = "#368BB1", mid = "#F9F766", high = "#E81D16", 
                         midpoint = 0.5, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), 
                         na.value = "#368BB1", name = str_wrap("Probability", width = 24)) +
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
  ggsave(paste0(output_dir, '/', nipah_type, '_prediction_map.pdf'), device = "pdf", width = 10, height = 10)
  write_rds(prediction_map, paste0(output_dir, '/', nipah_type, '_prediction_map.rds'))
}

# Run models
run_model('nipah_human')
run_model('nipah_zoonotic')

# Stop parallel computing
stopCluster(cl)