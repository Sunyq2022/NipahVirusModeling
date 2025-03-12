# Nipah Virus Spatial Modeling and Visualization

This repository contains R scripts for modeling and evaluating the spatial distribution of Nipah virus outbreaks.

## Repository Structure
├── Data/
│   ├── sasea_map.rds              # Processed South Asia and Southeast Asia map
│   ├── ten_segment_line.rds       # Processed ten-segment line map
│   ├── nipah_occurrences.rds      # Processed Nipah virus occurrence points
│   ├── env_vars_processed.rds     # Processed environmental variables grid
│   ├── bat_pathogens_sasea.rds    # Processed bat pathogen data
│   ├── chikv_dengue_sasea.rds     # Processed Chikungunya and Dengue data
│   ├── bat_deng_chikv.rds         # Combined bat, Chikungunya, and Dengue data
├── Outputs/                       # Output directory (not tracked)
├── model_and_visualize.R          # Script to model and generate prediction maps
├── evaluate_and_visualize.R       # Script for model evaluation and visualization
└── README.md                      # This documentation