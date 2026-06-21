# =============================================================================
# meta_analysis_helpers.r
# =============================================================================
# Purpose: Shared helper functions for meta-analysis scripts
# =============================================================================

# =============================================================================
# FUNCTION: prepare_forest_data
# =============================================================================
# Prepare data for forest plot with hierarchical grouping
# 
# Parameters:
#   data: Data frame with effect sizes and metadata
#   effect_col: Name of the effect size column (e.g., "fisher_z" or "kappa")
#   var_col: Name of the variance column (e.g., "var_z" or "var_kappa")
#   rve_model: RVE model object (from robu)
#   convert_to_r: Should effect sizes be converted back to r? (TRUE for Fisher's z)
#
# Returns:
#   Data frame sorted and formatted for forest plot with header rows
# =============================================================================

prepare_forest_data <- function(data, effect_col, var_col, rve_model, convert_to_r = TRUE) {
  library(dplyr)
  
  # Add row index before any transformations for matching weights from robu model
  data <- data %>% mutate(.row_idx = row_number())
  
  # Extract weights from robu model and match by row index
  # The robu model's data.full preserves the original row order of the input data
  # We use effect size values to match weights back to the correct rows
  # This handles cases where robu() may have reordered data internally
  
  # Get weights from model - these are in the order of rve_model$data.full
  model_weights <- rve_model$data.full$r.weights
  model_effect_sizes <- rve_model$data.full$effect.size
  
  # Match weights by effect size values (unique identifier for each observation)
  # Round to handle floating point comparison issues
  data_effect_sizes <- data[[effect_col]]
  
  # Create a lookup: for each row in data, find matching row in model data
  matched_weights <- sapply(seq_len(nrow(data)), function(i) {
    # Find the matching row in model data by effect size
    matches <- which(abs(model_effect_sizes - data_effect_sizes[i]) < 1e-10)
    if (length(matches) == 1) {
      return(model_weights[matches])
    } else if (length(matches) > 1) {
      # If multiple matches (rare), use the first unused one
      return(model_weights[matches[1]])
    } else {
      # Fallback to inverse-variance weight
      return(1 / data[[var_col]][i])
    }
  })
  
  # Normalize weights to sum to 100
  matched_weights <- matched_weights / sum(matched_weights) * 100
  
  # Prepare base data
  data_forest <- data %>%
    mutate(
      es_display = if (convert_to_r) tanh(.data[[effect_col]]) else .data[[effect_col]],
      ci_lower = if (convert_to_r) {
        tanh(.data[[effect_col]] - 1.96 * sqrt(.data[[var_col]]))
      } else {
        .data[[effect_col]] - 1.96 * sqrt(.data[[var_col]])
      },
      ci_upper = if (convert_to_r) {
        tanh(.data[[effect_col]] + 1.96 * sqrt(.data[[var_col]]))
      } else {
        .data[[effect_col]] + 1.96 * sqrt(.data[[var_col]])
      },
      weight = matched_weights,
      study = gsub("^([^,]+),.*?(\\(\\d{4}\\))", "\\1 et al. \\2", author_year),
      sample = sample_subsample
    ) %>%
    group_by(study_group) %>%
    mutate(study_mean_es = mean(es_display)) %>%
    ungroup() %>%
    arrange(study_mean_es, sample) %>%
    mutate(ID = row_number())
  
  # Add header rows for studies with multiple samples
  for (sg in unique(data_forest$study_group)) {
    data.grp <- data_forest[data_forest$study_group == sg, ]
    if (nrow(data.grp) > 1) {
      minID <- min(data.grp$ID)
      study_name <- data.grp$study[1]
      # Move study-level symbols (*/^) onto the sample labels for multi-sample
      # studies so the annotation is consistently shown on the sample, not the
      # study header (single-sample studies keep the symbol on the author).
      study_symbols <- regmatches(study_name, regexpr("[*^]+", study_name))
      if (length(study_symbols) == 0) study_symbols <- ""
      study_name <- trimws(gsub("[*^]+", "", study_name))
      # Indent sub-sample labels (append symbols to each sample)
      indx <- data_forest$study_group == sg
      data_forest$study[indx] <- paste0("        ", data_forest$sample[indx], study_symbols)
      # Add header row
      header_row <- data.frame(
        study_group = sg,
        study = paste0("    ", study_name),
        ID = minID - 0.5,
        n_total = sum(data.grp$n_total),
        weight = sum(data.grp$weight),
        es_display = NA,
        ci_lower = NA,
        ci_upper = NA,
        study_mean_es = data.grp$study_mean_es[1],
        stringsAsFactors = FALSE
      )
      data_forest <- bind_rows(data_forest, header_row)
    } else {
      # Single-sample studies
      indx <- data_forest$study_group == sg
      data_forest$study[indx] <- paste0("    ", data_forest$study[indx])
    }
  }
  
  # Re-sort by ID
  data_forest <- data_forest %>% arrange(ID)
  
  return(data_forest)
}


# =============================================================================
# FUNCTION: create_forest_plot
# =============================================================================
# Create forest plot with hierarchical study grouping
#
# Parameters:
#   data_forest: Prepared data from prepare_forest_data()
#   output_file: Path for output TIFF file
#   main_title: Main title for the plot
#   xlab: X-axis label (e.g., "Pearson's r" or "Cohen's kappa")
#   xlim: X-axis limits (default c(-5, 2.5))
#   alim: Axis limits for forest plot (default c(-1, 1))
# =============================================================================

create_forest_plot <- function(data_forest, output_file, main_title, xlab,
                                xlim = c(-5, 2.5), alim = c(-1, 1)) {
  library(metafor)
  
  # Source forest.col functions
  source("Utilities/forest.col.r")
  
  # Split header rows vs data rows
  n_total_rows <- nrow(data_forest)
  studyHeadingLocations <- (1:n_total_rows)[is.na(data_forest$es_display)]
  sampleLocations <- (1:n_total_rows)[!is.na(data_forest$es_display)]
  studyHeadings <- data_forest$study[is.na(data_forest$es_display)]
  data_for_plot <- data_forest[!is.na(data_forest$es_display), ]
  
  # Set up plot parameters
  cex <- 1.5
  ylim <- c(-1, n_total_rows + 3)
  
  # Prepare additional info columns
  ilab <- cbind(data_for_plot$n_total, format(round(data_for_plot$weight, 2), trim = FALSE))
  ilab.xpos <- c(-1.8, -1.2)
  ilab.pos <- 2
  anno.pos <- 2.5
  
  # Calculate row positions
  row_positions_all <- n_total_rows:1
  sampleRowPositions <- row_positions_all[sampleLocations]
  headerRowPositions <- row_positions_all[studyHeadingLocations]
  
  # Calculate dimensions
  n_studies <- nrow(data_forest)
  plot_height <- 3 + 30 * n_studies / 75
  
  # Create TIFF file
  tiff(
    filename = output_file,
    width = 18,
    height = plot_height,
    units = "cm",
    res = 600,
    compression = "lzw",
    pointsize = 6
  )
  
  # Plot forest
  forest.col(
    x = data_for_plot$es_display,
    vi = (data_for_plot$ci_upper - data_for_plot$ci_lower)^2 / (2 * 1.96)^2,
    xlim = xlim,
    alim = alim,
    ylim = ylim,
    digits = 2,
    xlab = xlab,
    main = main_title,
    slab = data_for_plot$study,
    ilab = ilab,
    ilab.xpos = ilab.xpos,
    ilab.pos = ilab.pos,
    rows = sampleRowPositions,
    efac = 0.4,
    cex = cex,
    refline = 0,
    anno.pos = anno.pos
  )
  
  # Add column labels
  text(
    x = c(ilab.xpos, anno.pos),
    y = n_total_rows + 2,
    labels = c("N", "Weight", paste0(xlab, " [95% CI]")),
    cex = cex,
    pos = ilab.pos,
    font = 2
  )
  
  # Add header rows (study names only, no N/Weight)
  if (length(studyHeadingLocations) > 0) {
    text(
      x = xlim[1],
      y = headerRowPositions,
      labels = studyHeadings,
      cex = cex,
      pos = 4,
      font = 1
    )
  }
  
  dev.off()
  
  cat(sprintf("  - Saved: %s\n", output_file))
}


# =============================================================================
# FUNCTION: add_pooled_effect_to_plot
# =============================================================================
# Add pooled effect summary to an existing forest plot
#
# Parameters:
#   rve_model: RVE model object
#   n_total_rows: Total number of rows in forest plot
#   convert_to_r: Should effect be converted back to r? (TRUE for Fisher's z)
#   cex: Character expansion factor
#   anno.pos: Position for annotation
#   ilab.xpos: Positions for N and Weight columns
# =============================================================================

add_pooled_effect_to_plot <- function(rve_model, n_total_rows, convert_to_r = TRUE,
                                       cex = 1.5, anno.pos = 2.5, ilab.xpos = c(-1.8, -1.2)) {
  source("Utilities/forest.col.r")
  
  # Calculate pooled effect
  if (convert_to_r) {
    pooled_es <- tanh(rve_model$reg_table$b.r)
    pooled_ci_lower <- tanh(rve_model$reg_table$CI.L)
    pooled_ci_upper <- tanh(rve_model$reg_table$CI.U)
  } else {
    pooled_es <- rve_model$reg_table$b.r
    pooled_ci_lower <- rve_model$reg_table$CI.L
    pooled_ci_upper <- rve_model$reg_table$CI.U
  }
  
  # Add pooled effect polygon
  addpoly.col(
    x = pooled_es,
    ci.lb = pooled_ci_lower,
    ci.ub = pooled_ci_upper,
    rows = -1,
    cex = cex,
    mlab = "Overall RE Model",
    efac = 0.4,
    digits = 2,
    font = 2,
    anno.pos = anno.pos
  )
  
  # Add total N and weight
  total_n <- sum(rve_model$data.full$N, na.rm = TRUE)
  text(
    x = c(ilab.xpos[1], ilab.xpos[2]),
    y = -1,
    labels = c(total_n, "100.00"),
    cex = cex,
    pos = 2,
    font = 2
  )
}
