# Shared calculation and rendering for matched cross-cohort contrasts.
# dplyr/ggplot2 resolve column names through non-standard evaluation.
# nolint start: object_usage_linter.
harmonize_trait <- function(trait) {
  trait <- sub("TF_perc$", "F_perc", trait)
  dplyr::recode(trait, FG = "GLU", PCr_NMR = "SCR")
}

matched_cohort_contrasts <- function(result_table) {
  models <- result_table |>
    tidyr::unnest(modres) |>
    dplyr::filter(
      !(cohort == "ALSPAC" & grepl("F_perc$", trait) &
          !grepl("TF_perc$", trait))
    ) |>
    dplyr::mutate(trait = harmonize_trait(trait)) |>
    dplyr::filter(trait != "WHR")

  # Do not pool sex-specific quartiles or recruitment/imaging visits.
  ukb_models <- models |>
    dplyr::filter(cohort == "UKB") |>
    dplyr::select(sex, trait, ukb_visit = reftime,
                  ukb_model = mod1, ukb_summary = bmsum)
  if (anyDuplicated(ukb_models[c("sex", "trait")])) {
    stop("Multiple UKB models per sex/trait; select a visit explicitly.")
  }
  replication_models <- models |>
    dplyr::filter(cohort %in% c("ALSPAC", "RS"))
  unmatched <- dplyr::anti_join(
    replication_models, ukb_models, by = c("sex", "trait")
  )
  if (nrow(unmatched) > 0) {
    warning("Replication models without UKB matches: ", paste(
      paste(unmatched$cohort, unmatched$sex, unmatched$trait, sep = "/"),
      collapse = ", "
    ))
  }
  pairs <- dplyr::inner_join(
    replication_models, ukb_models, by = c("sex", "trait")
  )
  if (nrow(pairs) == 0) {
    stop("No replication and UKB models share sex and harmonized trait.")
  }

  result <- purrr::pmap_dfr(
    pairs[c("cohort", "sex", "trait", "reftime", "mod1", "bmsum",
            "ukb_visit", "ukb_model", "ukb_summary")],
    function(cohort, sex, trait, reftime, mod1, bmsum,
             ukb_visit, ukb_model, ukb_summary) {
      bmi_values <- as.numeric(bmsum$bmi[2:4])
      slope_values <- as.numeric(bmsum$Slope[c(2, 4)])
      if (any(!is.finite(c(bmi_values, slope_values))) ||
            slope_values[2] <= slope_values[1]) {
        stop("Invalid BMI/slope quartiles for ", cohort, "/", sex)
      }
      purrr::map2_dfr(
        list(mod1, ukb_model), c(cohort, "UKB"),
        function(model, estimate_cohort) {
          # Models use an identity link: contrasts are in biomarker units.
          if (model$family$link != "identity") {
            stop("Matched contrasts require identity-link models.")
          }
          # Keep other covariates at each model's reference values. Age is
          # cohort-specific; these are not age-standardized comparisons.
          new_data <- as.data.frame(lapply(
            model$var.summary,
            function(x) if (is.numeric(x)) x[2] else x[1]
          ))
          model_summary <- if (estimate_cohort == "UKB") {
            ukb_summary
          } else {
            bmsum
          }
          new_data$age <- model_summary$age[3]
          new_data <- new_data[rep(1, 3), , drop = FALSE]
          new_data$bmi <- bmi_values
          new_data$Slope <- slope_values[1]
          low_matrix <- predict(model, newdata = new_data, type = "lpmatrix")
          new_data$Slope <- slope_values[2]
          high_matrix <- predict(model, newdata = new_data, type = "lpmatrix")
          # Covariance of both predictions is essential for the paired CI.
          contrast <- high_matrix - low_matrix
          estimate <- as.vector(contrast %*% coef(model))
          std_error <- sqrt(rowSums((contrast %*% model$Vp) * contrast))
          tibble::tibble(
            reference_cohort = cohort, estimate_cohort, sex, trait,
            reference_visit = reftime,
            model_visit = if (estimate_cohort == "UKB") ukb_visit else reftime,
            bmi_percentile = c("p25", "p50", "p75"), bmi = bmi_values,
            slope_p25 = slope_values[1], slope_p75 = slope_values[2],
            age = model_summary$age[3], estimate, std_error,
            conf_low = estimate - qnorm(0.975) * std_error,
            conf_high = estimate + qnorm(0.975) * std_error,
            # Only marginal central-range support is available in bmsum.
            # This flag does not establish joint BMI/slope support.
            within_ukb_central_range =
              bmi_values >= ukb_summary$bmi[1] &
              bmi_values <= ukb_summary$bmi[5] &
              slope_values[1] >= ukb_summary$Slope[1] &
              slope_values[2] <= ukb_summary$Slope[5]
          )
        }
      )
    }
  )

  # ALSPAC android fat percentage was provided inverted (as lean percentage:
  # lean% = 100 - fat%). The ALSPAC model therefore predicts lean%, so its
  # contrast is a difference in lean%, not fat%. Because the models use an
  # identity link, converting to fat% is an affine transform with slope -1:
  # the estimate flips sign with unchanged magnitude and standard error, and
  # the confidence interval mirrors around zero. The paired UKB estimates are
  # already on the fat% scale and are left unchanged.
  inverted <- result$reference_cohort == "ALSPAC" &
    result$estimate_cohort == "ALSPAC" &
    result$trait == "AndF_perc"
  result$estimate[inverted] <- -result$estimate[inverted]
  lo <- result$conf_low[inverted]
  hi <- result$conf_high[inverted]
  result$conf_low[inverted] <- -hi
  result$conf_high[inverted] <- -lo
  result
}

matched_cohort_plot <- function(contrasts, reference) {
  biomarker_groups <- tibble::tibble(
    trait = c(
      "ALT", "AST", "GGT", "HDL", "LDL", "TG", "SBP", "DBP", "GLU",
      "CRP", "SCR", "BF_perc", "AndF_perc", "GynF_perc",
      "ArmsF_perc", "LegsF_perc", "TrunkF_perc"
    ),
    group = rep(
      c("Liver", "Lipids", "BP", "Glucose", "Inflam", "Kidney",
        "Imaging"),
      c(3, 3, 2, 1, 1, 1, 6)
    ),
    block = rep(c(1, 2), c(9, 8))
  )
  plot_data <- contrasts |>
    dplyr::filter(reference_cohort == reference) |>
    dplyr::inner_join(biomarker_groups, by = "trait") |>
    dplyr::mutate(
      group = factor(group, unique(biomarker_groups$group)),
      trait = factor(trait, biomarker_groups$trait),
      estimate_cohort = factor(estimate_cohort, c(reference, "UKB")),
      # Discrete spacing keeps paired estimates aligned; labels are actual
      # BMI values, not percentile names. Predictions remain unrounded.
      bmi_label = factor(bmi, levels = sort(unique(bmi)))
    )
  slope_values <- plot_data |>
    dplyr::distinct(sex, slope_p25, slope_p75) |>
    dplyr::arrange(sex)
  slope_caption <- paste0(
    if (reference == "RS") "Rotterdam" else reference,
    " slope (kg/m²/year), p25 / p75: ",
    paste(sprintf(
      "%s %.3f / %.3f", slope_values$sex,
      slope_values$slope_p25, slope_values$slope_p75
    ), collapse = "; "),
    ".\nUKB contrasts use these same slope values."
  )
  panels <- purrr::map(sort(unique(plot_data$block)), function(block_value) {
    panel_data <- dplyr::filter(plot_data, block == block_value)
    ggplot2::ggplot(panel_data, ggplot2::aes(bmi_label, estimate)) +
      ggplot2::geom_hline(yintercept = 0, lty = "11", linewidth = 0.2) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = conf_low, ymax = conf_high,
                     group = estimate_cohort),
        linewidth = 0.3, position = ggplot2::position_dodge(width = 0.25)
      ) +
      ggplot2::geom_point(
        ggplot2::aes(color = estimate_cohort, group = estimate_cohort),
        size = 1.2, position = ggplot2::position_dodge(width = 0.25)
      ) +
      # Share x scales down each sex column; label only the bottom edge.
      ggh4x::facet_nested(
        group + trait ~ sex, scales = "free", axes = "margins",
        labeller = ggplot2::labeller(
          trait = function(x) sub("_perc$", "%", x)
        ),
        strip = ggh4x::strip_nested(size = "variable")
      ) +
      ggplot2::scale_x_discrete(labels = function(x) {
        format(round(as.numeric(x), 2), trim = TRUE, nsmall = 2)
      }) +
      ggplot2::scale_color_manual(
        values = stats::setNames(c("#F8766D", "#00BFC4"),
                                 c(reference, "UKB")),
        labels = c(if (reference == "RS") "Rotterdam" else reference, "UKB")
      ) +
      ggplot2::labs(
        x = expression("BMI (kg/" * m^2 * ")"),
        y = if (block_value == 1) {
          "Difference between\n75th - 25th slope percentiles"
        } else {
          NULL
        },
        color = NULL
      ) +
      ggplot2::theme_bw()
  })
  patchwork::wrap_plots(panels, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      caption = slope_caption,
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(size = 7, hjust = 0)
      )
    ) &
    ggplot2::theme(
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 8),
      legend.key.width = grid::unit(1, "mm"),
      strip.text = ggplot2::element_text(
        size = 7, margin = ggplot2::margin(1, 1, 1, 1)
      ),
      panel.spacing = grid::unit(0.5, "mm"),
      axis.text.x = ggplot2::element_text(size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.title = ggplot2::element_text(size = 8),
      axis.ticks = ggplot2::element_line(linewidth = 0.2)
    )
}
# nolint end: object_usage_linter.
