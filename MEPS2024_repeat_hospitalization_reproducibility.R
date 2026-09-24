# Reproducibility script
# Selected Five-Condition Chronic Disease Count and Repeat Hospitalization
# among U.S. Adults With Inpatient Hospital Use: MEPS 2024
#
# Input:
#   MEPS HC-256 2024 Full-Year Consolidated public-use file (.dta)
#
# This script is cleaned from the original working analysis.
# It preserves the final analysis specifications used in the manuscript.
# It does not install packages, clear the workspace, or use a machine-specific path.

library(haven)
library(dplyr)
library(survey)
library(ggplot2)
library(scales)

options(survey.lonely.psu = "adjust")

# 1. Load HC-256
data_file <- "h256.dta"

if (!file.exists(data_file)) {
  stop(
    "HC-256 data file not found. Place the MEPS 2024 HC-256 .dta file ",
    "in the working directory and update data_file if its filename differs."
  )
}

h256 <- read_dta(data_file)

required_vars <- c(
  "DUPERSID", "AGE24X", "SEX", "RACETHX", "POVCAT24", "INSCOV24",
  "IPDIS24", "HIBPDX", "DIABDX_M18", "CHDDX", "ARTHDX", "ASTHDX",
  "PERWT24F", "VARSTR", "VARPSU"
)

missing_vars <- setdiff(required_vars, names(h256))
if (length(missing_vars) > 0) {
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
}

# 2. Construct analysis variables among all positive-weight adults
meps_adults <- h256 %>%
  filter(AGE24X >= 18, PERWT24F > 0) %>%
  mutate(
    hypertension = case_when(HIBPDX == 1 ~ 1, HIBPDX == 2 ~ 0, TRUE ~ NA_real_),
    diabetes = case_when(DIABDX_M18 == 1 ~ 1, DIABDX_M18 == 2 ~ 0, TRUE ~ NA_real_),
    chd = case_when(CHDDX == 1 ~ 1, CHDDX == 2 ~ 0, TRUE ~ NA_real_),
    arthritis = case_when(ARTHDX == 1 ~ 1, ARTHDX == 2 ~ 0, TRUE ~ NA_real_),
    asthma = case_when(ASTHDX == 1 ~ 1, ASTHDX == 2 ~ 0, TRUE ~ NA_real_),
    complete_conditions = complete.cases(hypertension, diabetes, chd, arthritis, asthma),
    analytic_domain = IPDIS24 >= 1 & complete_conditions,
    repeat_hosp = if_else(analytic_domain, as.numeric(IPDIS24 >= 2), NA_real_),
    condition_count = if_else(
      complete_conditions,
      hypertension + diabetes + chd + arthritis + asthma,
      NA_real_
    ),
    condition_cat = factor(
      case_when(
        condition_count == 0 ~ "0",
        condition_count == 1 ~ "1",
        condition_count == 2 ~ "2",
        condition_count == 3 ~ "3",
        condition_count >= 4 ~ "4-5",
        TRUE ~ NA_character_
      ),
      levels = c("0", "1", "2", "3", "4-5")
    ),
    multiple_selected = case_when(
      condition_count >= 2 ~ 1,
      condition_count <= 1 ~ 0,
      TRUE ~ NA_real_
    ),
    age10 = AGE24X / 10,
    sex_f = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    race_f = factor(
      RACETHX,
      levels = c(2, 1, 3, 4, 5),
      labels = c(
        "Non-Hispanic White", "Hispanic", "Non-Hispanic Black",
        "Non-Hispanic Asian", "Non-Hispanic Other/Multiple"
      )
    ),
    poverty_f = factor(
      POVCAT24,
      levels = c(5, 1, 2, 3, 4),
      labels = c("High income", "Poor/Negative", "Near poor", "Low income", "Middle income")
    ),
    insurance_f = factor(
      INSCOV24,
      levels = c(1, 2, 3),
      labels = c("Any private", "Public only", "Uninsured")
    )
  )

# 3. Verify analytic sample
analytic_n <- sum(meps_adults$analytic_domain)
stopifnot(analytic_n == 1247)

condition_counts <- table(meps_adults$condition_count[meps_adults$analytic_domain])
stopifnot(
  identical(as.integer(condition_counts), c(226L, 318L, 337L, 236L, 112L, 18L))
)

outcome_counts <- table(meps_adults$repeat_hosp[meps_adults$analytic_domain])
stopifnot(identical(as.integer(outcome_counts), c(963L, 284L)))

# 4. MEPS complex survey design and hospitalized analytic domain
meps_design <- svydesign(
  ids = ~VARPSU,
  strata = ~VARSTR,
  weights = ~PERWT24F,
  data = meps_adults,
  nest = TRUE
)

analytic_design <- subset(meps_design, analytic_domain)

cat("Analytic N:", nrow(analytic_design$variables), "\n")
cat("Survey design degrees of freedom:", degf(analytic_design), "\n")

# 5. Descriptive estimates
overall_repeat <- svymean(~repeat_hosp, design = analytic_design, na.rm = TRUE)
print(overall_repeat)
print(confint(overall_repeat))

weighted_by_count <- svyby(
  ~repeat_hosp,
  ~condition_cat,
  design = analytic_design,
  FUN = svymean,
  vartype = "ci",
  na.rm = TRUE
)
print(weighted_by_count)

age_est <- svymean(~AGE24X, analytic_design, na.rm = TRUE)
sex_est <- svymean(~factor(SEX), analytic_design, na.rm = TRUE)
race_est <- svymean(~factor(RACETHX), analytic_design, na.rm = TRUE)
poverty_est <- svymean(~factor(POVCAT24), analytic_design, na.rm = TRUE)
insurance_est <- svymean(~factor(INSCOV24), analytic_design, na.rm = TRUE)
count_est <- svymean(~factor(condition_count), analytic_design, na.rm = TRUE)
repeat_est <- svymean(~factor(repeat_hosp), analytic_design, na.rm = TRUE)

# 6. Primary continuous-count model
primary_model <- svyglm(
  repeat_hosp ~ condition_count + age10 + sex_f + race_f + poverty_f + insurance_f,
  design = analytic_design,
  family = quasibinomial()
)

print(summary(primary_model))

primary_results <- data.frame(
  term = names(coef(primary_model)),
  beta = coef(primary_model),
  SE = sqrt(diag(vcov(primary_model))),
  OR = exp(coef(primary_model)),
  CI_lower = exp(confint(primary_model)[, 1]),
  CI_upper = exp(confint(primary_model)[, 2]),
  p_value = summary(primary_model)$coefficients[, 4],
  row.names = NULL
)

print(primary_results)
write.csv(primary_results, "Table2_primary_model_full.csv", row.names = FALSE)

# 7. Categorical exposure analysis and overall Wald test
categorical_model <- svyglm(
  repeat_hosp ~ condition_cat + age10 + sex_f + race_f + poverty_f + insurance_f,
  design = analytic_design,
  family = quasibinomial()
)

print(summary(categorical_model))
print(exp(cbind(OR = coef(categorical_model), confint(categorical_model))))

categorical_wald <- regTermTest(categorical_model, ~condition_cat)
print(categorical_wald)

# 8. Binary sensitivity analysis: >=2 versus 0-1 selected conditions
binary_model <- svyglm(
  repeat_hosp ~ multiple_selected + age10 + sex_f + race_f + poverty_f + insurance_f,
  design = analytic_design,
  family = quasibinomial()
)

print(summary(binary_model))
print(exp(cbind(OR = coef(binary_model), confint(binary_model))))

# 9. Quadratic functional-form assessment
quadratic_model <- svyglm(
  repeat_hosp ~ condition_count + I(condition_count^2) +
    age10 + sex_f + race_f + poverty_f + insurance_f,
  design = analytic_design,
  family = quasibinomial()
)

quadratic_test <- regTermTest(quadratic_model, ~I(condition_count^2))
print(quadratic_test)

# 10. Rao-Scott likelihood-ratio comparison:
# continuous 0-5 count versus fully categorical 0-5 count
analytic_design$variables$count_linear <- analytic_design$variables$condition_count
analytic_design$variables$count_factor <- factor(
  analytic_design$variables$condition_count,
  levels = c(0, 1, 2, 3, 4, 5)
)

linear_shape_model <- svyglm(
  repeat_hosp ~ count_linear + age10 + sex_f + race_f + poverty_f + insurance_f,
  design = analytic_design,
  family = quasibinomial()
)

full_shape_model <- svyglm(
  repeat_hosp ~ count_factor + age10 + sex_f + race_f + poverty_f + insurance_f,
  design = analytic_design,
  family = quasibinomial()
)

rao_scott_test <- anova(linear_shape_model, full_shape_model, method = "LRT")
print(rao_scott_test)

# 11. Disease-combination count used in the limitations discussion
combo_check <- meps_adults %>%
  filter(analytic_domain) %>%
  mutate(
    disease_combination = paste0(
      ifelse(hypertension == 1, "HTN+", ""),
      ifelse(diabetes == 1, "DM+", ""),
      ifelse(chd == 1, "CHD+", ""),
      ifelse(arthritis == 1, "ART+", ""),
      ifelse(asthma == 1, "AST+", "")
    ),
    disease_combination = ifelse(
      disease_combination == "",
      "None",
      sub("\\+$", "", disease_combination)
    )
  ) %>%
  count(disease_combination, sort = TRUE)

print(combo_check, n = Inf)
cat("Observed disease combinations:", nrow(combo_check), "\n")
stopifnot(nrow(combo_check) == 31)

write.csv(combo_check, "disease_combination_counts.csv", row.names = FALSE)

# 12. Figure 1: survey-weighted prevalence by selected condition count
figure1_data <- weighted_by_count %>%
  mutate(
    condition_cat = factor(
      condition_cat,
      levels = c("0", "1", "2", "3", "4-5"),
      labels = c("0", "1", "2", "3", "4-5")
    ),
    prevalence = repeat_hosp
  )

fig1 <- ggplot(figure1_data, aes(x = condition_cat, y = prevalence)) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), width = 0.12, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = percent(prevalence, accuracy = 0.1)), vjust = -1.1, size = 3.8) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.52),
    breaks = seq(0, 0.5, 0.1)
  ) +
  labs(
    title = "Repeat Hospitalization by Selected Five-Condition Count",
    subtitle = "U.S. adults with inpatient hospital use, MEPS 2024",
    x = "Number of selected diagnosed chronic conditions",
    y = "Repeat hospitalization (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(
  "Figure1_Repeat_Hospitalization_MEPS2024.png",
  plot = fig1, width = 8, height = 5, dpi = 300
)

# 13. Figure 2: primary adjusted model
figure2_data <- primary_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    variable = c(
      "Selected condition count (per 1 condition)",
      "Age (per 10 years)",
      "Female vs Male",
      "Hispanic vs NH White",
      "NH Black vs NH White",
      "NH Asian vs NH White",
      "NH Other/Multiple vs NH White",
      "Poor/Negative vs High income",
      "Near poor vs High income",
      "Low income vs High income",
      "Middle income vs High income",
      "Public only vs Any private",
      "Uninsured vs Any private"
    ),
    label = sprintf("%.2f (%.2f-%.2f)", OR, CI_lower, CI_upper)
  )

figure2_data$variable <- factor(figure2_data$variable, levels = rev(figure2_data$variable))

fig2 <- ggplot(figure2_data, aes(x = OR, y = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.15, linewidth = 0.6) +
  geom_point(size = 2.7) +
  geom_text(aes(x = 5.4, label = label), hjust = 0, size = 3.2) +
  scale_x_log10(limits = c(0.1, 8), breaks = c(0.1, 0.3, 1, 3)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Adjusted Associations With Repeat Hospitalization",
    subtitle = "U.S. adults with inpatient hospital use, MEPS 2024",
    x = "Adjusted odds ratio (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 135, 10, 10)
  )

ggsave(
  "Figure2_Adjusted_OR_MEPS2024.png",
  plot = fig2, width = 10, height = 6, dpi = 300
)

# 14. Figure 3: categorical exposure analysis
cat_coef <- coef(categorical_model)
cat_ci <- confint(categorical_model)
cat_terms <- grep("^condition_cat", names(cat_coef), value = TRUE)

figure3_data <- data.frame(
  term = cat_terms,
  OR = exp(cat_coef[cat_terms]),
  lower = exp(cat_ci[cat_terms, 1]),
  upper = exp(cat_ci[cat_terms, 2]),
  row.names = NULL
)

figure3_data$comparison <- c(
  "1 condition vs 0",
  "2 conditions vs 0",
  "3 conditions vs 0",
  "4-5 conditions vs 0"
)

figure3_data$label <- sprintf(
  "%.2f (%.2f-%.2f)",
  figure3_data$OR,
  figure3_data$lower,
  figure3_data$upper
)

figure3_data$comparison <- factor(
  figure3_data$comparison,
  levels = rev(figure3_data$comparison)
)

fig3 <- ggplot(figure3_data, aes(x = OR, y = comparison)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.12, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(x = 7.0, label = label), hjust = 0, size = 3.5) +
  scale_x_log10(limits = c(0.5, 11), breaks = c(0.5, 1, 2, 3, 5)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Selected Five-Condition Count and Repeat Hospitalization",
    subtitle = "Reference: 0 conditions | Overall Wald test: F(4, 138) = 3.28, p = .013",
    x = "Adjusted odds ratio (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 125, 10, 10)
  )

ggsave(
  "Figure3_Categorical_Condition_Count_MEPS2024.png",
  plot = fig3, width = 9.5, height = 5, dpi = 300
)

# 15. Save key derived outputs
write.csv(
  weighted_by_count,
  "weighted_repeat_hospitalization_by_count.csv",
  row.names = FALSE
)

capture.output(
  list(
    analytic_n = nrow(analytic_design$variables),
    survey_design_df = degf(analytic_design),
    overall_repeat = overall_repeat,
    overall_repeat_ci = confint(overall_repeat),
    primary_model = summary(primary_model),
    categorical_model = summary(categorical_model),
    categorical_wald = categorical_wald,
    binary_model = summary(binary_model),
    quadratic_test = quadratic_test,
    rao_scott_test = rao_scott_test
  ),
  file = "MEPS2024_analysis_output.txt"
)

cat("\nReproducibility analysis completed successfully.\n")
