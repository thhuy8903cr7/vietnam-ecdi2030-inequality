# ============================================================
# MICS6 ECDI2030 - ANALYSIS FOR THE JOURNAL OF PUBLIC HEALTH MANUSCRIPT
# ============================================================
# Run from the MICS6_ECDI2030_Analysis project root.
# The derived MICS6 input is obtained under UNICEF MICS data-access procedures;
# participant-level data are not part of the supplementary code package.
if (dir.exists("Rlib")) {
  .libPaths(c(normalizePath("Rlib"), .libPaths()))
}

required_packages <- c(
  "dplyr", "srvyr", "survey", "gtsummary", "ggplot2", "sf", "grid",
  "patchwork", "officer", "flextable", "ICEHmeasures", "healthequal",
  "marginaleffects", "rineq"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the analysis: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressWarnings(suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
}))

options(survey.lonely.psu = "adjust")
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. DATA IMPORT, CLEANING, AND VARIABLE RECODING
# ============================================================

mics_fill <- read.csv("mics6_data_final.csv")
# Select variables required for the analysis
data <- mics_fill %>%
  select(HH1, HL4, CAGE, ethnicity, cinsurance, HL12, HL13, HL16, HL17,
         cdisability, CA1, CA14, CA16, CA17, melevel,
         WB4, MSTATUS, MT9, MT10, MT1, MT2, MT3, LS1,
         DV1A, DV1B, DV1C, DV1D, DV1E, HC11, HC12, HC13,
         HH48, HH51, windex5, HH6, HH7, UB6, EC1,
         EC2A, EC2B, EC2C, stratum, chweight,
         num_range("EC", 21:40),
         starts_with("EC5"))
# Create the ECDI outcome variable
data <- data %>%
  mutate(EC27 = ifelse(EC26 %in% c("NO"), "NO", EC27))
data <- data %>%
  # STEP 1: Recode EC21 to EC40 into numeric scores
  mutate(
    # 1. EC21-EC38: YES = 1; NO/DK = 0; NO RESPONSE = missing
    across(EC21:EC38, ~ case_when(
      trimws(toupper(.)) == "YES" ~ 1,
      trimws(toupper(.)) %in% c("NO", "DK") ~ 0,
      trimws(toupper(.)) == "NO RESPONSE" | is.na(.) ~ NA_real_,
      TRUE ~ NA_real_
    )),
    
    # 2. Handle EC39 separately (internalizing behavior; reverse coded)
    EC39 = case_when(
      trimws(toupper(EC39)) %in% c("NEVER", "A FEW TIMES A YEAR", "MONTHLY") ~ 1,
      trimws(toupper(EC39)) %in% c("WEEKLY", "DAILY", "DK") ~ 0,
      trimws(toupper(EC39)) == "NO RESPONSE" | is.na(EC39) ~ NA_real_,
      TRUE ~ NA_real_
    ),
    
    # 3. Handle EC40 separately (externalizing behavior; reverse coded)
    EC40 = case_when(
      trimws(toupper(EC40)) %in% c("NOT AT ALL", "LESS") ~ 1,
      trimws(toupper(EC40)) %in% c("THE SAME", "MORE", "A LOT MORE", "DK") ~ 0,
      trimws(toupper(EC40)) == "NO RESPONSE" | is.na(EC40) ~ NA_real_,
      TRUE ~ NA_real_
    )
  ) %>%
  
  # STEPS 2 & 3: Calculate the total score; any NO RESPONSE is treated as missing
  mutate(
    # na.rm = FALSE ensures that any missing component makes the total score missing
    ecdi_total = rowSums(select(., EC21:EC40), na.rm = FALSE)
  ) %>%
  
  # STEP 4: Create the age-specific ECDI outcome variable
  mutate(
    # ECDI2030 coding for this analysis: 1 = on-track, 0 = not on-track.
    ecdi = case_when(
      is.na(CAGE) | is.na(ecdi_total) ~ NA_real_,
      
      CAGE >= 24 & CAGE <= 29 ~ ifelse(ecdi_total >= 7, 1, 0),
      CAGE >= 30 & CAGE <= 35 ~ ifelse(ecdi_total >= 9, 1, 0),
      CAGE >= 36 & CAGE <= 41 ~ ifelse(ecdi_total >= 11, 1, 0),
      CAGE >= 42 & CAGE <= 47 ~ ifelse(ecdi_total >= 13, 1, 0),
      CAGE >= 48 & CAGE <= 59 ~ ifelse(ecdi_total >= 15, 1, 0),
      
      TRUE ~ NA_real_ # Applies to children outside the 24-59 month age range
    )
  )
# Create component variables
## Number of stimulation activities by mother and father
data <- data %>%
  mutate(
    # STEP 1: Calculate total stimulation scores for mother and father
    # Use %in% to count MOTHER/FATHER without overwriting the original text values
    mother_stim_score = (EC5AA %in% "MOTHER") + (EC5BA %in% "MOTHER") + (EC5CA %in% "MOTHER") +
      (EC5DA %in% "MOTHER") + (EC5EA %in% "MOTHER") + (EC5FA %in% "MOTHER"),
    
    father_stim_score = (EC5AB %in% "FATHER") + (EC5BB %in% "FATHER") + (EC5CB %in% "FATHER") +
      (EC5DB %in% "FATHER") + (EC5EB %in% "FATHER") + (EC5FB %in% "FATHER"),
    
    # STEP 2: Classify maternal stimulation
    # Prioritize DK/NO RESPONSE before classifying the score
    mother_stim = case_when(
      if_any(c(EC5AA, EC5BA, EC5CA, EC5DA, EC5EA, EC5FA), ~ . %in% c("NO RESPONSE", "DK", "NO RESPONSE/DK")) ~ NA_real_,
      mother_stim_score >= 4 ~ 1,
      TRUE ~ 0
    ),
    mother_stim = factor(mother_stim, levels = c(1, 0)),
    # STEP 3: Classify paternal stimulation
    father_stim = case_when(
      if_any(c(EC5AB, EC5BB, EC5CB, EC5DB, EC5EB, EC5FB), ~ . %in% c("NO RESPONSE", "DK", "NO RESPONSE/DK")) ~ NA_real_,
      father_stim_score >= 4 ~ 1,
      TRUE ~ 0
    ),
    father_stim = factor(father_stim, levels = c(1, 0))
  )
## Toys
data <- data %>%
  mutate(
    # 1. Count toy types; only exact YES values contribute to the score
    total_toy = (EC2A %in% "YES") + (EC2B %in% "YES") + (EC2C %in% "YES"),
    # 2. Classify availability of at least two toy types
    toy = case_when(
      # Priority rule: if any item is DK/NO RESPONSE, classify as NO RESPONSE/DK
      if_any(c(EC2A, EC2B, EC2C), ~ . %in% c("NO RESPONSE", "DK", "NO RESPONSE/DK")) ~ NA_real_,
      # Two or more toy types = adequate (1)
      total_toy >= 2 ~ 1,
      
      # All other valid cases = not adequate (0)
      TRUE ~ 0
    ),
    toy = factor(toy, levels = c(1, 0))
  )
## Books
data <- data %>%
  mutate(
    book = case_when(
      # Priority rule: identify refusal/don't know responses first
      EC1 %in% c("NO RESPONSE", "DK", "NO RESPONSE/DK") ~ NA_real_,
      # Three or more books = adequate (1)
      EC1 %in% c("THREE", "FOUR", "FIVE", "SIX", "SEVEN", "EIGHT", "NINE", "TEN OR MORE BOOKS") ~ 1,
      # This recode assigns all remaining responses, including blanks, to 0.
      TRUE ~ 0
    ),
    book = factor(book, levels = c(1, 0))
  )
# Functional disability
data <- data %>%
  mutate(
    cdisability = factor(
      cdisability,
      levels = c("Has no functional difficulty", "Has functional difficulty")
    )
  )
# Maternal education variable
data <- data %>%
  mutate(
    melevel = factor(
      melevel,
      levels = c(
        "Pre-primary or none",
        "Primary",
        "Lower Secondary",
        "Upper Secondary",
        "Vocational High school",
        "University/College/Higher"
      )
    )
  )
# Ethnicity variable
data <- data %>%
  mutate(
    ethnicity = factor(
      ethnicity,
      levels = c("Kinh and Hoa", "Tay, Thai, Muong, Nung", "Khmer", "Mong", "Other/Missing")
    )
  )
# Area variable
data <- data %>%
  mutate(
    HH6 = factor(HH6, levels = c("URBAN", "RURAL"))
  )
# Geographic region variable
data <- data %>%
  mutate(
    HH7 = factor(
      HH7,
      levels = c(
        "RED RIVER DELTA",
        "NORTHERN MIDLANDS AND MOUNTAIN",
        "NORTH CENTRAL AND CENTRAL COASTAL",
        "CENTRAL HIGHLANDS",
        "SOUTH EAST",
        "MEKONG RIVER DELTA"
      )
    )
  )
# Early childhood education attendance
data <- data %>%
  mutate(
    UB6 = factor(UB6, levels = c("YES", "NO"))
  )
# BHYT
data <- data %>%
  mutate(
    cinsurance = factor(cinsurance, levels = c("With insurance", "Without insurance"))
  )
# Create the Recent_Illness variable
data <- data %>%
  mutate(
    recent_illness = ifelse(CA1 == "YES" | CA14 == "YES" | CA16 == "YES" | CA17 == "YES", 1, 0),
    recent_illness = factor(recent_illness, levels = c(0, 1))
  )
## Create parental co-residence variable
data <- data %>%
  mutate(
    # Apply the HL12/HL16 skip pattern before deriving parental co-residence.
    HL13 = ifelse(trimws(HL12) == "NO", "NO",
                  ifelse(trimws(HL12) == "DK", "DK", HL13)),
    
    HL17 = ifelse(trimws(HL16) == "NO", "NO",
                  ifelse(trimws(HL16) == "DK", "DK", HL17))
  )
data <- data %>%
  mutate(
    family = case_when(
      # STEP 1: First identify missing, DK, no response, or blank cases
      if_any(c(HL13, HL17), ~ trimws(.) %in% c("NO RESPONSE", "DK", "NO RESPONSE/DK", "")) ~ NA_real_,
      # STEP 2: Classify into four groups; trimws() avoids whitespace problems
      # Group 1: Living with both mother and father
      trimws(HL13) %in% "YES" & trimws(HL17) %in% "YES" ~ 1,
      # Group 2: Living with mother only
      trimws(HL13) %in% "YES" & trimws(HL17) %in% "NO" ~ 2,
      # Group 3: Living with father only
      trimws(HL13) %in% "NO" & trimws(HL17) %in% "YES" ~ 3,
      # Group 4: Living with neither parent
      trimws(HL13) %in% "NO" & trimws(HL17) %in% "NO" ~ 4,
      
      # Final safety rule: any remaining unexpected case is set to missing
      TRUE ~ NA_real_
    ),
    family = factor(family, levels = c(1, 2, 3, 4))
  )
## Household digital environment
data <- data %>%
  mutate(
    digital_enviroment = case_when(
      trimws(HC11) %in% c("", "9", "98", "DK", "NO RESPONSE") |
        trimws(HC12) %in% c("", "9", "98", "DK", "NO RESPONSE") |
        trimws(HC13) %in% c("", "9", "98", "DK", "NO RESPONSE") ~ NA_real_,
      
      # Logical TRUE/FALSE statements are summed directly to create a numeric count
      TRUE ~ (trimws(HC11) == "YES") +
        (trimws(HC12) == "YES") +
        (trimws(HC13) == "YES")
    ),
    digital_enviroment = factor(digital_enviroment, levels = c(3, 2, 1, 0))
  )
## Process HH51
data <- data %>%
  mutate(
    # Step 1: Convert to numeric; values greater than 2 are grouped as >2
    HH51_group = ifelse(as.numeric(HH51) > 2, ">2", HH51),
    
    # Step 2: Convert to factor and set 1 as the reference category
    HH51_group = factor(HH51_group, levels = c("1", "2", ">2"))
  )
## Process HH48
data <- data %>%
  mutate(
    HH48_num = as.numeric(as.character(HH48)),
    HH48_group = case_when(
      HH48_num >= 1 & HH48_num <= 4 ~ "1 - 4",
      HH48_num >= 5 & HH48_num <= 6 ~ "5 - 6",
      HH48_num > 6 ~ ">6",
      TRUE ~ NA_character_ # Handle missing or invalid cases if present
    ),
    HH48_group = factor(HH48_group, levels = c("1 - 4", "5 - 6", ">6"))
  )
## Process windex5
data <- data %>%
  mutate(windex5 = factor(windex5, levels = c("Richest", "Rich", "Middle", "Poor", "Poorest")))
## Build the descriptive sample and complete-case analytic sample
data_final <- data %>%
  mutate(
    CAGE_group = case_when(
      CAGE >= 24 & CAGE <= 29 ~ "24-29 months",
      CAGE >= 30 & CAGE <= 35 ~ "30-35 months",
      CAGE >= 36 & CAGE <= 41 ~ "36-41 months",
      CAGE >= 42 & CAGE <= 47 ~ "42-47 months",
      CAGE >= 48 & CAGE <= 59 ~ "48-59 months",
      TRUE ~ NA_character_
    ),
    melevel_group = case_when(
      melevel %in% c("Pre-primary or none", "Primary") ~ "Primary or below",
      melevel %in% c("Lower Secondary", "Upper Secondary") ~ "Secondary",
      melevel %in% c("Vocational High school", "University/College/Higher") ~ "Higher",
      TRUE ~ NA_character_
    ),
    digital_enviroment = case_when(
      trimws(HC11) %in% c("", "9", "98", "DK", "NO RESPONSE") |
        trimws(HC12) %in% c("", "9", "98", "DK", "NO RESPONSE") |
        trimws(HC13) %in% c("", "9", "98", "DK", "NO RESPONSE") ~ NA_real_,
      trimws(HC13) == "NO" ~ 0,
      trimws(HC13) == "YES" & trimws(HC12) == "YES" & trimws(HC11) == "YES" ~ 2,
      trimws(HC13) == "YES" & trimws(HC12) == "YES" & trimws(HC11) == "NO" ~ 1,
      TRUE ~ NA_real_
    ),
    family = case_when(
      if_any(c(HL13, HL17), ~ trimws(.) %in% c("NO RESPONSE", "DK", "NO RESPONSE/DK", "")) ~ NA_real_,
      trimws(HL13) == "YES" & trimws(HL17) == "YES" ~ 1,
      TRUE ~ 2
    ),
    ethnicity = ifelse(as.character(ethnicity) == "Kinh and Hoa", "Kinh and Hoa", "Other"),
    CAGE_group = factor(CAGE_group, levels = c("24-29 months", "30-35 months", "36-41 months", "42-47 months", "48-59 months")),
    mother_stim = factor(mother_stim, levels = c(0, 1)),
    father_stim = factor(father_stim, levels = c(0, 1)),
    toy = factor(toy, levels = c(0, 1)),
    book = factor(book, levels = c(0, 1)),
    family = factor(family, levels = c(2, 1)),
    digital_enviroment = factor(digital_enviroment, levels = c(0, 1, 2)),
    melevel_group = factor(melevel_group, levels = c("Primary or below", "Secondary", "Higher")),
    melevel_num = as.numeric(melevel_group),
    HH51_group = factor(HH51_group, levels = c("1", "2", ">2")),
    HH48_group = factor(HH48_group, levels = c("1 - 4", "5 - 6", ">6")),
    ethnicity = factor(ethnicity, levels = c("Kinh and Hoa", "Other")),
    HH6 = factor(HH6, levels = c("URBAN", "RURAL")),
    HH7 = factor(HH7, levels = c("NORTHERN MIDLANDS AND MOUNTAIN", "CENTRAL HIGHLANDS", "MEKONG RIVER DELTA", "NORTH CENTRAL AND CENTRAL COASTAL", "SOUTH EAST", "RED RIVER DELTA")),
    UB6 = factor(UB6, levels = c("NO", "YES")),
    cinsurance = factor(cinsurance, levels = c("Without insurance", "With insurance")),
    cdisability = factor(cdisability, levels = c("Has functional difficulty", "Has no functional difficulty")),
    recent_illness = factor(recent_illness, levels = c(0, 1)),
    windex5 = factor(windex5, levels = c("Poorest", "Poor", "Middle", "Rich", "Richest")),
    windex5_num = as.numeric(windex5)
  ) %>%
  filter(!is.na(ecdi)) %>%
  select(HH1, HL4, CAGE, CAGE_group, ethnicity, cinsurance, family,
         cdisability, melevel, melevel_group, melevel_num, recent_illness,
         HH48_group, HH51_group, windex5, windex5_num, digital_enviroment,
         HH6, HH7, UB6, book, toy, mother_stim, father_stim,
         stratum, chweight, ecdi)

data_model <- data_final %>%
  select(ecdi, windex5, windex5_num, CAGE_group, HL4, melevel_group, melevel_num,
         HH6, HH7, ethnicity, family, recent_illness, digital_enviroment,
         HH51_group, HH48_group, CAGE, book, toy, mother_stim, father_stim,
         HH1, stratum, chweight) %>%
  na.omit() %>%
  mutate(
    edu_secondary = ifelse(melevel_group == "Secondary", 1, 0),
    edu_higher = ifelse(melevel_group == "Higher", 1, 0),
    hh51_2child = ifelse(HH51_group == "2", 1, 0),
    hh51_gt2child = ifelse(HH51_group == ">2", 1, 0),
    digital_1 = ifelse(digital_enviroment == "1", 1, 0),
    digital_2 = ifelse(digital_enviroment == "2", 1, 0),
    reg_central_highlands = ifelse(HH7 == "CENTRAL HIGHLANDS", 1, 0),
    reg_mekong_delta = ifelse(HH7 == "MEKONG RIVER DELTA", 1, 0),
    reg_central_coast = ifelse(HH7 == "NORTH CENTRAL AND CENTRAL COASTAL", 1, 0),
    reg_south_east = ifelse(HH7 == "SOUTH EAST", 1, 0),
    reg_red_river_delta = ifelse(HH7 == "RED RIVER DELTA", 1, 0),
    hh48_5_6 = ifelse(HH48_group == "5 - 6", 1, 0),
    hh48_gt6 = ifelse(HH48_group == ">6", 1, 0)
  )

rownames(data_model) <- NULL

# ============================================================
# 3. SURVEY DESIGN OBJECTS
# ============================================================

mics6_design <- data_model %>%
  as_survey_design(ids = HH1,
                   strata = stratum,
                   weights = chweight,
                   nest = TRUE)

mics6_design_des <- data_final %>%
  as_survey_design(ids = HH1,
                   strata = stratum,
                   weights = chweight,
                   nest = TRUE)

# ============================================================
# 4. TABLE HELPER FUNCTIONS
# ============================================================

create_table_block <- function(design_data, vars) {
  t_desc <- design_data %>%
    tbl_svysummary(
      include = all_of(vars),
      statistic = list(all_categorical() ~ "{n_unweighted} ({p}%)"),
      type = list(all_dichotomous() ~ "categorical"),
      missing = "ifany",
      missing_text = "Missing",
      missing_stat = "{N_miss_unweighted}"
    ) %>%
    modify_header(stat_0 ~ "**N (weighted %)**")
  
  ci_df <- lapply(vars, function(v) {
    design_data %>%
      group_by(.data[[v]]) %>%
      summarise(
        prev = survey_mean(ecdi, vartype = "ci", proportion = TRUE, na.rm = TRUE)
      ) %>%
      filter(!is.na(.data[[v]])) %>%
      mutate(
        variable = v,
        label = as.character(.data[[v]]),
        ci_str = sprintf("(%.0f%%, %.0f%%)", prev_low * 100, prev_upp * 100)
      ) %>%
      select(variable, label, ci_str)
  }) %>%
    bind_rows()
  
  t_prev <- design_data %>%
    tbl_svysummary(
      by = ecdi,
      percent = "row",
      include = all_of(vars),
      statistic = list(all_categorical() ~ "{p}%"),
      type = list(all_dichotomous() ~ "categorical"),
      missing = "no"
    ) %>%
    add_p(pvalue_fun = ~style_pvalue(.x, digits = 3)) %>%
    modify_table_body(~ .x %>%
                        left_join(ci_df, by = c("variable", "label")) %>%
                        mutate(stat_2 = ifelse(!is.na(ci_str), paste(stat_2, ci_str), stat_2))
    ) %>%
    modify_header(stat_2 ~ "**Prevalence (95% CI)**") %>%
    modify_column_hide(columns = stat_1)
  
  tbl_merge(tbls = list(t_desc, t_prev), tab_spanner = FALSE, quiet = TRUE)
}

vars_child <- c("CAGE_group", "HL4", "ethnicity", "cinsurance", "cdisability", "recent_illness")
tbl_child <- create_table_block(mics6_design_des, vars_child)

tbl_ub6 <- subset(mics6_design_des, CAGE >= 36) %>%
  create_table_block(c("UB6"))

vars_parent <- c("melevel_group", "family", "mother_stim", "father_stim")
tbl_parent <- create_table_block(mics6_design_des, vars_parent)

vars_hh <- c("windex5", "HH48_group", "HH51_group", "digital_enviroment", "book", "toy")
tbl_hh <- create_table_block(mics6_design_des, vars_hh)

vars_comm <- c("HH6", "HH7")
tbl_comm <- create_table_block(mics6_design_des, vars_comm)

# ============================================================
# 5. MAIN RESULT TABLES
# ============================================================

# Main Table 1: Baseline characteristics
final_table_1 <- tbl_stack(
  tbls = list(tbl_child, tbl_ub6, tbl_parent, tbl_hh, tbl_comm),
  group_header = c(
    "Child characteristics",
    "Early childhood education (>=36 months)",
    "Parental / caregiver characteristics",
    "Household characteristics",
    "Community characteristics"
  ),
  quiet = TRUE
) %>%
  bold_labels() %>%
  italicize_levels()

# ============================================================
# 6. ANALYTICAL TABLES
# ============================================================
advanced_table_warning_option <- getOption("warn")
options(warn = -1)

# 6.1 Supplementary measurement and inequality tables

# Supplementary Table S1: ECDI2030 measurement domains and age-specific thresholds
appendix_table_1_thresholds <- data.frame(
  Domain = c(
    "Learning",
    "Psychosocial well-being",
    "Health",
    "Overall ECDI2030 age-specific threshold"
  ),
  Items = c(
    "EC21-EC28",
    "EC29-EC36, EC39-EC40",
    "EC37-EC38",
    "Total number of age-appropriate milestones achieved"
  ),
  Criterion = c(
    "Domain items coded as achieved/not achieved according to ECDI2030 scoring guidance.",
    "Domain items coded as achieved/not achieved; selected behavior items are reverse coded.",
    "Health-related items coded as achieved/not achieved.",
    "On-track status is assigned using the age-specific total-score cut-off."
  ),
  `Age-specific threshold` = c(
    "Included in total ECDI2030 score",
    "Included in total ECDI2030 score",
    "Included in total ECDI2030 score",
    "24-29 months: >=7; 30-35 months: >=9; 36-41 months: >=11; 42-47 months: >=13; 48-59 months: >=15"
  ),
  check.names = FALSE
)

format_inequality_result <- function(obj) {
  est <- as.numeric(obj$estimate)
  se <- as.numeric(obj$se)
  null_val <- ifelse(!is.null(obj$measure) && obj$measure %in% c("r", "rii"), 1, 0)
  p_val <- 2 * pnorm(-abs((est - null_val) / se))
  sprintf(
    "%.3f (%.3f, %.3f); %s",
    est,
    as.numeric(obj$lowerci),
    as.numeric(obj$upperci),
    ifelse(p_val < 0.001, "p<0.001", sprintf("p=%.3f", p_val))
  )
}

estimate_inequality_indices <- function(data_input, design_input, strat_var) {
  regions <- c("All Nation", levels(data_input$HH7))
  results <- data.frame()
  
  for (reg in regions) {
    sub_data <- if (reg == "All Nation") data_input else subset(data_input, HH7 == reg)
    sub_design <- if (reg == "All Nation") design_input else subset(design_input, HH7 == reg)
    
    mean_ecdi <- svymean(~ecdi, sub_design, na.rm = TRUE)
    subgroup_est <- svyby(as.formula("~ecdi"), as.formula(paste0("~", strat_var)), sub_design, svymean, na.rm = TRUE)
    subgroup_est$pop <- as.numeric(svytable(as.formula(paste0("~", strat_var)), sub_design))
    subgroup_order <- subgroup_est[[strat_var]]
    
    temp <- data.frame(
      Region = reg,
      `ECDI prevalence (%)` = round(as.numeric(mean_ecdi) * 100, 3),
      `Absolute difference (D)` = format_inequality_result(d(
        est = subgroup_est$ecdi,
        se = subgroup_est$se,
        favourable_indicator = 1,
        ordered_dimension = 1,
        subgroup_order = subgroup_order
      )),
      `Relative ratio (R)` = format_inequality_result(r(
        est = subgroup_est$ecdi,
        se = subgroup_est$se,
        favourable_indicator = 1,
        ordered_dimension = 1,
        subgroup_order = subgroup_order
      )),
      `Slope index of inequality (SII)` = format_inequality_result(sii(
        est = sub_data$ecdi,
        weight = sub_data$chweight,
        psu = sub_data$HH1,
        strata = sub_data$stratum,
        subgroup_order = sub_data[[strat_var]],
        favourable_indicator = 1,
        force = TRUE
      )),
      `Erreygers Concentration Index (ECI)` = format_inequality_result(suppressMessages(suppressWarnings(rci(
        est = sub_data$ecdi,
        weight = sub_data$chweight,
        psu = sub_data$HH1,
        strata = sub_data$stratum,
        subgroup_order = sub_data[[strat_var]],
        favourable_indicator = 1,
        force = TRUE,
        method = "erreygers",
        lmin = 0,
        lmax = 1
      )))),
      `Population attributable risk (PAR)` = format_inequality_result(parisk(
        est = subgroup_est$ecdi,
        pop = subgroup_est$pop,
        favourable_indicator = 1,
        ordered_dimension = 1,
        subgroup_order = subgroup_order,
        scaleval = 1,
        force = TRUE
      )),
      check.names = FALSE
    )
    
    results <- rbind(results, temp)
  }
  
  results
}

wealth_inequality_indices <- suppressWarnings(suppressMessages(
  estimate_inequality_indices(data_model, mics6_design, "windex5_num")
)) %>%
  mutate(Stratifier = "Household wealth")

education_inequality_indices <- suppressWarnings(suppressMessages(
  estimate_inequality_indices(data_model, mics6_design, "melevel_num")
)) %>%
  mutate(Stratifier = "Mother's education")

# Supplementary Table S2: inequality indices by wealth and mother's education
appendix_table_2_inequality <- bind_rows(
  wealth_inequality_indices,
  education_inequality_indices
) %>%
  select(Stratifier, everything())

# 6.2 Main analytical result tables

# Main Table 2: Multivariable logistic regression, children aged 24-59 months
model_full_visual <- svyglm(
  design = mics6_design,
  formula = ecdi ~ CAGE + HL4 + ethnicity + family + recent_illness +
    melevel_group + mother_stim + father_stim +
    windex5 + HH48_group + HH51_group + toy + book +
    digital_enviroment + HH6 + HH7,
  family = quasibinomial()
)

table_univariate_full <- mics6_design %>%
  tbl_uvregression(
    method = svyglm,
    y = ecdi,
    method.args = list(family = quasibinomial()),
    exponentiate = TRUE,
    hide_n = TRUE,
    include = c(CAGE, HL4, ethnicity, family, recent_illness, melevel_group,
                mother_stim, father_stim, windex5, HH48_group, HH51_group,
                toy, book, digital_enviroment, HH6, HH7),
    pvalue_fun = ~style_pvalue(.x, digits = 3)
  ) %>%
  modify_table_styling(
    columns = estimate,
    rows = reference_row %in% TRUE,
    missing_symbol = "1.00 (Ref)"
  ) %>%
  modify_header(
    estimate ~ "**Crude OR (95% CI)**",
    p.value ~ "**p-value**"
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

table_visual_full <- tbl_regression(
  model_full_visual,
  exponentiate = TRUE,
  pvalue_fun = ~style_pvalue(.x, digits = 3)
) %>%
  modify_table_styling(
    columns = estimate,
    rows = reference_row %in% TRUE,
    missing_symbol = "1.00 (Ref)"
  ) %>%
  modify_header(
    estimate ~ "**Adjusted OR (95% CI)**",
    p.value ~ "**p-value**"
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

main_table_logistic_24_59 <- tbl_merge(
  tbls = list(table_univariate_full, table_visual_full),
  tab_spanner = c("**Univariate analysis**", "**Multivariable analysis**"),
  quiet = TRUE
) %>%
  modify_header(label = "**Characteristics**") %>%
  modify_caption(sprintf("**Multivariable logistic regression for ECDI2030 status among children aged 24-59 months in Vietnam (N = %s)**", format(nrow(data_model), big.mark = ",")))

# Main Table 3: Wealth-related decomposition, children aged 24-59 months
model_full_decomposition <- svyglm(
  design = mics6_design,
  formula = ecdi ~ CAGE + HL4 + ethnicity + family + recent_illness +
    edu_secondary + edu_higher + hh51_gt2child + hh51_2child +
    digital_1 + digital_2 + reg_red_river_delta + reg_central_coast +
    reg_central_highlands + reg_south_east + reg_mekong_delta + hh48_gt6 + hh48_5_6 +
    book + toy + mother_stim + father_stim + HH6,
  family = quasibinomial()
)

create_decomposition_table <- function(model, data_input, ranker) {
  # Average marginal effects must use the same survey weights as the model.
  model_weights <- as.numeric(weights(model))
  mfx <- avg_slopes(model, wts = model_weights)
  ame_df <- as.data.frame(mfx)
  x_matrix <- model.matrix(model)[, -1, drop = FALSE]
  
  ame_df$matrix_name <- ifelse(
    paste0(ame_df$term, sub(" - .*", "", ame_df$contrast)) %in% colnames(x_matrix),
    paste0(ame_df$term, sub(" - .*", "", ame_df$contrast)),
    ame_df$term
  )
  
  ame_vector <- ame_df$estimate[match(colnames(x_matrix), ame_df$matrix_name)]
  names(ame_vector) <- colnames(x_matrix)

  if (anyNA(ame_vector)) {
    stop(
      "Average marginal effects could not be matched to all model-matrix columns: ",
      paste(names(ame_vector)[is.na(ame_vector)], collapse = ", "),
      call. = FALSE
    )
  }
  
  # rineq decomposes the ordinary concentration index as
  #   C_y = sum[(AME_k * mean(x_k) / mean(y)) * C_k] + residual.
  # For a binary outcome bounded by 0 and 1, multiplying every component by
  # 4 * mean(y) converts this identity to the Erreygers scale:
  #   E_y = 4 * sum[AME_k * mean(x_k) * C_k] + 4 * GCI_residual.
  # Therefore, determinants use ordinary concentration indices (citype = "CI").
  # Applying citype = "CIc" to every determinant would incorrectly add a
  # determinant-specific correction factor.
  decomp_result <- decomposition(
    outcome = data_input$ecdi,
    betas = ame_vector,
    mm = x_matrix,
    ranker = ranker,
    wt = data_input$chweight,
    correction = FALSE,
    citype = "CI"
  )
  
  outcome_mean <- weighted.mean(
    data_input$ecdi,
    w = data_input$chweight,
    na.rm = TRUE
  )
  erreygers_factor <- 4 * outcome_mean
  overall_eci <- erreygers_factor * decomp_result$overall_ci$concentration_index

  ci_limits <- t(decomp_result$confints)
  determinant_names <- names(decomp_result$ci_contribution)

  result <- data.frame(
    Determinant = determinant_names,
    `Average marginal effect` = c(NA_real_, decomp_result$betas),
    `Weighted mean` = c(NA_real_, decomp_result$averages),
    Elasticity = c(0, decomp_result$elasticities),
    `Concentration Index` = c(NA_real_, decomp_result$partial_cis),
    `CI lower 95%` = c(NA_real_, ci_limits[, 1]),
    `CI upper 95%` = c(NA_real_, ci_limits[, 2]),
    `Contribution (Abs)` = erreygers_factor * decomp_result$ci_contribution,
    `Contribution (%)` = 100 *
      (erreygers_factor * decomp_result$ci_contribution) / overall_eci,
    check.names = FALSE,
    row.names = NULL
  )

  result <- result %>%
    mutate(
      across(
        c(`Average marginal effect`, `Weighted mean`, Elasticity,
          `Concentration Index`, `CI lower 95%`, `CI upper 95%`),
        ~round(.x, 4)
      ),
      `Contribution (Abs)` = round(`Contribution (Abs)`, 6),
      `Contribution (%)` = round(`Contribution (%)`, 4)
    )

  attr(result, "outcome_mean") <- outcome_mean
  attr(result, "overall_eci") <- overall_eci
  result
}

main_table_decomposition_24_59 <- create_decomposition_table(
  model_full_decomposition,
  data_model,
  data_model$windex5_num
)

write.csv(
  main_table_decomposition_24_59,
  "outputs/Table_3_Erreygers_decomposition_24_59.csv",
  row.names = FALSE
)

# 6.3 Supplementary analytical tables

# Supplementary Tables S3 and S4: children aged 36-59 months
data_model_ub6 <- data_final %>%
  filter(CAGE >= 36 & CAGE <= 59) %>%
  select(ecdi, windex5, windex5_num, CAGE, HL4, melevel_group, melevel_num,
         HH6, HH7, ethnicity, family, recent_illness, digital_enviroment,
         HH51_group, HH48_group, book, toy, mother_stim, father_stim,
         UB6, HH1, stratum, chweight) %>%
  na.omit() %>%
  droplevels() %>%
  mutate(
    edu_secondary = ifelse(melevel_group == "Secondary", 1, 0),
    edu_higher = ifelse(melevel_group == "Higher", 1, 0),
    hh51_2child = ifelse(HH51_group == "2", 1, 0),
    hh51_gt2child = ifelse(HH51_group == ">2", 1, 0),
    digital_1 = ifelse(digital_enviroment == "1", 1, 0),
    digital_2 = ifelse(digital_enviroment == "2", 1, 0),
    reg_red_river_delta = ifelse(HH7 == "RED RIVER DELTA", 1, 0),
    reg_central_coast = ifelse(HH7 == "NORTH CENTRAL AND CENTRAL COASTAL", 1, 0),
    reg_central_highlands = ifelse(HH7 == "CENTRAL HIGHLANDS", 1, 0),
    reg_south_east = ifelse(HH7 == "SOUTH EAST", 1, 0),
    reg_mekong_delta = ifelse(HH7 == "MEKONG RIVER DELTA", 1, 0),
    ub6_yes = ifelse(UB6 == "YES", 1, 0),
    hh48_5_6 = ifelse(HH48_group == "5 - 6", 1, 0),
    hh48_gt6 = ifelse(HH48_group == ">6", 1, 0)
  )

rownames(data_model_ub6) <- NULL

# Compare observed sample sizes with those currently reported in the manuscript.
# A mismatch is reported, but observations are never removed merely to match text.
observed_sample_n <- c(
  descriptive = nrow(data_final),
  analysis = nrow(data_model),
  age_36_59 = nrow(data_model_ub6)
)
expected_sample_n <- c(descriptive = 2691L, analysis = 2667L, age_36_59 = 1859L)
if (!identical(observed_sample_n, expected_sample_n)) {
  warning(
    "Sample sizes differ from the manuscript: observed ",
    paste(names(observed_sample_n), observed_sample_n, sep = "=", collapse = ", "),
    "; expected ",
    paste(names(expected_sample_n), expected_sample_n, sep = "=", collapse = ", "),
    call. = FALSE
  )
}

mics6_design_ub6 <- data_model_ub6 %>%
  as_survey_design(ids = HH1, strata = stratum, weights = chweight, nest = TRUE)

model_ub6_visual <- suppressWarnings(svyglm(
  design = mics6_design_ub6,
  formula = ecdi ~ CAGE + HL4 + ethnicity + family + recent_illness +
    melevel_group + mother_stim + father_stim +
    windex5 + HH48_group + HH51_group + toy + book +
    digital_enviroment + HH6 + HH7 + UB6,
  family = quasibinomial()
))

table_univariate_ub6 <- mics6_design_ub6 %>%
  tbl_uvregression(
    method = svyglm,
    y = ecdi,
    method.args = list(family = quasibinomial()),
    exponentiate = TRUE,
    hide_n = TRUE,
    include = c(CAGE, HL4, ethnicity, family, recent_illness, melevel_group,
                mother_stim, father_stim, windex5, HH48_group, HH51_group,
                toy, book, digital_enviroment, HH6, HH7, UB6),
    pvalue_fun = ~style_pvalue(.x, digits = 3)
  ) %>%
  modify_table_styling(
    columns = estimate,
    rows = reference_row %in% TRUE,
    missing_symbol = "1.00 (Ref)"
  ) %>%
  modify_header(
    estimate ~ "**Crude OR (95% CI)**",
    p.value ~ "**p-value**"
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

table_visual_ub6 <- tbl_regression(
  model_ub6_visual,
  exponentiate = TRUE,
  pvalue_fun = ~style_pvalue(.x, digits = 3)
) %>%
  modify_table_styling(
    columns = estimate,
    rows = reference_row %in% TRUE,
    missing_symbol = "1.00 (Ref)"
  ) %>%
  modify_header(
    estimate ~ "**Adjusted OR (95% CI)**",
    p.value ~ "**p-value**"
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

appendix_table_3_logistic_36_59 <- tbl_merge(
  tbls = list(table_univariate_ub6, table_visual_ub6),
  tab_spanner = c("**Univariate analysis**", "**Multivariable analysis**"),
  quiet = TRUE
) %>%
  modify_header(label = "**Characteristics**") %>%
  modify_caption(sprintf("**Multivariable logistic regression for ECDI2030 status among children aged 36-59 months in Vietnam (N = %s)**", format(nrow(data_model_ub6), big.mark = ",")))

formula_ub6 <- ecdi ~ CAGE + HL4 + ethnicity + family + recent_illness +
  edu_secondary + edu_higher + hh51_2child + hh51_gt2child +
  digital_1 + digital_2 + reg_red_river_delta + reg_central_coast +
  hh48_gt6 + hh48_5_6 + reg_central_highlands + reg_south_east +
  reg_mekong_delta + book + toy + mother_stim + father_stim + HH6 + ub6_yes

model_start_ub6 <- suppressWarnings(svyglm(design = mics6_design_ub6, formula = formula_ub6, family = quasibinomial()))
model_full_ub6 <- suppressWarnings(svyglm(
  design = mics6_design_ub6,
  formula = formula_ub6,
  family = quasibinomial(),
  start = coef(model_start_ub6)
))

appendix_table_4_decomposition_36_59 <- create_decomposition_table(
  model_full_ub6,
  data_model_ub6,
  data_model_ub6$windex5_num
)

write.csv(
  appendix_table_4_decomposition_36_59,
  "outputs/Table_S4_Erreygers_decomposition_36_59.csv",
  row.names = FALSE
)

options(warn = advanced_table_warning_option)

# ============================================================
# 7. FIGURES
# ============================================================

# ============================================================
# Main Figure 1: On-track ECDI2030 prevalence by MICS6 region
# ============================================================
normalize_admin_name <- function(x) {
  x <- enc2utf8(x)
  if (requireNamespace("stringi", quietly = TRUE)) {
    x <- stringi::stri_trans_general(x, "Latin-ASCII")
  } else {
    x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  }
  toupper(gsub("[^A-Z]", "", x))
}

assign_mics_region <- function(province_name) {
  province_name <- enc2utf8(province_name)
  
  province_region_lookup <- c(
    "Hà Giang" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Cao Bằng" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Bắc Kạn" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Tuyên Quang" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Lào Cai" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Yên Bái" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Thái Nguyên" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Lạng Sơn" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Bắc Giang" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Phú Thọ" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Điện Biên" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Lai Châu" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Sơn La" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Hoà Bình" = "NORTHERN MIDLANDS AND MOUNTAIN",
    "Hà Nội" = "RED RIVER DELTA",
    "Vĩnh Phúc" = "RED RIVER DELTA",
    "Bắc Ninh" = "RED RIVER DELTA",
    "Quảng Ninh" = "RED RIVER DELTA",
    "Hải Dương" = "RED RIVER DELTA",
    "Hải Phòng" = "RED RIVER DELTA",
    "Hưng Yên" = "RED RIVER DELTA",
    "Thái Bình" = "RED RIVER DELTA",
    "Hà Nam" = "RED RIVER DELTA",
    "Nam Định" = "RED RIVER DELTA",
    "Ninh Bình" = "RED RIVER DELTA",
    "Thanh Hóa" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Nghệ An" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Hà Tĩnh" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Quảng Bình" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Quảng Trị" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Thừa Thiên Huế" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Đà Nẵng" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Quảng Nam" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Quảng Ngãi" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Bình Định" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Phú Yên" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Khánh Hòa" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Ninh Thuận" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Bình Thuận" = "NORTH CENTRAL AND CENTRAL COASTAL",
    "Kon Tum" = "CENTRAL HIGHLANDS",
    "Gia Lai" = "CENTRAL HIGHLANDS",
    "Đắk Lắk" = "CENTRAL HIGHLANDS",
    "Đắk Nông" = "CENTRAL HIGHLANDS",
    "Lâm Đồng" = "CENTRAL HIGHLANDS",
    "Bình Phước" = "SOUTH EAST",
    "Tây Ninh" = "SOUTH EAST",
    "Bình Dương" = "SOUTH EAST",
    "Đồng Nai" = "SOUTH EAST",
    "Bà Rịa - Vũng Tàu" = "SOUTH EAST",
    "Hồ Chí Minh" = "SOUTH EAST",
    "Long An" = "MEKONG RIVER DELTA",
    "Tiền Giang" = "MEKONG RIVER DELTA",
    "Bến Tre" = "MEKONG RIVER DELTA",
    "Trà Vinh" = "MEKONG RIVER DELTA",
    "Vĩnh Long" = "MEKONG RIVER DELTA",
    "Đồng Tháp" = "MEKONG RIVER DELTA",
    "An Giang" = "MEKONG RIVER DELTA",
    "Kiên Giang" = "MEKONG RIVER DELTA",
    "Cần Thơ" = "MEKONG RIVER DELTA",
    "Hậu Giang" = "MEKONG RIVER DELTA",
    "Sóc Trăng" = "MEKONG RIVER DELTA",
    "Bạc Liêu" = "MEKONG RIVER DELTA",
    "Cà Mau" = "MEKONG RIVER DELTA"
  )
  
  unname(province_region_lookup[province_name])
}

prev_region_svy <- suppressWarnings(survey::svyby(
  ~ecdi, ~HH7, mics6_design,
  survey::svymean,
  na.rm = TRUE,
  vartype = c("se", "ci"),
  drop.empty.groups = FALSE
))

prev_region <- as.data.frame(prev_region_svy) %>%
  transmute(
    region_mics = as.character(HH7),
    prevalence = ecdi * 100,
    se = se * 100,
    ci_low = ci_l * 100,
    ci_high = ci_u * 100
  ) %>%
  arrange(region_mics)

prevalence_check <- data_model %>%
  group_by(HH7) %>%
  summarise(
    n_unweighted = n(),
    prevalence_weighted_manual = weighted.mean(ecdi, chweight, na.rm = TRUE) * 100,
    prevalence_unweighted = mean(ecdi, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  rename(region_mics = HH7) %>%
  mutate(region_mics = as.character(region_mics)) %>%
  arrange(region_mics)

prev_region_check <- prev_region %>%
  left_join(prevalence_check, by = "region_mics") %>%
  mutate(absolute_difference = abs(prevalence - prevalence_weighted_manual))

stopifnot(max(prev_region_check$absolute_difference, na.rm = TRUE) < 1e-8)
write.csv(prev_region_check, "outputs/Figure_1_ECDI2030_region_estimates.csv", row.names = FALSE)

shape_dir <- "."
province_shp <- file.path(shape_dir, "gadm41_VNM_1_UTM.shp")
admin3_shp <- file.path(shape_dir, "gadm41_VNM_3_UTM.shp")

missing_shape_files <- province_shp[!file.exists(province_shp)]
if (length(missing_shape_files) > 0) {
  stop(
    "The province shapefile is required for Figure 1 but was not found: ",
    paste(missing_shape_files, collapse = ", "),
    call. = FALSE
  )
}

vn <- sf::st_read(province_shp, quiet = TRUE) %>%
  sf::st_as_sf() %>%
  sf::st_make_valid() %>%
  sf::st_transform(4326) %>%
  mutate(region_mics = assign_mics_region(NAME_1))

unmatched_provinces <- vn %>%
  filter(is.na(region_mics)) %>%
  distinct(NAME_1)

if (nrow(unmatched_provinces) > 0) {
  stop(
    "Some provinces could not be matched to MICS6 regions: ",
    paste(unmatched_provinces$NAME_1, collapse = ", "),
    call. = FALSE
  )
}

if (file.exists(admin3_shp)) {
  nearshore_islands <- sf::st_read(admin3_shp, quiet = TRUE) %>%
    sf::st_as_sf() %>%
    sf::st_make_valid() %>%
    sf::st_transform(4326) %>%
    filter(ENGTYPE_3 == "Island") %>%
    select(NAME_1, NAME_2, NAME_3, VARNAME_3, TYPE_3, ENGTYPE_3, geometry)
  
  nearshore_island_points <- if (nrow(nearshore_islands) > 0) {
    suppressWarnings(
      nearshore_islands %>%
        sf::st_transform(3857) %>%
        sf::st_point_on_surface() %>%
        sf::st_transform(4326)
    )
  } else {
    nearshore_islands
  }
} else {
  message("Admin-3 island shapefile not found. Using fallback point coordinates for Bach Long Vi and Con Co.")
  nearshore_islands <- NULL
  nearshore_island_points <- data.frame(
    NAME_1 = c("Hai Phong", "Quang Tri"),
    NAME_2 = c("Bach Long Vi", "Con Co"),
    NAME_3 = c("Bach Long Vi", "Con Co"),
    lon = c(107.73, 107.34),
    lat = c(20.13, 17.16)
  ) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
}

vn_province_map <- vn %>%
  left_join(prev_region, by = "region_mics")

vn_region_map <- vn_province_map %>%
  group_by(region_mics) %>%
  summarise(
    prevalence = first(prevalence),
    geometry = sf::st_union(geometry),
    .groups = "drop"
  )

region_label_points <- vn_region_map %>%
  sf::st_transform(3857) %>%
  mutate(
    prevalence_label = sprintf("%.1f%%", prevalence),
    geometry = sf::st_point_on_surface(geometry)
  ) %>%
  sf::st_transform(4326)

region_labels <- bind_cols(
  sf::st_drop_geometry(region_label_points),
  as.data.frame(sf::st_coordinates(region_label_points))
)

set.seed(123)
archipelago_points <- bind_rows(
  data.frame(group = "Hoang Sa", lon = rnorm(18, 112.3, 0.35), lat = rnorm(18, 16.5, 0.30)),
  data.frame(group = "Truong Sa", lon = rnorm(28, 114.2, 0.85), lat = rnorm(28, 10.5, 0.65))
)

archipelago_labels <- data.frame(
  label = c("Hoang Sa", "Truong Sa"),
  lon = c(112.55, 114.75),
  lat = c(16.4, 10.3)
)

map_prevalence_ecdi_final <- ggplot() +
  geom_sf(data = vn_province_map, aes(fill = prevalence), color = "white", linewidth = 0.12) +
  geom_sf(data = vn_region_map, inherit.aes = FALSE, fill = NA, color = "white", linewidth = 0.38) +
  {if (!is.null(nearshore_islands)) geom_sf(data = nearshore_islands, inherit.aes = FALSE, fill = "grey35", color = "white", linewidth = 0.15)} +
  geom_sf(data = nearshore_island_points, inherit.aes = FALSE, shape = 21, size = 1.1, stroke = 0.2, fill = "grey35", color = "white") +
  geom_label(
    data = region_labels,
    aes(x = X, y = Y, label = prevalence_label),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold",
    color = "#1E2C32",
    fill = "white",
    alpha = 0.86,
    label.padding = unit(0.14, "lines"),
    label.r = unit(0.06, "lines"),
    linewidth = 0
  ) +
  geom_point(data = archipelago_points, aes(x = lon, y = lat), inherit.aes = FALSE, color = "grey35", size = 0.45, alpha = 0.85) +
  geom_text(data = archipelago_labels, aes(x = lon, y = lat, label = label), inherit.aes = FALSE, size = 3, fontface = "bold", hjust = 0, color = "grey25") +
  annotate("text", x = 112.7, y = 22.4, label = "VIET NAM", size = 8, fontface = "bold", color = "#2B6C8A") +
  annotate("text", x = 112.4, y = 13.1, label = "East Sea", size = 4, angle = 310, fontface = "bold", color = "grey45") +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    limits = c(55, 90),
    breaks = c(60, 70, 80, 90),
    name = "ECDI2030 (%)",
    na.value = "#F0F0F0"
  ) +
  coord_sf(xlim = c(101.8, 117), ylim = c(7.8, 24), expand = FALSE) +
  labs(
    title = "Prevalence of on-track ECDI2030",
    subtitle = "Children aged 24-59 months, Viet Nam MICS6"
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey35"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.4, "cm"),
    plot.background = element_rect(fill = "#F8FBFF", color = NA),
    panel.background = element_rect(fill = "#F8FBFF", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"

ggsave(
  "outputs/Figure_1_On_track_ECDI2030_publication.png",
  map_prevalence_ecdi_final,
  width = 7,
  height = 8.5,
  dpi = 1000,
  device = png_device
)

ggsave(
  "outputs/Figure_1_On_track_ECDI2030_publication.pdf",
  map_prevalence_ecdi_final,
  width = 7,
  height = 8.5,
  device = cairo_pdf
)

# ============================================================
# Main Figure 2: prevalence patterns, SII, and ECI by wealth and education
# ============================================================

extract_estimate_ci <- function(x) {
  matched <- regexec("([-0-9.]+) \\(([-0-9.]+), ([-0-9.]+)\\)", x)
  parts <- regmatches(x, matched)
  data.frame(
    estimate = as.numeric(vapply(parts, function(z) if (length(z) >= 4) z[2] else NA_character_, character(1))),
    lower = as.numeric(vapply(parts, function(z) if (length(z) >= 4) z[3] else NA_character_, character(1))),
    upper = as.numeric(vapply(parts, function(z) if (length(z) >= 4) z[4] else NA_character_, character(1)))
  )
}

format_region_label <- function(x) {
  recode(
    x,
    "National" = "National level",
    "All Nation" = "National level",
    "NORTHERN MIDLANDS AND MOUNTAIN" = "Northern Midlands & Mountain",
    "NORTH CENTRAL AND CENTRAL COASTAL" = "North Central & Central Coastal",
    "CENTRAL HIGHLANDS" = "Central Highlands",
    "MEKONG RIVER DELTA" = "Mekong River Delta",
    "RED RIVER DELTA" = "Red River Delta",
    "SOUTH EAST" = "South East",
    .default = x
  )
}

format_region_factor <- function(x) {
  factor(
    format_region_label(x),
    levels = rev(c(
      "National level", "Red River Delta", "Northern Midlands & Mountain",
      "North Central & Central Coastal", "Central Highlands",
      "South East", "Mekong River Delta"
    ))
  )
}

prepare_inequality_plot <- function(data, measure_col) {
  data %>%
    transmute(Region, estimate_text = .data[[measure_col]]) %>%
    bind_cols(extract_estimate_ci(.$estimate_text)) %>%
    mutate(Region = format_region_factor(Region))
}

plot_forest_sii <- function(data, title_text) {
  ggplot(data, aes(y = Region, x = estimate)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.55) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18, linewidth = 0.65) +
    geom_point(shape = 15, size = 2.4) +
    scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2)) +
    theme_classic(base_family = "sans") +
    labs(title = title_text, x = "Slope Index of Inequality (SII)", y = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 8.5),
      axis.text.x = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, size = 8.5, face = "bold")
    )
}

plot_forest_eci <- function(data, title_text) {
  ggplot(data, aes(y = Region, x = estimate)) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.55) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18, linewidth = 0.65, color = "darkblue") +
    geom_point(shape = 17, size = 2.4, color = "darkblue") +
    scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2)) +
    theme_classic(base_family = "sans") +
    labs(title = title_text, x = "Erreygers Concentration Index (ECI)", y = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 8.5),
      axis.text.x = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, size = 8.5, face = "bold")
    )
}

plot_data_wealth <- suppressWarnings(bind_rows(
  survey::svyby(~ecdi, ~windex5, mics6_design, survey::svymean, na.rm = TRUE) %>%
    mutate(Region = "National"),
  survey::svyby(~ecdi, ~HH7 + windex5, mics6_design, survey::svymean, na.rm = TRUE) %>%
    rename(Region = HH7)
)) %>%
  mutate(Region = format_region_factor(Region))

plot_data_melevel <- suppressWarnings(bind_rows(
  survey::svyby(~ecdi, ~melevel_group, mics6_design, survey::svymean, na.rm = TRUE) %>%
    mutate(Region = "National"),
  survey::svyby(~ecdi, ~HH7 + melevel_group, mics6_design, survey::svymean, na.rm = TRUE) %>%
    rename(Region = HH7)
)) %>%
  mutate(
    Region = format_region_factor(Region),
    melevel_group = factor(melevel_group, levels = c("Primary or below", "Secondary", "Higher"))
  )

wealth_sii_data <- prepare_inequality_plot(wealth_inequality_indices, "Slope index of inequality (SII)")
wealth_eci_data <- prepare_inequality_plot(wealth_inequality_indices, "Erreygers Concentration Index (ECI)")
education_sii_data <- prepare_inequality_plot(education_inequality_indices, "Slope index of inequality (SII)")
education_eci_data <- prepare_inequality_plot(education_inequality_indices, "Erreygers Concentration Index (ECI)")

p_equi_wealth <- ggplot(plot_data_wealth, aes(x = ecdi * 100, y = Region, color = windex5)) +
  geom_line(aes(group = Region), color = "grey70", linewidth = 0.75) +
  geom_point(size = 2.4, stroke = 0.35) +
  scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, 10)) +
  scale_color_manual(
    name = "Wealth Quintile:",
    values = c("red", "darkorange", "forestgreen", "dodgerblue", "navy"),
    labels = c("Q1 (Poorest)", "Q2", "Q3", "Q4", "Q5 (Richest)")
  ) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5)) +
  theme_classic(base_family = "sans") +
  labs(title = "A. Household socioeconomic status", x = "Prevalence (%)", y = NULL) +
  theme(
    axis.text.y = element_text(size = 8.5),
    axis.title.x = element_text(size = 8.5),
    axis.text.x = element_text(size = 7.5),
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold")
  )

p_equi_edu <- ggplot(plot_data_melevel, aes(x = ecdi * 100, y = Region, color = melevel_group)) +
  geom_line(aes(group = Region), color = "grey70", linewidth = 0.75) +
  geom_point(size = 2.4, stroke = 0.35) +
  scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, 10)) +
  scale_color_manual(
    name = "Mother's Education:",
    values = c("deeppink", "slategray", "darkcyan"),
    labels = c("Primary or below", "Secondary", "Higher")
  ) +
  guides(color = guide_legend(nrow = 1, title.position = "left", title.vjust = 0.5)) +
  theme_classic(base_family = "sans") +
  labs(title = "D. Maternal educational status", x = "Prevalence (%)", y = NULL) +
  theme(
    axis.text.y = element_text(size = 8.5),
    axis.title.x = element_text(size = 8.5),
    axis.text.x = element_text(size = 7.5),
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold")
  )

figure_2_inequality <- (
  (p_equi_wealth +
     plot_forest_sii(wealth_sii_data, "B. Wealth-related absolute\ninequality (SII)") +
     plot_forest_eci(wealth_eci_data, "C. Wealth-related socioeconomic\nconcentration (ECI)") +
     plot_layout(widths = c(3, 1.15, 1.15))) /
    (p_equi_edu +
       plot_forest_sii(education_sii_data, "E. Education-related absolute\ninequality (SII)") +
       plot_forest_eci(education_eci_data, "F. Education-related socioeconomic\nconcentration (ECI)") +
       plot_layout(widths = c(3, 1.15, 1.15)))
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.5),
    legend.margin = margin(t = 0, b = 0)
  )

ggsave(
  "outputs/Figure_2_Socioeconomic_inequalities_ECDI2030.png",
  figure_2_inequality,
  width = 15,
  height = 10,
  dpi = 1000,
  device = png_device,
  limitsize = FALSE
)

ggsave(
  "outputs/Figure_2_Socioeconomic_inequalities_ECDI2030.pdf",
  figure_2_inequality,
  width = 15,
  height = 10,
  device = cairo_pdf,
  limitsize = FALSE
)


# ============================================================
# 8. ANALYSIS REPORT EXPORT (not the submission manuscript)
# ============================================================

output_dir <- "outputs"
output_docx <- file.path(output_dir, "MICS6_ECDI2030_analysis_results.docx")

export_warning_option <- getOption("warn")
options(warn = -1)

overall_ecdi <- survey::svymean(~ecdi, mics6_design, na.rm = TRUE)
ecdi_distribution <- as.data.frame(table(data_final$ecdi, useNA = "ifany"))
names(ecdi_distribution) <- c("ECDI", "n")

region_prevalence <- survey::svyby(~ecdi, ~HH7, mics6_design, survey::svymean, na.rm = TRUE) %>%
  as.data.frame() %>%
  transmute(
    Region = as.character(HH7),
    `Weighted prevalence (%)` = round(ecdi * 100, 1),
    `SE (%)` = round(se * 100, 1)
  )

doc <- read_docx()
doc <- body_add_par(doc, "MICS6 ECDI2030 - Reproduced Analysis Results", style = "heading 1")
doc <- body_add_par(doc, format(Sys.time(), "Export date: %Y-%m-%d %H:%M"), style = "Normal")
doc <- body_add_par(doc, "Note: ECDI2030 was coded as 1 = on-track and 0 = not on-track. The percentages below use MICS6 complex survey weights.", style = "Normal")

doc <- body_add_par(doc, "Quality Checks", style = "heading 1")
doc <- body_add_par(doc, "ECDI Variable Check", style = "heading 2")
doc <- body_add_par(
  doc,
  sprintf(
    "The overall survey-weighted prevalence of on-track ECDI2030 was %.1f%% (SE %.1f%%).",
    as.numeric(overall_ecdi) * 100,
    SE(overall_ecdi) * 100
  ),
  style = "Normal"
)
doc <- body_add_flextable(doc, flextable(ecdi_distribution) %>% autofit())

doc <- body_add_par(doc, "ECDI2030 Prevalence by Region", style = "heading 2")
doc <- body_add_par(doc, "Survey-weighted regional prevalence and standard errors in the complete-case analytic sample.", style = "Normal")
doc <- body_add_flextable(doc, flextable(region_prevalence) %>% autofit())

doc <- body_add_par(doc, "Main Result Tables", style = "heading 1")
doc <- body_add_par(doc, "Table 1. Baseline Characteristics", style = "heading 2")
doc <- body_add_par(doc, "Counts are unweighted; percentages, prevalence estimates, confidence intervals, and p-values account for the complex survey design.", style = "Normal")
doc <- body_add_flextable(doc, as_flex_table(final_table_1) %>% autofit())

doc <- body_add_par(doc, sprintf("Table 2. Multivariable logistic regression for ECDI2030 status among children aged 24-59 months in Vietnam (N = %s)", format(nrow(data_model), big.mark = ",")), style = "heading 2")
doc <- body_add_flextable(doc, as_flex_table(main_table_logistic_24_59) %>% autofit())

doc <- body_add_par(doc, sprintf("Table 3. Decomposition of wealth-related inequality in ECDI2030 status among children aged 24-59 months in Vietnam (N = %s)", format(nrow(data_model), big.mark = ",")), style = "heading 2")
doc <- body_add_flextable(doc, flextable(main_table_decomposition_24_59) %>% autofit())

doc <- body_add_par(doc, "Supplementary Tables", style = "heading 1")
doc <- body_add_par(doc, "Supplementary Table S1. Measurement domains and age-specific thresholds of the ECDI2030", style = "heading 2")
doc <- body_add_flextable(doc, flextable(appendix_table_1_thresholds) %>% autofit())

doc <- body_add_par(doc, sprintf("Supplementary Table S2. Inequality indices in the Early Childhood Development Index (ECDI2030) by household wealth and mother's education in Vietnam (N = %s)", format(nrow(data_model), big.mark = ",")), style = "heading 2")
doc <- body_add_flextable(doc, flextable(appendix_table_2_inequality) %>% autofit())

doc <- body_add_par(doc, sprintf("Supplementary Table S3. Multivariable logistic regression for ECDI2030 status among children aged 36-59 months in Vietnam (N = %s)", format(nrow(data_model_ub6), big.mark = ",")), style = "heading 2")
doc <- body_add_flextable(doc, as_flex_table(appendix_table_3_logistic_36_59) %>% autofit())

doc <- body_add_par(doc, sprintf("Supplementary Table S4. Decomposition of wealth-related inequality in ECDI2030 status among children aged 36-59 months in Vietnam (N = %s)", format(nrow(data_model_ub6), big.mark = ",")), style = "heading 2")
doc <- body_add_flextable(doc, flextable(appendix_table_4_decomposition_36_59) %>% autofit())

doc <- body_add_par(doc, "Main Figures", style = "heading 1")
if (file.exists("outputs/Figure_1_On_track_ECDI2030_publication.png")) {
  doc <- body_add_par(doc, "Figure 1. Regional prevalence of on-track ECDI2030", style = "heading 2")
  doc <- body_add_par(doc, "This figure uses survey-weighted ECDI2030 prevalence across the six MICS6 regions.", style = "Normal")
  doc <- body_add_img(doc, src = "outputs/Figure_1_On_track_ECDI2030_publication.png", width = 5.8, height = 7.0)
}

if (file.exists("outputs/Figure_2_Socioeconomic_inequalities_ECDI2030.png")) {
  doc <- body_add_par(doc, sprintf("Figure 2. Socioeconomic inequalities in ECDI2030 by household wealth and maternal education (N = %s)", format(nrow(data_model), big.mark = ",")), style = "heading 2")
  doc <- body_add_img(doc, src = "outputs/Figure_2_Socioeconomic_inequalities_ECDI2030.png", width = 7.2, height = 4.7)
}

invisible(print(doc, target = output_docx))
options(warn = export_warning_option)
