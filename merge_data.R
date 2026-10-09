# Load required libraries
library(haven)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# --- IMPORT DATA ---
df_ch <- read_sav("data_raw/ch.sav")
df_wm <- read_sav("data_raw/wm.sav")
df_hl <- read_sav("data_raw/hl.sav")
df_hh <- read_sav("data_raw/hh.sav")

# 🔴 CRITICAL STEP TO PREVENT MISSING DATA:
# Force R to extract all labels as text (e.g., "YES", "NO") to match pandas.read_spss behavior
convert_to_labels <- function(df) {
  df %>% mutate(across(where(is.labelled), as_factor))
}

df_ch <- convert_to_labels(df_ch)
df_wm <- convert_to_labels(df_wm)
df_hl <- convert_to_labels(df_hl)
df_hh <- convert_to_labels(df_hh)

# Ensure CAGE is always a numeric variable before filtering
df_ch$CAGE <- suppressWarnings(as.numeric(as.character(df_ch$CAGE)))

# --- STEP 1: Filter target age group (24-59 months) ---
ch_5y <- df_ch %>%
  filter(CAGE >= 24 & CAGE <= 59)

# --- STEP 2: Standardize identification keys ---
# Exactly simulates the command: pd.to_numeric().fillna(-1).astype(int).astype(str).str.strip()
standardize_keys <- function(df, keys) {
  df %>%
    mutate(across(any_of(keys), ~ {
      val <- as.character(.) 
      val <- suppressWarnings(as.numeric(val)) # Convert to numeric, errors coerced to NA
      val <- replace_na(val, -1)               # Fill NAs with -1
      str_trim(as.character(as.integer(val)))  # Remove decimals, convert back to a clean string
    }))
}

keys_to_fix <- c("HH1", "HH2", "LN", "HL1", "UF3", "UF4")

ch_5y <- standardize_keys(ch_5y, keys_to_fix)
df_wm <- standardize_keys(df_wm, keys_to_fix)
df_hl <- standardize_keys(df_hl, keys_to_fix)
df_hh <- standardize_keys(df_hh, keys_to_fix)

# --- STEP 3: Merge CH + WM (Extract Mother/Caretaker information) ---
df1 <- ch_5y %>%
  left_join(df_wm, 
            by = c("HH1" = "HH1", "HH2" = "HH2", "UF4" = "LN"), 
            suffix = c("", "_wm"))

# --- STEP 4: Merge with HL (Extract Child demographic information) ---
df2 <- df1 %>%
  left_join(df_hl, 
            by = c("HH1" = "HH1", "HH2" = "HH2", "LN" = "HL1"), 
            suffix = c("", "_hl"))

# --- STEP 5: Merge with HH (Extract Household information) ---
df_final <- df2 %>%
  left_join(df_hh, 
            by = c("HH1", "HH2"), 
            suffix = c("", "_hh"))

# --- CHECK & EXPORT FILE ---
print(paste("✅ Successfully created df_final. Number of columns:", ncol(df_final), "- Number of rows:", nrow(df_final)))
write.csv(df_final, file = "data_raw/mics6_data_final.csv")