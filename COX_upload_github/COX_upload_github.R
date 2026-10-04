################################################################################
# Required upstream objects (must exist before sourcing this file):
#   - final_data_merged_2 : contains 高峰号, 接种时间, 疫苗, days_from_enrollment
#   - data_59701          : analysis cohort with baseline covariates, outcomes,
#                           and death time
#   - plot_data_dedup     : long-format disease onset records
#                           (高峰号, days_from_enrollment, disease)
#   - data0_67572         : additional variables, e.g., continuous age 年龄
#
# Note: Due to data management regulations, the original individual-level data 
# cannot be provided.
################################################################################

# ------------------------------------------------------------------------------
# 0. Packages and global parameters
# ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(lubridate)
library(survival)
library(purrr)
library(broom)
library(tibble)
library(cmprsk)
library(data.table)

lag_days        <- 14L
protection_days <- 365L

# ------------------------------------------------------------------------------
# 1. Extract vaccination records
# ------------------------------------------------------------------------------
vaccine_data_deva <- final_data_merged_2 %>%
  select(高峰号, 接种时间, 疫苗, days_from_enrollment) %>%
  rename(
    vaccine_date = 接种时间,
    vaccine_type = 疫苗
  ) %>%
  mutate(
    vaccine_date   = as.Date(vaccine_date),
    vaccine_status = 1
  ) %>%
  arrange(高峰号, vaccine_date)

# ------------------------------------------------------------------------------
# 2. Helper function: merge overlapping influenza protection intervals
# ------------------------------------------------------------------------------
merge_intervals <- function(df) {
  if (nrow(df) == 0) return(df)
  
  df <- df %>% arrange(vacc_start)
  
  flu_df   <- df %>% filter(vaccine_type == "流感")
  other_df <- df %>% filter(vaccine_type != "流感")
  
  if (nrow(flu_df) == 0) return(df)
  
  merged <- list()
  current_start <- flu_df$vacc_start[1]
  current_end   <- flu_df$vacc_end[1]
  
  if (nrow(flu_df) == 1) {
    merged[[1]] <- data.frame(
      高峰号       = flu_df$高峰号[1],
      vaccine_type = "流感",
      vacc_start   = current_start,
      vacc_end     = current_end
    )
  } else {
    for (i in 2:nrow(flu_df)) {
      next_start <- flu_df$vacc_start[i]
      next_end   <- flu_df$vacc_end[i]
      
      if (next_start <= current_end) {
        current_end <- max(current_end, next_end)
      } else {
        merged[[length(merged) + 1]] <- data.frame(
          高峰号       = flu_df$高峰号[1],
          vaccine_type = "流感",
          vacc_start   = current_start,
          vacc_end     = current_end
        )
        current_start <- next_start
        current_end   <- next_end
      }
    }
    merged[[length(merged) + 1]] <- data.frame(
      高峰号       = flu_df$高峰号[1],
      vaccine_type = "流感",
      vacc_start   = current_start,
      vacc_end     = current_end
    )
  }
  
  bind_rows(merged) %>% bind_rows(other_df)
}

# ------------------------------------------------------------------------------
# 3. Analyses restricted to influenza and zoster vaccines
# ------------------------------------------------------------------------------

## 3.1 Dementia outcome: influenza and herpes zoster only
baseline_data_deva <- data_59701 %>%
  select(
    高峰号,
    调查日期, dementia_censor_date,
    dementia_event,
    是否患高血压, 是否患冠心病, 是否患脑卒中, 是否患糖尿病,
    是否患甲状腺疾病,
    是否患COPD, 是否患慢性肾病, 吸烟, 喝酒, 喝茶,
    每周体锻, 焦虑抑郁_new,
    visit_category, 年龄分组,
    文化程度_new, 实际睡眠时间_new, 入组年份, 地区, 性别, BMI分类,
    是否患帕金森病, 亲属患帕金森病,
    亲属患阿尔茨海默病,
    是否患抑郁症, 亲属患抑郁症
  ) %>%
  rename(
    enrollment_date = 调查日期,
    censor_date     = dementia_censor_date,
    event           = dementia_event
  ) %>%
  mutate(
    enrollment_date = as.Date(enrollment_date),
    censor_date     = as.Date(censor_date)
  ) %>%
  mutate(across(c(是否患高血压:性别, BMI分类), as.factor))

start_data_deva <- baseline_data_deva %>%
  select(高峰号, enrollment_date, censor_date, event, everything()) %>%
  mutate(
    start_time       = 0,
    stop_time        = as.numeric(difftime(censor_date, enrollment_date, units = "days")),
    enrollment_month = month(enrollment_date),
    enrollment_year  = year(enrollment_date)
  ) %>%
  distinct(高峰号, .keep_all = TRUE) %>%
  filter(stop_time > 0)

followup_intervals_deva <- vaccine_data_deva %>%
  left_join(start_data_deva %>% select(高峰号, stop_time), by = "高峰号") %>%
  mutate(
    vacc_start = days_from_enrollment + lag_days,
    vacc_end   = case_when(
      vaccine_type == "流感" ~ days_from_enrollment + lag_days + protection_days,
      vaccine_type == "带疱" ~ stop_time,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(vacc_start < stop_time) %>%
  select(高峰号, vacc_start, vacc_end, vaccine_type)

# Keep only the first herpes zoster record per subject
herpes_rows <- followup_intervals_deva %>%
  filter(vaccine_type == "带疱") %>%
  group_by(高峰号) %>%
  slice_min(vacc_start, n = 1, with_ties = FALSE) %>%
  ungroup()
other_rows <- followup_intervals_deva %>%
  filter(vaccine_type != "带疱")
followup_intervals_deva <- bind_rows(herpes_rows, other_rows)

# Check for events occurring shortly after the first vaccination
first_vaccine <- vaccine_data_deva %>%
  filter(days_from_enrollment >= 0) %>%
  group_by(高峰号) %>%
  summarise(first_vaccine_day = min(days_from_enrollment), .groups = "drop")

early_event_check <- start_data_deva %>%
  left_join(first_vaccine, by = "高峰号") %>%
  mutate(
    days_to_event = stop_time,
    event_within_14d_post_vaccine = if_else(
      !is.na(first_vaccine_day) &
        event == 1 &
        days_to_event >= first_vaccine_day &
        days_to_event <= first_vaccine_day + 30,
      1, 0
    )
  )
table(early_event_check$event_within_14d_post_vaccine)

# Merge individual protection windows
merged_intervals_deva <- followup_intervals_deva %>%
  group_split(高峰号) %>%
  map_dfr(merge_intervals)

vaccine_on_events_deva <- merged_intervals_deva %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_start,
    vaccine_current = 1
  )

vaccine_off_events_deva <- merged_intervals_deva %>%
  filter(vaccine_type == "流感") %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_end,
    vaccine_current = 0
  )

vaccine_events_deva <- bind_rows(
  vaccine_on_events_deva,
  vaccine_off_events_deva
) %>% arrange(高峰号, event_time)

# tmerge: build time-dependent dataset
td_data_deva <- tmerge(
  data1 = start_data_deva,
  data2 = start_data_deva,
  id    = 高峰号,
  event = event(stop_time, event)
)

vaccine_wide <- vaccine_events_deva %>%
  pivot_wider(
    id_cols      = c(高峰号, event_time),
    names_from   = vaccine_type,
    values_from  = vaccine_current,
    values_fill  = 0,
    names_prefix = "vaccine_"
  )

td_data_deva <- tmerge(
  data1 = td_data_deva,
  data2 = vaccine_wide,
  id    = 高峰号,
  vaccine_flu    = tdc(event_time, vaccine_流感),
  vaccine_zoster = tdc(event_time, vaccine_带疱)
)

# Fill missing vaccine indicators with 0
td_data_deva <- td_data_deva %>%
  mutate(
    vaccine_flu    = ifelse(is.na(vaccine_flu), 0, vaccine_flu),
    vaccine_zoster = ifelse(is.na(vaccine_zoster), 0, vaccine_zoster)
  )

cox_model_deva <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_zoster +
    vaccine_flu +
    是否患高血压 + 是否患冠心病 + 是否患脑卒中 + 是否患糖尿病 +
    是否患甲状腺疾病 + 是否患COPD + 是否患慢性肾病 +
    吸烟 + 喝酒 + 喝茶 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    地区 + 性别 + BMI分类 + 年龄分组 +
    是否患帕金森病 +
    亲属患帕金森病 +
    亲属患阿尔茨海默病 +
    是否患抑郁症 +
    亲属患抑郁症,
  data = td_data_deva
)
summary(cox_model_deva)

## 3.2 Alzheimer's disease (AD) outcome
baseline_data_ADva <- data_59701 %>%
  select(
    高峰号,
    调查日期, AD_censor_date,
    AD_event,
    是否患高血压, 是否患冠心病, 是否患脑卒中, 是否患糖尿病,
    是否患甲状腺疾病,
    是否患COPD, 是否患慢性肾病, 吸烟, 喝酒, 喝茶,
    每周体锻, 焦虑抑郁_new,
    visit_category, 年龄分组,
    文化程度_new, 实际睡眠时间_new, 入组年份, 地区, 性别, BMI分类,
    是否患帕金森病, 亲属患帕金森病,
    亲属患阿尔茨海默病,
    是否患抑郁症, 亲属患抑郁症
  ) %>%
  rename(
    enrollment_date = 调查日期,
    censor_date     = AD_censor_date,
    event           = AD_event
  ) %>%
  mutate(
    enrollment_date = as.Date(enrollment_date),
    censor_date     = as.Date(censor_date)
  ) %>%
  mutate(across(c(是否患高血压:性别, BMI分类), as.factor))

start_data_ADva <- baseline_data_ADva %>%
  select(高峰号, enrollment_date, censor_date, event, everything()) %>%
  mutate(
    start_time       = 0,
    stop_time        = as.numeric(difftime(censor_date, enrollment_date, units = "days")),
    enrollment_month = month(enrollment_date),
    enrollment_year  = year(enrollment_date)
  ) %>%
  distinct(高峰号, .keep_all = TRUE) %>%
  filter(stop_time > 0)

followup_intervals_ADva <- vaccine_data_deva %>%
  left_join(start_data_ADva %>% select(高峰号, stop_time), by = "高峰号") %>%
  mutate(
    vacc_start = days_from_enrollment + lag_days,
    vacc_end   = case_when(
      vaccine_type == "流感" ~ days_from_enrollment + lag_days + protection_days,
      vaccine_type == "带疱" ~ stop_time,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(vacc_start < stop_time) %>%
  select(高峰号, vacc_start, vacc_end, vaccine_type)

herpes_rows <- followup_intervals_ADva %>%
  filter(vaccine_type == "带疱") %>%
  group_by(高峰号) %>%
  slice_min(vacc_start, n = 1, with_ties = FALSE) %>%
  ungroup()
other_rows <- followup_intervals_ADva %>%
  filter(vaccine_type != "带疱")
followup_intervals_ADva <- bind_rows(herpes_rows, other_rows)

merged_intervals_ADva <- followup_intervals_ADva %>%
  group_split(高峰号) %>%
  map_dfr(merge_intervals)

vaccine_on_events_ADva <- merged_intervals_ADva %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_start,
    vaccine_current = 1
  )

vaccine_off_events_ADva <- merged_intervals_ADva %>%
  filter(vaccine_type == "流感") %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_end,
    vaccine_current = 0
  )

vaccine_events_ADva <- bind_rows(
  vaccine_on_events_ADva,
  vaccine_off_events_ADva
) %>% arrange(高峰号, event_time)

td_data_ADva <- tmerge(
  data1 = start_data_ADva,
  data2 = start_data_ADva,
  id    = 高峰号,
  event = event(stop_time, event)
)

vaccine_wide <- vaccine_events_ADva %>%
  pivot_wider(
    id_cols      = c(高峰号, event_time),
    names_from   = vaccine_type,
    values_from  = vaccine_current,
    values_fill  = 0,
    names_prefix = "vaccine_"
  )

td_data_ADva <- tmerge(
  data1 = td_data_ADva,
  data2 = vaccine_wide,
  id    = 高峰号,
  vaccine_flu    = tdc(event_time, vaccine_流感),
  vaccine_zoster = tdc(event_time, vaccine_带疱)
)

td_data_ADva <- td_data_ADva %>%
  mutate(
    vaccine_flu    = ifelse(is.na(vaccine_flu), 0, vaccine_flu),
    vaccine_zoster = ifelse(is.na(vaccine_zoster), 0, vaccine_zoster)
  )

cox_model_ADva <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_flu +
    vaccine_zoster +
    是否患高血压 + 是否患冠心病 + 是否患脑卒中 + 是否患糖尿病 +
    是否患甲状腺疾病 + 是否患COPD + 是否患慢性肾病 +
    吸烟 + 喝酒 + 喝茶 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    地区 + 性别 + BMI分类 + 年龄分组,
  data = td_data_ADva
)
summary(cox_model_ADva)

# ------------------------------------------------------------------------------
# 4. Analyses restricted to influenza, herpes zoster and pneumococcal vaccines
# ------------------------------------------------------------------------------

## 4.1 Dementia outcome: influenza, herpes zoster and pneumococcal
baseline_data_deva <- data_59701 %>%
  select(
    高峰号,
    调查日期, dementia_censor_date,
    dementia_event,
    是否患高血压, 是否患冠心病, 是否患脑卒中, 是否患糖尿病,
    是否患甲状腺疾病,
    是否患COPD, 是否患慢性肾病, 吸烟, 喝酒, 喝茶,
    每周体锻, 焦虑抑郁_new,
    visit_category, 年龄分组,
    文化程度_new, 实际睡眠时间_new, 入组年份, 地区, 性别, BMI分类,
    是否患帕金森病, 亲属患帕金森病,
    亲属患阿尔茨海默病,
    是否患抑郁症, 亲属患抑郁症,
    肺炎疫苗接种
  ) %>%
  rename(
    enrollment_date = 调查日期,
    censor_date     = dementia_censor_date,
    event           = dementia_event
  ) %>%
  mutate(
    enrollment_date = as.Date(enrollment_date),
    censor_date     = as.Date(censor_date)
  ) %>%
  mutate(across(c(是否患高血压:性别, BMI分类), as.factor))

start_data_deva <- baseline_data_deva %>%
  select(高峰号, enrollment_date, censor_date, event, everything()) %>%
  mutate(
    start_time       = 0,
    stop_time        = as.numeric(difftime(censor_date, enrollment_date, units = "days")),
    enrollment_month = month(enrollment_date),
    enrollment_year  = year(enrollment_date)
  ) %>%
  distinct(高峰号, .keep_all = TRUE) %>%
  filter(stop_time > 0)

followup_intervals_deva <- vaccine_data_deva %>%
  left_join(start_data_deva %>% select(高峰号, stop_time), by = "高峰号") %>%
  mutate(
    vacc_start = days_from_enrollment + lag_days,
    vacc_end   = case_when(
      vaccine_type == "流感" ~ days_from_enrollment + lag_days + protection_days,
      vaccine_type == "带疱" ~ stop_time,
      vaccine_type == "肺炎" ~ stop_time,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(vacc_start < stop_time) %>%
  select(高峰号, vacc_start, vacc_end, vaccine_type)

herpes_rows <- followup_intervals_deva %>%
  filter(vaccine_type == "带疱") %>%
  group_by(高峰号) %>%
  slice_min(vacc_start, n = 1, with_ties = FALSE) %>%
  ungroup()
other_rows <- followup_intervals_deva %>%
  filter(vaccine_type != "带疱")
followup_intervals_deva <- bind_rows(herpes_rows, other_rows)

# Add baseline pneumococcal vaccination (from enrollment)
pneumonia_data <- subset(start_data_deva, `肺炎疫苗接种` == 1)
pneumonia_data <- pneumonia_data[, c("高峰号", "肺炎疫苗接种", "stop_time")]
names(pneumonia_data) <- c("高峰号", "vaccine_type", "vacc_end")
new_rows <- data.frame(
  高峰号       = pneumonia_data$高峰号,
  vaccine_type = "肺炎",
  vacc_start   = 0,
  vacc_end     = pneumonia_data$vacc_end
)
followup_intervals_deva2 <- rbind(followup_intervals_deva, new_rows)

merged_intervals_deva2 <- followup_intervals_deva2 %>%
  group_split(高峰号) %>%
  map_dfr(merge_intervals)

vaccine_on_events_deva2 <- merged_intervals_deva2 %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_start,
    vaccine_current = 1
  )

vaccine_off_events_deva2 <- merged_intervals_deva2 %>%
  filter(vaccine_type == "流感") %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_end,
    vaccine_current = 0
  )

vaccine_events_deva2 <- bind_rows(
  vaccine_on_events_deva2,
  vaccine_off_events_deva2
) %>% arrange(高峰号, event_time)

td_data_deva <- tmerge(
  data1 = start_data_deva,
  data2 = start_data_deva,
  id    = 高峰号,
  event = event(stop_time, event)
)

vaccine_wide2 <- vaccine_events_deva2 %>%
  pivot_wider(
    id_cols      = c(高峰号, event_time),
    names_from   = vaccine_type,
    values_from  = vaccine_current,
    values_fill  = 0,
    names_prefix = "vaccine_"
  )

td_data_deva2 <- tmerge(
  data1 = td_data_deva,
  data2 = vaccine_wide2,
  id    = 高峰号,
  vaccine_flu    = tdc(event_time, vaccine_流感),
  vaccine_zoster = tdc(event_time, vaccine_带疱),
  vaccine_pneu   = tdc(event_time, vaccine_肺炎)
)

# Fill missing indicators with 0, and keep pneumococcal/herpes zoster status = 1
# once vaccinated (never reverts to 0)
td_data_deva2 <- td_data_deva2 %>%
  mutate(
    vaccine_flu    = ifelse(is.na(vaccine_flu), 0, vaccine_flu),
    vaccine_zoster = ifelse(is.na(vaccine_zoster), 0, vaccine_zoster),
    vaccine_pneu   = ifelse(is.na(vaccine_pneu), 0, vaccine_pneu)
  ) %>%
  group_by(高峰号) %>%
  arrange(tstart) %>%
  mutate(
    vaccine_pneu   = ifelse(cummax(vaccine_pneu) == 1, 1, vaccine_pneu),
    vaccine_zoster = ifelse(cummax(vaccine_zoster) == 1, 1, vaccine_zoster)
  ) %>%
  ungroup() %>%
  arrange(高峰号, tstart)

cox_model_deva2 <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    是否患高血压 + 是否患糖尿病 +
    是否患甲状腺疾病 + 是否患COPD + 是否患慢性肾病 +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_deva2
)
summary(cox_model_deva2)

## 4.2 AD outcome: influenza, herpes zoster and pneumococcal
baseline_data_ADva <- data_59701 %>%
  select(
    高峰号,
    调查日期, AD_censor_date,
    AD_event,
    是否患高血压, 是否患冠心病, 是否患脑卒中, 是否患糖尿病,
    是否患甲状腺疾病,
    是否患COPD, 是否患慢性肾病, 吸烟, 喝酒, 喝茶,
    每周体锻, 焦虑抑郁_new,
    visit_category, 年龄分组,
    文化程度_new, 实际睡眠时间_new, 入组年份, 地区, 性别, BMI分类,
    是否患帕金森病, 亲属患帕金森病,
    亲属患阿尔茨海默病,
    是否患抑郁症, 亲属患抑郁症,
    肺炎疫苗接种
  ) %>%
  rename(
    enrollment_date = 调查日期,
    censor_date     = AD_censor_date,
    event           = AD_event
  ) %>%
  mutate(
    enrollment_date = as.Date(enrollment_date),
    censor_date     = as.Date(censor_date)
  ) %>%
  mutate(across(c(是否患高血压:性别, BMI分类), as.factor))

start_data_ADva <- baseline_data_ADva %>%
  select(高峰号, enrollment_date, censor_date, event, everything()) %>%
  mutate(
    start_time       = 0,
    stop_time        = as.numeric(difftime(censor_date, enrollment_date, units = "days")),
    enrollment_month = month(enrollment_date),
    enrollment_year  = year(enrollment_date)
  ) %>%
  distinct(高峰号, .keep_all = TRUE) %>%
  filter(stop_time > 0)

followup_intervals_ADva <- vaccine_data_deva %>%
  left_join(start_data_ADva %>% select(高峰号, stop_time), by = "高峰号") %>%
  mutate(
    vacc_start = days_from_enrollment + lag_days,
    vacc_end   = case_when(
      vaccine_type == "流感" ~ days_from_enrollment + lag_days + protection_days,
      vaccine_type == "带疱" ~ stop_time,
      vaccine_type == "肺炎" ~ stop_time,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(vacc_start < stop_time) %>%
  select(高峰号, vacc_start, vacc_end, vaccine_type)

herpes_rows <- followup_intervals_ADva %>%
  filter(vaccine_type == "带疱") %>%
  group_by(高峰号) %>%
  slice_min(vacc_start, n = 1, with_ties = FALSE) %>%
  ungroup()
other_rows <- followup_intervals_ADva %>%
  filter(vaccine_type != "带疱")
followup_intervals_ADva <- bind_rows(herpes_rows, other_rows)

pneumonia_data <- subset(start_data_ADva, `肺炎疫苗接种` == 1)
pneumonia_data <- pneumonia_data[, c("高峰号", "肺炎疫苗接种", "stop_time")]
names(pneumonia_data) <- c("高峰号", "vaccine_type", "vacc_end")
new_rows <- data.frame(
  高峰号       = pneumonia_data$高峰号,
  vaccine_type = "肺炎",
  vacc_start   = 0,
  vacc_end     = pneumonia_data$vacc_end
)
followup_intervals_ADva2 <- rbind(followup_intervals_ADva, new_rows)

merged_intervals_ADva2 <- followup_intervals_ADva2 %>%
  group_split(高峰号) %>%
  map_dfr(merge_intervals)

vaccine_on_events_ADva2 <- merged_intervals_ADva2 %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_start,
    vaccine_current = 1
  )

vaccine_off_events_ADva2 <- merged_intervals_ADva2 %>%
  filter(vaccine_type == "流感") %>%
  transmute(
    高峰号,
    vaccine_type,
    event_time     = vacc_end,
    vaccine_current = 0
  )

vaccine_events_ADva2 <- bind_rows(
  vaccine_on_events_ADva2,
  vaccine_off_events_ADva2
) %>% arrange(高峰号, event_time)

td_data_ADva2 <- tmerge(
  data1 = start_data_ADva,
  data2 = start_data_ADva,
  id    = 高峰号,
  event = event(stop_time, event)
)

vaccine_wide2 <- vaccine_events_ADva2 %>%
  pivot_wider(
    id_cols      = c(高峰号, event_time),
    names_from   = vaccine_type,
    values_from  = vaccine_current,
    values_fill  = 0,
    names_prefix = "vaccine_"
  )

td_data_ADva2 <- tmerge(
  data1 = td_data_ADva2,
  data2 = vaccine_wide2,
  id    = 高峰号,
  vaccine_flu    = tdc(event_time, vaccine_流感),
  vaccine_zoster = tdc(event_time, vaccine_带疱),
  vaccine_pneu   = tdc(event_time, vaccine_肺炎)
)

td_data_ADva2 <- td_data_ADva2 %>%
  mutate(
    vaccine_flu    = ifelse(is.na(vaccine_flu), 0, vaccine_flu),
    vaccine_zoster = ifelse(is.na(vaccine_zoster), 0, vaccine_zoster),
    vaccine_pneu   = ifelse(is.na(vaccine_pneu), 0, vaccine_pneu)
  ) %>%
  group_by(高峰号) %>%
  arrange(tstart) %>%
  mutate(
    vaccine_pneu   = ifelse(cummax(vaccine_pneu) == 1, 1, vaccine_pneu),
    vaccine_zoster = ifelse(cummax(vaccine_zoster) == 1, 1, vaccine_zoster)
  ) %>%
  ungroup() %>%
  arrange(高峰号, tstart)

cox_model_ADva2 <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    是否患高血压 + 是否患糖尿病 +
    是否患甲状腺疾病 + 是否患COPD + 是否患慢性肾病 +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_ADva2
)
summary(cox_model_ADva2)

# ------------------------------------------------------------------------------
# 5. Vaccines plus time-dependent diseases
# ------------------------------------------------------------------------------

## 5.1 Dementia outcome: vaccines plus time-dependent diseases
baseline_data_deva <- data_59701 %>%
  select(
    高峰号,
    调查日期, dementia_censor_date,
    dementia_event,
    是否患高血压, 是否患冠心病, 是否患脑卒中, 是否患糖尿病,
    是否患甲状腺疾病,
    是否患COPD, 是否患慢性肾病, 吸烟, 喝酒, 喝茶,
    每周体锻, 焦虑抑郁_new,
    visit_category, 年龄分组,
    文化程度_new, 实际睡眠时间_new, 入组年份, 地区, 性别, BMI分类,
    是否患帕金森病, 亲属患帕金森病,
    亲属患阿尔茨海默病,
    是否患抑郁症, 亲属患抑郁症,
    肺炎疫苗接种
  ) %>%
  rename(
    enrollment_date = 调查日期,
    censor_date     = dementia_censor_date,
    event           = dementia_event
  ) %>%
  mutate(
    enrollment_date = as.Date(enrollment_date),
    censor_date     = as.Date(censor_date)
  ) %>%
  mutate(across(c(是否患高血压:性别, BMI分类), as.factor))

start_data_deva <- baseline_data_deva %>%
  select(高峰号, enrollment_date, censor_date, event, everything()) %>%
  mutate(
    start_time       = 0,
    stop_time        = as.numeric(difftime(censor_date, enrollment_date, units = "days")),
    enrollment_month = month(enrollment_date),
    enrollment_year  = year(enrollment_date)
  ) %>%
  distinct(高峰号, .keep_all = TRUE) %>%
  filter(stop_time > 0)

# Incident diseases during follow-up (from long-format plot_data_dedup)
followup_intervals_dedi <- plot_data_dedup %>%
  left_join(start_data_deva %>% select(高峰号, stop_time), by = "高峰号") %>%
  mutate(
    disease_start = days_from_enrollment,
    disease_end   = stop_time
  ) %>%
  select(高峰号, disease_start, disease_end, disease)

disease_columns <- c(
  "是否患高血压", "是否患冠心病", "是否患脑卒中", "是否患糖尿病",
  "是否患甲状腺疾病", "是否患COPD", "是否患慢性肾病",
  "是否患帕金森病", "是否患抑郁症"
)

start_data_deva[disease_columns] <- lapply(
  start_data_deva[disease_columns],
  function(x) as.numeric(as.character(x))
)

# Baseline prevalent diseases (convert wide to long, keep disease = 1)
disease_baseline <- start_data_deva %>%
  pivot_longer(
    cols      = all_of(disease_columns),
    names_to  = "disease",
    values_to = "value"
  ) %>%
  mutate(disease = case_when(
    disease == "是否患高血压"   ~ "data_hypertension",
    disease == "是否患冠心病"   ~ "data_CHD",
    disease == "是否患脑卒中"   ~ "data_stroke",
    disease == "是否患糖尿病"   ~ "data_DM",
    disease == "是否患甲状腺疾病" ~ "data_jiakang",
    disease == "是否患COPD"     ~ "data_COPD",
    disease == "是否患慢性肾病" ~ "data_CKD",
    disease == "是否患帕金森病" ~ "data_parkin",
    disease == "是否患抑郁症"   ~ "data_depress",
    TRUE ~ disease
  )) %>%
  filter(value == 1) %>%
  select(高峰号, disease, value, stop_time)

disease_baseline <- disease_baseline %>%
  mutate(disease_start = 0) %>%
  rename(disease_end = stop_time) %>%
  select(-value)

followup_intervals_dedi2 <- rbind(followup_intervals_dedi, disease_baseline)

# Keep the earliest onset for each disease per subject (deduplicate)
followup_intervals_dedi2 <- followup_intervals_dedi2 %>%
  group_by(高峰号, disease) %>%
  slice_min(disease_start, n = 1) %>%
  ungroup() %>%
  arrange(高峰号, disease)

disease_events_dedi <- followup_intervals_dedi2 %>%
  transmute(
    高峰号,
    disease,
    event_time      = disease_start,
    disease_current = 1
  )

# Combine with vaccine events (renamed as diseases)
vaccine_events_dedi <- vaccine_events_deva2 %>%
  rename(
    disease         = vaccine_type,
    disease_current = vaccine_current
  )

# Pneumococcal and herpes zoster: keep only the earliest vaccination per subject
pneumonia_zoster_earliest <- vaccine_events_dedi %>%
  filter(disease %in% c("肺炎", "带疱")) %>%
  group_by(高峰号, disease) %>%
  slice_min(event_time, n = 1, with_ties = FALSE) %>%
  ungroup()
other_diseases <- vaccine_events_dedi %>%
  filter(!disease %in% c("肺炎", "带疱"))
vaccine_events_dedi <- bind_rows(pneumonia_zoster_earliest, other_diseases)

disease_events_dedi <- full_join(vaccine_events_dedi, disease_events_dedi, by = NULL)

# Wide format with LOCF so that statuses do not interfere across intervals
disease_wide <- disease_events_dedi %>%
  pivot_wider(
    id_cols      = c(高峰号, event_time),
    names_from   = disease,
    values_from  = disease_current,
    values_fill  = NA
  ) %>%
  arrange(高峰号, event_time)

disease_wide <- disease_wide %>%
  group_by(高峰号) %>%
  arrange(event_time) %>%
  tidyr::fill(everything(), .direction = "down") %>%
  ungroup() %>%
  mutate(across(-c(高峰号, event_time), ~ replace_na(.x, 0)))

td_data_deva <- tmerge(
  data1 = start_data_deva,
  data2 = start_data_deva,
  id    = 高峰号,
  event = event(stop_time, event)
)

td_data_dedi <- tmerge(
  data1 = td_data_deva,
  data2 = disease_wide,
  id    = 高峰号,
  vaccine_flu        = tdc(event_time, 流感),
  vaccine_zoster     = tdc(event_time, 带疱),
  vaccine_pneu       = tdc(event_time, 肺炎),
  data_hypertension  = tdc(event_time, data_hypertension),
  data_CHD           = tdc(event_time, data_CHD),
  data_stroke        = tdc(event_time, data_stroke),
  data_DM            = tdc(event_time, data_DM),
  data_jiakang       = tdc(event_time, data_jiakang),
  data_COPD          = tdc(event_time, data_COPD),
  data_CKD           = tdc(event_time, data_CKD),
  data_parkin        = tdc(event_time, data_parkin),
  data_depress       = tdc(event_time, data_depress)
)

# Fill missing with 0 and keep status = 1 once present (never reverts to 0)
td_data_dedi <- td_data_dedi %>%
  mutate(across(
    c(vaccine_flu, vaccine_zoster, vaccine_pneu,
      data_hypertension, data_CHD, data_stroke, data_DM,
      data_jiakang, data_COPD, data_CKD, data_parkin, data_depress),
    ~ ifelse(is.na(.), 0, .)
  )) %>%
  group_by(高峰号) %>%
  arrange(tstart) %>%
  mutate(
    vaccine_zoster     = ifelse(cummax(vaccine_zoster) == 1, 1, vaccine_zoster),
    vaccine_pneu       = ifelse(cummax(vaccine_pneu) == 1, 1, vaccine_pneu),
    data_hypertension  = ifelse(cummax(data_hypertension) == 1, 1, data_hypertension),
    data_CHD           = ifelse(cummax(data_CHD) == 1, 1, data_CHD),
    data_stroke        = ifelse(cummax(data_stroke) == 1, 1, data_stroke),
    data_DM            = ifelse(cummax(data_DM) == 1, 1, data_DM),
    data_jiakang       = ifelse(cummax(data_jiakang) == 1, 1, data_jiakang),
    data_COPD          = ifelse(cummax(data_COPD) == 1, 1, data_COPD),
    data_CKD           = ifelse(cummax(data_CKD) == 1, 1, data_CKD),
    data_parkin        = ifelse(cummax(data_parkin) == 1, 1, data_parkin),
    data_depress       = ifelse(cummax(data_depress) == 1, 1, data_depress)
  ) %>%
  ungroup() %>%
  arrange(高峰号, tstart)

# Remove the original Chinese baseline disease columns
td_data_dedi <- td_data_dedi %>%
  select(-all_of(disease_columns))

# NOTE: The disease model below is reported as THREE Cox models:
#   - Base model (primary analysis): without CHD, stroke, Parkinson's disease,
#     depression, and family history of PD/AD/depression;
#   - Sensitivity analysis 1: base model additionally adjusted for CHD and stroke;
#   - Sensitivity analysis 2: base model additionally adjusted for Parkinson's
#     disease, depression, and family history of PD/AD/depression.

# Base model (primary analysis)
cox_model_dedi_base <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    data_hypertension + data_DM +
    data_jiakang + data_COPD + data_CKD +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_dedi
)
summary(cox_model_dedi_base)

# Sensitivity analysis 1: additionally adjusted for CHD and stroke
cox_model_dedi_sens_cvd <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    data_hypertension + data_DM +
    data_jiakang + data_COPD + data_CKD +
    data_CHD + data_stroke +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_dedi
)
summary(cox_model_dedi_sens_cvd)

# Sensitivity analysis 2: additionally adjusted for Parkinson's disease,
# depression, and family history of PD/AD/depression
cox_model_dedi_sens_neuro <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    data_hypertension + data_DM +
    data_jiakang + data_COPD + data_CKD +
    data_parkin + 亲属患帕金森病 + 亲属患阿尔茨海默病 + data_depress + 亲属患抑郁症 +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_dedi
)
summary(cox_model_dedi_sens_neuro)

## 5.2 AD outcome: vaccines plus time-dependent diseases
baseline_data_ADva <- data_59701 %>%
  select(
    高峰号,
    调查日期, AD_censor_date,
    AD_event,
    是否患高血压, 是否患冠心病, 是否患脑卒中, 是否患糖尿病,
    是否患甲状腺疾病,
    是否患COPD, 是否患慢性肾病, 吸烟, 喝酒, 喝茶,
    每周体锻, 焦虑抑郁_new,
    visit_category, 年龄分组,
    文化程度_new, 实际睡眠时间_new, 入组年份, 地区, 性别, BMI分类,
    是否患帕金森病, 亲属患帕金森病,
    亲属患阿尔茨海默病,
    是否患抑郁症, 亲属患抑郁症,
    肺炎疫苗接种
  ) %>%
  rename(
    enrollment_date = 调查日期,
    censor_date     = AD_censor_date,
    event           = AD_event
  ) %>%
  mutate(
    enrollment_date = as.Date(enrollment_date),
    censor_date     = as.Date(censor_date)
  ) %>%
  mutate(across(c(是否患高血压:性别, BMI分类), as.factor))

start_data_ADva <- baseline_data_ADva %>%
  select(高峰号, enrollment_date, censor_date, event, everything()) %>%
  mutate(
    start_time       = 0,
    stop_time        = as.numeric(difftime(censor_date, enrollment_date, units = "days")),
    enrollment_month = month(enrollment_date),
    enrollment_year  = year(enrollment_date)
  ) %>%
  distinct(高峰号, .keep_all = TRUE) %>%
  filter(stop_time > 0)

followup_intervals_ADdi <- plot_data_dedup %>%
  left_join(start_data_ADva %>% select(高峰号, stop_time), by = "高峰号") %>%
  mutate(
    disease_start = days_from_enrollment,
    disease_end   = stop_time
  ) %>%
  select(高峰号, disease_start, disease_end, disease)

disease_columns <- c(
  "是否患高血压", "是否患冠心病", "是否患脑卒中", "是否患糖尿病",
  "是否患甲状腺疾病", "是否患COPD", "是否患慢性肾病",
  "是否患帕金森病", "是否患抑郁症"
)

start_data_ADva[disease_columns] <- lapply(
  start_data_ADva[disease_columns],
  function(x) as.numeric(as.character(x))
)

disease_baseline <- start_data_ADva %>%
  pivot_longer(
    cols      = all_of(disease_columns),
    names_to  = "disease",
    values_to = "value"
  ) %>%
  mutate(disease = case_when(
    disease == "是否患高血压"   ~ "data_hypertension",
    disease == "是否患冠心病"   ~ "data_CHD",
    disease == "是否患脑卒中"   ~ "data_stroke",
    disease == "是否患糖尿病"   ~ "data_DM",
    disease == "是否患甲状腺疾病" ~ "data_jiakang",
    disease == "是否患COPD"     ~ "data_COPD",
    disease == "是否患慢性肾病" ~ "data_CKD",
    disease == "是否患帕金森病" ~ "data_parkin",
    disease == "是否患抑郁症"   ~ "data_depress",
    TRUE ~ disease
  )) %>%
  filter(value == 1) %>%
  select(高峰号, disease, value, stop_time)

disease_baseline <- disease_baseline %>%
  mutate(disease_start = 0) %>%
  rename(disease_end = stop_time) %>%
  select(-value)

followup_intervals_ADdi2 <- rbind(followup_intervals_ADdi, disease_baseline)

followup_intervals_ADdi2 <- followup_intervals_ADdi2 %>%
  group_by(高峰号, disease) %>%
  slice_min(disease_start, n = 1) %>%
  ungroup() %>%
  arrange(高峰号, disease)

disease_events_ADdi <- followup_intervals_ADdi2 %>%
  transmute(
    高峰号,
    disease,
    event_time      = disease_start,
    disease_current = 1
  )

vaccine_events_ADdi <- vaccine_events_ADva2 %>%
  rename(
    disease         = vaccine_type,
    disease_current = vaccine_current
  )

pneumonia_zoster_earliest <- vaccine_events_ADdi %>%
  filter(disease %in% c("肺炎", "带疱")) %>%
  group_by(高峰号, disease) %>%
  slice_min(event_time, n = 1, with_ties = FALSE) %>%
  ungroup()
other_diseases <- vaccine_events_ADdi %>%
  filter(!disease %in% c("肺炎", "带疱"))
vaccine_events_ADdi <- bind_rows(pneumonia_zoster_earliest, other_diseases)

disease_events_ADdi <- full_join(vaccine_events_ADdi, disease_events_ADdi, by = NULL)

disease_wide <- disease_events_ADdi %>%
  pivot_wider(
    id_cols      = c(高峰号, event_time),
    names_from   = disease,
    values_from  = disease_current,
    values_fill  = NA
  ) %>%
  arrange(高峰号, event_time)

disease_wide <- disease_wide %>%
  group_by(高峰号) %>%
  arrange(event_time) %>%
  tidyr::fill(everything(), .direction = "down") %>%
  ungroup() %>%
  mutate(across(-c(高峰号, event_time), ~ replace_na(.x, 0)))

td_data_ADva <- tmerge(
  data1 = start_data_ADva,
  data2 = start_data_ADva,
  id    = 高峰号,
  event = event(stop_time, event)
)

td_data_ADdi <- tmerge(
  data1 = td_data_ADva,
  data2 = disease_wide,
  id    = 高峰号,
  vaccine_flu        = tdc(event_time, 流感),
  vaccine_zoster     = tdc(event_time, 带疱),
  vaccine_pneu       = tdc(event_time, 肺炎),
  data_hypertension  = tdc(event_time, data_hypertension),
  data_CHD           = tdc(event_time, data_CHD),
  data_stroke        = tdc(event_time, data_stroke),
  data_DM            = tdc(event_time, data_DM),
  data_jiakang       = tdc(event_time, data_jiakang),
  data_COPD          = tdc(event_time, data_COPD),
  data_CKD           = tdc(event_time, data_CKD),
  data_parkin        = tdc(event_time, data_parkin),
  data_depress       = tdc(event_time, data_depress)
)

td_data_ADdi <- td_data_ADdi %>%
  mutate(across(
    c(vaccine_flu, vaccine_zoster, vaccine_pneu,
      data_hypertension, data_CHD, data_stroke, data_DM,
      data_jiakang, data_COPD, data_CKD, data_parkin, data_depress),
    ~ ifelse(is.na(.), 0, .)
  )) %>%
  group_by(高峰号) %>%
  arrange(tstart) %>%
  mutate(
    vaccine_zoster     = ifelse(cummax(vaccine_zoster) == 1, 1, vaccine_zoster),
    vaccine_pneu       = ifelse(cummax(vaccine_pneu) == 1, 1, vaccine_pneu),
    data_hypertension  = ifelse(cummax(data_hypertension) == 1, 1, data_hypertension),
    data_CHD           = ifelse(cummax(data_CHD) == 1, 1, data_CHD),
    data_stroke        = ifelse(cummax(data_stroke) == 1, 1, data_stroke),
    data_DM            = ifelse(cummax(data_DM) == 1, 1, data_DM),
    data_jiakang       = ifelse(cummax(data_jiakang) == 1, 1, data_jiakang),
    data_COPD          = ifelse(cummax(data_COPD) == 1, 1, data_COPD),
    data_CKD           = ifelse(cummax(data_CKD) == 1, 1, data_CKD),
    data_parkin        = ifelse(cummax(data_parkin) == 1, 1, data_parkin),
    data_depress       = ifelse(cummax(data_depress) == 1, 1, data_depress)
  ) %>%
  ungroup() %>%
  arrange(高峰号, tstart)

td_data_ADdi <- td_data_ADdi %>%
  select(-all_of(disease_columns))

# NOTE: As for the dementia model, the AD disease model is reported as
# THREE Cox models (base + two sensitivity analyses; see 5.1 for the design).

# Base model (primary analysis)
cox_model_ADdi_base <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    data_hypertension + data_DM +
    data_jiakang + data_COPD + data_CKD +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_ADdi
)
summary(cox_model_ADdi_base)

# Sensitivity analysis 1: additionally adjusted for CHD and stroke
cox_model_ADdi_sens_cvd <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    data_hypertension + data_DM +
    data_jiakang + data_COPD + data_CKD +
    data_CHD + data_stroke +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_ADdi
)
summary(cox_model_ADdi_sens_cvd)

# Sensitivity analysis 2: additionally adjusted for Parkinson's disease,
# depression, and family history of PD/AD/depression
cox_model_ADdi_sens_neuro <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    data_hypertension + data_DM +
    data_jiakang + data_COPD + data_CKD +
    data_parkin + 亲属患帕金森病 + 亲属患阿尔茨海默病 + data_depress + 亲属患抑郁症 +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_ADdi
)
summary(cox_model_ADdi_sens_neuro)

# ------------------------------------------------------------------------------
# 6. Stratified analyses
# ------------------------------------------------------------------------------

## 6.1 Stratification by four age groups (40-49, 50-59, 60-69, 70-79)
age_groups <- c("40-49岁", "50-59岁", "60-69岁", "70-79岁")
results_list <- list()

for (age_grp in age_groups) {
  subset_data <- td_data_ADdi %>%   # replace with td_data_dedi as needed
    filter(`年龄分组` == age_grp)
  
  cox_model_subtype <- coxph(
    Surv(tstart, tstop, event) ~
      vaccine_pneu +
      vaccine_flu +
      vaccine_zoster +
      data_hypertension + data_DM +
      data_jiakang + data_COPD + data_CKD +
      吸烟 + 喝酒 +
      每周体锻 + 焦虑抑郁_new +
      visit_category +
      文化程度_new + 实际睡眠时间_new +
      性别 + BMI分类 + 年龄分组,
    data = subset_data
  )
  
  model_summary <- summary(cox_model_subtype)
  coef_table    <- model_summary$coefficients
  conf_int      <- model_summary$conf.int
  
  result_df <- data.frame(
    年龄段       = age_grp,
    变量         = rownames(coef_table),
    HR           = exp(coef_table[, "coef"]),
    HR_lower_95  = conf_int[, "lower .95"],
    HR_upper_95  = conf_int[, "upper .95"],
    p_value      = coef_table[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  results_list[[age_grp]] <- result_df
}
final_results <- bind_rows(results_list)

## 6.2 Stratification by age >= 60 / < 60
td_data_ADdi <- td_data_ADdi %>%
  left_join(data0_67572 %>% select(高峰号, 年龄), by = "高峰号")
td_data_dedi <- td_data_dedi %>%
  left_join(data0_67572 %>% select(高峰号, 年龄), by = "高峰号")

td_data_ADdi <- td_data_ADdi %>%
  mutate(年龄分组1 = ifelse(年龄 >= 60, "≥60岁", "<60岁"))
td_data_dedi <- td_data_dedi %>%
  mutate(年龄分组1 = ifelse(年龄 >= 60, "≥60岁", "<60岁"))

age_groups <- c("<60岁", "≥60岁")
results_list <- list()

for (age_grp in age_groups) {
  subset_data <- td_data_ADdi %>%   # replace with td_data_dedi as needed
    filter(年龄分组1 == age_grp)
  
  if (nrow(subset_data) == 0) {
    warning(paste("年龄段", age_grp, "无数据，跳过该组。"))
    next
  }
  
  cox_model_subtype <- coxph(
    Surv(tstart, tstop, event) ~
      vaccine_pneu +
      vaccine_flu +
      vaccine_zoster +
      data_hypertension + data_DM +
      data_jiakang + data_COPD + data_CKD +
      吸烟 + 喝酒 +
      每周体锻 + 焦虑抑郁_new +
      visit_category +
      文化程度_new + 实际睡眠时间_new +
      性别 + BMI分类,
    data = subset_data
  )
  
  model_summary <- summary(cox_model_subtype)
  coef_table    <- model_summary$coefficients
  conf_int      <- model_summary$conf.int
  
  result_df <- data.frame(
    年龄段       = age_grp,
    变量         = rownames(coef_table),
    HR           = conf_int[, "exp(coef)"],
    HR_lower_95  = conf_int[, "lower .95"],
    HR_upper_95  = conf_int[, "upper .95"],
    p_value      = coef_table[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  results_list[[age_grp]] <- result_df
}
final_resultsde60 <- bind_rows(results_list)

## 6.3 Stratification by anxiety/depression status (焦虑抑郁_new) among >= 60
subset_data <- td_data_dedi %>%
  filter(年龄分组1 == "≥60岁")
subset_data_AD <- td_data_ADdi %>%
  filter(年龄分组1 == "≥60岁")

results_list <- list()
性别_groups <- unique(subset_data$焦虑抑郁_new)

for (hyp_grp in 性别_groups) {
  sub_data <- subset_data %>% filter(焦虑抑郁_new == hyp_grp)
  
  cox_model <- coxph(
    Surv(tstart, tstop, event) ~
      vaccine_pneu + vaccine_flu + vaccine_zoster +
      年龄分组 +
      data_jiakang + data_COPD + data_CKD + data_hypertension +
      visit_category + data_DM + 吸烟 +
      文化程度_new + 实际睡眠时间_new + BMI分类 + 性别 + 每周体锻 +
      喝酒,
    data = sub_data
  )
  
  model_summary <- summary(cox_model)
  coef_table    <- model_summary$coefficients
  conf_int      <- model_summary$conf.int
  
  result_df <- data.frame(
    高血压状态   = ifelse(hyp_grp == "有", "是", "否"),
    变量         = rownames(coef_table),
    HR           = conf_int[, "exp(coef)"],
    HR_lower_95  = conf_int[, "lower .95"],
    HR_upper_95  = conf_int[, "upper .95"],
    p_value      = coef_table[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  results_list[[as.character(hyp_grp)]] <- result_df
}
final_results <- bind_rows(results_list)

# ------------------------------------------------------------------------------
# 7. Competing-risk analysis for death
# ------------------------------------------------------------------------------
result_deathde <- subset(data_59701, death_time == dementia_censor_date,
                         select = c(高峰号, death_time, dementia_censor_date))
result_deathad <- subset(data_59701, death_time == AD_censor_date,
                         select = 高峰号)

td_data_ADdi_death <- td_data_ADdi %>%
  group_by(高峰号) %>%
  mutate(
    event = ifelse(高峰号 %in% result_deathad$高峰号 & event == 0 & tstop == max(tstop), 2, event)
  ) %>%
  ungroup() %>%
  arrange(高峰号) %>%
  mutate(event_death = ifelse(event == 2, 1, 0))

fit_death <- coxph(
  Surv(tstart, tstop, event_death) ~
    vaccine_pneu +
    vaccine_flu +
    vaccine_zoster +
    data_hypertension + data_DM +
    data_jiakang + data_COPD + data_CKD +
    吸烟 + 喝酒 +
    每周体锻 + 焦虑抑郁_new +
    visit_category +
    文化程度_new + 实际睡眠时间_new +
    性别 + BMI分类 + 年龄分组,
  data = td_data_ADdi_death
)
summary(fit_death)

# ------------------------------------------------------------------------------
# 8. Interaction models
# ------------------------------------------------------------------------------

## 8.1 Interaction between pneumococcal vaccination and age group
fit_sex <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu * 年龄分组1 +
    vaccine_flu +
    vaccine_zoster +
    BMI分类 +
    data_jiakang + data_COPD + data_CKD +
    data_hypertension + data_DM + 吸烟 + 喝酒 + 焦虑抑郁_new +
    文化程度_new + 实际睡眠时间_new + 每周体锻 + visit_category +
    性别,
  data = td_data_dedi
)
summary(fit_sex)

fit_AD <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu * 年龄分组1 +
    vaccine_flu +
    vaccine_zoster +
    BMI分类 +
    data_jiakang + data_COPD + data_CKD +
    data_hypertension + data_DM + 吸烟 + 喝酒 + 焦虑抑郁_new +
    文化程度_new + 实际睡眠时间_new + 每周体锻 + visit_category +
    性别,
  data = td_data_ADdi
)
summary(fit_AD)

## 8.2 Interaction between pneumococcal and influenza vaccination
fit_sex <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu * vaccine_flu +
    vaccine_zoster +
    BMI分类 +
    data_jiakang + data_COPD + data_CKD +
    data_hypertension + data_DM + 吸烟 + 喝酒 + 焦虑抑郁_new +
    文化程度_new + 实际睡眠时间_new + 每周体锻 + visit_category +
    性别 + 年龄分组,
  data = subset_data
)
summary(fit_sex)

fit_AD <- coxph(
  Surv(tstart, tstop, event) ~
    vaccine_pneu * vaccine_flu +
    vaccine_zoster +
    BMI分类 +
    data_jiakang + data_COPD + data_CKD +
    data_hypertension + data_DM + 吸烟 + 喝酒 + 焦虑抑郁_new +
    文化程度_new + 实际睡眠时间_new + 每周体锻 + visit_category +
    性别 + 年龄分组,
  data = subset_data_AD
)
summary(fit_AD)

## 8.3 Unified output for vaccine x BMI interactions
run_model <- function(interaction_var) {
  formula_str <- paste0(
    "Surv(tstart, tstop, event) ~ ",
    ifelse(interaction_var == "vaccine_pneu",
           "vaccine_pneu*BMI分类 + vaccine_flu + vaccine_zoster",
           ifelse(interaction_var == "vaccine_flu",
                  "vaccine_pneu + vaccine_flu*BMI分类 + vaccine_zoster",
                  "vaccine_pneu + vaccine_flu + vaccine_zoster*BMI分类")),
    " + visit_category + data_jiakang + data_COPD + data_CKD +",
    "性别 + 年龄分组 + data_hypertension + data_DM + 喝酒 +",
    " 文化程度_new + 实际睡眠时间_new + 吸烟 + 焦虑抑郁_new + 每周体锻"
  )
  
  fit <- coxph(as.formula(formula_str), data = subset_data_AD)
  return(fit)
}

fit_pneu   <- run_model("vaccine_pneu")
fit_flu    <- run_model("vaccine_flu")
fit_zoster <- run_model("vaccine_zoster")

extract_interaction_p <- function(fit, interaction_term) {
  coef_summary <- summary(fit)$coefficients
  row_name <- grep(interaction_term, rownames(coef_summary), value = TRUE)
  if (length(row_name) == 0) return(NA)
  return(coef_summary[row_name, "Pr(>|z|)"])
}

result_table <- data.table(
  交互项 = c("vaccine_pneu × BMI分类",
          "vaccine_flu × BMI分类",
          "vaccine_zoster × BMI分类"),
  P值 = c(
    extract_interaction_p(fit_pneu, "vaccine_pneu:BMI分类"),
    extract_interaction_p(fit_flu, "vaccine_flu:BMI分类"),
    extract_interaction_p(fit_zoster, "vaccine_zoster:BMI分类")
  )
)

result_table[, 显著性 := fifelse(P值 < 0.001, "***",
                              fifelse(P值 < 0.01, "**",
                                      fifelse(P值 < 0.05, "*", "")))]

print(result_table)
summary(fit_zoster)

# ------------------------------------------------------------------------------
# 9. Session information
# ------------------------------------------------------------------------------
sessionInfo()
