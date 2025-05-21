#' Function to fit generalised additive model to Covid-19 cases/hospitalisation data
#' Adapted by KD (May 2022) from sc2growth.R downloaded from 
#' https://gist.github.com/emvolz/3d858e03f3780b953067901121e3c7c0 on 4 May 2022
#' @param cases_df Dataframe containing case/hospitalisation data with columns: Date,cases,time,wday.
#' @param GAM_smooth_function Code for smoothing function in generalised 
#'    additive model from mgcv package. See help(smooth.terms) for descriptions.
#'    Options are: "tp","ts","ds","cr","cs","sos","ps","cp","re","gp" or "so". 
#' @param deg_free_k Degrees of freedom used for time element in GAM model fit.
#'    See help(choose.k) for further information.
#' @return m Generalised additive model (GAM)

gam_fitting <- function( cases_df, GAM_smooth_function = "cc", deg_free_k = 30 ) 
{
  library(lubridate)
  library(mgcv)
  library(nlme)
  
  # Fit generalised additive model (GAM)
  m <- mgcv::gam( cases ~ s( time, bs=GAM_smooth_function, k = deg_free_k) + s( wday, bs="cc", k = 7), data = cases_df ) #family = quasipoisson,
  #m <- mgcv::gam( cases ~ s( time, bs=GAM_smooth_function, k = deg_free_k ), family = quasipoisson, data = cases_df )
  
  # Plot graphs to check quality of model fit
  par(mfrow=c(2,2))
  #gam.check(m)
  #summary(m)
  return(m) # Return from function
}



library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(mgcv)

Sys.setenv(LANGUAGE = "en")
# Try one of these:
Sys.setlocale("LC_TIME", "English_United States.1252")
Sys.setlocale("LC_ALL", "English_United States.1252")


# Read the COG-UK metadata file
COGUK <- read.csv("E:/COG-UK/cog_global_2021-07-25_reduced.csv", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

# Read the WHO global COVID-19 data file
WHO <- read.csv("E:/COG-UK/WHO-COVID-19-global-data.csv", 
                header = TRUE, 
                stringsAsFactors = FALSE)

# For the COGUK dataset, filter where the 'country' column equals "UK"
COGUK_UK <- subset(COGUK, country == "UK")

# For the WHO dataset, filter where the 'Country' column equals "United Kingdom"
WHO_UK <- subset(WHO, Country == "United Kingdom of Great Britain and Northern Ireland")


# Convert the sample_date column to a Date type
COGUK_UK <- COGUK_UK %>%
  mutate(sample_date = as.Date(sample_date))


# Summarise by day:
COGUK_daily_summary <- COGUK_UK %>%
  group_by(sample_date) %>%
  summarise(
    total_seq = n(),
    b177_seq = sum(grepl("^(B\\.1\\.177|Z|Y|W|U|V|AA)", lineage, perl = TRUE) & lineage != "Unassigned")
  )

# View the result
head(COGUK_daily_summary)


# Pivot the data to long format for plotting
COGUK_daily_long <- COGUK_daily_summary %>%
  pivot_longer(
    cols = c(total_seq, b177_seq),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(Category = recode(Category,
                           total_seq = "Total Sequences",
                           b177_seq = "B.1.177 Sequences"))

ggplot(COGUK_daily_long, aes(x = sample_date, y = Count, color = Category)) +
  # Thinner lines
  geom_line(size = 0.5) +
  labs(
    title = "Daily Sequenced Cases and B.1.177 Cases Over Time",
    x = "Date",
    y = "Count",
    color = "Category"
  ) +
  # Adjust how dates are displayed on the x-axis
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  # Use a minimal theme and rotate x-axis labels
  theme_minimal(base_family = "sans") +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


####### Fitting into Weeks and merge COGUK and WHO #######

WHO_weekly <- WHO_UK %>%
  rename(Date = Date_reported, cases = New_cases) %>%  # rename for consistency
  mutate(Date = as.Date(Date),
         week = floor_date(Date, unit = "week", week_start = 1)) %>%
  group_by(week) %>%
  summarise(WHO_total_cases = sum(cases, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(time = as.numeric(week))

# 2. Aggregate COG-UK Data by Week
COGUK_weekly <- COGUK_UK %>%
  mutate(sample_date = as.Date(sample_date),
         week = floor_date(sample_date, unit = "week", week_start = 1)) %>%
  group_by(week) %>%
  summarise(
    total_seq = n(),
    b177_seq = sum(grepl("^(B\\.1\\.177|Z|Y|W|U|V|AA)", lineage, perl = TRUE) & lineage != "Unassigned")
  ) %>%
  ungroup() %>%
  mutate(time = as.numeric(week),
         p_b177 = b177_seq / total_seq)

# 3. Select correct period for smoothing

# Find first week with B.1.177 sequences
first_b177_week <- COGUK_weekly %>% 
  filter(b177_seq > 0) %>% 
  slice_min(week) %>% 
  pull(week)

# Filter both datasets to overlapping period
WHO_filtered <- WHO_weekly %>% 
  filter(week >= first_b177_week, week <= as.Date("2021-07-19"))

COGUK_filtered <- COGUK_weekly %>% 
  filter(week >= first_b177_week, week <= as.Date("2021-07-19"))

# 4. Merge weekly information

merged_weekly <- inner_join(
  WHO_filtered %>% select(week, WHO_total_cases),
  COGUK_filtered %>% select(week, total_seq, b177_seq),
  by = "week"
) %>% 
  mutate(
    time = as.numeric(week - min(week)),  # Reset time index from 0
    p_b177 = ifelse(total_seq > 0, b177_seq/total_seq, NA)
  )

# 5.Fit Robust GAMs
# Fit beta regression for proportions (avoids <0 values)
p_b177_gam <- gam(
  p_b177 ~ s(time, bs = "tp", k = 30),
  data = merged_weekly,
  family = betar(),
  method = "REML"
)

# Fit quasi-poisson model for sequencing effort
gam_seq <- gam(
  total_seq ~ s(time, bs = "tp", k = 30) + offset(log(WHO_total_cases)),
  data = merged_weekly,
  family = quasipoisson(link = "log")
)

# Calculate sequencing proportion (sequences per case)
merged_weekly <- merged_weekly %>%
  mutate(
    seq_prop = exp(predict(gam_seq, type = "link") - log(WHO_total_cases))
  )

merged_weekly <- merged_weekly %>%
  mutate(
    p_b177_adj = (b177_seq + 0.5) / (total_seq + 1)  # Avoids 0/0 or 0/N
  )

# Refit beta regression
gam_beta <- gam(
  p_b177_adj ~ s(time, bs = "tp", k = 30),
  data = merged_weekly,
  family = betar(link = "logit"),
  method = "REML"
)

merged_weekly <- merged_weekly %>%
  mutate(
    p_b177_pred = predict(gam_beta, type = "response"),
    b177_cases = total_seq*p_b177_pred/seq_prop
  )



merged_weekly <- merged_weekly %>%
  filter(
    total_seq > 0,
    seq_prop <= 1  # Sequencing can't exceed 100% of cases
  )

ggplot(merged_weekly, aes(x = week)) +
  # 1) WHO total cases as a bar (fill)
  geom_col(aes(y = WHO_total_cases, fill = "WHO total cases"), width = 5) +
  # 2) Estimated B.1.177 as a line (color)
  geom_line(aes(y = b177_cases, color = "Estimated B.1.177 Cases"), linewidth = 1) +
  # 3) Sequenced B.1.177 as points (color)
  geom_point(aes(y = b177_seq, color = "Sequenced B.1.177 Cases"), size = 2) +
  
  # Separate legends: fill scale for the bar, color scale for line/points
  scale_fill_manual(
    name = "WHO Data", 
    values = c("WHO total cases" = "gray90")
  ) +
  scale_color_manual(
    name = "B.1.177 Data",
    values = c("Estimated B.1.177 Cases" = "red",
               "Sequenced B.1.177 Cases" = "black")
  ) +
  
  labs(
    title = "Corrected B.1.177 Case Estimates",
    x = "Week",
    y = "Cases"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





par(mfrow=c(2,2))
gam.check(gam_beta)
summary(gam_beta)

gam.check(gam_seq)
summary(gam_seq)

gam.check(p_b177_gam)
summary(p_b177_gam)



####### Estimate with daily cases #######

UKHSA <- read.csv("C:/Users/xusun/Desktop/Phylogenetic Project/Codes/UKHSA.csv")

UKHSA <- UKHSA %>%
  arrange(as.Date(date))



####### B,1,177 case before smoothing #######

### 1. Process UKHSA Daily Data
UKHSA_daily <- UKHSA %>%
  mutate(date = as.Date(date)) %>%        # ensure date is Date type
  arrange(date) %>%
  filter(date >= as.Date("2020-03-23")) %>%  # filter out dates before 2020-03-23
  rename(daily_cases = newCasesBySpecimenDate) %>%
  mutate(time = as.numeric(date - min(date)))  # time in days from first date (post-filter)

### 2. Process COG‑UK Daily Data
COGUK_daily <- COGUK_UK %>%
  mutate(sample_date = as.Date(sample_date)) %>%
  filter(sample_date >= as.Date("2020-03-23")) %>%  # filter out dates before 2020-03-23
  group_by(sample_date) %>%
  summarise(
    total_seq = n(),
    b177_seq = sum(grepl("^(B\\.1\\.177|Z|Y|W|U|V|AA)", lineage, perl = TRUE) & 
                     lineage != "Unassigned")
  ) %>%
  ungroup() %>%
  rename(date = sample_date) %>%
  mutate(time = as.numeric(date - min(date)),  # reset time index starting at filtered first date
         p_b177 = ifelse(total_seq > 0, b177_seq / total_seq, NA))

### 3. Merge Daily Datasets (using the full date range, post-filter)
merged_daily <- inner_join(
  UKHSA_daily %>% select(date, daily_cases),
  COGUK_daily %>% select(date, total_seq, b177_seq, p_b177),
  by = "date"
) %>% 
  mutate(time = as.numeric(date - min(date)))  # reset time index starting at first filtered date

### 4. Adjust p_b177 and Compute b177 Cases
merged_daily <- merged_daily %>%
  mutate(
    # Adjust p_b177 to avoid exact 0 or 1 values.
    p_b177_adj = ifelse(p_b177 == 0, 0.001,
                        ifelse(p_b177 == 1, 0.999, p_b177)),
    # Compute raw estimated b177 cases using adjusted proportion.
    b177_cases_adj = daily_cases * p_b177_adj
  )

# For smoothing, ensure that b177_cases_adj is strictly positive.
merged_daily <- merged_daily %>%
  mutate(b177_cases_adj_pos = ifelse(b177_cases_adj <= 0, 0.001, b177_cases_adj))

### 5. Fit a GAM to Smooth the Adjusted b177 Cases (if you are using GAM)
gam_b177 <- gam(
  b177_cases_adj_pos ~ s(time, bs = "tp", k = 40),
  data = merged_daily,
  family = Gamma(link = "log"),
  method = "REML"
)


gam.check(gam_b177)


merged_daily <- merged_daily %>%
  mutate(b177_cases_est = predict(gam_b177, newdata = merged_daily, type = "response"))


####### Try a poisson model as well #######

gam_b177_poisson <- gam(
  round(b177_cases_adj_pos) ~ s(time, bs = "tp", k = 40),
  data = merged_daily,
  family = poisson(link = "log"),
  method = "REML"
)

gam.check(gam_b177_poisson)
merged_daily <- merged_daily %>%
  mutate(b177_cases_est_poisson = predict(gam_b177_poisson, newdata = merged_daily, type = "response"))

# Constrain the estimated values so they do not exceed daily_cases
merged_daily <- merged_daily %>%
  mutate(b177_cases_est_poisson = pmin(b177_cases_est_poisson, daily_cases))


# Constrain the smoothed estimate so it does not exceed daily_cases.
merged_daily <- merged_daily %>%
  mutate(b177_cases_est = pmin(b177_cases_est, daily_cases))

### 6. Plot the Daily Estimates
smoothing_plot <- ggplot(merged_daily, aes(x = date)) +
  # UKHSA daily cases as bars
  geom_col(aes(y = daily_cases, fill = "UKHSA daily cases"), width = 1) +
  # Smoothed estimated B.1.177 cases as a red line
  geom_line(aes(y = b177_cases_est, color = "Smoothed Estimated B.1.177 Cases"), size = 1) +
  # Unsmoothed B.1.177 cases as a blue line
  geom_line(aes(y = b177_cases_adj, color = "Unsmoothed B.1.177 Cases"), size = 0.5) +
  # Raw sequenced B.1.177 cases as black points
  geom_point(aes(y = b177_seq, color = "Sequenced B.1.177 Cases"), size = 0.5) +
  scale_fill_manual(
    name = "UKHSA Data", 
    values = c("UKHSA daily cases" = "gray90")
  ) +
  scale_color_manual(
    name = "COG-UK Data",
    values = c("Smoothed Estimated B.1.177 Cases" = "red",
               "Unsmoothed B.1.177 Cases" = "blue",
               "Sequenced B.1.177 Cases" = "black")
  ) +
  labs(
    title = "B.1.177 Case Estimates via Smoothing (Daily)",
    x = "Date",
    y = "Cases"
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(smoothing_plot)


head(merged_daily)



# Create the poisson plot

poisson_plot <- ggplot(merged_daily, aes(x = date)) +
  # UKHSA daily cases as bars
  geom_col(aes(y = daily_cases, fill = "UKHSA daily cases"), width = 1) +
  # Poisson smoothed estimated B.1.177 cases as an orange (or your chosen color) line
  geom_line(aes(y = b177_cases_est_poisson, color = "Smoothed Poisson Estimate"), size = 1) +
  # Unsmoothed B.1.177 cases as a blue line
  geom_line(aes(y = b177_cases_adj, color = "Unsmoothed B.1.177 Cases"), size = 0.5) +
  # Raw sequenced B.1.177 cases as black points
  geom_point(aes(y = b177_seq, color = "Sequenced B.1.177 Cases"), size = 0.5) +
  scale_fill_manual(
    name = "UKHSA Data", 
    values = c("UKHSA daily cases" = "gray90")
  ) +
  scale_color_manual(
    name = "COG-UK Data",
    values = c("Smoothed Poisson Estimate" = "orange",
               "Unsmoothed B.1.177 Cases" = "blue",
               "Sequenced B.1.177 Cases" = "black")
  ) +
  labs(
#    title = "B.1.177 Case Estimates via Poisson Smoothing (Daily)",
    x = "Date",
    y = "Cases"
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(poisson_plot)




####### Create a cumulative infection plot #######
# Create cumulative columns
merged_daily <- merged_daily %>%
  arrange(date) %>%
  mutate(
    cum_unsmoothed = cumsum(b177_cases_adj),
    cum_smoothed   = cumsum(b177_cases_est)
  )

# Plot cumulative unsmoothed and smoothed B.1.177 cases
cumulative_plot <- ggplot(merged_daily, aes(x = date)) +
  geom_line(aes(y = cum_unsmoothed, color = "Cumulative Unsmooth B.1.177 Cases"), size = 1) +
  geom_line(aes(y = cum_smoothed, color = "Cumulative Smoothed B.1.177 Cases"), size = 1) +
  labs(
    title = "Cumulative B.1.177 Case Estimates (Proxy for I Compartment)",
    x = "Date",
    y = "Cumulative Cases"
  ) +
  scale_color_manual(
    name = "Cumulative Data",
    values = c("Cumulative Unsmooth B.1.177 Cases" = "blue",
               "Cumulative Smoothed B.1.177 Cases" = "red")
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(cumulative_plot)

# Plot a poisson one as well

merged_daily <- merged_daily %>%
  arrange(date) %>%
  mutate(
    cum_smoothed_poisson = cumsum(b177_cases_est_poisson)
  )

# Plot cumulative unsmoothed and Poisson-smoothed B.1.177 cases
cumulative_plot_poisson <- ggplot(merged_daily, aes(x = date)) +
  geom_line(aes(y = cum_unsmoothed, color = "Cumulative Unsmooth B.1.177 Cases"), size = 1) +
  geom_line(aes(y = cum_smoothed_poisson, color = "Cumulative Poisson-Smoothed B.1.177 Cases"), size = 1) +
  labs(
    title = "Cumulative B.1.177 Case Estimates (Proxy for I Compartment) - Poisson",
    x = "Date",
    y = "Cumulative Cases"
  ) +
  scale_color_manual(
    name = "Cumulative Data",
    values = c("Cumulative Unsmooth B.1.177 Cases" = "blue",
               "Cumulative Poisson-Smoothed B.1.177 Cases" = "orange")
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(cumulative_plot_poisson)

####### Finite Differences Approach #######
offset <- 0.001
merged_daily <- merged_daily %>%
  arrange(time) %>%
  mutate(b177_cases_offset = b177_cases_est + offset)

### Finite Differences Approach (for reference) ###

# Compute finite differences using a central difference scheme:
n <- nrow(merged_daily)
dI_dt <- numeric(n)

for(i in 1:n) {
  if(i == 1) {
    dI_dt[i] <- (merged_daily$b177_cases_offset[i+1] - merged_daily$b177_cases_offset[i]) /
      (merged_daily$time[i+1] - merged_daily$time[i])
  } else if(i == n) {
    dI_dt[i] <- (merged_daily$b177_cases_offset[i] - merged_daily$b177_cases_offset[i-1]) /
      (merged_daily$time[i] - merged_daily$time[i-1])
  } else {
    dI_dt[i] <- (merged_daily$b177_cases_offset[i+1] - merged_daily$b177_cases_offset[i-1]) /
      (merged_daily$time[i+1] - merged_daily$time[i-1])
  }
}

merged_daily <- merged_daily %>%
  mutate(deriv = dI_dt)

# Set the recovery rate gamma to 1/7
gamma <- 1/7

# Calculate the instantaneous growth rate r and then beta:
merged_daily <- merged_daily %>%
  mutate(r = deriv / b177_cases_offset,
         beta_est_fd = r + gamma)

print(merged_daily %>% select(time, b177_cases_est, b177_cases_offset, deriv, r, beta_est_fd))


### Smoothing Spline Approach on Log-Transformed Data using log1p ###
# Use log1p transformation for stability (i.e., log(1 + I)) which is less volatile for I near zero.
merged_daily <- merged_daily %>%
  mutate(logI = log1p(b177_cases_offset))

# Fit a smoothing spline to the log1p-transformed infection counts
spline_fit <- smooth.spline(merged_daily$time, merged_daily$logI, spar = 0.7)

# Predict the first derivative (d(log1p(I))/dt) and second derivative if needed
first_deriv <- predict(spline_fit, merged_daily$time, deriv = 1)$y
second_deriv <- predict(spline_fit, merged_daily$time, deriv = 2)$y  # may be useful for diagnostics

# Under the assumption S/N ~ 1, recall that for the original model:
#    d(logI)/dt = β - γ.
# With the log1p transformation, we have:
#    d/dt log1p(I) = I'/(1 + I).
# For very small I, (1+I) ~ 1, so the derivative is similar to that of log(I).
# Thus, as an approximation we set:
merged_daily <- merged_daily %>%
  mutate(d_logI = first_deriv,      # instantaneous growth rate r(t) from the transformed data
         d2_logI = second_deriv,     # curvature (not used directly to obtain beta)
         beta_est_spline = d_logI + gamma)  # estimated beta

# Ensure 'date' is a Date class for plotting:
merged_daily <- merged_daily %>% mutate(date = as.Date(date))

# Quick plot of the estimated beta (spline approach) over time with month ticks on the x-axis:
ggplot(merged_daily, aes(x = date, y = beta_est_spline)) +
  geom_line() +
  geom_hline(yintercept = 1/7, color = "red", linetype = "dashed") +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  labs(title = "Estimated Transmission Rate (beta) Over Time",
       x = "Date", y = expression(beta(t))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



####### Try with cumulative smoothed case for estimating beta #######

# Use the pure cumulative smoothed values and compute its log
merged_daily <- merged_daily %>%
  arrange(time) %>%
  mutate(logCum = log(cum_smoothed))

# Fit a smoothing spline to the log-transformed cumulative data
spline_fit_cum <- smooth.spline(merged_daily$time, merged_daily$logCum, spar = 0.7)

# Predict the first derivative (d(logCum)/dt) and second derivative (for diagnostics)
first_deriv_cum <- predict(spline_fit_cum, merged_daily$time, deriv = 1)$y
second_deriv_cum <- predict(spline_fit_cum, merged_daily$time, deriv = 2)$y

# Set gamma (recovery rate) to 1/7
gamma <- 1/7

# Under the assumption S/N ~ 1, we have: d(logCum)/dt ≈ β - γ,
# so we can estimate β as: beta_est_cum = d(logCum)/dt + γ.
merged_daily <- merged_daily %>%
  mutate(d_logCum = first_deriv_cum,
         d2_logCum = second_deriv_cum,  # For diagnostic use
         beta_est_cum = d_logCum + gamma)

# Ensure 'date' is a Date class
merged_daily <- merged_daily %>% mutate(date = as.Date(date))

# Plot the estimated β (from cumulative data) with monthly x-axis ticks,
# and add a red dashed horizontal line at y = 1/7 (the recovery rate)
ggplot(merged_daily, aes(x = date, y = beta_est_cum)) +
  geom_line() +
  geom_hline(yintercept = 1/7, color = "red", linetype = "dashed") +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  labs(title = "Estimated Transmission Rate (β) Over Time Using Cumulative Data",
       x = "Date", y = expression(beta(t))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





####### Try with Cumulative Poisson-Smoothed Cases for Estimating β #######


# 1. Compute the log of the cumulative Poisson-smoothed cases
merged_daily <- merged_daily %>%
  mutate(logCum_poisson = log(cum_smoothed_poisson))

# 2. Fit a smoothing spline to the log-transformed cumulative Poisson-smoothed data
spline_fit_cum_poisson <- smooth.spline(merged_daily$time, merged_daily$logCum_poisson, spar = 0.7)

# 3. Predict the first derivative (and optionally the second derivative for diagnostics)
first_deriv_cum_poisson <- predict(spline_fit_cum_poisson, merged_daily$time, deriv = 1)$y
second_deriv_cum_poisson <- predict(spline_fit_cum_poisson, merged_daily$time, deriv = 2)$y  # optional

# 4. Set the recovery rate γ (fixed at 1/7)
gamma <- 1/7

# 5. Estimate β(t) using the derivative of the log cumulative cases:
#    d/dt(logCum_poisson) ≈ β(t) - γ, so β(t) = d/dt(logCum_poisson) + γ.
merged_daily <- merged_daily %>%
  mutate(d_logCum_poisson = first_deriv_cum_poisson,
         d2_logCum_poisson = second_deriv_cum_poisson,  # For diagnostic purposes
         beta_est_cum_poisson = d_logCum_poisson + gamma)

# 6. Ensure 'date' is a Date object
merged_daily <- merged_daily %>% mutate(date = as.Date(date))

# 7. Plot the estimated β(t) using Poisson cumulative smoothing:
poisson_cum_beta_plot <- ggplot(merged_daily, aes(x = date, y = beta_est_cum_poisson)) +
  geom_line() +
  geom_hline(yintercept = 1/7, color = "red", linetype = "dashed") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(#title = "Estimated Transmission Rate (β) via Poisson Cumulative Smoothing",
       x = "Date", y = expression(beta(t))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(poisson_cum_beta_plot)






####### Try of anchor time works #######

# Get the anchor (the first time point)
anchor_time <- min(merged_daily$time)
anchor_value <- merged_daily$logCum[merged_daily$time == anchor_time]

# Create a data frame for the anchor
anchor_df <- data.frame(time = anchor_time, logCum = anchor_value)

# Combine the anchor with the original data
gam_data <- bind_rows(merged_daily %>% select(time, logCum), anchor_df)

# Create a weight vector: use weight 1 for original data and a very high weight for the anchor.
weights <- c(rep(1, nrow(merged_daily)), 1000)

# Fit a gam that smoothes logCum as a function of time, now anchored at the beginning.
gam_fit <- gam(logCum ~ s(time, bs = "tp"), 
               data = gam_data, 
               weights = weights, 
               family = gaussian(), 
               method = "REML")

# Predict the fitted values over the original data's time points
merged_daily <- merged_daily %>%
  mutate(fitted_logCum = predict(gam_fit, newdata = data.frame(time = time)))

# Now compute the first derivative using finite differences on the fitted spline;
# Alternatively, you could use the mgcv::predict.gam(..., type="lpmatrix") approach or backsolve the derivative.
# For simplicity, we'll use a simple finite-difference on the fitted values:
n <- nrow(merged_daily)
d_logCum_fd <- numeric(n)
time_vals <- merged_daily$time

for(i in 1:n) {
  if(i == 1){
    d_logCum_fd[i] <- (merged_daily$fitted_logCum[i+1] - merged_daily$fitted_logCum[i]) /
      (time_vals[i+1] - time_vals[i])
  } else if(i == n){
    d_logCum_fd[i] <- (merged_daily$fitted_logCum[i] - merged_daily$fitted_logCum[i-1]) /
      (time_vals[i] - time_vals[i-1])
  } else {
    d_logCum_fd[i] <- (merged_daily$fitted_logCum[i+1] - merged_daily$fitted_logCum[i-1]) /
      (time_vals[i+1] - time_vals[i-1])
  }
}

# Set gamma (recovery rate) to 1/7
gamma <- 1/7

# Compute estimated beta: under S/N ~ 1, d(logCum)/dt ~ β - γ, so:
merged_daily <- merged_daily %>%
  mutate(d_logCum = d_logCum_fd,
         beta_est_cum = d_logCum + gamma)

# Plot the estimated beta (from cumulative data) with monthly breaks on the x-axis 
# and add a red dashed horizontal line at y=1/7 (the recovery rate)
merged_daily <- merged_daily %>% mutate(date = as.Date(date))

ggplot(merged_daily, aes(x = date, y = beta_est_cum)) +
  geom_line() +
  geom_hline(yintercept = 1/7, color = "red", linetype = "dashed") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Estimated Transmission Rate (β) Over Time Using Cumulative Data (Anchored)",
       x = "Date", y = expression(beta(t))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
