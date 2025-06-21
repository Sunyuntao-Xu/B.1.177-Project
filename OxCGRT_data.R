# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# Define file paths
file_2020 <- "C:/Users/xusun/Desktop/Phylogenetic Project/OxfordTracker/OxCGRT_GBR_differentiated_withnotes_2020.csv"
file_2021 <- "C:/Users/xusun/Desktop/Phylogenetic Project/OxfordTracker/OxCGRT_GBR_differentiated_withnotes_2021.csv"

# Load both CSVs
oxcgrt_2020 <- read_csv(file_2020)
oxcgrt_2021 <- read_csv(file_2021)

# Merge them row-wise
oxcgrt_merged <- bind_rows(oxcgrt_2020, oxcgrt_2021)

# Optional: convert date column if needed
# Assuming there's a "Date" column in yyyymmdd format
oxcgrt_merged <- oxcgrt_merged %>%
  mutate(Date = as.Date(as.character(Date), format = "%Y%m%d"))

# Select only desired columns
oxcgrt_filtered <- oxcgrt_merged %>%
  select(
    CountryName,
    `RegionName`,
    Date,
    C8NV_International_travel_controls = `C8NV_International travel controls`,
    C8V_International_travel_controls = `C8V_International travel controls`,
    C8EV_International_travel_controls = `C8EV_International travel controls`,
    C8_Notes
  )


# Replace NA in Region Name with "UK (national)"
oxcgrt_filtered$RegionName[is.na(oxcgrt_filtered$RegionName)] <- "UK (national)"

# Convert to factor for easier use
oxcgrt_filtered$RegionName <- as.factor(oxcgrt_filtered$RegionName)

# Summarize
table(oxcgrt_filtered$RegionName)


# Define regions to plot
regions <- c("UK (national)", "England", "Scotland", "Wales", "Northern Ireland")

df_plot <- oxcgrt_filtered %>%
  filter(RegionName %in% regions)

# Create the faceted plot
ggplot(df_plot, aes(x = Date, y = C8EV_International_travel_controls)) +
  geom_step(direction = "hv", color = "steelblue") +
  scale_y_continuous(
    breaks = 0:4,
    labels = c(
      "0 - No restrictions",
      "1 - Screening",
      "2 - Quarantine",
      "3 - Ban (some regions)",
      "4 - Total closure"
    ),
    limits = c(0, 4)
  ) +
  facet_wrap(~ RegionName, ncol = 1, scales = "free_x") +
  labs(
    title = "C8EV: International Travel Controls Across UK Regions",
    x = "Date",
    y = "Restriction Level"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Filter UK national data
uk_nat <- oxcgrt_filtered %>%
  filter(RegionName == "UK (national)") %>%
  arrange(Date) %>%
  select(Date, C8EV_International_travel_controls)

# Detect changes in C8EV
uk_nat_changes <- uk_nat %>%
  mutate(Change = C8EV_International_travel_controls != lag(C8EV_International_travel_controls)) %>%
  filter(Change == TRUE | is.na(lag(C8EV_International_travel_controls))) %>%
  select(Date, C8EV_International_travel_controls)

# View the change points
print(uk_nat_changes)
