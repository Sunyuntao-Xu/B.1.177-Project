library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(readr)
library(tidyr)

# Load the metadata.tsv fast
GISAID_metadata <- fread("E:/metadata_tsv_2025_02_09.tar/metadata.tsv", sep="\t", header=TRUE, showProgress=TRUE)


head(GISAID_metadata)

# Keep only "complete" sequences
GISAID_metadata <- GISAID_metadata[`Is complete?` == TRUE]

str(GISAID_metadata)
table(GISAID_metadata$`Is low coverage?`, useNA = "always")  # Show value distribution

#     TRUE     <NA> 
#  1218098 15543968 


# Remove sequences explicitly marked as low coverage
GISAID_metadata <- GISAID_metadata %>%
  filter(`Is low coverage?` != TRUE | is.na(`Is low coverage?`))

# Verify that filtering was successful
table(GISAID_metadata$`Is low coverage?`, useNA = "always")

####### Extract B.1.177 from the metadata #######

# Ensure column name is correctly referenced and exists in the dataset
B.1.177_GISAID <- GISAID_metadata %>%
  filter(
    grepl("^(B\\.1\\.177|Z|Y|W|U|V|AA)", `Pango lineage`, perl = TRUE) &  # Matches B.1.177, sub-lineages, and aliases
      `Pango lineage` != "Unassigned"  # Exclude "Unassigned" entries
  )

# Leave 170,934

####### Change date format #######
B.1.177_GISAID$`Collection date` <- as.Date(B.1.177_GISAID$`Collection date`, format="%Y-%m-%d")


# Remove rows where Collection Date is NA
B.1.177_GISAID <- B.1.177_GISAID %>%
  filter(!is.na(`Collection date`))

# Leave 165,888


####### Seperate region for analysis #######

# Ensure Location column is in character format
B.1.177_GISAID$Location <- as.character(B.1.177_GISAID$Location)

# Split Location into three new columns: Continent, Country, and Region
B.1.177_GISAID <- B.1.177_GISAID %>%
  separate(Location, into = c("Continent", "Country", "Region"), sep = " / ", extra = "drop", fill = "right")

B.1.177_GISAID <- B.1.177_GISAID %>%
  mutate(Region = ifelse(is.na(Region) | Region == Country, NA, Region)) %>%  # Keep valid region values, set others to NA
  mutate(Country = case_when(
    Country %in% c("Hong Kong", "Taiwan") ~ "China",
    Country %in% c("Canary Islands") ~ "Spain",
    Country %in% c("Guadeloupe", "Mayotte", "Reunion", "French Guiana") ~ "France",
    Country %in% c("Bermuda", "British Virgin Islands", "Gibraltar") ~ "United Kingdom",
    Country %in% c("Curacao", "Aruba", "Bonaire") ~ "Netherlands",
    Country %in% c("Kosovo") ~ "Serbia",
    TRUE ~ Country
  )) %>%
  mutate(Region = case_when(
    Country == "China" & is.na(Region) ~ "Hong Kong / Taiwan",  # Assign a meaningful region for merged areas
    Country == "Spain" & is.na(Region) ~ "Canary Islands",
    Country == "France" & is.na(Region) ~ "Overseas Territories",
    Country == "United Kingdom" & is.na(Region) ~ "Overseas Territories",
    Country == "Netherlands" & is.na(Region) ~ "Dutch Caribbean",
    Country == "Serbia" & is.na(Region) ~ "Kosovo",
    TRUE ~ Region  # Keep original region if it exists
  ))

####### Load the earliest tip from UShER #######

# Define the file path (same location where you saved it)
input_file <- "C:/Users/xusun/Desktop/Phylogenetic Project/Codes/Earliest_Tips_Metadata.csv"

# Load the CSV file into a dataframe
earliest_tips_metadata <- read.csv(input_file, stringsAsFactors = FALSE)

# Display the dataframe
library(ace_tools)
ace_tools::display_dataframe_to_user(name = "Loaded Earliest Tips Metadata", dataframe = earliest_tips_metadata)

# Check the first few rows to confirm it loaded correctly
head(earliest_tips_metadata)


####### Create column with identical name format to align #######
# Create a new 'strain' column without the 'hCoV-19/' prefix
B.1.177_GISAID <- B.1.177_GISAID %>%
  mutate(strain = gsub("^hCoV-19/", "", `Virus name`))


# Ensure both dataframes have a cleaned 'strain' column
common_strains <- intersect(earliest_tips_metadata$strain, B.1.177_GISAID$strain)

# Extract rows where strains match in earliest_tips_metadata
matched_earliest <- earliest_tips_metadata %>%
  filter(strain %in% common_strains)




# Extract rows where strains match in B.1.177_GISAID
matched_gisaid <- B.1.177_GISAID %>%
  filter(strain %in% common_strains)

# Merge the matched data for easier comparison
cross_matched_data <- left_join(matched_earliest, matched_gisaid, by = "strain")

####### Extract the unmatched for alignment reference #######
# Find unmatched strains
unmatched_earliest <- earliest_tips_metadata %>%
  filter(!strain %in% common_strains)

# Extract the 'genbank_accession' column from unmatched entries
unmatched_accessions <- unmatched_earliest %>%
  select(strain, genbank_accession)


####### GISAID dataset information check #######
# Check sequence length distribution

B.1.177_GISAID$`Sequence length` <- as.numeric(B.1.177_GISAID$`Sequence length`)

summary(B.1.177_GISAID$`Sequence length`)

hist(B.1.177_GISAID$`Sequence length`, breaks = 50, col = "blue",
     main = "Distribution of Genome Lengths",
     xlab = "Genome Length", ylab = "Frequency")


# Check country distribution against time
# Get unique number of countries to set color palette dynamically
num_countries <- length(unique(B.1.177_GISAID$Country))
color_palette <- rainbow(num_countries)

# Plot distribution of B.1.177 through time
B.1.177_plot <- ggplot(B.1.177_GISAID, aes(x = `Collection date`, fill = Country)) +
  geom_bar(position = "stack", alpha = 0.7) +  # Stacked bars for country distribution
  scale_fill_manual(values = color_palette, name = "Country") +  # Dynamically assigned colors
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +  # Format x-axis
  labs(title = "B.1.177 Distribution Over Time", x = "Sample Time", y = "Sequence Count") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(B.1.177_plot)

B.1.177_GISAID$month <- format(B.1.177_GISAID$`Collection date`, "%Y-%m")

B.1.177_GISAID <- B.1.177_GISAID %>%
  filter(month >= "2020-01" & month <= "2021-12")

# down to 165,868

####### Extract Europe samples #######
B.1.177_GISAID_Europe <- B.1.177_GISAID %>%
  filter(Continent == "Europe")

# down to 164,759

table(B.1.177_GISAID_Europe$Country)


####### Extract names for fasta alignment analysis #######
# Extract just the sequence names from the 'Virus name' column
# Replace spaces with underscores in the 'Virus name' column

GISAID_selected_sequences <- B.1.177_GISAID_Europe$`Virus name`


# Save to a text file for later filtering
write_lines(GISAID_selected_sequences, "E:/GISAID_selected_sequences.txt", sep = "\n")


####### Exclude countries with too small samples #######

# Define exclusion threshold
exclusion_threshold <- 20  # Countries with <20 samples will be excluded

# Identify countries to exclude
excluded_countries <- names(which(table(B.1.177_GISAID_Europe$Country) < exclusion_threshold))

# Filter dataset to remove these countries
B.1.177_GISAID_Europe <- B.1.177_GISAID_Europe %>%
  filter(!Country %in% excluded_countries)

# down to 164,720

# Use a better color palette from RColorBrewer (easier to distinguish than rainbow)
num_countries_europe <- length(unique(B.1.177_GISAID_Europe$Country))

# Choose an appropriate palette dynamically (Set3 is good for categorical variables)
color_palette_europe <- colorRampPalette(brewer.pal(12, "Set3"))(num_countries_europe)

B.1.177_GISAID_Europe <- B.1.177_GISAID_Europe %>%
  mutate(Week = floor_date(`Collection date`, unit = "week"))  # Adjust to 2-week bins

# Improve the plot with readable colors, labels, and formatting
B.1.177_Europe_plot <- ggplot(B.1.177_GISAID_Europe, aes(x = `Collection date`, fill = Country)) +
  geom_bar(position = "stack", alpha = 1, color = NA, linewidth = 0.0001) +  # FIX: Replace size with linewidth
  scale_fill_manual(values = color_palette_europe, name = "Country") +  
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m") +  
  labs(title = "B.1.177 Distribution Over Time in Europe", x = "Sample Time", y = "Sequence Count") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Display the improved plot
print(B.1.177_Europe_plot)



# Plot in week

B.1.177_Europe_plot_month <- ggplot(B.1.177_GISAID_Europe, aes(x = Week, fill = Country)) +
  geom_bar(position = "stack", alpha = 0.85, color = NA) +  # Remove black outline
  scale_fill_manual(values = color_palette_europe, name = "Country") +  
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +  # Set x-axis to months
  labs(title = "B.1.177 Distribution Over Time in Europe", x = "Sample Time (Monthly Binned)", y = "Sequence Count") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(B.1.177_Europe_plot_month)




####### Downsample for tree running and alignment #######
# Set target sample size
target_samples <- 135  # Between 120-150

# Check unique country counts
num_countries <- length(unique(B.1.177_GISAID_Europe$Country))

# Stratified sampling: Ensure samples come from different countries and months
set.seed(1234)  # For reproducibility
sampled_data <- B.1.177_GISAID_Europe %>%
  group_by(Country, month = format(`Collection date`, "%Y-%m")) %>%  # Stratify by country & month
  slice_sample(n = max(1, round(target_samples / num_countries))) %>%  # Proportional per country
  ungroup() %>%
  sample_n(min(n(), target_samples))  # Final adjustment to hit target range


# Print sample count breakdown
print(table(sampled_data$Country))
print(paste("Total selected:", nrow(sampled_data)))


####### Load Masking Sites #######
cleaning_site <- read.csv("C:/Users/xusun/Desktop/Phylogenetic Project/Codes/cleaning_site_masking.csv")

# Keep only rows with FILTER == "mask"
mask_only <- cleaning_site[cleaning_site$FILTER == "mask", ]

# Minus 54 for the row
# Expand ranges like "1-55" into a sequence of numbers
expand_pos <- function(pos_str) {
  if (grepl("-", pos_str)) {
    range_vals <- as.integer(unlist(strsplit(pos_str, "-")))
    return(seq(range_vals[1], range_vals[2]))
  } else {
    return(as.integer(pos_str))
  }
}

# Apply the function to every row of POS and unnest
mask_expanded <- mask_only %>%
  rowwise() %>%
  mutate(POS_expanded = list(expand_pos(POS))) %>%
  unnest(POS_expanded) %>%
  ungroup()

# Subtract 54 from all positions
mask_expanded <- mask_expanded %>%
  mutate(POS_final = POS_expanded - 54)

# Keep only relevant columns, if needed
mask_final <- mask_expanded %>%
  select(POS_final, REF, ALT, FILTER, EXC, GENE, AA_POS, AA_REF, AA_ALT)

# Filter out negative positions
mask_final <- mask_expanded %>%
  mutate(POS_final = POS_expanded - 54) %>%
  filter(POS_final > 0) %>%
  select(POS_final, REF, ALT, FILTER, EXC, GENE, AA_POS, AA_REF, AA_ALT)

# Manually add the new sites
extra_sites <- data.frame(
  POS_final = c(1095, 5575, 6797, 7274, 13893, 28041, 29308),
  REF = c("G", "G", "A", "G", "A", "A", "C"),
  ALT = c("T", NA, NA, NA, "T", NA, NA),
  FILTER = "mask",
  EXC = NA,
  GENE = NA,
  AA_POS = NA,
  AA_REF = NA,
  AA_ALT = NA
)

# Bind to the existing mask_final
mask_final_combined <- bind_rows(mask_final, extra_sites)

mask_final_combined <- mask_final_combined %>%
  arrange(POS_final)

writeLines(paste(mask_final_combined$POS_final, collapse = ","), "E:/positions_comma.txt")




# Optional: Save it to a new CSV file
write.csv(mask_only, "C:/Users/xusun/Desktop/Phylogenetic Project/Codes/cleaning_site_mask_only.csv", row.names = FALSE)



####### Create a compared list of wanted B.1.177 #######
exclude_ids <- readLines("E:/B.1.177_exclude_ids.txt", encoding = "UTF-8")
exclude_ids <- trimws(exclude_ids)
exclude_ids <- exclude_ids[exclude_ids != ""]


# Replace all spaces with underscores to match the FASTA headers
GISAID_selected_sequences_nospace <- gsub(" ", "_", GISAID_selected_sequences)

desired_ids <- setdiff(GISAID_selected_sequences_nospace, exclude_ids)


writeLines(desired_ids, "E:/B.1.177_desired_ids.txt")

####### Identify identical sequences #######
# Load duplicated cluster info from seqkit
identical_clusters <- read_delim("E:/identical.detail.txt", delim = "\t", col_names = FALSE)


# Clean and explode cluster members
identical_clusters_df <- identical_clusters %>%
  rename(ClusterSize = X1, IDs = X2) %>%
  separate_rows(IDs, sep = ",\\s*") %>%
  rename(ID = IDs)

# Join identical clusters to metadata
cluster_meta <- identical_clusters_df %>%
  left_join(B.1.177_GISAID_Europe, by = c("ID" = "Virus name"))

# Group by sequence identity + country + date
redundant_ids <- cluster_meta %>%
  group_by(ClusterSize, Country, `Collection date`) %>%  # column name might be "Collection.date"
  filter(n() > 1) %>%
  slice(-1) %>%  # keep one per group, discard the rest
  ungroup()

# Write list of IDs to exclude
write_lines(redundant_ids$ID, "E:/redundant_same_country_date_ids.txt")

