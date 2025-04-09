Sys.setenv(LANGUAGE = "en")

install.packages("devtools")
install.packages('bbmle')
install.packages('rjson')
install.packages('ggplot2')
install.packages('optim')
devtools::install_github("emvolz-phylodynamics/phydynR")
install_github( 'emvolz/treedater')
install.packages('migrations')
install.packages("future")
install.packages('future.apply')
install.packages("doParallel")
install.packages("foreach")

library(phydynR)
library(ape)
library(devtools)
library(rjson)
library(bbmle)
library(ggplot2)
library(treedater)
library(dplyr)
library(data.table)
library(future)
library(lubridate)
library(future.apply)
library(doParallel)
library(foreach)

# library(migrations) (Not working)

plan(multisession, workers = 14)  # Use 8 CPU cores (or set to the number of cores you want)


####### Estimating transmission with SIR model #######


# Obtain the phylogenetic tree
tree <- read.tree(system.file("extdata/sirModel0.nwk", package = "phydynR"))
plot(tree) # A tree with 75 samples (cases)

# Obtain the already created epidemiological data
epidata <- rjson::fromJSON(file=system.file("extdata/sirModel0.json", package = "phydynR"))

file.show( system.file("extdata/sirModel0.xml", package = "phydynR")) 


# Create the paramters for SIR compartment
parms_truth <- list( beta = 0.00020002, # Transmission rate
                     gamma = 1,         # Recovery rate
                     S0 = 9999,         # Initial susceptible population
                     t0 = 0 )           # Time of start

# Create sample time for 75 cases, which all sampled on day 12, and add time as name to each tip
sampleTimes <- rep(12, 75)
names(sampleTimes) <- tree$tip.label

# Create dated tree 'bdt'
bdt <- DatedTree( phylo = tree, 
                  sampleTimes = sampleTimes)
bdt

plot(bdt)

# Define birth and death events of infectious population
births <- c( I = "beta * S * I" )
deaths <- c( I = "gamma * I" )

# Define change of susceptible population
nonDemeDynamics <- c(S = "-beta * S * I")


# initial number of I and S
x0 <- c(I = 1, S = unname(parms_truth$S0))

# initial t0 (time of origin of the process)
t0 <- bdt$maxSampleTime - max(bdt$heights) - 1

# Define the change of demographic, prepare to integrate into with tree information
dm <- build.demographic.process(births = births,
                                nonDemeDynamics = nonDemeDynamics,
                                deaths = deaths,
                                parameterNames = names(parms_truth),
                                rcpp = TRUE,
                                sde = FALSE)


# Plot the demographic graph based on the defined model
show.demographic.process(demo.model = dm,
                         theta = list(beta = parms_truth$beta, gamma = parms_truth$gamma, S0 = parms_truth$S0, t0 = parms_truth$t0),
                         x0 = c(I = 1, S = parms_truth$S0),
                         t0 = 0,
                         t1 = 30)


# Calculate likelihood of observing such tree under the defined demographic changes
# Result obtain likelihood (how well model explain the tree)
print(system.time(print(phydynR::colik(tree = bdt,
                                       theta = parms_truth,
                                       demographic.process.model = dm,
                                       x0 = x0,
                                       t0 = t0,
                                       res = 1000,
                                       integrationMethod = "rk4")
)))

likelihood_result <- colik(tree = bdt,
                           theta = parms_truth,
                           demographic.process.model = dm,
                           x0 = x0,
                           t0 = t0,
                           res = 1000,
                           integrationMethod = "rk4")


# Define a function that optimise transmission rate and infectious population
# Use mll (maximum log-likelihood) for negative log-likelihood to optimise result

obj_fun <- function(lnbeta, lnI0){
  
  beta <- exp(lnbeta)
  I0 <- exp(lnI0)
  parms <- parms_truth
  parms$beta <- beta
  x0 <- c(I = unname(I0), S = unname(parms$S0) )
  
  mll <- -phydynR::colik(tree = bdt,
                         theta = parms_truth,
                         demographic.process.model = dm,
                         x0 = x0,
                         t0 = t0,
                         res = 1000,
                         integrationMethod = "rk4")
  
  print(paste(mll, beta, I0))
  mll
}


# Find the optimised beta and I value with maximum likelihood estimation
fit <- mle2(obj_fun,
            start = list(lnbeta = log(parms_truth$beta), lnI0 = log(1)),
            method = "Nelder-Mead",
            optimizer = "optim",
            control = list(trace=6, reltol=1e-8))




AIC(fit)                    # Measure model quality (lower the better)
# [1] -145.7974
logLik(fit)                 # Log-likelihood of the fitted model
# 'log Lik.' 74.89871 (df=2)
coef(fit)                   # Estimated parameter values (log-scales)
# lnbeta       lnI0 
# -8.4748155  0.1351695
exp(coef(fit))              # Convert log parameters to original scales
# lnbeta         lnI0 
# 0.0002086577 1.1447308446 

# Use this to see how is it different (how bias is estimate)
exp(coef(fit)["lnbeta"]) - parms_truth$beta
# lnbeta 
# 8.637689e-06



# Extraction of fitted parameters of transmission rate and infectious population
beta <- exp(coef(fit)["lnbeta"]) # Exponentiate transmission rate
# 0.0002086577 
I0 <- exp(coef(fit)["lnI0"])     # Exponentiate initial infectious population
# 1.144731

# Update the model parameter
parms <- parms_truth
parms$beta <- beta

# Initial condition setting
x0 <- c(I = unname(I0), S = unname(parms$S0) )
#        I           S 
# 1.144731 9999.000000 


# Simulate fitted model over:
o <- dm(parms,                     # Parameter with estimated transmission rate
        x0,                        # Initial condition
        t0,                        # Initial time
        t1 = bdt$maxSampleTime,    # End time matching maximum sapling time from the phylogenetic tree
        res = 1e3,                 # Specify resolution of simulation, number of time steps
        integrationMethod='rk4')   # Runge-Kutta 4 method used for numerical integration
o <- o[[5]]                        # Extract the fifth element from output, simulated number of infected overtime

# Simulate the real model over
otruth <- dm(parms_truth,
             x0,
             t0, 
             t1 = bdt$maxSampleTime,
             res = 1e3, 
             integrationMethod='rk4')
otruth <- otruth[[5]]

# Preparing observed data plot
rdata <- data.frame(time = epidata$t, I = epidata$I)

ggplot(rdata, aes(x = time, y = I)) +
  geom_point() +                                                # Actual number of infection
  geom_line(data = o, aes(x = time, y = I), col='blue') +       # SIR prediction combining both
  geom_line(data = otruth, aes(x = time, y = I), col='red') +   # Estimated number of infection
  theme_bw() +
  theme(element_text(size = 12)) +
  xlab("Time") +
  ylab("Number of infected")


profbeta <- profile(
  fit,
  which = "lnbeta",  # Specify the parameter to profile
  alpha = 0.05,
  std.err = 1,
  trace = TRUE,
  method = "Brent",  # Use Brent optimization method for reliability
  tol.newmin = 1
)


load( system.file("extdata/sirModel0-profbeta.RData", package = "phydynR"))


c( exp( confint( profbeta ) ), TrueVal = parms_truth$beta )
#>        2.5 %       97.5 %      TrueVal 
#> 0.0001901721 0.0002282366 0.0002000200


plot(profbeta)
abline( v = log( parms_truth$beta) , col = "red")








####### Prepare the tree file and the 120 tips tree as a try for the B.1.177 #######


tree_metadata_test <- read.csv('C:/Users/xusun/Desktop/Phylogenetic Project/Tree files/5k_treedater_metadata.csv')
head(tree_metadata_test)

B.1.177_tr_test <- read.tree('C:/Users/xusun/Desktop/Phylogenetic Project/Tree files/B.1.177_5k_treedater_alias.nwk')




# Ensure at least one sample per country is retained
set.seed(1234)  # For reproducibility
representative_samples <- tree_metadata_test %>%
  group_by(new_continent) %>%        # Group by continent (UK seperate)
  slice_sample(n = 1)                # Select one sample per continent

# Remove these representative samples from the pool to avoid double sampling
remaining_samples <- tree_metadata %>%
  filter(!tip_label %in% representative_samples$tip_label)

# Calculate total target size excluding representative samples
target_size <- 120 - nrow(representative_samples)

# Determine sampling proportion for remaining samples
sampling_proportion <- target_size / nrow(remaining_samples)

# Subsample the remaining data proportionally by month and country
set.seed(1234)
subsampled_remaining <- remaining_samples %>%
  group_by(month, new_continent) %>%
  sample_frac(size = sampling_proportion, replace = FALSE)



# Combine representative and subsampled data
final_sampled_data <- bind_rows(representative_samples, subsampled_remaining)

# Extract sequence IDs for pruning the tree
sampled_ids <- final_sampled_data$tip_label


# Prune the tree based on sampled sequence IDs
pruned_tree <- drop.tip(B.1.177_tr, setdiff(B.1.177_tr$tip.label, sampled_ids))
unrooted_tree <- unroot(pruned_tree)

# Plot the pruned tree to confirm the result
plot(pruned_tree, no.margin = TRUE, cex = 0.5)

# final_sampled_data is the subset of metadata corresponding to the tree tips
# Check that tip labels are unique:
final_sampled_data <- final_sampled_data %>% distinct(tip_label, .keep_all = TRUE)


####### Align sample date to trees #######


# Extract tips and states
tip_labels <- pruned_tree$tip.label

# Assign states based on whether samples are from the UK
states <- ifelse(tree_metadata$is_uk == "Y", "I", "src") # "I" for UK, "src" for others

# Align states with tip labels
sampleStates <- data.frame(
  tip_label = tree_metadata$tip_label,
  state = states
)


sampleStates <- sampleStates[match(tip_labels, sampleStates$tip_label), ]

# Create matrix for infection status (UK infected vs. source population)
I <- src <- rep(0, length(tip_labels))
I[sampleStates$state == "I"] <- 1
src[sampleStates$state == "src"] <- 1
sampleStates_matrix <- cbind(I, src)
rownames(sampleStates_matrix) <- tip_labels

# Extract sampleTimes from final_sampled_data and set names accordingly
sampleTimes <- final_sampled_data$year_fractional
names(sampleTimes) <- final_sampled_data$tip_label

# Reorder sampleTimes to match pruned_tree$tip.label
sampleTimes <- sampleTimes[pruned_tree$tip.label]



# Create DatedTree
dated_tree <- DatedTree(
  phylo = pruned_tree,
  sampleTimes = sampleTimes,
  sampleStates = sampleStates_matrix,
  minEdgeLength = 0.0001,
  tol = 1e-6
)



# Test the result
print(dated_tree)
head(sampleStates_matrix)
head(sampleTimes)

plot(dated_tree, no.margin = TRUE, cex = 0.5)



####### Now add in the compartment model #######

####### Revised SARS-CoV-2 model with importation from Europe #######



# Define parameters
demes <- c('I_UK', 'src')      # Infectious in UK and source population
nonDemes <- c('S_UK', 'R_UK')  # Susceptible and recovered in the UK

parms <- list( 
  beta = 0.5,         # Transmission rate within UK
  gamma = 0.1,        # Recovery rate
  mu = 0.01,          # Importation rate 
  r = 0.01,           # Define a growth rate for src
  N_UK = 67e6,        # UK Population
  S0_UK = 67e6 - 100, # Initial S population UK
  N0_src = 5e7        # Initial source population
)


# Initial state of the compartments
x0 <- c(
  I_UK = 100,         # UK infected
  S_UK = parms$S0_UK,
  R_UK = 0,
  src = parms$N0_src  # External reservoir initial value
)


# Create matrix for birth and migration

births <- matrix('0', nrow = 2, ncol = 2, dimnames = list(demes, demes))

# I_UK generates new infections within the UK:
births['I_UK', 'I_UK'] <- 'parms$beta * S_UK * I_UK / parms$N_UK'

# Let src grow exponentially on its own:
births['src', 'src'] <- 'parms$r * src'


migrations <- matrix('0', nrow = 2, ncol = 2, dimnames = list(demes, demes))
migrations['src', 'I_UK'] <- 'parms$mu * src'  # importation from src to I_UK


deaths <- c(
  I_UK = 'parms$gamma * I_UK',
  src  = '0' # explicitly set src death rate to zero
)

nonDemeDynamics <- c(
  S_UK = '-parms$beta * S_UK * I_UK / parms$N_UK',  # S_UK decreases
  R_UK = 'parms$gamma * I_UK'                       # R_UK increases
  # Removed 'src' since it's already a deme
)


# Create demographic model
dm <- build.demographic.process(
  births = births,
  deaths = deaths,
  migrations = migrations,
  nonDemeDynamics = nonDemeDynamics,
  parameterNames = names(parms),
  rcpp = FALSE,
  sde = FALSE
)



# Validate demographic process
show.demographic.process(
  dm,
  theta = parms,
  x0 = x0,
  t0 = 0,
  t1 = 100,
  res = 100
)



# Calculate initial likelihood
likelihood <- colik(
  tree = dated_tree,
  theta = parms,
  demographic.process.model = dm,
  x0 = x0,
  t0 = 0,
  res = 1000,
  integrationMethod = "rk4"
)
print(paste("Initial log-likelihood:", likelihood))


# Define objective function
obj_fun <- function(lnbeta, lngamma, lnmu){
  theta <- list(
    beta = exp(lnbeta),
    gamma = exp(lngamma),
    mu = exp(lnmu),
    r = parms$r,
    N_UK = parms$N_UK,
    S0_UK = parms$S0_UK,
    N0_src = parms$N0_src
  )
  
  -colik(
    tree = dated_tree,
    theta = theta,
    demographic.process.model = dm,
    x0 = x0,
    t0 = 0,
    res = 1000
  )
}


# Run optimization
fit <- mle2(
  obj_fun,
  start = list(lnbeta = log(0.5), 
               lngamma = log(0.1),
               lnmu = log(0.01)),
  method = "Nelder-Mead",
  control = list(maxit = 500)
)


AIC(fit)                    # Measure model quality (lower the better)
# [1] 2232.061
logLik(fit)                 # Log-likelihood of the fitted model
# 'log Lik.' -1113.03 (df=3)
coef(fit)                   # Estimated parameter values (log-scales)
#     lnbeta    lngamma       lnmu 
# -1.6795856  0.1405608 -4.0889808
exp(coef(fit))              # Convert log parameters to original scales
#    lnbeta   lngamma      lnmu 
# 0.1864512 1.1509190 0.0167563 

# Use this to see how is it different (how bias is estimate)
exp(coef(fit)["lnbeta"]) - parms$beta
#     lnbeta 
# -0.3135488




####### Extract metadata from UShER #######

UShER_metadata <- read.csv("C:/Users/xusun/Desktop/Phylogenetic Project/Tree files/public-latest.metadata.tsv")
head(UShER_metadata)

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)

# Step 1: Replace \t with | to standardize delimiters
raw_metadata_clean <- UShER_metadata %>%
  mutate(across(everything(), ~ str_replace_all(., "\t", "|")))

# Step 2: Verify the number of columns after splitting
# Split the first few rows to check consistency in the number of fields
split_check <- str_split_fixed(raw_metadata_clean[[1]], "\\|", n = Inf)
print(head(split_check))

cleaned_metadata <- as.data.frame(split_check)

colnames(cleaned_metadata) <- c("strain", "genbank_accession", "date", "genbank_accession_2", 
                                "date_2", "country", "host", "completeness", "length", 
                                "Nextstrain_clade", "pangolin_lineage", 
                                "Nextstrain_clade_usher", "pango_lineage_usher")

# Remove duplicated columns
cleaned_metadata <- cleaned_metadata %>%
  select(-c(genbank_accession_2, date_2))


# Confirm changes
head(cleaned_metadata)

# Convert the 'date' column to proper date format, handling missing values
cleaned_metadata <- cleaned_metadata %>%
  mutate(date = as.Date(date, format="%Y-%m-%d")) %>%  # Convert valid date strings
  filter(!is.na(date)) %>%  # Remove missing values ("?")
  arrange(date)  # Sort by earliest date


B.1.177_UShER <- cleaned_metadata %>%
  filter(
    grepl("^(B\\.1\\.177|Z|Y|W|U|V|AA)", pangolin_lineage) &  # Matches only these lineages
      pangolin_lineage != "Unassigned" )

B.1.177_UShER <- B.1.177_UShER %>%
  filter(date != "?" & !is.na(date) & date != "")

# Identify different formats in the dataset
B.1.177_UShER <- B.1.177_UShER %>%
  mutate(
    date = case_when(
      grepl("^\\d{4}-\\d{2}-\\d{2}$", date) ~ date,  # Already in YYYY-MM-DD format
      grepl("^\\d{4}-\\d{2}$", date) ~ paste0(date, "-15"),  # YYYY-MM → YYYY-MM-01
      grepl("^\\d{4}$", date) ~ paste0(date, "-06-30"),  # YYYY → YYYY-01-01
      TRUE ~ NA_character_  # Invalid dates
    )
  )

B.1.177_UShER$date <- as.Date(B.1.177_UShER$date, format="%Y-%m-%d")

B.1.177_UShER <- B.1.177_UShER %>%
  mutate(country = case_when(
    country %in% c("Hong Kong", "Taiwan") ~ "China",
    TRUE ~ country
  ))

B.1.177_UShER <- B.1.177_UShER %>%
  mutate(country = case_when(
    country %in% c("Northern Ireland", "Northern_Ireland") ~ "Northern Ireland",
    TRUE ~ country
  ))

B.1.177_UShER$month <- format(B.1.177_UShER$date, "%Y-%m")


# Check sequence length distribution

B.1.177_UShER$length <- as.numeric(B.1.177_UShER$length)

summary(B.1.177_UShER$length)

hist(B.1.177_UShER$length, breaks = 50, col = "blue",
     main = "Distribution of Genome Lengths",
     xlab = "Genome Length", ylab = "Frequency")

# Filter out extreme values for sequence length
B.1.177_UShER <- B.1.177_UShER %>%
  filter(length >= 29500 & length <= 30000)

# Verify removal
summary(B.1.177_UShER$length)


# Convert month to Date format for accurate filtering
# Filtering out tips outside the desired period
B.1.177_UShER <- B.1.177_UShER %>%
  filter(month >= "2020-01" & month <= "2021-10")



####### Plot distribution of B.1.177 through time #######
B.1.177_plot <- ggplot(B.1.177_UShER, aes(x = date, fill = country)) +
  geom_histogram(bins = 500, position = "stack", alpha = 0.7) +
  scale_fill_manual(values = rainbow(41), name = "Country") +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  labs(title = "B", x = "Sample Time", y = "Sequence Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(B.1.177_plot)

####### Load the tree for UShER #######

UShER_tree <- read.tree("C:/Users/xusun/Desktop/Phylogenetic Project/Tree files/public-2025-02-08.all.nwk")
head(UShER_tree)

# Extract only the accession number (second field in the "|" delimited string)
UShER_tree$tip.label <- sapply(strsplit(UShER_tree$tip.label, "\\|"), `[`, 1)

# Check if the format matches metadata
head(UShER_tree$tip.label)

# Find the ones that are missing
unmatched_tips <- setdiff(B.1.177_UShER$strain, UShER_tree$tip.label)

# Print a few examples
head(unmatched_tips)


# Use the cleaned list
UShER_B.1.177_tree <- keep.tip(UShER_tree, intersect(B.1.177_UShER$strain, UShER_tree$tip.label))

UShER_B.1.177_tree_plot <- plot(UShER_B.1.177_tree, no.mar=F, cex = .2 )

####### Change the unit of branch length into substitutions per site #######

UShER_B.1.177_tree$edge.length <- UShER_B.1.177_tree$edge.length/29903


UShER_B.1.177_node_depth <- node.depth.edgelength(UShER_B.1.177_tree)

# Extract tip labels and their corresponding depths
tip_depths <- UShER_B.1.177_node_depth[1:length(UShER_B.1.177_tree$tip.label)]
names(tip_depths) <- UShER_B.1.177_tree$tip.label

# Define a threshold to remove outliers (e.g., top 5% longest branches)
threshold <- quantile(tip_depths, 0.95)

# Identify tips with very high depth (long branches)
long_tips_to_remove <- names(tip_depths[tip_depths > threshold])

UShER_B.1.177_tree <- drop.tip(UShER_B.1.177_tree, long_tips_to_remove)

print(UShER_B.1.177_tree_plot)


####### Now select proportionally by month and country #######

# Count the number of samples per month and country
samples_per_group <- B.1.177_UShER %>%
  group_by(month, country) %>%
  tally()


total_sample <- nrow(B.1.177_UShER)
target_size <- 5000
sampling_proportion <- target_size/total_sample

# Calculate the number of samples to take from each month and country
samples_per_group <- samples_per_group %>%
  mutate(samples_to_take = round(n * sampling_proportion))

samples_per_group$samples_to_take <- ifelse(samples_per_group$samples_to_take == 0, yes=1, no=samples_per_group$samples_to_take)


set.seed(1234) # to make sure will sample the same seqs if re-running script
sampled_data <- B.1.177_UShER %>%
  group_by(month, country) %>%
  group_map(~ sample_n(.x, size = min(nrow(.x), samples_per_group$samples_to_take[samples_per_group$month == .y$month & samples_per_group$country == .y$country])))

# Flatten the sampled_data list and extract sequence IDs
sampled_ids <- unlist(lapply(sampled_data, function(B.1.177_UShER) B.1.177_UShER$strain))

sampled_id_df <- B.1.177_UShER[B.1.177_UShER$strain %in% sampled_ids,]
sts <- decimal_date(as.Date(sampled_id_df$date))
names(sts) <- sampled_id_df$strain

# Assuming you have the correct tree object: B.1.177_tr_clean_rooted
# Prune the tree based on the sampled sequence IDs (check your tree object's name)
pruned_tree <- drop.tip(UShER_B.1.177_tree, setdiff(UShER_B.1.177_tree$tip.label, sampled_ids))
plot( pruned_tree , no.mar=F, cex = .2 )



# Keep only the dates for tips present in the pruned tree
sts_filtered <- sts[names(sts) %in% pruned_tree$tip.label]

# Ensure the names are in the same order as the tree tips
sts_filtered <- sts_filtered[match(pruned_tree$tip.label, names(sts_filtered))]


UShER_B.1.177_tree_rtt <- rtt(t=pruned_tree, tip.dates = sts_filtered, objective = "rms")



plot(UShER_B.1.177_tree_rtt, show.tip.label = F)
pruned_tree_unrooted <- unroot(UShER_B.1.177_tree_rtt)
plot(pruned_tree_unrooted, no.mar=F, cex = .2 )

valid_tips <- intersect(pruned_tree_unrooted$tip.label, names(sts_filtered))
length(valid_tips)  # Should match tree tip count
summary(sts_filtered)  # Should be reasonable years
length(pruned_tree_unrooted$tip.label)

####### Set parallel running for the core #######
num_cores <- detectCores() - 1  # Use all but one core to avoid freezing system
cl <- makeCluster(14)  # Create cluster
registerDoParallel(cl)  # Register parallel backend


####### Treedater running first round #######
timetr <- dater(pruned_tree_unrooted, sts_filtered, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0002, 0.0010), clock = 'strict')


plot( timetr , no.mar=F, cex = .2 )
axisPhylo(root.time = timetr$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage (before filtering outliers)")

rootToTipRegressionPlot( timetr )
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage Prior to Outlier Removal')

outliers <- outlierTips( timetr , alpha = 0.05)
hist(outliers$q)


####### Treedater running second round #######
pruned_tree_refined <- ape::drop.tip( pruned_tree_unrooted, rownames(outliers[outliers$q < 0.05,]) ) 

# Run again with refined pruned tree
timetr_refined <- dater(pruned_tree_refined, sts_filtered, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0004, 0.0008), clock = 'strict')

plot(timetr_refined, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_refined$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage (after filtering outliers)")

rtt_plot_timetr <- rootToTipRegressionPlot(timetr_refined)
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage After Outlier Removal')

# Manually define the additional singleton outliers
manual_singletons <- c("Germany/18139036/2021", "Slovakia/RNA/2021", "Germany/18061124/2021")

# Identify outliers using Treedater (already refined tree)
outliers_refined <- outlierTips(timetr_refined, alpha = 0.05)

# Combine automatic and manual outliers
final_outliers <- unique(c(rownames(outliers_refined[outliers_refined$q < 0.05, ]), manual_singletons))

####### Treedater running third round #######
pruned_tree_final <- ape::drop.tip( pruned_tree_refined, final_outliers) 

# Finally run with dater for refined tree
timetr_final <- dater(pruned_tree_final, sts_filtered, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0003, 0.0008), clock = 'strict')

plot(timetr_final, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_final$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage (after secondary outliers filtering)")

rtt_plot_timetr_final <- rootToTipRegressionPlot(timetr_final)
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage After Secondary Outlier Removal')

# Manually define additional singelton outliers again 
manual_singletons_2 <- c('USA/MA-MGH-00883/2020', 'Switzerland/BS-UHB-42514567/2020', 'England/CAMC-BBDF91/2020', 'EGY/CUNCI-HGC11I008/2021')

# Identify outliers using Treedater (already refined tree)
outliers_refined_2 <- outlierTips(timetr_final, alpha = 0.05)

# Combine automatic and manual outliers
final_outliers_2 <- unique(c(rownames(outliers_refined_2[outliers_refined_2$q < 0.05, ]), manual_singletons_2))


####### Treedater tunning forth round #######
pruned_tree_final_2 <- ape::drop.tip( pruned_tree_final, final_outliers_2) 

# Finally run with dater for refined tree excluding singleton
timetr_final_2 <- dater(pruned_tree_final_2, sts_filtered, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0003, 0.0008), clock = 'strict')

date_decimal(2019.85746687627)

plot(timetr_final_2, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_final$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage (after tertiary outliers filtering)")

rtt_plot_timetr_final_2 <- rootToTipRegressionPlot(timetr_final_2)
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage After Tertiary Outlier Removal')



####### Extract date of 20 earliest tips #######
# Extract tip labels and their assigned dates from the dated tree
tip_dates <- data.frame(
  tip = names(timetr_final_2$sts),  # Extract tip names
  date = timetr_final_2$sts         # Extract the dated values
)

# Sort tips by their assigned date
tip_dates_sorted <- tip_dates[order(tip_dates$date), ]

# Select the 20 earliest tips
earliest_tips <- head(tip_dates_sorted, 20)

# Convert decimal dates to standard date format
earliest_tips <- earliest_tips %>%
  mutate(date_formatted = as.Date(date_decimal(date)))

# Print them for review
print(earliest_tips)


####### Try to refine by including the earliest tip from Wuhan #######

B.1.177_UShER_alternative <- cleaned_metadata %>%
  filter(
    (grepl("^(B\\.1\\.177|Z|Y|W|U|V|AA)", pangolin_lineage) & pangolin_lineage != "Unassigned") |
      (strain == "Wuhan/IPBCAMS-WH-01/2019")  # Ensure Wuhan strain is included
  )


B.1.177_UShER_alternative <- B.1.177_UShER_alternative %>%
  mutate(country = case_when(
    country %in% c("Hong Kong", "Taiwan") ~ "China",
    TRUE ~ country
  ))

B.1.177_UShER_alternative <- B.1.177_UShER_alternative %>%
  mutate(country = case_when(
    country %in% c("Northern Ireland", "Northern_Ireland") ~ "Northern Ireland",
    TRUE ~ country
  ))

B.1.177_UShER_alternative$month <- format(B.1.177_UShER_alternative$date, "%Y-%m")


# Check sequence length distribution

B.1.177_UShER_alternative$length <- as.numeric(B.1.177_UShER_alternative$length)

summary(B.1.177_UShER_alternative$length)

hist(B.1.177_UShER$length, breaks = 50, col = "blue",
     main = "Distribution of Genome Lengths",
     xlab = "Genome Length", ylab = "Frequency")

# Filter out extreme values for sequence length
B.1.177_UShER_alternative <- B.1.177_UShER_alternative %>%
  filter(length >= 29500 & length <= 30000)

# Verify removal
summary(B.1.177_UShER_alternative$length)


# Convert month to Date format for accurate filtering
# Filtering out tips outside the desired period
B.1.177_UShER_alternative <- B.1.177_UShER_alternative %>%
  filter(month >= "2019-12" & month <= "2021-10")

# Use the cleaned list
UShER_B.1.177_tree_alternative <- keep.tip(UShER_tree, intersect(B.1.177_UShER_alternative$strain, UShER_tree$tip.label))


UShER_B.1.177_tree_alternative$edge.length <- UShER_B.1.177_tree_alternative$edge.length/29903


UShER_B.1.177_node_depth_alternative <- node.depth.edgelength(UShER_B.1.177_tree_alternative)

# Extract tip labels and their corresponding depths
tip_depths_alternative <- UShER_B.1.177_node_depth_alternative[1:length(UShER_B.1.177_tree_alternative$tip.label)]
names(tip_depths_alternative) <- UShER_B.1.177_tree_alternative$tip.label

# Define a threshold to remove outliers (e.g., top 5% longest branches)
threshold <- quantile(tip_depths_alternative, 0.95)

# Identify tips with very high depth (long branches)
long_tips_to_remove_alternative <- names(tip_depths_alternative[tip_depths_alternative > threshold])

UShER_B.1.177_tree_alternative <- drop.tip(UShER_B.1.177_tree_alternative, long_tips_to_remove_alternative)



# Count the number of samples per month and country
samples_per_group_alternative <- B.1.177_UShER_alternative %>%
  group_by(month, country) %>%
  tally()


total_sample <- nrow(B.1.177_UShER)
target_size <- 6000
sampling_proportion <- target_size/total_sample

# Calculate the number of samples to take from each month and country
samples_per_group_alternative <- samples_per_group_alternative %>%
  mutate(samples_to_take = round(n * sampling_proportion))

samples_per_group_alternative$samples_to_take <- ifelse(samples_per_group_alternative$samples_to_take == 0, yes=1, no=samples_per_group_alternative$samples_to_take)


set.seed(1234) # to make sure will sample the same seqs if re-running script
sampled_data_alternative <- B.1.177_UShER_alternative %>%
  group_by(month, country) %>%
  group_map(~ sample_n(.x, size = min(nrow(.x), samples_per_group_alternative$samples_to_take[samples_per_group_alternative$month == .y$month & samples_per_group_alternative$country == .y$country])))

# Flatten the sampled_data list and extract sequence IDs
sampled_ids_alternative <- unlist(lapply(sampled_data_alternative, function(B.1.177_UShER_alternative) B.1.177_UShER_alternative$strain))

sampled_id_df_alternative <- B.1.177_UShER_alternative[B.1.177_UShER_alternative$strain %in% sampled_ids_alternative,]
sts_alternative <- decimal_date(as.Date(sampled_id_df_alternative$date))
names(sts_alternative) <- sampled_id_df_alternative$strain

# Assuming you have the correct tree object: B.1.177_tr_clean_rooted
# Prune the tree based on the sampled sequence IDs (check your tree object's name)
pruned_tree_alternative <- drop.tip(UShER_B.1.177_tree_alternative, setdiff(UShER_B.1.177_tree_alternative$tip.label, sampled_ids_alternative))
plot( pruned_tree , no.mar=F, cex = .2 )

# Check if "Wuhan/IPBCAMS-WH-01/2019" is in the tree
"Wuhan/IPBCAMS-WH-01/2019" %in% pruned_tree_alternative$tip.label


length(pruned_tree_alternative$tip.label)
length(sts_alternative)

setdiff(pruned_tree_alternative$tip.label, names(sts_alternative))

sts_filtered_alternative <- sts_alternative[names(sts_alternative) %in% pruned_tree_alternative$tip.label]

sts_filtered_alternative <- sts_filtered_alternative[match(pruned_tree_alternative$tip.label, names(sts_filtered_alternative))]

length(sts_filtered_alternative)

UShER_B.1.177_tree_alternative_rtt <- rtt(t=pruned_tree_alternative, tip.dates = sts_filtered_alternative, objective = "rms")



plot(UShER_B.1.177_tree_alternative_rtt, show.tip.label = F)
pruned_tree_unrooted_alternative <- unroot(UShER_B.1.177_tree_alternative_rtt)
plot(pruned_tree_unrooted_alternative, no.mar=F, cex = .2 )

####### Treedater running first round including Wuhan #######

timetr_alternative <- dater(pruned_tree_unrooted_alternative, sts_filtered_alternative, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0002, 0.0010), clock = 'strict')

plot(timetr_alternative, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_alternative$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage including Wuhan-Hu-1 (before outliers filtering)")


rtt_plot_timetr_alternative <- rootToTipRegressionPlot(timetr_alternative)
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage including Wuhan-Hu-1 Before Outlier Removal')

manual_singletons_alternative <- c('USA/MA-MGH-00883/2020', 'England/MILK-B8FC28/2020', 'EGY/CUNCI-HGC11I008/2021')

outliers_alternative <- outlierTips(timetr_alternative, alpha = 0.05)

# Combine automatic and manual outliers
outliers_alternative_set <- unique(c(rownames(outliers_alternative[outliers_alternative$q < 0.05, ]), manual_singletons_alternative))

####### Treedater running second round including Wuhan #######

pruned_tree_unrooted_alternative_2 <- ape::drop.tip( pruned_tree_unrooted_alternative, outliers_alternative_set) 

timetr_alternative_2 <- dater(pruned_tree_unrooted_alternative_2, sts_filtered_alternative, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0003, 0.0010), clock = 'strict')

plot(timetr_alternative_2, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_alternative_2$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage including Wuhan-Hu-1 (after outliers filtering)")

rtt_plot_timetr_alternative_2 <- rootToTipRegressionPlot(timetr_alternative_2)
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage including Wuhan-Hu-1 After Outlier Removal')


manual_singletons_alternative_2 <- c('England/LIVE-1D7247/2020', 'Switzerland/BS-UHB-42505189/2020', 'Switzerland/BL−UHB−42515623/2020',
                                     'Northern Ireland/NIRE-23B389/2021', 'Switzerland/ZH-ETHZ-461408/2021')

outliers_alternative_2 <- outlierTips(timetr_alternative_2, alpha = 0.05)

outliers_alternative_set_2 <- unique(c(rownames(outliers_alternative_2[outliers_alternative_2$q < 0.05, ]), manual_singletons_alternative_2))

####### Treedater running third round including Wuhan #######

pruned_tree_unrooted_alternative_3 <- ape::drop.tip( pruned_tree_unrooted_alternative_2, outliers_alternative_set_2) 

timetr_alternative_3 <- dater(pruned_tree_unrooted_alternative_3, sts_filtered_alternative, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0003, 0.0010), clock = 'strict')

plot(timetr_alternative_3, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_alternative_3$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage including Wuhan-Hu-1 (after secondary outliers filtering)")

rtt_plot_timetr_alternative_3 <- rootToTipRegressionPlot(timetr_alternative_3)
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage including Wuhan-Hu-1 After Secondary Outlier Removal')


manual_singletons_alternative_3 <- c('Switzerland/BL-UHB-42515623/2020', 'Northern_Ireland/NIRE-23B389/2021', 'England/LIVE-1DD77F/2020', 'England/LIVE-1DD5A2/2020',
                                     'England/LIVE-1DD7AC/2020', 'TUN/TUN-2020-4874/2020')

outliers_alternative_3 <- outlierTips(timetr_alternative_3, alpha = 0.05)

outliers_alternative_set_3 <- unique(c(rownames(outliers_alternative_3[outliers_alternative_3$q < 0.05, ]), manual_singletons_alternative_3))

####### Treedater running forth round including Wuhan #######

pruned_tree_unrooted_alternative_4 <- unroot(ape::drop.tip( pruned_tree_unrooted_alternative_3, outliers_alternative_set_3))

timetr_alternative_4 <- dater(pruned_tree_unrooted_alternative_4, sts_filtered_alternative, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0003, 0.0010), clock = 'strict')

plot(timetr_alternative_4, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_alternative_4$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage including Wuhan-Hu-1 (after tertiary outliers filtering)")

rtt_plot_timetr_alternative_4 <- rootToTipRegressionPlot(timetr_alternative_4)
title(main = 'Temporal Evolution of Evolutionary Distance for 6k Sequenced B.1.177 Lineage including Wuhan-Hu-1 After Tertiary Outlier Removal')


####### Extract earliest 30 tips for alternative #######

# Extract tip labels and their assigned dates from the dated tree
tip_dates_alternative <- data.frame(
  tip = names(timetr_alternative_4$sts),  # Extract tip names
  date = timetr_alternative_4$sts         # Extract the dated values
)

# Sort tips by their assigned date
tip_dates_sorted_alternative <- tip_dates_alternative[order(tip_dates_alternative$date), ]

# Select the 20 earliest tips
earliest_tips_alternative <- head(tip_dates_sorted_alternative, 30)

# Convert decimal dates to standard date format
earliest_tips_alternative <- earliest_tips_alternative %>%
  mutate(date_formatted = as.Date(date_decimal(date)))

# Print them for review
print(earliest_tips_alternative)

manual_singletons_alternative_4 <- c('Germany/IMS-10213-CVDP-C139F839-D699-483D-9462-54EC34CEFD86/2021', 'England/LIVE-E24BDA/2021', 'England/PHEC-U307U9FC/2021',
                                     'DNK/144-Day1/2021', 'FRA/IHUCOVID-033024_Nova1/2021', 'England/PHEP-011898/2021', 'FRA/21062240386/2021', 'England/PHEP-011898/2021', 
                                     'England/LOND-128B79A/2021', 'England/LOND-128B299/2021', 'ESP/220-Day1/2021')



pruned_tree_alternative_5 <- drop.tip(pruned_tree_unrooted_alternative_4, manual_singletons_alternative_4)


timetr_alternative_5 <- dater(pruned_tree_alternative_5, sts_filtered_alternative, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0003, 0.0010), clock = 'strict')



####### Extract the full information of earliest tip from metadata #######
# Subset UShER B.1.177 metadata to retain only the earliest tips
earliest_tips_metadata <- B.1.177_UShER %>%
  filter(strain %in% earliest_tips_alternative$tip)

# Ensure column names match before merging
colnames(earliest_tips_alternative)[colnames(earliest_tips_alternative) == "tip"] <- "strain"

# Merge 'date_formatted' as 'inferred_date' into earliest_tips_metadata
earliest_tips_metadata <- earliest_tips_metadata %>%
  left_join(earliest_tips_alternative %>% select(strain, date_formatted), by = "strain") %>%
  rename(inferred_date = date_formatted)  # Rename to 'inferred_date'



####### Export the earliest tip dataframe as CSV #######

# Define the file path
output_file <- "C:/Users/xusun/Desktop/Phylogenetic Project/Codes/Earliest_Tips_Metadata.csv"

# Export dataframe as CSV
write.csv(earliest_tips_metadata, file = output_file, row.names = FALSE)







####### Try to too the tree based on the Wuhan tip #######
pruned_tree_rooted_alternative <- root(pruned_tree_alternative, outgroup = "Wuhan/IPBCAMS-WH-01/2019", resolve.root = TRUE)

plot(pruned_tree_rooted_alternative, no.mar =F, cex= .2)

timetr_alternative_rooted <- dater(pruned_tree_rooted_alternative, sts_filtered_alternative, s=29903, minblen=1/365, quiet=F, omega0=0.0008, meanRateLimits=c(0.0003, 0.0010), clock = 'strict')


plot(timetr_alternative_rooted, no.mar=F, cex = .2)
axisPhylo(root.time = timetr_alternative_rooted$timeOfMRCA, backward = F)
title(main="Randomly sampled 6000 sequences tree for B.1.177 lineage including Wuhan-Hu-1 rooted before outleir filtering")



###### See how many B.1.177 excluded #######





B.1.177_UShER_check <- cleaned_metadata %>%
  filter(
    grepl("^(B\\.1\\.177|Z|Y|W|U|V|AA)", pangolin_lineage) &  # Matches only these lineages
      pangolin_lineage != "Unassigned" )


####### Run a trial with 5k and 6k tree #######

# Step 1: Remove Duplicates from B.1.177_UShER
B.1.177_UShER_unique <- B.1.177_UShER_alternative %>%
  distinct(strain, .keep_all = TRUE)  # Keep only unique strain names

####### First extract information #######
B.1.177_UShER_pruned <- B.1.177_UShER_pruned %>%
  filter(strain != "Wuhan/IPBCAMS-WH-01/2019")

B.1.177_UShER_pruned <- B.1.177_UShER_unique %>%
  filter(strain %in% timetr_alternative_4$tip.label) 

missing_tips <- setdiff(timetr_alternative_4$tip.label, B.1.177_UShER_pruned$strain)
extra_tips <- setdiff(B.1.177_UShER_pruned$strain, timetr_alternative_4$tip.label)

# Create 'state' column based on UK vs. non-UK classification
B.1.177_UShER_pruned <- B.1.177_UShER_pruned %>%
  mutate(state = ifelse(country %in% c("England", "Scotland", "Wales", "Northern Ireland", "United Kingdom"), "I", "src"))

# Create state dataframe for phydynR

B.1.177_state_6k <- B.1.177_UShER_pruned$state

# Align states with tip labels
sampleStates <- data.frame(
  tip_label = timetr_alternative_4$tip.label,
  state = B.1.177_state_6k
)

sampleStates <- sampleStates[match(timetr_alternative_4$tip.label, sampleStates$tip_label), ]

# Extract sampleTimes from final_sampled_data and set names accordingly
sts_data <- data.frame(
  strain = names(timetr_alternative_4$sts),  # Extract strain names from the tree
  sts = as.numeric(timetr_alternative_4$sts)  # Extract numerical sts values
)


# Step 2: Ensure exact row alignment by matching strain names
B.1.177_UShER_pruned <- B.1.177_UShER_pruned %>%
  inner_join(sts_data, by = "strain") %>%  # Merge based on strain name
  arrange(match(strain, sts_data$strain))  # Ensure same order as sts_data

# Create matrix for infection status (UK infected vs. source population)
I <- src <- rep(0, length(timetr_alternative_4$tip.label))
I[sampleStates$state == "I"] <- 1
src[sampleStates$state == "src"] <- 1
sampleStates_matrix <- cbind(I, src)
rownames(sampleStates_matrix) <- timetr_alternative_4$tip.label

sampleTimes <- B.1.177_UShER_pruned$sts
names(sampleTimes) <- B.1.177_UShER_pruned$strain

# Reorder sampleTimes to match pruned_tree$tip.label
sampleTimes <- sampleTimes[timetr_alternative_4$tip.label]

dated_tree <- DatedTree(
  phylo = timetr_alternative_4,
  sampleTimes = sampleTimes,
  sampleStates = sampleStates_matrix,
  minEdgeLength = 0.0001,
  tol = 1e-6
)




# Calculate initial likelihood
likelihood <- colik(
  tree = dated_tree,
  theta = parms,
  demographic.process.model = dm,
  x0 = x0,
  t0 = 0,
  res = 1000,
  integrationMethod = "rk4"
)
print(paste("Initial log-likelihood:", likelihood))

# [1] "Initial log-likelihood: -150842.980517727"

# Define objective function
obj_fun <- function(lnbeta, lngamma, lnmu){
  theta <- list(
    beta = exp(lnbeta),
    gamma = exp(lngamma),
    mu = exp(lnmu),
    r = parms$r,
    N_UK = parms$N_UK,
    S0_UK = parms$S0_UK,
    N0_src = parms$N0_src
  )
  
  -colik(
    tree = dated_tree,
    theta = theta,
    demographic.process.model = dm,
    x0 = x0,
    t0 = 0,
    res = 1000
  )
}


# Run optimization
fit <- mle2(
  obj_fun,
  start = list(lnbeta = log(0.5), 
               lngamma = log(0.1),
               lnmu = log(0.01)),
  method = "Nelder-Mead",
  control = list(maxit = 500)
)



AIC(fit)                    # Measure model quality (lower the better)
# [1] 139981.9
logLik(fit)                 # Log-likelihood of the fitted model
# 'log Lik.' -69987.95 (df=3)
coef(fit)                   # Estimated parameter values (log-scales)
#     lnbeta   lngamma      lnmu 
#  -7.763712 13.983749 -4.199340
exp(coef(fit))              # Convert log parameters to original scales
#       lnbeta      lngamma         lnmu 
# 4.248765e-04 1.183219e+06 1.500548e-02 

# Use this to see how is it different (how bias is estimate)
exp(coef(fit)["lnbeta"]) - parms$beta
#     lnbeta 
# -0.4995751
