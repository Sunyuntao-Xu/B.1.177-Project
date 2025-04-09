library(dplyr)
library(phydynR)

Sys.setenv(LANGUAGE = "en")

####### Input B.1.177 tree for a try #######

tree_metadata_test <- read.csv('C:/Users/xusun/Desktop/Phylogenetic Project/Tree files/5k_treedater_metadata.csv')
head(tree_metadata_test)

B.1.177_tr_test <- read.tree('C:/Users/xusun/Desktop/Phylogenetic Project/Tree files/B.1.177_5k_treedater_alias.nwk')

# First extract tree and ensure correct format

# tree needed (contain tip labels, time states)
# metadata needed (contain tip labels, time, states)

####### Align sample date to trees #######

# Step 1: Get the final set of tip labels from pruned tree
tip_labels <- B.1.177_tr_test$tip.label
sampled_ids <- tree_metadata_test$tip_label
pruned_tree <- drop.tip(B.1.177_tr_test, setdiff(tip_labels, sampled_ids))
tip_labels <- pruned_tree$tip.label  # Update tip_labels to match pruned_tree

# Step 2: Create sample state dataframe, match and filter to tips
states <- ifelse(tree_metadata_test$is_uk == "Y", "I_UK", "src")

sampleStates_df <- data.frame(
  tip_label = tree_metadata_test$tip_label,
  state = states,
  stringsAsFactors = FALSE
) %>%
  filter(tip_label %in% tip_labels) %>%
  arrange(match(tip_label, tip_labels))

# Step 3: Create matrix with proper columns named as demes
I <- as.numeric(sampleStates_df$state == "I_UK")
src <- as.numeric(sampleStates_df$state == "src")
sampleStates_matrix <- cbind(I_UK = I, src = src)
rownames(sampleStates_matrix) <- sampleStates_df$tip_label

# Final check: tip label order must match
stopifnot(all(rownames(sampleStates_matrix) == tip_labels))

# Step 4: Align sampleTimes to match tree
sampleTimes <- tree_metadata_test$year_fractional
names(sampleTimes) <- tree_metadata_test$tip_label
sampleTimes <- sampleTimes[tip_labels]

# Step 5: Build dated tree
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



####### Define Model with Time-Varying beta(t) #######

# Time points for beta
betaTimes <- c(2020, 2020.5, 2021, 2022)
betaNames <- paste0("beta", betaTimes)

# Parameters
parms <- list(
  gamma = 0.1,
  mu = 0.01,
  r = 0.01,
  N_UK = 67e6,
  S0_UK = 67e6 - 100,
  N0_src = 5e7,
  
  beta2020 = 2,
  beta2020.5 = 5,
  beta2021 = 3,
  beta2022 = 1
)

# Define time-varying beta
parms$beta.t <- function(t, p) {
  approx(betaTimes, unlist(p[betaNames]), xout = t, rule = 2)$y
}

# Initial state
x0 <- c(
  I_UK = 100,
  S_UK = parms$S0_UK,
  R_UK = 0,
  src = parms$N0_src
)

# Demes
demes <- c('I_UK', 'src')
nonDemes <- c('S_UK', 'R_UK')

# Births matrix
births <- matrix('0', nrow = 2, ncol = 2, dimnames = list(demes, demes))
births['I_UK', 'I_UK'] <- 'parms$beta.t(t, parms) * S_UK * I_UK / parms$N_UK'
births['src', 'src'] <- 'parms$r * src'

# Migration matrix
migrations <- matrix('0', nrow = 2, ncol = 2, dimnames = list(demes, demes))
migrations['src', 'I_UK'] <- 'parms$mu * src'

# Deaths
deaths <- c(
  I_UK = 'parms$gamma * I_UK',
  src  = '0'
)

# Non-deme dynamics
nonDemeDynamics <- c(
  S_UK = '-parms$beta.t(t, parms) * S_UK * I_UK / parms$N_UK',
  R_UK = 'parms$gamma * I_UK'
)

# Build model
dm <- build.demographic.process(
  births = births,
  deaths = deaths,
  migrations = migrations,
  nonDemeDynamics = nonDemeDynamics,
  parameterNames = names(parms),
  rcpp = FALSE,
  sde = FALSE
)

# Visualize the model dynamics over time (e.g., 2020 to 2022.5)
show.demographic.process(
  dm,
  theta = parms,
  x0 = x0,
  t0 = 2020,
  t1 = 2022.5,
  res = 500
)


# Likelihood calculation
colik(
  tree = dated_tree,
  theta = parms,
  demographic.process.model = dm,
  x0 = x0,
  t0 = 0,
  res = 1000
)

# Objective function for MLE
obj_fun <- function(lnbeta2020, lnbeta2020.5, lnbeta2021, lnbeta2022, lngamma, lnmu){
  theta <- list(
    beta2020 = exp(lnbeta2020),
    beta2020.5 = exp(lnbeta2020.5),
    beta2021 = exp(lnbeta2021),
    beta2022 = exp(lnbeta2022),
    gamma = exp(lngamma),
    mu = exp(lnmu),
    r = parms$r,
    N_UK = parms$N_UK,
    S0_UK = parms$S0_UK,
    N0_src = parms$N0_src
  )
  
  theta$beta.t <- parms$beta.t
  
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
  start = list(
    lnbeta2020 = log(0.5),
    lnbeta2020.5 = log(0.3),
    lnbeta2021 = log(0.2),
    lnbeta2022 = log(0.1),
    lngamma = log(0.1),
    lnmu = log(0.01)
  ),
  method = "Nelder-Mead",
  control = list(maxit = 1000)
)

# Results
AIC(fit)
logLik(fit)
coef(fit)
exp(coef(fit))

