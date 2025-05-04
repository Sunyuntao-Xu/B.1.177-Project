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
# Define time points for beta(t)
betaTimes <- c(2020.4, 2020.6, 2020.75, 2020.95, 2021.15, 2021.2)
betaNames <- paste0("beta", betaTimes)

# Define parameters
parms <- list(
  gamma = 0.14,           # Recovery rate, e.g., 1/7 per day
  mu = 0.01,              # Natural death rate
  N0_src = 5e7,           # Initial number for Europe, if needed
  
  beta2020.4 = 1,      # Approximate beta in early June 2020
  beta2020.6 = 3,      # Approximate beta in late July 2020
  beta2020.75 = 5,     # Peak beta in September 2020
  beta2020.95 = 2,     # Beta in December 2020
  beta2021.15 = 1,     # Beta in February 2021
  beta2021.2 = 0.1,          # Beta after Febuary 2021
  m = 0.01             # migration rate 
)

# Define the time-varying beta function (force of infection)
parms$beta.t <- function(t, p){
  approx(betaTimes, unlist(p[betaNames]), xout = t, rule=2)$y
}

# Set initial state: Only IR compartments, no S.
x0 <- c(
  I_UK = 100,
  I_src = 1000,
  R_UK = 0,
  R_src = 0
)

# Define demes (infectious compartments) and nonDemes (recovered compartments)
demes <- c("I_UK", "I_src")
nonDemes <- c("R_UK", "R_src")

# Births matrix: representing the generation of new infections
# For both regions, new infections occur at rate beta(t)*I.
births <- matrix("0", nrow=2, ncol=2, dimnames=list(demes, demes))
births["I_UK", "I_UK"] <- "parms$beta.t(t, parms)*I_UK"
births["I_src", "I_src"] <- "parms$beta.t(t, parms)*I_src"

# We omit migration terms in the ODEs since balanced migration cancels out.
migrations <- matrix("0", nrow = 2, ncol = 2, dimnames = list(demes, demes))
migrations["I_UK", "I_src"] <- "parms$m * I_UK"   # Migration from Europe to UK
migrations["I_src", "I_UK"] <- "parms$m * I_UK"   # Migration from UK to Europe


# Deaths represent losses from the infectious compartments due to recovery and natural death.
deaths <- c(
  I_UK = "parms$mu * I_UK",
  I_src = "parms$mu * I_src"
)

# Non-deme dynamics: Flow into recovered compartments.
nonDemeDynamics <- c(
  R_UK = "parms$gamma*I_UK - parms$mu*R_UK",
  R_src = "parms$gamma*I_src - parms$mu*R_src"
)

# Build the demographic process with the simplified IR model.
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
  t1 = 2022,
  res = 500
)


obj_fun <- function(lnbeta2020.4, lnbeta2020.6, lnbeta2020.75, 
                    lnbeta2020.95, lnbeta2021.15, lnbeta2021.2,
                    lnmu, lnm) {  # Add ln(m) here
  theta <- list(
    beta2020.4    = exp(lnbeta2020.4),
    beta2020.6    = exp(lnbeta2020.6),
    beta2020.75   = exp(lnbeta2020.75),
    beta2020.95   = exp(lnbeta2020.95),
    beta2021.15   = exp(lnbeta2021.15),
    beta2021.2    = exp(lnbeta2021.2),
    gamma         = 0.14,
    mu            = exp(lnmu),
    m             = exp(lnm),     # Add migration rate here
    N0_src        = parms$N0_src
  )
  
  theta$beta.t <- parms$beta.t
  
  ll <- tryCatch(
    colik(
      tree = dated_tree,
      theta = theta,
      demographic.process.model = dm,
      x0 = x0,
      t0 = 2020,
      res = 1000
    ),
    error = function(e) {
      message("Error in colik: ", e$message)
      return(NA)
    }
  )
  
  if(is.na(ll) || ll == -Inf) {
    return(1e12)
  } else {
    return(-ll)
  }
}

# Run optimization with updated starting values (note: beta2021.2 initial guess is log(1e-6))
fit <- mle2(
  obj_fun,
  start = list(
    lnbeta2020.4  = log(1),
    lnbeta2020.6  = log(3),
    lnbeta2020.75 = log(5),
    lnbeta2020.95 = log(2),
    lnbeta2021.15 = log(1),
    lnbeta2021.2  = log(0.1),
    lnmu          = log(0.01),
    lnm           = log(0.01)   # initial guess for migration rate
  ),
  method = "Nelder-Mead",
  control = list(maxit = 3000)
)

print(elapsed_time)



# Check the results
AIC(fit)
# [1] 55368.19

logLik(fit)
# 'log Lik.' -27676.09 (df=8)

coefs <- coef(fit)
print(coefs)
# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2          lnmu 
#    -1.389644      2.787364      1.752361      2.678468      2.825441    -10.883130     -6.420074 

#        lnm 
# -0.5024209  


print(exp(coefs))

# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2          lnmu 
# 2.491640e-01  1.623816e+01  5.768204e+00  1.456276e+01  1.686838e+01  1.877226e-05  1.628535e-03 

#          lnm 
# 6.050646e-01
