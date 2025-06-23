install.packages('optimParallel')

library(dplyr)
library(phydynR)
library(optimParallel)

Sys.setenv(LANGUAGE = "en")

cl <- makeCluster(detectCores() /2)  # Use all but one core
setDefaultCluster(cl)

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
states <- ifelse(tree_metadata_test$is_uk == "Y", "I_UK", "I_src")  

sampleStates_df <- data.frame(
  tip_label = tree_metadata_test$tip_label,
  state = states,
  stringsAsFactors = FALSE
) %>%
  filter(tip_label %in% tip_labels) %>%
  arrange(match(tip_label, tip_labels))

# Step 3: Create matrix with proper columns named as demes
I <- as.numeric(sampleStates_df$state == "I_UK")
I_src <- as.numeric(sampleStates_df$state == "I_src")  
sampleStates_matrix <- cbind(I_UK = I, I_src = I_src)  
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

mTimes <- c(2020.0, 2020.44, 2020.98)  # e.g. Jan 1, June 8, Dec 24 in decimal
mNames <- c("m1", "m2", "m3")

# Define parameters
parms <- list(
  gamma = 365 / 7,           # Recovery rate, e.g., 1/7 per day, and in year
  
  beta2020.4 = 50,      # Approximate beta in early June 2020
  beta2020.6 = 100,      # Approximate beta in late July 2020
  beta2020.75 = 400,     # Peak beta in September 2020
  beta2020.95 = 300,     # Beta in December 2020
  beta2021.15 = 20,     # Beta in February 2021
  beta2021.2 = 10,          # Beta after Febuary 2021
  m1 = 0.1,             # migration rate
  m2 = 0.03,
  m3 = 0.01
)

# Define the time-varying beta function (force of infection)
parms$beta.t <- function(t, p){
  approx(betaTimes, unlist(p[betaNames]), xout = t, rule=2)$y
}

parms$m.t <- function(t, p) {
  approx(mTimes, unlist(p[mNames]), xout = t, rule = 2)$y
}


# Set initial state: Only IR compartments, no S.
x0 <- c(
  I_UK = 10,
  I_src = 100,
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
migrations["I_UK", "I_src"] <- "parms$m.t(t, parms) * I_UK"   # Migration from Europe to UK
migrations["I_src", "I_UK"] <- "parms$m.t(t, parms) * I_UK"   # Migration from UK to Europe


# Deaths represent losses from the infectious compartments due to recovery and natural death.
deaths <- c(
  I_UK = "parms$gamma * I_UK",
  I_src = "parms$gamma * I_src"
)


# Non-deme dynamics: Flow into recovered compartments.
nonDemeDynamics <- c(
  R_UK = "parms$gamma*I_UK",
  R_src = "parms$gamma*I_src"
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
                    lnm1, lnm2, lnm3,
                    lnI0_UK, lnI0_src) {
  
  # Add lnI0_UK and lnI0_src here
  theta <- list(
    beta2020.4    = exp(lnbeta2020.4),
    beta2020.6    = exp(lnbeta2020.6),
    beta2020.75   = exp(lnbeta2020.75),
    beta2020.95   = exp(lnbeta2020.95),
    beta2021.15   = exp(lnbeta2021.15),
    beta2021.2    = exp(lnbeta2021.2),
    gamma         = 365 / 7,
    m1            = exp(lnm1),
    m2            = exp(lnm2),
    m3            = exp(lnm3)
    )
  
  # Now set x0 dynamically, with estimated I_UK
  x0_dyn <- c(
    I_UK = exp(lnI0_UK),
    I_src = exp(lnI0_src),  # now estimated dynamically
    R_UK = 0,
    R_src = 0
  )
  
  
  theta$beta.t <- parms$beta.t
  theta$m.t <- parms$m.t
  
  ll <- tryCatch(
    colik(
      tree = dated_tree,
      theta = theta,
      demographic.process.model = dm,
      x0 = x0_dyn,        # use the dynamic x0 here
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
    lnbeta2020.4  = log(50),
    lnbeta2020.6  = log(100),
    lnbeta2020.75 = log(400),
    lnbeta2020.95 = log(300),
    lnbeta2021.15 = log(20),
    lnbeta2021.2  = log(10),
    lnm1          = log(0.1),
    lnm2          = log(0.03),
    lnm3          = log(0.01),
    lnI0_UK       = log(10),        
    lnI0_src      = log(100)
  ),
  method = "Nelder-Mead",
  control = list(maxit = 2000)
)




# Check the results
AIC(fit)
# [1] 55198.51

logLik(fit)
# 'log Lik.' -27590.26 (df=9)

coefs <- coef(fit)
print(coefs)
# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2          lnmu 
#   -0.8357754     2.9775500     1.6429063     2.7180352    -1.9776397    -0.1313246    -4.8527977 

#        lnm       lnI0_UK 
# -0.2589397     3.9504695   


print(exp(coefs))

# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2          lnmu 
#  0.433538174  19.639640352   5.170173678  15.150525446   0.138395502   0.876933106   0.007806507 

#          lnm       lnI0_UK 
# 0.771869567  51.959755905 



######## Both I_UK and I_src estimated without N0_src and N0_UK #######
AIC(fit)
# [1] 53040.1

logLik(fit)
# 'log Lik.' -26510.05 (df=10)

coefs <- coef(fit)
print(coefs)
# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2          lnmu 
#    2.2642873     3.0505157     1.8385908     3.3203062     8.1957259     3.7373889   -23.6530533 

#          lnm       lnI0_UK      lnI0_src 
#   -0.1104735    -0.8702661     0.2844055 

print(exp(coefs))
# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2          lnmu 
# 9.624263e+00  2.112624e+01  6.287671e+00  2.766882e+01  3.625422e+03  4.198821e+01  5.340839e-11 

#          lnm       lnI0_UK      lnI0_src 
# 8.954100e-01  4.188401e-01  1.328972e+00 



####### mu removed, I_UK an I_src estimated, adjusted to per year
AIC(fit)
# [1] 196301.9

logLik(fit)
# 'log Lik.' -98141.93 (df=9)

coefs <- coef(fit)
print(coefs)

# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2           lnm 
#     3.590588      3.441364      3.675585      4.285186     -7.394489     -2.082821     -6.706446 

#      lnI0_UK      lnI0_src 
#     2.768535      2.194466 

print(exp(coefs))

# lnbeta2020.4  lnbeta2020.6 lnbeta2020.75 lnbeta2020.95 lnbeta2021.15  lnbeta2021.2           lnm 
# 3.625540e+01  3.122953e+01  3.947175e+01  7.261605e+01  6.146308e-04  1.245783e-01  1.223003e-03 

#      lnI0_UK      lnI0_src 
# 1.593528e+01  8.975211e+00

