#------------------------------------------------------------------------------
# Multistate mark-recapture model for the Delta GSSMB study.
# Simulates capature histories, builds model, and then fits the model
# to simulated data. This uses the dDHMM distribution from nimbleEcology.
#
# R. Perry  4/24/2025
#------------------------------------------------------------------------------

library(nimble)
library(nimbleEcology)
library(abind)
library(MCMCvis)

# States:
# 1 = Sacramento R
# 2 = Sutter Sl
# 3 = Steamboat Sl
# 4 = Georgiana Sl
nstate = 4

# Set some true parameter values for filling matrices and 
# simulating capture histories.
#-- Survival
true_pars = list(
  S_sac1 = 0.9,
  S_sac2 = 0.8,
  S_sac3 = 0.7,
  S_sac4 = 0.7,
  S_sac5 = 0.6,
  phi_sut_min = 0.5,
  phi_sut_stm = 0.4,
  S_min = 0.7,
  S_stm1 = 0.9,
  S_stm2 = 0.8,
  S_geo1 = 0.5,
  S_geo2 = 0.4,
  #-- Detection
  p_sac1 = 0.9,
  p_sac2 = 0.8,
  p_sac3 = 0.7,
  p_sac4 = 0.7,
  p_sac5 = 0.9,
  p_sut = 0.8,
  p_min = 0.8,
  p_stm1 = 0.9,
  p_stm2 = 0.9,
  p_geo1 = 0.5,
  p_geo2 = 0.5,
  lambda = 0.8,
  #-- Routing
  psi_sut = 0.3,
  psi_stm = 0.3,
  psi_geo = 0.5)

# assign parameters in vector to scalars.
for(i in 1:length(true_pars)) assign(names(true_pars)[i], true_pars[[i]])

#-------------------------------------------------------------------------------
# Construct transition and detection matrices and populate with true values.
#-------------------------------------------------------------------------------
# Create a template for transition and detection arrays.
temp_mat = matrix(0, nrow = nstate + 1, nstate + 1)
# For transition matrix, if in death state, stay in death state with prob 1.
# For observation matrix, if true state is dead, probability of not observing is 1.
temp_mat[nstate+1, nstate+1] = 1

#---------
# Set true values into detection matrices and create detection array.
# Detection prob at Sac release = 1.
p_mat_rel = temp_mat
p_mat_rel[1, 1] <- 1

p_mat1 = temp_mat
p_mat1[1, ] <- c(p_sac1, 0, 0, 0, (1-p_sac1))

p_mat2 = temp_mat
p_mat2[1, ] <- c(p_sac2, 0, 0, 0, (1-p_sac2))
p_mat2[2, ] <- c(0, p_sut, 0, 0, (1-p_sut))
p_mat2[3, ] <- c(0, 0, p_stm1, 0, (1-p_stm1))

p_mat3 = temp_mat
p_mat3[1, ] <- c(p_sac3, 0, 0, 0, (1-p_sac3))
p_mat3[2, nstate+1] <- 1
p_mat3[3, nstate+1] <- 1
p_mat3[4, ] <- c(0, 0, 0, p_geo1, (1-p_geo1))

p_mat4 = temp_mat
p_mat4[1, ] <- c(p_sac4, 0, 0, 0, (1-p_sac4))
p_mat4[2, ] <- c(0, p_min, 0, 0, (1-p_min))
p_mat4[3, ] <- c(0, 0, p_stm2, 0, (1-p_stm2))
p_mat4[4, ] <- c(0, 0, 0, p_geo2, (1-p_geo2))

p_mat5 = temp_mat
p_mat5[1, ] <- c(p_sac5, 0, 0, 0, (1-p_sac5))

p_mat6 = temp_mat
p_mat6[1, 1] <- 1

p_arr <- abind(p_mat_rel, p_mat1, p_mat2, p_mat3, p_mat4, p_mat5, p_mat6, along = 3)
dim(p_arr) 

# Create transition matrices
tr_rel_1 = temp_mat
tr_rel_1[1, ] <- c(S_sac1, 0, 0, 0, 1-S_sac1)

tr_12 = temp_mat
tr_12[1, ] <- c(S_sac2*(1-psi_sut-psi_stm), S_sac2*psi_sut, S_sac2*psi_stm, 0, (1-S_sac2))
rowSums(tr_12)

tr_23 = temp_mat
tr_23[1, ] <- c(S_sac3*(1-psi_geo), 0, 0, S_sac3*psi_geo, (1-S_sac3))
tr_23[2, 2] <- 1
tr_23[3, 3] <- 1
rowSums(tr_23)

tr_34 = temp_mat
tr_34[1, ] <- c(S_sac4, 0, 0, 0, (1-S_sac4))
tr_34[2, ] <- c(0, phi_sut_min, phi_sut_stm, 0, (1-phi_sut_min-phi_sut_stm))
tr_34[3, ] <- c(0, 0, S_stm1, 0, (1-S_stm1))
tr_34[4, ] <- c(0, 0, 0, S_geo1, (1-S_geo1))
rowSums(tr_34)

tr_45 = temp_mat
tr_45[1, ] <- c(S_sac5, 0, 0, 0, (1-S_sac5))
tr_45[2, ] <- c(S_min, 0, 0, 0, (1-S_min))
tr_45[3, ] <- c(S_stm2, 0, 0, 0, (1-S_stm2))
tr_45[4, ] <- c(S_geo2, 0, 0, 0, (1-S_geo2))
rowSums(tr_45)

tr_56 = temp_mat
tr_56[1, ] <- c(lambda, 0, 0, 0, (1-lambda))
rowSums(tr_56)

tr_arr = abind(tr_rel_1, tr_12, tr_23, tr_34, tr_45, tr_56, along = 3)
dim(tr_arr)
p_arr
tr_arr

#------------------------------------------------------------------------------
# Now let's simulate capture histories using rDHMMo function.
# Let's set it  by release group, but add functionality to add indivdiual 
# covariates later.  Parameters above are the same for each release group, but
# the model will estimate them separately, by release group.
#------------------------------------------------------------------------------
# Vector of initial state probabilities. Release at Sac (state 1) with prob 1.
rel_state = c(1, 0, 0, 0, 0)

nfish_per_rel = 200
nrel = 5
ntot = nrel * nfish_per_rel
# Release group of each individual
rel_vec = rep(1:5, each = nfish_per_rel)

ch_mat = matrix(nrow = ntot, ncol = 7)

for(i in 1:ntot){
  ch_mat[i, ] = rDHMMo(n = 1,
                       init = rel_state,
                       probObs = p_arr,
                       probTrans = tr_arr,
                       len = 7,
                       checkRowSums = F)
}

# Capture histories for first few fish:
ch_mat[1:10, ]

# Check that likelihood function works.
dDHMMo(ch_mat[1, 1:7],
       init =  rel_state,
       probObs = p_arr[1:5, 1:5, 1:7],
       probTrans = tr_arr,
       len = 7,
       checkRowSums = F)

#-------------------------------------------------------------------------------
# Write the nimble model to estimate parameters from the simulated data.
# 
# Parameters are individual-specific to allow future flexibility, but for this
# simulation, parameters are assigned for each release group and then to each 
# indiviual.
#-------------------------------------------------------------------------------
# Nimble model
nimCode = nimbleCode({
  # Priors and parameters are release-specific.
  for(r in 1:nrel){
    # Priors
    S_sac1[r] ~ dbeta(1, 1)  #beta(1, 1) ~ uniform(0, 1)
    S_sac2[r] ~ dbeta(1, 1)
    S_sac3[r] ~ dbeta(1, 1)
    S_sac4[r] ~ dbeta(1, 1)
    S_sac5[r] ~ dbeta(1, 1)
    phi_sut_min[r] ~ dbeta(1, 1)
    phi_sut_stm[r] ~ dbeta(1, 1)
    S_min[r] ~ dbeta(1, 1)
    S_stm1[r] ~ dbeta(1, 1)
    S_stm2[r] ~ dbeta(1, 1)
    S_geo1[r] ~ dbeta(1, 1)
    S_geo2[r] ~ dbeta(1, 1)
    #-- Detection
    p_sac1[r] ~ dbeta(1, 1)
    p_sac2[r] ~ dbeta(1, 1)
    p_sac3[r] ~ dbeta(1, 1)
    p_sac4[r] ~ dbeta(1, 1)
    p_sac5[r] ~ dbeta(1, 1)
    p_sut[r] ~ dbeta(1, 1)
    p_min[r] ~ dbeta(1, 1)
    p_stm1[r] ~ dbeta(1, 1)
    p_stm2[r] ~ dbeta(1, 1)
    p_geo1[r] ~ dbeta(1, 1)
    p_geo2[r] ~ dbeta(1, 1)
    lambda[r] ~ dbeta(1, 1)
    #-- Routing
    psi_sut[r] ~ dbeta(1, 1)
    psi_stm[r] ~ dbeta(1, 1)
    psi_geo[r] ~ dbeta(1, 1)
    
    #----------------------------------------------------------------
    # Summaries of route-specific and total survival for each release.
    # Survival of fish that enter Sutter Slough is sum of phi terms.
    S_sut[r] <- phi_sut_min[r] + phi_sut_stm[r]
    
    # probability of taking Sac route and Geo route at both junctions.
    # Pr(sac at junction 1 then sac at junction 2)
    psi_sac_sac[r] <- (1-psi_sut[r] - psi_stm[r]) * (1-psi_geo[r]) 
    # Pr(sac at junction 1 then geo at junction 2)
    psi_sac_geo[r] <- (1-psi_sut[r] - psi_stm[r]) * psi_geo[r]
    
    # Survival probability for each route from Freeport to Chipps
    S_sac_rte[r] <- S_sac2[r]*S_sac3[r]*S_sac4[r]*S_sac5[r]
    S_sut_rte[r] <- S_sac2[r]*(phi_sut_min[r]*S_min[r] + phi_sut_stm[r] * S_stm2[r])
    S_stm_rte[r] <- S_sac2[r] * S_stm1[r] * S_stm2[r]
    S_geo_rte[r] <- S_sac2[r]*S_sac3[r]*S_geo1[r]*S_geo2[r]
    
    # Total through-delta survival from Freeport to Chipps
    S_delta[r] <- S_sac_rte[r] * psi_sac_sac[r] +
                  S_sut_rte[r] * psi_sut[r] +
                  S_stm_rte[r] * psi_stm[r] +
                  S_geo_rte[r] * psi_sac_geo[r]
  }

  #-- Now assigne release-specific parameters for each inividual
  for(i in 1:nfish){  
    #------ Define Detection Arrays -----
    p_arr[1, 1:5, 1, i] <- c(1, 0, 0, 0, 0)
    p_arr[2, 1:5, 1, i] <- rep(0, 5)
    p_arr[3, 1:5, 1, i] <- rep(0, 5)
    p_arr[4, 1:5, 1, i] <- rep(0, 5)
    p_arr[5, 1:5, 1, i] <- c(0, 0, 0, 0, 1)
    
    p_arr[1, 1:5, 2, i] <- c(p_sac1[rel[i]],0,0,0,(1-p_sac1[rel[i]]))
    p_arr[2, 1:5, 2, i] <- rep(0, 5)
    p_arr[3, 1:5, 2, i] <- rep(0, 5)
    p_arr[4, 1:5, 2, i] <- rep(0, 5)
    p_arr[5, 1:5, 2, i] <- c(0, 0, 0, 0, 1)
  
    p_arr[1, 1:5, 3, i] <- c(p_sac2[rel[i]],0,0,0,(1-p_sac2[rel[i]]))
    p_arr[2, 1:5, 3, i] <- c(0,p_sut[rel[i]],0,0,(1-p_sut[rel[i]]))
    p_arr[3, 1:5, 3, i] <- c(0,0,p_stm1[rel[i]],0,(1-p_stm1[rel[i]]))
    p_arr[4, 1:5, 3, i] <- rep(0, 5)
    p_arr[5, 1:5, 3, i] <- c(0, 0, 0, 0, 1)
    
    p_arr[1, 1:5, 4, i] <- c(p_sac3[rel[i]],0,0,0,(1-p_sac3[rel[i]]))
    p_arr[2, 1:5, 4, i] <- c(0,0,0,0,1)
    p_arr[3, 1:5, 4, i] <- c(0,0,0,0,1)
    p_arr[4, 1:5, 4, i] <- c(0,0,0,p_geo1[rel[i]],(1-p_geo1[rel[i]]))
    p_arr[5, 1:5, 4, i] <- c(0, 0, 0, 0, 1)
  
    p_arr[1, 1:5, 5, i] <- c(p_sac4[rel[i]],0,0,0,(1-p_sac4[rel[i]]))
    p_arr[2, 1:5, 5, i] <- c(0,p_min[rel[i]],0,0,(1-p_min[rel[i]]))
    p_arr[3, 1:5, 5, i] <- c(0,0,p_stm2[rel[i]],0,(1-p_stm2[rel[i]]))
    p_arr[4, 1:5, 5, i] <- c(0,0,0,p_geo2[rel[i]],(1-p_geo2[rel[i]]))
    p_arr[5, 1:5, 5, i] <- c(0, 0, 0, 0, 1)
  
    p_arr[1, 1:5, 6, i] <- c(p_sac5[rel[i]],0,0,0,(1-p_sac5[rel[i]]))
    p_arr[2, 1:5, 6, i] <- rep(0, 5)
    p_arr[3, 1:5, 6, i] <- rep(0, 5)
    p_arr[4, 1:5, 6, i] <- rep(0, 5)
    p_arr[5, 1:5, 6, i] <- c(0, 0, 0, 0, 1)
  
    p_arr[1, 1:5, 7, i] <- c(1, 0, 0, 0, 0)
    p_arr[2, 1:5, 7, i] <- rep(0, 5)
    p_arr[3, 1:5, 7, i] <- rep(0, 5)
    p_arr[4, 1:5, 7, i] <- rep(0, 5)
    p_arr[5, 1:5, 7, i] <- c(0, 0, 0, 0, 1)
    
    #------ Define Transition Arrays -----
    tr_arr[1, 1:5, 1, i] <- c(S_sac1[rel[i]],0,0,0,(1-S_sac1[rel[i]]))
    tr_arr[2, 1:5, 1, i] <- rep(0, 5)
    tr_arr[3, 1:5, 1, i] <- rep(0, 5)
    tr_arr[4, 1:5, 1, i] <- rep(0, 5)
    tr_arr[5, 1:5, 1, i] <- c(0, 0, 0, 0, 1)
    
    tr_arr[1, 1:5, 2, i] <- c(S_sac2[rel[i]]*(1-psi_sut[rel[i]]-psi_stm[rel[i]]), S_sac2[rel[i]]*psi_sut[rel[i]], S_sac2[rel[i]]*psi_stm[rel[i]], 0, (1-S_sac2[rel[i]]))
    tr_arr[2, 1:5, 2, i] <- rep(0, 5)
    tr_arr[3, 1:5, 2, i] <- rep(0, 5)
    tr_arr[4, 1:5, 2, i] <- rep(0, 5)
    tr_arr[5, 1:5, 2, i] <- c(0, 0, 0, 0, 1)
    
    tr_arr[1, 1:5, 3, i] <- c(S_sac3[rel[i]]*(1-psi_geo[rel[i]]), 0, 0, S_sac3[rel[i]]*psi_geo[rel[i]], (1-S_sac3[rel[i]]))
    tr_arr[2, 1:5, 3, i] <- c(0,1,0,0,0)
    tr_arr[3, 1:5, 3, i] <- c(0,0,1,0,0)
    tr_arr[4, 1:5, 3, i] <- rep(0, 5)
    tr_arr[5, 1:5, 3, i] <- c(0, 0, 0, 0, 1)
  
    tr_arr[1, 1:5, 4, i] <- c(S_sac4[rel[i]], 0, 0, 0, (1-S_sac4[rel[i]]))
    tr_arr[2, 1:5, 4, i] <- c(0,phi_sut_min[rel[i]], phi_sut_stm[rel[i]] ,0 ,(1-phi_sut_min[rel[i]] - phi_sut_stm[rel[i]]))
    tr_arr[3, 1:5, 4, i] <- c(0,0,S_stm1[rel[i]], 0, (1-S_stm1[rel[i]]))
    tr_arr[4, 1:5, 4, i] <- c(0,0,0,S_geo1[rel[i]],(1-S_geo1[rel[i]]))
    tr_arr[5, 1:5, 4, i] <- c(0, 0, 0, 0, 1)
  
    tr_arr[1, 1:5, 5, i] <- c(S_sac5[rel[i]],0,0,0,(1-S_sac5[rel[i]]))
    tr_arr[2, 1:5, 5, i] <- c(S_min[rel[i]],0,0,0,(1-S_min[rel[i]]))
    tr_arr[3, 1:5, 5, i] <- c(S_stm2[rel[i]],0,0,0,(1-S_stm2[rel[i]]))
    tr_arr[4, 1:5, 5, i] <- c(S_geo2[rel[i]],0,0,0,(1-S_geo2[rel[i]]))
    tr_arr[5, 1:5, 5, i] <- c(0, 0, 0, 0, 1)
  
    tr_arr[1, 1:5, 6, i] <- c(lambda[rel[i]],0,0,0,(1-lambda[rel[i]]))
    tr_arr[2, 1:5, 6, i] <- rep(0, 5)
    tr_arr[3, 1:5, 6, i] <- rep(0, 5)
    tr_arr[4, 1:5, 6, i] <- rep(0, 5)
    tr_arr[5, 1:5, 6, i] <- c(0, 0, 0, 0, 1)
    
    # Initial state probability vector for all fish. Released in state 1.
    # rel_vec <- c(1, 0, 0, 0, 0)
    # HMM likelihood from nimbleEcology package.
    ch_mat[i, 1:7] ~ dDHMMo(init = c(1, 0, 0, 0, 0)[1:5],
                         probObs = p_arr[1:5, 1:5, 1:7, i],
                         probTrans = tr_arr[1:5, 1:5, 1:6, i],
                         len = 7,
                         checkRowSums = 0)
  }
})

# Initial values.  Need only for routing and phi to ensure that initial values
# sum to < 1 across routes.
inits = list(phi_sut_min = rep(0.2, nrel),
             phi_sut_stm = rep(0.4, nrel),
             psi_sut = rep(0.1, nrel),
             psi_stm = rep(0.3, nrel),
             psi_geo = rep(0.2, nrel))

# Define model.
nimMod = nimbleModel(code = nimCode,
                     inits = inits,
                     data = list(ch_mat = ch_mat),
                     constants = list(nfish = ntot,
                                      nrel = nrel,
                                      rel = rel_vec))

# Configure MCMC
confMCMC = configureMCMC(nimMod, onlySlice = TRUE)

# Monitor parameters not included by default
confMCMC$addMonitors(c("S_sut", "psi_sac_sac", "psi_sac_geo", "S_sac_rte",
                       "S_sut_rte", "S_stm_rte", "S_geo_rte", "S_delta"))

# Build the mcmc, compile the model, and compile the mcmc and model together.
MCMC = buildMCMC(confMCMC)
CModel = compileNimble(nimMod)
CMCMC <- compileNimble(MCMC, project = CModel)

# Run the model!
mcmc_out = runMCMC(CMCMC, 
                niter = 2000,  # <- May need to increase for real data, look at Rhats (<1.1) and traceplots
                nchains = 3,
                nburnin = 1000, # <- Increase burnin it if needed (chains not stabilized in traceplots)
                thin = 1,
                samplesAsCodaMCMC = T)

# Examine trace plots for all monitored parameters.  Set pdf = T to write to file.
MCMCtrace(mcmc_out,
          #params = names(true_pars),
          pdf = F,
          ind = T, 
          iter = 10000,
          Rhat = T,
          n.eff = T,
          gvals = rep(true_pars, each = 5)) #<- Maggie, comment this line out for real data.
                                            # This plots the true values for the simulation,
                                            # which you won't have for real data.

# Table of parameter estimates and posterior summaries
MCMCsummary(mcmc_out,
            digits = 3)
