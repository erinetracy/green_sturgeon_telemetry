#------------------------------------------------------------------------------
# Multistate mark-recapture model for the Delta GSSMB study.
# Simulates capature histories, builds model, and then fits the model
# to simulated data. This uses the dDHMM distribution from nimbleEcology.
#
# R. Perry  4/24/2025
#------------------------------------------------------------------------------

# Simulate capture histories and fit a multistate model in Nimble using
# dDHMM distribution from nimbleEcology.
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

# Now let's simulate capture histories using rDHMMo function.
# Vector of initial state probabilities. Release at Sac (state 1) with prob 1.
rel_vec = c(1, 0, 0, 0, 0)

nfish = 1000

ch_mat = matrix(nrow = nfish, ncol = 7)

for(i in 1:nfish){
  ch_mat[i, ] = rDHMMo(n = 1,
                       init = rel_vec,
                       probObs = p_arr,
                       probTrans = tr_arr,
                       len = 7,
                       checkRowSums = F)
}

p_arr[,,4]
tr_arr[,,4]

dDHMMo(ch_mat[131, 1:7],
       init = rel_vec,
       probObs = p_arr,#[1:5, 1:5, 1:7],
       probTrans = tr_arr,
       len = 7,
       checkRowSums = F)


#-------------------------------------------------------------------------------
# Write the nimble model to estimate parameters from the simulated data.
#-------------------------------------------------------------------------------
# Nimble model
nimCode = nimbleCode({
  # Priors
  S_sac1 ~ dbeta(1, 1)  #beta(1, 1) ~ uniform(0, 1)
  S_sac2 ~ dbeta(1, 1)
  S_sac3 ~ dbeta(1, 1)
  S_sac4 ~ dbeta(1, 1)
  S_sac5 ~ dbeta(1, 1)
  phi_sut_min ~ dbeta(1, 1)
  phi_sut_stm ~ dbeta(1, 1)
  S_min ~ dbeta(1, 1)
  S_stm1 ~ dbeta(1, 1)
  S_stm2 ~ dbeta(1, 1)
  S_geo1 ~ dbeta(1, 1)
  S_geo2 ~ dbeta(1, 1)
  #-- Detection
  p_sac1 ~ dbeta(1, 1)
  p_sac2 ~ dbeta(1, 1)
  p_sac3 ~ dbeta(1, 1)
  p_sac4 ~ dbeta(1, 1)
  p_sac5 ~ dbeta(1, 1)
  p_sut ~ dbeta(1, 1)
  p_min ~ dbeta(1, 1)
  p_stm1 ~ dbeta(1, 1)
  p_stm2 ~ dbeta(1, 1)
  p_geo1 ~ dbeta(1, 1)
  p_geo2 ~ dbeta(1, 1)
  lambda ~ dbeta(1, 1)
  #-- Routing
  psi_sut ~ dbeta(1, 1)
  psi_stm ~ dbeta(1, 1)
  psi_geo ~ dbeta(1, 1)
  
  #------ Define Detection Arrays -----
  p_arr[1, 1:5, 1] <- c(1, 0, 0, 0, 0)
  p_arr[2, 1:5, 1] <- rep(0, 5)
  p_arr[3, 1:5, 1] <- rep(0, 5)
  p_arr[4, 1:5, 1] <- rep(0, 5)
  p_arr[5, 1:5, 1] <- c(0, 0, 0, 0, 1)
  
  p_arr[1, 1:5, 2] <- c(p_sac1,0,0,0,(1-p_sac1))
  p_arr[2, 1:5, 2] <- rep(0, 5)
  p_arr[3, 1:5, 2] <- rep(0, 5)
  p_arr[4, 1:5, 2] <- rep(0, 5)
  p_arr[5, 1:5, 2] <- c(0, 0, 0, 0, 1)

  p_arr[1, 1:5, 3] <- c(p_sac2,0,0,0,(1-p_sac2))
  p_arr[2, 1:5, 3] <- c(0,p_sut,0,0,(1-p_sut))
  p_arr[3, 1:5, 3] <- c(0,0,p_stm1,0,(1-p_stm1))
  p_arr[4, 1:5, 3] <- rep(0, 5)
  p_arr[5, 1:5, 3] <- c(0, 0, 0, 0, 1)
  
  p_arr[1, 1:5, 4] <- c(p_sac3,0,0,0,(1-p_sac3))
  p_arr[2, 1:5, 4] <- c(0,0,0,0,1)
  p_arr[3, 1:5, 4] <- c(0,0,0,0,1)
  p_arr[4, 1:5, 4] <- c(0,0,0,p_geo1,(1-p_geo1))
  p_arr[5, 1:5, 4] <- c(0, 0, 0, 0, 1)

  p_arr[1, 1:5, 5] <- c(p_sac4,0,0,0,(1-p_sac4))
  p_arr[2, 1:5, 5] <- c(0,p_min,0,0,(1-p_min))
  p_arr[3, 1:5, 5] <- c(0,0,p_stm2,0,(1-p_stm2))
  p_arr[4, 1:5, 5] <- c(0,0,0,p_geo2,(1-p_geo2))
  p_arr[5, 1:5, 5] <- c(0, 0, 0, 0, 1)

  p_arr[1, 1:5, 6] <- c(p_sac5,0,0,0,(1-p_sac5))
  p_arr[2, 1:5, 6] <- rep(0, 5)
  p_arr[3, 1:5, 6] <- rep(0, 5)
  p_arr[4, 1:5, 6] <- rep(0, 5)
  p_arr[5, 1:5, 6] <- c(0, 0, 0, 0, 1)

  p_arr[1, 1:5, 7] <- c(1, 0, 0, 0, 0)
  p_arr[2, 1:5, 7] <- rep(0, 5)
  p_arr[3, 1:5, 7] <- rep(0, 5)
  p_arr[4, 1:5, 7] <- rep(0, 5)
  p_arr[5, 1:5, 7] <- c(0, 0, 0, 0, 1)
  
  #------ Define Transition Arrays -----
  tr_arr[1, 1:5, 1] <- c(S_sac1,0,0,0,(1-S_sac1))
  tr_arr[2, 1:5, 1] <- rep(0, 5)
  tr_arr[3, 1:5, 1] <- rep(0, 5)
  tr_arr[4, 1:5, 1] <- rep(0, 5)
  tr_arr[5, 1:5, 1] <- c(0, 0, 0, 0, 1)
  
  tr_arr[1, 1:5, 2] <- c(S_sac2*(1-psi_sut-psi_stm), S_sac2*psi_sut, S_sac2*psi_stm, 0, (1-S_sac2))
  tr_arr[2, 1:5, 2] <- rep(0, 5)
  tr_arr[3, 1:5, 2] <- rep(0, 5)
  tr_arr[4, 1:5, 2] <- rep(0, 5)
  tr_arr[5, 1:5, 2] <- c(0, 0, 0, 0, 1)
  
  tr_arr[1, 1:5, 3] <- c(S_sac3*(1-psi_geo), 0, 0, S_sac3*psi_geo, (1-S_sac3))
  tr_arr[2, 1:5, 3] <- c(0,1,0,0,0)
  tr_arr[3, 1:5, 3] <- c(0,0,1,0,0)
  tr_arr[4, 1:5, 3] <- rep(0, 5)
  tr_arr[5, 1:5, 3] <- c(0, 0, 0, 0, 1)

  tr_arr[1, 1:5, 4] <- c(S_sac4, 0, 0, 0, (1-S_sac4))
  tr_arr[2, 1:5, 4] <- c(0,phi_sut_min,phi_sut_stm,0,(1-phi_sut_min-phi_sut_stm))
  tr_arr[3, 1:5, 4] <- c(0,0,S_stm1,0,(1-S_stm1))
  tr_arr[4, 1:5, 4] <- c(0,0,0,S_geo1,(1-S_geo1))
  tr_arr[5, 1:5, 4] <- c(0, 0, 0, 0, 1)

  tr_arr[1, 1:5, 5] <- c(S_sac5,0,0,0,(1-S_sac5))
  tr_arr[2, 1:5, 5] <- c(S_min,0,0,0,(1-S_min))
  tr_arr[3, 1:5, 5] <- c(S_stm2,0,0,0,(1-S_stm2))
  tr_arr[4, 1:5, 5] <- c(S_geo2,0,0,0,(1-S_geo2))
  tr_arr[5, 1:5, 5] <- c(0, 0, 0, 0, 1)

  tr_arr[1, 1:5, 6] <- c(lambda,0,0,0,(1-lambda))
  tr_arr[2, 1:5, 6] <- rep(0, 5)
  tr_arr[3, 1:5, 6] <- rep(0, 5)
  tr_arr[4, 1:5, 6] <- rep(0, 5)
  tr_arr[5, 1:5, 6] <- c(0, 0, 0, 0, 1)
  
  # Initial state probability vector for all fish. Released in state 1.
  #rel_vec <- c(1, 0, 0, 0, 0)
  
  for(i in 1:nfish){
    ch_mat[i, 1:7] ~ dDHMMo(init = c(1, 0, 0, 0, 0)[1:5],
                         probObs = p_arr[1:5, 1:5, 1:7],
                         probTrans = tr_arr[1:5, 1:5, 1:6],
                         len = 7,
                         checkRowSums = 0)
  }
})


# Initial values.  Need only for routing and phi to ensure that initial values
# sum to < 1 across routes.
inits = list(phi_sut_min = 0.2,
             phi_sut_stm = 0.4,
             psi_sut = 0.1,
             psi_stm = 0.3,
             psi_geo = 0.2)

nimMod = nimbleModel(code = nimCode,
                     inits = true_pars,
                     data = list(ch_mat = ch_mat),
                     constants = list(nfish = nfish))

confMCMC = configureMCMC(nimMod, onlySlice = TRUE)
MCMC = buildMCMC(confMCMC)
CModel = compileNimble(nimMod)
CMCMC <- compileNimble(MCMC, project = CModel)

mcmc_out = runMCMC(CMCMC, 
                niter = 2000,
                nchains = 3,
                nburnin = 1000,
                thin = 1,
                samplesAsCodaMCMC = T)

MCMCtrace(mcmc_out,
          params = names(true_pars),
          pdf = F,
          ind = T, 
          iter = 10000,
          Rhat = T,
          n.eff = T,
          gvals = true_pars)
