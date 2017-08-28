# Species parameters
# u_c : matrix of dimension D X R with niche optima for each species and environmental variable
# s_c : matrix of dimension D X R with niche breadth for each species and environmental variable
# u_e : matrix of dimension D X R with niche optima for each species and environmental variable
# s_e : matrix of dimension D X R with niche breadth for each species and environmental variable
# alpha : scalar, average dispersal distance
# m : scalar, immigration probability
# c_0 : vector of dimension R with colonization probability at 0 interaction
# e_0 : vector of dimension R with extinction probability at 0 interaction
# c_max: vector of dimension R maximal colonization probability at infinite positive interactions
# e_min: vector of dimension R maximal extinction probability at infinite positive interactions
# d_c: scalar, strength of the effect of interactions on the colonization probability
# d_e: scalar, strength of the effect of interactions on the extinction probability
# A: matrix of dimension R X R with the effect of each species (column) on other species (rows)

#########################
# Example 1
# Pure metapopulation dynamics, no sorting
#########################

# # Number of patches
N = 100

# # Number of environmental variables
D = 1

# # Regional species richness
R = 25

# # Effect of the environment on colonization
u_c = matrix(nr = D, nc = R)
u_c[1,] = seq(0.1,0.9,length=R)
u_c[2,] = rep(0.5, R)
s_c = matrix(Inf, nr = D, nc = R)

# # Effect of the environment on extinction
u_e = matrix(nr = D, nc = R)
u_e[1,] = rep(0.5, R)
u_e[2,] = rep(0.5, R)
s_e = matrix(Inf, nr = D, nc = R)

# # Mean dispersal
alpha = 0.3

# # Immigration
m = 0.001

# # Colonization function
c_0 = rep(0.8, R) # Colonization at 0 interactions
c_max = rep(0.8, R) # Colonization at max interactions

# # Extinction function
e_0 = rep(0.2, R) # Extinction at 0 interactions
e_min = rep(0.05, R) # Exinction at max interactions

# # Sensitivity to interactions
d_c = 0
d_e = 0

# # Interaction matrix
A = matrix(0, nr = R, nc = R)

# # Collect all parameters into a single list
pars = list(u_c = u_c, u_e = u_e, s_c = s_c, s_e = s_e, alpha = alpha, m = m, 
c_0 = c_0, e_0 = e_0, c_max = c_max, e_min = e_min, d_c = d_c, d_e = d_e, A = A)

save(pars, file = "metapop.Rds")

#########################
# Example 2
# Species sorting with dispersal limitations
#########################

# # Number of patches
N = 100

# # Number of environmental variables
D = 1

# # Regional species richness
R = 25

# # Effect of the environment on colonization
u_c = matrix(nr = D, nc = R)
u_c[1,] = seq(0.1,0.9,length=R)
u_c[2,] = rep(0.5, R)
s_c = matrix(0.5, nr = D, nc = R)

# # Effect of the environment on extinction
u_e = matrix(nr = D, nc = R)
u_e[1,] = rep(0.5, R)
u_e[2,] = rep(0.5, R)
s_e = matrix(Inf, nr = D, nc = R)

# # Mean dispersal
alpha = 0.3

# # Immigration
m = 0.001

# # Colonization function
c_0 = rep(0.8, R) # Colonization at 0 interactions
c_max = rep(0.8, R) # Colonization at max interactions

# # Extinction function
e_0 = rep(0.2, R) # Extinction at 0 interactions
e_min = rep(0.05, R) # Exinction at max interactions

# # Sensitivity to interactions
d_c = 0.8
d_e = 0

# # Interaction matrix
A = matrix(0, nr = R, nc = R)

# # Collect all parameters into a single list
pars = list(u_c = u_c, u_e = u_e, s_c = s_c, s_e = s_e, alpha = alpha, m = m, 
c_0 = c_0, e_0 = e_0, c_max = c_max, e_min = e_min, d_c = d_c, d_e = d_e, A = A)

save(pars, file = "sorting.Rds")

#########################
# Example 3
# Species sorting on colonization + dispersal
#########################

# Number of patches
N = 100

# Number of environmental variables
D = 1

# Regional species richness
R = 25

# Effect of the environment on colonization
u_c = matrix(nr = D, nc = R)
u_c[1,] = seq(0.1,0.9,length=R)
u_c[2,] = rep(0.5, R)
s_c = matrix(1, nr = D, nc = R)

# Effect of the environment on extinction
u_e = matrix(nr = D, nc = R)
u_e[1,] = rep(0.5, R)
u_e[2,] = rep(0.5, R)
s_e = matrix(Inf, nr = D, nc = R)

# Mean dispersal
alpha = 0.3

# Immigration
m = 0.001

# Colonization function
c_0 = rep(0.8, R) # Colonization at 0 interactions
c_max = rep(0.8, R) # Colonization at max interactions

# Extinction function
e_0 = rep(0.2, R) # Extinction at 0 interactions
e_min = rep(0.05, R) # Exinction at max interactions

# Sensitivity to interactions
d_c = 0.8
d_e = 0

# Interaction matrix
A = matrix(-1, nr = R, nc = R)
# 25% of links possible
A[runif(R*R,0,1) > 0.25] = 0
diag(A) = 0

# Collect all parameters into a single
pars = list(u_c = u_c, u_e = u_e, s_c = s_c, s_e = s_e, alpha = alpha, m = m, 
c_0 = c_0, e_0 = e_0, c_max = c_max, e_min = e_min, d_c = d_c, d_e = d_e, A = A)

save(pars, file = "sorting_interactions.Rds")

#########################
# Example 4
# Neutral dynamics
#########################

# Number of patches
N = 100

# Number of environmental variables
D = 1

# Regional species richness
R = 25

# Effect of the environment on colonization
u_c = matrix(nr = D, nc = R)
u_c[1,] = seq(0.1,0.9,length=R)
u_c[2,] = rep(0.5, R)
s_c = matrix(Inf, nr = D, nc = R)

# Effect of the environment on extinction
u_e = matrix(nr = D, nc = R)
u_e[1,] = rep(0.5, R)
u_e[2,] = rep(0.5, R)
s_e = matrix(Inf, nr = D, nc = R)

# Mean dispersal
alpha = 0.3

# Immigration
m = 0.001

# Colonization function
c_0 = rep(0.8, R) # Colonization at 0 interactions
c_max = rep(0.8, R) # Colonization at max interactions

# Extinction function
e_0 = rep(0.2, R) # Extinction at 0 interactions
e_min = rep(0.05, R) # Exinction at max interactions

# Sensitivity to interactions
d_c = 0.8
d_e = 0

# Interaction matrix
A = matrix(-1, nr = R, nc = R)

# Collect all parameters into a single list
pars = list(u_c = u_c, u_e = u_e, s_c = s_c, s_e = s_e, alpha = alpha, m = m, 
c_0 = c_0, e_0 = e_0, c_max = c_max, e_min = e_min, d_c = d_c, d_e = d_e, A = A)

save(pars, file = "neutral.Rds")

