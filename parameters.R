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


# The following example has competitive interactions, species sorting over two gradients and dispersal
N = 100
D = 2
R = 5

u_c = matrix(nr = D, nc = R)
u_c[1,] = seq(0.1,0.9,length=R)
u_c[2,] = rep(0.5, R)

u_e = matrix(nr = D, nc = R)
u_e[1,] = rep(0.5, R)
u_e[2,] = seq(0.9,0.1,length=R)

s_c = matrix(0.2, nr = D, nc = R)
s_e = matrix(0.2, nr = D, nc = R)

alpha = 0.3
m = 0.001

c_0 = rep(0.2, R)
e_0 = rep(0.1, R)

c_max = rep(0.4, R)
e_min = rep(0.05, R)

d_c = 0.5
d_e = 0.5

A = matrix(-1, nr = R, nc = R)
diag(A) = 0

pars = list(u_c = u_c, u_e = u_e, s_c = s_c, s_e = s_e, alpha = alpha, m = m, 
c_0 = c_0, e_0 = e_0, c_max = c_max, e_min = e_min, d_c = d_c, d_e = d_e, A = A)



