---
title: A general metacommunity toy model
author: Dominique Gravel
---

# Background

I decided to turn to another version of the toy model, that is more general than the model we used with Fournier. My motivations for not using the model based on differential equations derived in Fournier's paper:

- We want to have the possibility to have multiple species in a single spatial unit. It is more realistic to do so, and in particular it could handle cases such as host-parasite interactions, which a single individual per microsite can not. Further, I don't like the idea that an environmental condition is attached to an individual, I prefer to link it to a location (microsite). 
- In any case, it is nonetheless possible to tweak the model to have a single occupant per unit time, or just a few. 
- I wanted to avoir the extreme case of zero-sum dynamics. ZS generates negative covariance by default, a property that is not desirable. I wish the model could span all ranges of spatial covariances, from negative to positive. 
- Turning a system of differential equations into a stochastic simulation model is a tricky mathematical operation. Some approximation errors could arise by the discretization in space and time when coding the simulation model, and I would prefer to formulate the model to represent exactly what we simulate. 
- Investigation of a system of differential equations requires to solve the model at equilibrium, which seems problematic to some people in the group. We had many discussions over this, and while I believe it is not critical for what we'll do, I prefer to be cleaner and formulate the model in a way that does not require to find the analytical solution. 

Problems do arise as with simulate multiple species per parch. The occupancy of a species is a marginal probability, across all possible combinations of species presence-absence. Deriving the covariance between species requires to keep trakc of all of the conditional occurrence probabilities. There are two solutions to this problem:

- We build up a markov chain model in which we keep track all possible states (as We did in Cazelles et al. 2016)
- We run stochastic simulations

Since we are anyway interested in running spatially explicit simulations, we will adopt the stochastic simulation approach. We nonetheless look for a solution to the expected equilibrium occurrence probability of a species in a location, given the environment, the neighborhood and the community

# Model description
## Generalities about the metacommunity dynamics

We run metacommunity dynamics over a spatially heterogeneous landscape, with multiple environmental variables that could either be totally randomly distributed or spatially autocorrelated. Each location has a coordinate, and all possible coordinates are feasible. In other words, we are not restricting ourselves to a lattice or some other kind of regular spatial arrangement of spatial units. A spatial unit could be empty, occupied by a single species or several ones. It is up to us to determine if a location holds one or a few individuals (e.g. 5x5m cells in a forest) or entire populations (e.g. link in ponds). Spatial dynamics occurs as a result of colonization events, in both empty patches and patches that are occupied by other species, and extinction events. Species occurrence is a result of a dynamic balance between these events. Ecological interactions could impact both the colonization and the extinction probabilities. For instance, the presence of a competitor pre-empting a patch would reduce the colonization probability by another competitor. Alternatively, the presence of a predator in a patch could increase the extinction probability of its prey. Similarly, the environment could influence both the colonization and the extinction probabilities.  

# Definitions

We define $X_{i,z,t}$ a stochastic variable representing the occurrence of species $i$ at location $z$ and time $t$. $X_{i,z,t}$ takes a value of 1 when the species is present and a value of 0 when it is absent. Similarly, we define $\mathbf{Y}_{z,t} = (X_{1,z,t}, X_{2,z,t},..., X_{R,z,t})$ a vector containing the presence-absence of each species from the regional pool $R$.

## Colonization

We consider a discrete-time markovian process to represent the dynamics of presence-absence of all species. 

The challenge we face is to incorporate the effect of dispersal, environmental filtering and ecological interactions in such a way we could cover all four paradgims of Metacommunity 1.0, along with other types of spatial dynamics such as predator-prey interactions (Gravel et al. 2011a,b), priority effects (REF) or mutualistic interactions (Gilarranz2015). 
Followingly, a colonization event from time $t$ to $t+\Delta$ corresponds to:

$P(X_{i,z,t+\Delta t} = 1 | X_{i,z,t} = 0) = I_{i,z,t}S_{i,z,t}C_{i,z,t}$

Where $I_{i,z,t}$ is the amount of immigrants of species $i$ reaching location $z$ at time $t$, $S_{i,z,t}$ is the effect of environmental filtering on the probability of establishing a viable local population and $C_{i,z,t}$ is the effect of ecological interactions on the establishment probability. We note because that we represent a stochastic process, the product of these three functions has to be bounded between 0 and 1. We consequently define these quantities as:

$I_{i,z,t}  = \frac{\sum k(z,\omega)X_{i,y,t}}{\sum k(z,\omega)}$

which is a weighted average of species $i$ occurrence probability in the neighborhood of $z$. The function $k(z,\omega)$ is a dispersal kernel that depends on the location of the locality $z$ and the neighbour $\omega$. For convencience, we will consider an exponential function of the eucliedean distance between localities. We add to the kernel a low distance and neighborhoud -independent constant $m$ in order to account from immigration from outside the simulated metacommunity. This assumption is required to prevent total extinction by drift under perfectly neutral dynamics. 

The effect of the environment is given by a product of the establishment performance over all environmental variables $E_n$: 

$S_{i,z,t} = \prod f(E_{n,z}, \mu_{i,n}, \sigma_{i,n})$

In our simulations, we will consider for convencience that the function $f$ has a simple gaussian form, for all species and all environmental variables. The function is therefore bounded between 0 and 1.  

Finally we have to deal with all possible ecological interactions. We start by representing the interaction network by a community matrix $\mathbf{A}$ of $R$ species that we incorporate into the model. The elements $\alpha_{ij}$ of $\mathbf{A}$ quantify the effect of species $j$ on the dynamics of species $i$. When $\alpha_{ij}$ is negative, the colonization probability of species $i$ decreases and/or its extinction probability increases when $j$ is found locally. Inversely, when $\alpha_{ij}$ is positive, the colonization probability increases and/or the extinction probability decreases. To account for the cumulative effects of local interactions on transition probabilities, we make coloniza tion and extinction probabilities community dependent. As explained above, at a time $t$, the $\mathbf{Y}_{z,t}$ vector gives the local assemblages. We calculate the sum of interactions at any time and for each species as $v = \mathbf{A}_{z,t}\mathbf{Y}_{z,t}$. transpose operator). Our approach can be interpreted as a spatial analogue to the generalized Lotka–Volterra model because it takes into account the impact of the whole network of interactions on each species dynamics and can deal with any type of interaction. We now define the function $C_{i,z,t} = g(v_{i,z,t})$ representing the total effect of ecological interactions on the colonization probability. For convenience, we will use a sigmoid function, with $g$ ranging between $c_{min}$ at high negative interactions and $c_{max}$ at high positive interactions. $c_{max}$ should be interpreted as the maximal colonization probability when the environmental conditions are optimal and their is no dispersal limitations.   

## Extinction probability

The definition of the extinction probability follows exactly the same rules as above, except that extinction is independent of composition in the neighbourhood. We follow the same logic to define the effect of ecological interactions and of variation in the environment.  Consequently, we get the markovian process: 

$P(X_{i,z,t+\Delta t} = 0 | X_{i,z,t} = 1) = M_{i,z,t}E_{i,z,t}$

Where $M_{i,z,t}$ and $E_{i,z,t}$ are respectively the responses of the extinction probability to the local environment and to ecological interactions. The difference with the above functions for colonization is that the extinction probability has to be larger when interactions are negative and smaller when they are positive, and that the extinction rate should be minimal (instead of maximal) at environmental optimum. 


## Occurrence probability

The probability of observing a species at time $t + \Delta t$ in a location, given the environment, the neighborhood and the composition in that location:

$P(X_{i,z,t+\Delta t} = 1) = P(X_{i,z,t+\Delta t} = 1 | X_{i,z,t} = 0)P(X_{i,z,t} = 0) + P(X_{i,z,t+\Delta t} = 0 | X_{i,z,t} = 1)P(X_{i,z,t} = 1)$

It is convenient to derive what is the equilibrium occurrence probability for a given location. Equilibrium means the occurrence probability is the same at two moments in time. In other words:

    P(X_{i,z,t} = 1) = P(X_{i,z,t+\Delta t} = 1)

For simplifaction of the notation, we'll define the equilibrium occurrence probability as $\overline{p}_{iz}$. After some re-arrangement of the equation for occurrence probability, we find that: 

$\frac{\overline{p}_{iz}}{(1-\overline{p}_{iz})} = \frac{I_{i,z,t}S_{i,z,t}C_{i,z,t}}{M_{i,z,t}E_{i,z,t}}$

Taking the log of this equation, we get a log-linear model of occurrence probability: 

$logit(\overline{p}_{iz}) = log(I_{i,z}) + log(S_{i,z}) + log(C_{i,z}) - log(M_{i,z}) - log(E_{i,z})$

It is important to note that this derivation requires that the occurrence probability is locally stationary, but it does not require that the occurrence probabilities across the entire landscape are themselves stationary. 

# Interpretations

It is possible to turn the switch on and off so that we could represent different archetypes. 

- collection of independent metapopulations (Hanski's approach)
- neutral dynamics (no species sorting, strong negative interactions)
- species sorting (no dispersal limitation, effect of the environment, strong negative interactions)
- mass effect + dispersal limitations
- competition-colonization
- priority effect + species sorting
- predator-prey dynamics
- mutualistic dynamics

**Table 1** Example of parameterization for typical examples

**Table 2** Expectations for E, S and covariances for archetypes




