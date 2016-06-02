# Instructions to run the simulation model


## Scripts

The following scripts are provided in the package to run the simulation model:

**main.R** : the main script running the model

**functions.R**: functions that are necessary to compute the colonization and extinction

**landscape.R**: functions that are necessary to generate landscape properties (spatial location of spatial units & environmental characteristics)

**parameters.R**: definition of the parameters

**example.R**: example of a single simulation run


## Playing with different scenarios

It is straightforward to run different simulation scenarios, with as many species as the user wants. The time to run the model scales exponentially with the number of patches that are simulated (see the figure in the depot). I haven't done the check, but I suspect that it scales more or less linearly with the number of species in the regional pool. 

Here are examples of how to tweak parameters to run different scenarios:

### Independent metapops

- No species-specific response to the environment
- All coefficients of the interaction matrix are set to 0
- All species do have the same colonization and extinction probabilities

### Neutral

- Same as the independent metapops, but very strong negative interactions
- All interactions do have the same strength
- Both the colonization and the extinction probabilities do respond to ecological interactions 

### Species sorting

- Environment is heterogeneous
- Species specific response to the environment
- A single gradient is possible and it could either be the colonization or the extinction probabilities that do vary with the environment
- Be aware that the more environmental gradients there is, the more difficult it will be to persist. The reason is the multiplication of gaussian responses. Each response is a probability, and so the multiplication of numbers smaller than 1 makes the final probability very low. It is recommended to increase the niche breadth if multiple gradients are imposed. 
- Species sorting will happen even if species are not interacting. However, adding reciprocal negative interactions will strengthen the species sorting. 

### Competition - colonization trade-off

- Environment is uniform and no specific response to it
- The matrix of interactions is asymetric: good competitors do have a negative effect on weak competitors, but not the inverse. 
- However, the max colonization rate is larger for the weak competitors
- The Tilman (1994) formulation could be recovered if there is a strong negative impact of the strong competitor on the weak competitor on both the colonization and the extinction probabilitlies

###  Predator prey interactions

- It is easy to implement predator-prey interactions with the addition of positive and negative interaction coefficients in the matrix **A**. 
- The extreme case of predator dependence on preys that was simulated in Gravel et al. (2011) can be recovered with high d values (indicating that the effect of interactions saturates after the presence of a single prey). Further, in this study there was no effect of the predator on the prey. Further, ecological interactions influenced both the colonization and the extinction probabilities. 