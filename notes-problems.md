# PROBLEMS

- Effet de l'environnement pour les simulations sans effet de l'environnement
    + Variance expliquée seulement; valeurs ne font pas de sens (e.g. le R2 = 0.04 pour Compet, ainsi que la variation est de 0.45 et 0.55 pour seulement effet environnement)
    
- Autocorrélation spatiale limitée
    + Données non corrélées ?
    + Modèle pas assez souple ?
    -- pourtant le R2adj d'une partition de variance avec RDA est entre 0.3 et 0.4 pour les modèles neutre et metapop (tel qu'attendu)

- Effet de co-distribution est quasi inexistant pour le modèle de compétition 
    + Peut être noyé dans les autres interactions qui sont inexistantes

- Pratiquement tous les scénarios ont la même signature, sauf le species sorting

- Compet
    + On ne parvient pas à trouver le signal des interactions sur la co-distribution, même si on utiliser seulement l'effet aléatoire 

- Neutral
    + effet d'environnement
    + effet de co-distribution for dans covariance, mais pas dans l'effet aléatoire

- Sorting
    --> cas parfait

- Metapo 
    --> 


# TODO

Case study: simulations COMPET. Il y a une structure dans ces données que des variables latentes devraient parvenir à révéler. Ce n'est pas le cas en ce moment. 

HYP - 

1. L'application du code HMSC est inappropriée

2. Il n'y a pas assez de signal dans les données
    -- essayer avec un cas à 3 espèces ?


3. Un autre package pourrait trouver du signal

