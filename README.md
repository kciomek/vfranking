# vfranking

The software is a package for [R environment](http://www.r-project.org "R project") and allows to analyse problem of constructing ranking on the set of alternatives (objects or actions) evaluated on multiple monotonic criteria (attributes). The ranking have to be consistent with provided preferences which can be expressed as pair-wise comparisons of alternatives (e.g., alternative _a_ is preferred over _b_, i.e., user requires _a_ to attain better position in ranking than _b_). The concept of additive value functions was used to model the preferences and finally to work out the recommendation - partial or complete ranking or stochastic characteristics related to alternatives and their possible positions in the ranking.

The model can be build with linear, piece-wise linear and general marginal value functions. Criteria can be required to be maximized (gain criteria) or minimized (cost criteria). The package allows to obtain the following results:

  * single value function calculated by maximizing epsilon (UTA method),
  * possible and necessary preference relation (Robust Ordinal Regression),
  * pairwise winning indices,
  * rank acceptability indices.

### Installation

To get the latest version, download it directly from GitHub with _devtools_:

    # install.packages("devtools")
    library(devtools)
    install_github("kciomek/vfranking")

### Tutorial

We have 4 objects evaluated on 2 criteria that have to be maximized (gain criteria - 'g'). We assume linear value functions (2 characteristic points per criterion). And one pair-wise comparison was provided as preference information: "_a<sub>2</sub>_ is strongly preferred over _a<sub>1</sub>_".

    > performances <- matrix(c(5, 2, 3, 7, 0.5, 0.9, 0.5, 0.4), ncol = 2)
    > performances
          [,1] [,2]
     [1,]    5  0.5
     [2,]    2  0.9
     [3,]    3  0.5
     [4,]    7  0.4

    > strongPreference <- rbind(c(1,2))

    > problem <- buildProblem(performances, c('g', 'g'), c(2, 2), strongPreference = strongPreference)
    > model <- buildModel(problem)

With _model_, it is possible to calculate results. To get a ranking imposed by a single value function obtained according to UTA method use function _maxEpsilonSolution_:

    > func <- maxEpsilonSolution(model)
    > func$ranks # positions in the ranking
    [1] 2 3 4 1
    
    > ranksToRanking(func$ranks)
    [1] "4-1-2-3"

To check preference relations between pair of alternatives with robustness analysis use:

    > isPreferred(model, 3, 2, necessarily=FALSE) # is a3 possibly preferred over a2?
    [1] TRUE
    > isPreferred(model, 2, 3, necessarily=FALSE) # is a2 possibly preferred over a3?
    [1] TRUE
    > isPreferred(model, 2, 4, necessarily=TRUE) # is a2 necessarily preferred over a4?
    [1] FALSE
    > isPreferred(model, 4, 2, necessarily=TRUE) # is a4 necessarily preferred over a2?
    [1] TRUE


To calculate stochastic results, it is needed to sample the parameters space:

    > samples <- sampleParameters(model, numberOfSamples=10000)
    > nrow(samples)
    [1] 10000

Then, you can calculate pairwise winning indices and ranking acceptability indices:

    > pwi(model, samples)
         [,1]   [,2]   [,3] [,4]
    [1,]    0 1.0000 1.0000    0
    [2,]    0 0.0000 0.5356    0
    [3,]    0 0.4644 0.0000    0
    [4,]    1 1.0000 1.0000    0
    
    > rai(model, samples)
         [,1] [,2]   [,3]   [,4]
    [1,]    0    1 0.0000 0.0000
    [2,]    0    0 0.5356 0.4644
    [3,]    0    0 0.4644 0.5356
    [4,]    1    0 0.0000 0.0000

