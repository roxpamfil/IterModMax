# IterModMax version 1.0

This code determines "optimal" resolution and interlayer coupling parameters for multilayer modularity maximisation, as explained in the following reference:
A. R. Pamfil, S. D. Howison, R. Lambiotte, and M. A. Porter. [Relating modularity maximization and stochastic block models in multilayer networks](https://arxiv.org/abs/1804.01964). arXiv preprint arXiv:1804.01964, 2018. 

## Dependencies

This code uses GenLouvain2.1, available at https://github.com/GenLouvain/GenLouvain. You need to add this code to the Matlab path for IterModMax to be able to use it.

## How to run

There are three main functions, which implement the iterative modularity maximisation method for different types of multilayer structure. These functions are `it_mod_max_temporal.m`, `it_mod_max_multilevel.m`, and `it_mod_max_multiplex.m`.

For temporal and multiplex networks, the mandatory inputs are a cell array `A` such that each entry is the adjacency matrix for the corresponding layer, and initial guesses for the parameters `gamma0` and `omega0`. Sample usage:

```matlab
[gamma, omega, ~, S, Q, converged] = it_mod_max_temporal(A, gamma0, omega0)
```

For more details, and for a description of optional inputs, see the documentation inside each Matlab file. Some examples are provided in the `Examples/` subdirectory.

