# IterModMax version 1.0

This code determines "optimal" resolution and interlayer coupling parameters for multilayer modularity maximisation.

## Dependencies

This code uses GenLouvain2.1, available at https://github.com/GenLouvain/GenLouvain.
For this code to be able to use it, place it in the relative path ../GenLouvain2.1.

## How to run

The main function for temporal networks is `it_mod_max_temporal.m`. It takes as mandatory inputs a cell array `A` such that each entry is the adjacency matrix for the corresponding layer, and initial guesses for the parameters `gamma0` and `omega0`.
For more details, and for a description of optional inputs, see the documentation inside the code.

Sample usage:

```matlab
[gamma, omega, beta, S, Q, converged] = it_mod_max_temporal(A, gamma0, omega0)
```
