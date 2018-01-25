# sublinear-Li-Stephens
Performs an arithmetically equivalent computation to the forward algorithm for the Li Stephens haplotype hidden Markov model, but in average-case sub-O(nk) time where n = number of genetic sites, k = number of reference haplotypes versus O(nk^2) for the standard forward algorithm

Can also do forward-backward with time complexity < O(nk) but larger than for forward algorithm, where the difference in time complexity rises monotonically with the number of forward-backward matrix entries needed.

## tests
Tests implemented using Catch (https://github.com/philsquared/Catch)
