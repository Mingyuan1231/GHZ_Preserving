# GHZPreserving.jl usage guide

GHZPreserving.jl simulates GHZ-preserving Clifford operations on binary stabilizer-sign labels. This guide describes the current Julia implementation accompanying [GHZ-Preserving Gates and Optimized Distillation Circuits](https://arxiv.org/abs/2510.25854), by Mingyuan Wang, Guus Avis, and Stefan Krastanov.

For installation and paper/software citations, see the [repository README](https://github.com/Mingyuan1231/GHZ_Preserving#readme). The software archive is available at [DOI: 10.5281/zenodo.17470505](https://doi.org/10.5281/zenodo.17470505).

## State representation

```julia
using GHZPreserving

n = 3  # qubits per GHZ state
m = 2  # number of GHZ copies held in the register
state = GHZState(n, m)

@assert state.qubit_num == n
@assert state.ghz_num == m
@assert length(state.phases) == n * m
```

`state.phases` is a `BitVector` with one block of `n` stabilizer-sign bits per GHZ copy. `false` represents a positive stabilizer sign and `true` a negative sign. `GHZState(n, m)` initializes every sign to positive.

One value of `GHZState` represents a product of GHZ-basis states, rather than an entire mixed-state probability distribution. Sampling many such values represents a GHZ-diagonal input ensemble:

```julia
using Random

Random.seed!(42)
sample = rand(GHZState, 3, 2, 0.9)
```

For each copy, this constructor selects the all-positive target label with probability `0.9`; otherwise it samples uniformly from the non-target labels. A single sample has a definite basis label. The probability parameter describes the ensemble, not the fidelity of every individual sample.

## Homogeneous and bilocal gates

An H-group operation applies the same two-qubit gate at every node, between two selected GHZ copies. A B-group operation applies its local gate at an adjacent pair of nodes in the implementation's node ordering.

```julia
state = GHZState(3, 2)
apply!(state, Hgroup{3}(2, 1, 2))     # homogeneous CNOT, copy 1 -> copy 2
apply!(state, Bgroup{3}(4, 1, 2, 1))  # local CZ at nodes 1 and 2

@assert all(!, state.phases)
```

| Constructor | Arguments |
| --- | --- |
| `Hgroup{n}(g, a, b)` | H-group gate index `g`; GHZ copy indices `a` and `b` |
| `Bgroup{n}(g, a, b, j)` | B-group gate index `g`; copies `a` and `b`; node pair `j, j+1`, with `1 <= j < n` |
| `PauliGroup{n}(g, a, j)` | Pauli index `g`; GHZ copy `a`; node `j` |

The H-group indices are `1: SWAP`, `2: CNOT12`, `3: INV_CNOT21`, `4: INV_CNOT12`, `5: identity`, and `6: CNOT21`. The `INV_CNOT` names follow the source-code convention.

The B-group indices are `1: Phase1`, `2: CZ+Phase1`, `3: identity`, `4: CZ`, `5: Phase2`, `6: CZ+Phase2`, `7: BothPhase`, and `8: CZ+BothPhase`. Pauli indices are `1: X`, `2: Y`, `3: Z`, and `4: identity`.

Use valid, distinct GHZ copy indices for two-copy operations. The `{n}` parameter denotes the number of qubits per GHZ state, not the number of stored copies. Detailed gate construction is in [src/GHZPreserving.jl](https://github.com/Mingyuan1231/GHZ_Preserving/blob/main/src/GHZPreserving.jl).

## Measurements and Pauli errors

`GHZMeasure(n, basis, copy)` uses `basis=1` for an X-type syndrome and `basis=3` for a Z-type syndrome. Use these two bases; the measurement implementation does not provide a Y-basis path.

```julia
state = GHZState(3, 2)
state, syndrome = measure!(state, GHZMeasure(3, 3, 2))
@assert syndrome == [false, false]
```

The returned values are GHZ stabilizer-syndrome bits: one bit for X and `n-1` bits for Z. `measure!` resets the corresponding phase bits in place; it does not remove a copy from the register or itself implement a complete distillation protocol.

Independent Pauli errors can be sampled with:

```julia
state = GHZState(3, 2)
apply!(state, PauliNoiseOp(3, 1, 0.01, 0.01, 0.01))
```

Here the three probabilities are for X, Y, and Z on each qubit of copy 1. Use nonnegative probabilities whose sum is at most one. The remaining probability represents no error.

The package also contains `NoisyGHZMeasureNoisyReset` and status-based postselection routines. Check their X/Z-specific probability conventions in the source before using them in a numerical study; do not infer those conventions from the parameter name alone.

## Circuit-search scripts and parameter conventions

The genetic-search implementation is in [src/genetic_alg.jl](https://github.com/Mingyuan1231/GHZ_Preserving/blob/main/src/genetic_alg.jl), with a related search script in [src/CNOT_only.jl](https://github.com/Mingyuan1231/GHZ_Preserving/blob/main/src/CNOT_only.jl). They are not included automatically by `using GHZPreserving`.

| Search-script field | Meaning | Paper notation |
| --- | --- | --- |
| `n` | Total raw GHZ-state count | `N` |
| `q` | Qubits per GHZ state | `n` |
| `k` | Retained output count | `K` |
| `r` | Simultaneously available registers per node | `R` |
| `f_in` | Probability of the target GHZ-basis label in sampled inputs | Input fidelity |

In particular, the search scripts' `n` differs from the simulator constructor's `n`. In `noisify`, the gate parameter `p2` enters each Pauli-error probability as `(1-p2)/3`; it should not be read directly as the paper's gate-error probability. Inspect the measurement parameter separately because its conventions depend on the routine and basis.

The scripts contain top-level calls to `run!`, large Monte Carlo evaluations, parameter sweeps, plotting, and writes to result files. They also include exploratory sections that rely on interactive state or additional packages. Read and adapt the selected experiment before executing it; running an entire script is not a minimal smoke test or a guaranteed reproduction command.

Saved results include:

- [circuit_data_HF.txt](https://github.com/Mingyuan1231/GHZ_Preserving/blob/main/circuit_data_HF.txt): serialized `Individual`-style records with circuit operations and estimated performance.
- [examples/testing data/circuit_data.txt](https://github.com/Mingyuan1231/GHZ_Preserving/blob/main/examples/testing%20data/circuit_data.txt): additional stored circuit records.
- [examples/testing data/testing data.jl](https://github.com/Mingyuan1231/GHZ_Preserving/blob/main/examples/testing%20data/testing%20data.jl): plotting code with tabulated example outputs; it imports `Plots`, which is not a direct dependency in the root `Project.toml`.

The text outputs are Julia-style records, not a documented language-independent interchange format. Record the commit, parameter choices, sampling budget, and random-number setup when conducting a new search. Stochastic runs need not recover identical circuits.

## Cached permutations and scaling

The first H- or B-group construction for a given `n` generates basis states and permutation tables. Later constructions reuse module-level caches. Each permutation on two GHZ copies contains `2^(2n)` output labels; the cache stores multiple gate permutations.

The table access itself is constant time. The current implementation additionally copies slices of `BitVector`s, concatenates labels, converts between bits and integers, and writes the output bits. These steps depend on label length. Separate initialization, memory, and complete gate-application costs when measuring performance. This code does not implement a general non-Pauli density-matrix simulation or an end-to-end network timing model.

## Tests and documentation

The [test directory](https://github.com/Mingyuan1231/GHZ_Preserving/tree/main/test) contains gate, state, measurement, noise, and QuantumClifford-comparison checks. Several checks use substantial stochastic sampling; the small examples above are a more focused starting point for checking package loading and ideal gate operations.

To build this guide with the repository environment:

```sh
julia --project=. docs/make.jl
```

Documenter writes the HTML output to `docs/build/`. The Markdown source can also be read directly on GitHub.
