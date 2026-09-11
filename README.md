# GHZPreserving.jl: GHZ Entanglement Distillation and Circuit Optimization

[![Paper](https://img.shields.io/badge/arXiv-2510.25854-b31b1b.svg)](https://arxiv.org/abs/2510.25854)
[![Software DOI](https://zenodo.org/badge/812791172.svg)](https://doi.org/10.5281/zenodo.17470505)

**GHZPreserving.jl** is a Julia research package for simulating GHZ-preserving Clifford operations and exploring **GHZ entanglement distillation**, also called **entanglement purification**. It represents Greenberger-Horne-Zeilinger (GHZ) basis states with binary stabilizer-sign labels and applies gates through cached permutations. The repository also contains genetic circuit-search scripts and saved circuit results for studying noisy operations and limited quantum-register capacity.

## Paper and software

This repository accompanies **[GHZ-Preserving Gates and Optimized Distillation Circuits](https://arxiv.org/abs/2510.25854)** by **Mingyuan Wang, Guus Avis, and Stefan Krastanov** (2025).

- **Research paper:** [arXiv:2510.25854](https://arxiv.org/abs/2510.25854) · [HTML full text](https://arxiv.org/html/2510.25854v1) · [paper DOI](https://doi.org/10.48550/arXiv.2510.25854)
- **Archived software:** [Zenodo DOI: 10.5281/zenodo.17470505](https://doi.org/10.5281/zenodo.17470505)
- **Usage guide:** [state labels, gate indices, measurements, and implementation notes](docs/src/index.md)

The paper develops a classification of GHZ-preserving operations using homogeneous (H), bilocal (B), and Pauli operations. It uses the resulting simulation framework to optimize GHZ-distillation circuits and compare them with pumping and nested protocols. See the paper for the mathematical results, benchmark conditions, and graph-state extensions.

## What is in this repository?

| Task | Implementation or material |
| --- | --- |
| Represent products of GHZ-basis states using binary phase labels | `GHZState` in [src/GHZPreserving.jl](src/GHZPreserving.jl) |
| Apply homogeneous, bilocal, and Pauli Clifford operations | `Hgroup`, `Bgroup`, `PauliGroup`, and `GHZGate` |
| Sample GHZ-diagonal inputs and Pauli errors | `rand(GHZState, ...)` and `PauliNoiseOp` |
| Work with X/Z measurement syndromes and measurement/reset routines | `GHZMeasure` and `NoisyGHZMeasureNoisyReset` |
| Explore genetic circuit optimization and register constraints | [src/genetic_alg.jl](src/genetic_alg.jl) and [src/CNOT_only.jl](src/CNOT_only.jl) |
| Inspect stored candidate circuits and numerical outputs | [circuit_data_HF.txt](circuit_data_HF.txt) and [examples/testing data](examples/testing%20data) |
| Compare operations with QuantumClifford representations | [test/](test/) and `toQCcircuit` |

## Install from the repository

The checked-in `Manifest.toml` was generated with **Julia 1.11.3** and records dependency versions. Start from the repository environment:

```sh
git clone https://github.com/Mingyuan1231/GHZ_Preserving.git
cd GHZ_Preserving
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=.
```

The repository name is `GHZ_Preserving`; the Julia package is loaded with `using GHZPreserving`.

## Minimal gate-simulation example

Run this in Julia after activating the repository environment:

```julia
using GHZPreserving

# Two copies of a three-qubit GHZ state; all stabilizer signs start positive.
state = GHZState(3, 2)

# H-group index 2: each node applies CNOT from copy 1 to copy 2.
apply!(state, Hgroup{3}(2, 1, 2))

# B-group index 4: nodes 1 and 2 each apply CZ between their local copies.
apply!(state, Bgroup{3}(4, 1, 2, 1))

@assert length(state.phases) == 6
@assert all(!, state.phases)
```

This example illustrates ideal GHZ-preserving gates. It does not perform noisy distillation or estimate output fidelity. The [usage guide](docs/src/index.md) explains the labels and constructors in more detail.

## Circuit optimization and saved results

The genetic-search scripts explore sequences of H-group gates, B-group gates, and measurements. Their parameters distinguish the number of raw GHZ states, qubits per GHZ state, retained outputs, and available registers. Candidate evaluation uses Monte Carlo sampling to estimate acceptance probability and conditional output fidelity.

**The optimization files are research scripts, not a one-command reproduction interface.** They contain top-level searches, parameter sweeps, plotting, and output-writing code. Inspect and adapt the relevant sections before running them; `using GHZPreserving` loads the simulator without launching those searches. The [usage guide](docs/src/index.md#circuit-search-scripts-and-parameter-conventions) identifies the parameter conventions and saved-result files.

## Simulation scope and computational cost

- A `GHZState` stores a product of GHZ-basis states. GHZ-diagonal noise is sampled as an ensemble of these labels; the package is not a general non-Pauli density-matrix simulator.
- H- and B-group constructors generate and cache permutation tables. Each two-copy permutation has `2^(2n)` entries for `n` qubits per GHZ state, so initialization and memory grow with `n`.
- A table lookup is constant time, but the current `apply!` implementation also slices, concatenates, encodes, decodes, and writes binary labels. The complete implementation should not be described as an unconditional constant-time update independent of `n`.
- The scripts study finite circuits. They do not simulate end-to-end quantum-network scheduling, storage times, or entanglement-generation latency.

## Citation

For the classification, simulation method, and distillation results, cite the paper:

```bibtex
@misc{wang2025ghzpreserving,
  title         = {{GHZ}-Preserving Gates and Optimized Distillation Circuits},
  author        = {Wang, Mingyuan and Avis, Guus and Krastanov, Stefan},
  year          = {2025},
  eprint        = {2510.25854},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
  doi           = {10.48550/arXiv.2510.25854},
  url           = {https://arxiv.org/abs/2510.25854}
}
```

If you use the code, also cite the [archived software release](https://doi.org/10.5281/zenodo.17470505) and record the Git commit used for your calculations. Machine-readable citation metadata is provided in [CITATION.cff](CITATION.cff).
