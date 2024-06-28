# QuantumSE.jl

A prototype tool for symbolic execution of quantum programs (QSE) with symbolic stabilizer states.

By symbolizing the phases of stabilizer states, **symbolic stabilizer states**  enable us to use symbolic expressions
to characterize the possible adversarial errors in quantum error correction (QEC) programs.
In this way, QSE with symbolic stabilizer states facilitate the efficient analysis of QEC programs.

See [the arXiv paper](https://arxiv.org/abs/2311.11313) for more details.

## Evaluations

We have evaluated QuantumSE.jl in [the arXiv paper](https://arxiv.org/abs/2311.11313).

Clone this repo and cd to `./QuantumSE.jl/example`.

```bash
git clone https://github.com/njuwfang/QuantumSE.jl.git && cd ./QuantumSE.jl/example
```

## Installation

To fully utililize the capabilities of `QuantumSE` package, install [Julia 1.10+](https://julialang.org/downloads/) and `QuantumSE` and `Bitwuzla` packages respectively by following the instructions below. Its recommended to familiarize yourself with [Julia's REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), especially "julia", "pkg" and "shell" modes; and [Julia's environments](https://docs.julialang.org/en/v1/manual/code-loading/#Environments)

### Clone QuantumSE repo

Clone this repo and cd to `./QuantumSE.jl/`.
```bash
git clone https://github.com/njuwfang/QuantumSE.jl.git && cd ./QuantumSE.jl
```


### SMT Solver(Manual Installation)(Ubuntu 22.04, 20.04)

QuantumSE.jl uses [Bitwuzla](https://github.com/bitwuzla/bitwuzla) as the default solver.

1. Install required dependencies:
    
    ```bash
    sudo apt-get install python3 python3-pip pkg-config m4 libgmp-dev
    ```

    ```bash
    pip install meson ninja
    ```
2. Clone our forked repo of Bitwuzla in `example` folder and cd to it.

    > NOTE: If using ubuntu in WSL, open Julia REPL under root privilege and execute these commands under shell mode of Julia REPL.

    ```bash
    cd example/ && git clone https://github.com/njuwfang/bitwuzla-for-QuantumSE.git && cd bitwuzla-for-QuantumSE
    ```
3. Configure it with `--kissat`
    
    ```bash
    ./configure.py --kissat
    ```
4. cd to `./build` and build
    
    ```bash
    cd ./build
    ninja
    ```
5. Install `bitwuzla`.

    ```bash
    ninja install
    ```


### Activate QuantumSE package

In julia's REPL pkg mode in `QuantumSE.jl` folder, execute:

```bash
(@v1.10) pkg> activate .
(QuantumSE) pkg> instantiate
```


### Finding Bugs in QEC Programs

> NOTE: If running ubuntu in WSL, after installing `QuantumSE` and `Bitwuzla` in julia REPL, run `julia> include("example\TannerCode.jl")`.

1. Repetition codes:
    ```bash
    julia RepetitionCode.jl
    ```

2. Toric codes:
    ```bash
    julia ToricCode.jl
    ```
3. Quantum Tanner codes:
    ```
    julia TannerCode.jl
    ```
The performance results are stored in `.dat` files.

### Comparing [symQV](https://github.com/fabianbauermarquart/symQV)

1. Clone symQV's repo and install its dependencies.
    ```bash
    cd example && git clone https://github.com/fabianbauermarquart/symQV.git && ./symQV/install.sh
    ```
2. Run scripts.
    ```bash
    python3 comp_symqv.py
    ```

    ```bash
    julia comp_quantumse.jl
    ```

The performance results are stored in `.dat` files.

### Comparing [Stim](https://github.com/quantumlib/Stim)

1. Install Stim (as of 7/11/2023, the latest stable version is 1.12.0) with `pip`.
    ```bash
    pip install stim==1.12.0
    ```
2. Generate benchmark files (in .stim format).
    ```bash
    julia generate_randomcircuits.jl
    ```
3. Benchmark Stim and QuantumSE.jl.
    ```bash
    python3 benchmark_stim.py
    ```

    ```bash
    julia benchmark_quantumse.jl
    ```

    The benchmark results are stored in `.dat` files.
