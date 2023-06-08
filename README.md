# QuantumCircuits.jl

**QuantumCircuits.jl** is a numerical framework written in Julia for simulating quantum evolution, with a focus on superconducting circuits, and builds on [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl).

## Usage
To use and test the package while still in development,
1. Clone the package
2. Navigate to the location of the repository on your computer
3. Open a julia REPL session by typing `julia` in the terminal
4. Go into `Pkg` mode and activate / resolve the environment `] activate .`, `resolve`

### Running tests
#### From the package manager:
Run tests: `] test`

#### From the REPL
```include("src/QuantumCircuits.jl")
using .QuantumCircuits
using .QuantumCircuits.Tests
```

### Running examples
#### From the REPL
activate environment using `] activate .`
Then `include("examples/monitored_rabi.jl")` etc.
