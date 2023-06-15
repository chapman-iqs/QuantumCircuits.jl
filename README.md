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
```
include("src/QuantumCircuits.jl")
using .QuantumCircuits
using .QuantumCircuits.Tests

# to run tests as testsets (comment out the test functions you don't want)
runtests(; functions = [
                         test_integration,
                         test_positive_trajectory,
                         test_single_timestep,
                         test_lindblad, 
                         test_timedelay
                        ])

# to run tests individually, as functions (you'll have to go into sourcecode to see what's available)
# for example:
returns_solution(; solve=rouchon)
lindblad_timedep(0.5, 0.2; solve=bayesian)
```

### Running examples
#### From the REPL
activate environment using `] activate .`
Then `include("examples/monitored_rabi.jl")` etc.
