# QC-notebooks

Interactive notebooks based on QuantumCircuits.jl focusing on circuit QED and simulation of superconducting quantum devices.

## Usage
You will need the latest versions of Julia and Pluto installed.

Clone the repository on your computer. Then, you can follow these steps to run a notebook:
* open a Terminal window in the `QC-notebooks` directory
* remove the Manifest.toml file: `rm Manifest.toml`
* start a Julia REPL in the same folder: `julia`
* activate the environment: `] activate .`
* resolve: `resolve`
* escape out of packages < ctrl C >
* `using Pluto`
* `Pluto.run()`

A browser window will open and you will see the Pluto interface. Then, you can open a notebook by entering its path, e.g. `notebooks/excitation-feedback.jl`
