# Thermodynamic modelling of mRNA solubility
This repository contains the data, parameters and code related to the _Thermodynamic Modeling of mRNA with the Addition of Precipitants_ publication. 

![mRNA_solubility](/assets/solubility_fluc_water_salts.png)

The repository is divided into the following folders:
* `data`: Contains all experimental mRNA solubilities obtained within the work for both the FLUC and COVID sequences.
* `parameters`: Contains all SAFT-$\gamma$ Mie parameters that were regressed within the work.
* `examples`: Contains example code to obtain the solubility of an arbitrary mRNA sequence in a mixed solvent with added salt.

## Installing Julia

To download and install Julia, please go to
https://julialang.org/install

On most systems, it should be as simple as running:
```bash
curl -fsSL https://install.julialang.org | sh
```

## Installing dependencies
If you make your way to this repo and open the REPL (run julia in the command line), you can install all the dependencies by running:
```julia
julia> ] activate .

julia> instantiate
```
This will install Clapeyron and Plotting tools.