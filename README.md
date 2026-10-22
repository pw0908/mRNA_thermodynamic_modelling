# Thermodynamic modelling of mRNA solubility
This repository contains the data, parameters and code related to the _Thermodynamic Modeling of mRNA with the Addition of Precipitants_ publication. 

<img src="/assets/solubility_fluc_water_salts.png" alt="mRNA_solubility" width="400"/>

The repository is divided into the following folders:
* `data`: Contains all experimental mRNA solubilities obtained within the work: the FLuc and COVID sequences at 4 °C and 25 °C (in water, water+ethanol and water+PEG6k with added NaCl, NaAc or NH4Ac), and sequences A, B and C at 25 °C (in ~12% w/w PEG6k with added NaCl, NaAc or NH4Ac).
* `parameters`: Contains all SAFT-$\gamma$ Mie parameters that were regressed within the work.
* `examples`: Contains example code to obtain the solubility of an arbitrary mRNA sequence in a mixed solvent with added salt, and to refit the solid-phase parameters to new solubility data.

## Sequence information
The sequences studied in this work are proprietary. As SAFT-γ Mie accounts for sequence only through nucleotide composition, the number of each nucleotide (given below, and in Table SI.11 of the Supporting Information) is all that is needed to reproduce the calculations. N1-methyl-pseudouridine is used in place of uridine (U) in all sequences.

| Sequence | A | C | G | U (N1-mΨ) | Total |
|---|---|---|---|---|---|
| FLuc  | 526  | 691  | 587  | 308 | 2112 |
| COVID | 1102 | 1315 | 1060 | 801 | 4278 |
| A     | 502  | 532  | 537  | 415 | 1986 |
| B     | 524  | 539  | 480  | 398 | 1941 |
| C     | 502  | 569  | 535  | 389 | 1995 |

The solid-phase parameters (per nucleotide) used in the publication are given in Tables SI.2 (FLuc, COVID and all predictions) and SI.3 (sequences A, B and C in PEG6k). The default parameters in `examples/method.jl` are those of Table SI.2.

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
This will install Clapeyron, Metaheuristics and the other dependencies.

## Examples
Both examples are run from the `examples` folder:
```bash
cd examples
julia --project=.. mrna_solubility.jl
julia --project=.. fit_solid_parameters.jl
```
* `mrna_solubility.jl`: computes the solubility of the FLuc sequence in water with added NaCl. The mRNA is specified only through its nucleotide composition `(nA, nC, nG, nU)`; `build_mRNA_model` in `method.jl` also accepts a full sequence string.
* `fit_solid_parameters.jl`: fits the four solid-phase parameters (per nucleotide: ΔG, ΔH, number of added cations and number of water molecules in the precipitate) to experimental solubility data using Metaheuristics.jl, with all SAFT-γ Mie group parameters held fixed. The example uses the FLuc/NaCl data in `examples/fluc_NaCl.csv` (columns: temperature `T` in K, salt mass fraction `ws`, mRNA solubility `out_c` in μg/mL). To fit a different construct or salt, change `sequence`, `salt` and the data file at the top of the script. The fitted parameters can then be passed to `build_mRNA_model(sequence, solvents, salts; solid=θ)`.

Refitting is recommended when the solvent environment differs substantially from the training conditions (e.g., high PEG fractions), or for constructs with highly stable secondary or tertiary structure. For a new, largely disordered single-stranded mRNA under similar conditions, the default parameters can be used directly.