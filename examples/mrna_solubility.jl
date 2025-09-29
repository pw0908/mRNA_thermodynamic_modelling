using Clapeyron
include("method.jl")

# Specify temperature and pressure
T = 298.15  # Temperature in Kelvin
p = 1.01325e5  # Pressure in Pascals

# Specify the number of nucleobases in the mRNA sequence
nA = 526
nC = 691
nG = 587
nU = 308

model = build_mRNA_model((nA,nC,nG,nU), ["water"], ["NaCl"])

n = [1.0, 1e-30, 1e-30, 0.01, 0.01]  # Molar amounts of each component
z0 = n ./ sum(n)  # Mole fractions

z = sle_solubility_mrna(model, p, T, z0; solute=["Na.mRNA"], x0=[-7.])  # Calculate molar solubility

v = volume(model, p, T, z)  # Calculate molar volume
cs = z[3]*Clapeyron.mw(model.fluid)[3]/v # Convert to concentration in μg/mL