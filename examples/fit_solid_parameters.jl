using Clapeyron, CSV, DataFrames, Metaheuristics
include("method.jl")

# Fit the solid-phase parameters of the mRNA precipitate to experimental solubility data.
# The SAFT-γ Mie group parameters are held fixed; only the four solid-phase parameters
# (all per nucleotide) are fitted:
#   θ[1] = ΔG_mRNA† [kJ/mol], θ[2] = ΔH_mRNA† [kJ/mol],
#   θ[3] = m_c+ (added cations in the complex), θ[4] = m_H2O (water in the complex)

# Specify the number of nucleobases in the mRNA sequence (here, FLuc).
# A sequence string, e.g. "AUGGC...", can be given instead of the tuple.
sequence = (526, 691, 587, 308)  # (A, C, G, U)
solvents = ["water"]
salt = "NaCl"

# Experimental data: temperature T [K], salt mass fraction ws, and mRNA solubility out_c [μg/mL]
data = CSV.read("fluc_NaCl.csv", DataFrame; header=3)

p = 1e5

function solubility(model, T, ws)
    Mw_salt = sum(Clapeyron.mw(model.fluid)[end-1:end])
    m = ws/Mw_salt*1e3/(1-ws) # salt molality
    z0 = [1e3/Clapeyron.mw(model.fluid)[1], 1e-30, 1e-30, m, m]
    z0 ./= sum(z0)
    z = sle_solubility_mrna(model, p, T, z0; solute=[model.solid.components[2]], x0=[-7.])
    v = volume(model.fluid, p, T, z)
    return z[2]*Clapeyron.mw(model.fluid)[2]/v # μg/mL
end

# Weighted least-squares objective on relative deviations (Eq. 3 of the main text, equal weights)
function objective(θ)
    model = build_mRNA_model(sequence, solvents, [salt]; solid=θ)
    f = 0.0
    for r in eachrow(data)
        c = try
            solubility(model, r.T, r.ws)
        catch
            return 1e10 # failed solubility calculation: penalize this parameter set
        end
        f += ((r.out_c - c)/r.out_c)^2
    end
    return f
end

#           ΔG     ΔH    m_c+  m_H2O
lower  = [ 5.0,   5.0,  0.0,  0.0]
upper  = [30.0,  50.0,  1.0, 20.0]
bounds = [lower upper]'

logger(st) = println("Iteration ", st.iteration, ": f = ", st.best_sol.f, ", θ = ", st.best_sol.x)

method = ECA(; options=Options(iterations=100, seed=1))
results = Metaheuristics.optimize(objective, bounds, method; logger=logger)
θ = minimizer(results)

println("\nFitted solid-phase parameters (per nucleotide):")
println("ΔG    = ", θ[1], " kJ/mol")
println("ΔH    = ", θ[2], " kJ/mol")
println("m_c+  = ", θ[3])
println("m_H2O = ", θ[4])

# Compare the fitted model with the data
model = build_mRNA_model(sequence, solvents, [salt]; solid=θ)
c_model = [solubility(model, r.T, r.ws) for r in eachrow(data)]
AAD = 100*sum(abs.((data.out_c .- c_model)./data.out_c))/nrow(data)
println("AAD% = ", AAD)
