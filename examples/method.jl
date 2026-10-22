using NLSolvers, ForwardDiff, Clapeyron
import Clapeyron: ESElectrolyteModel, find_water_indx, Solvers

# `solid` optionally overrides the solid-phase parameters, given per nucleotide as
# (ΔG [kJ/mol], ΔH [kJ/mol], m_c+ (added cations), m_H2O). If `nothing`, the fitted
# values from the paper (Table SI.2) are used.
function build_mRNA_model(sequence, solvents=["water"], salts=["NaCl"]; solid=nothing)
    if typeof(sequence) <: String
        nA = count(==('A'),sequence)
        nC = count(==('C'),sequence)
        nG = count(==('G'),sequence)
        nU = count(==('U'),sequence)
    else
        nA, nC, nG, nU = sequence
    end
    

    # From these, specify the number of groups in the mRNA sequence
    naCH = nA*1+nC*2+nU*1
    naN2 = nA*2
    naCNH2 = nA*1 + nC*1 + nG*1
    nafC = nA*2 + nG*2
    naCH_5 = nA*1 + nG*1
    naN_5 = nA*1 + nG*1
    ncN = nA*1 + nG*1 + nC*1 + nU*1
    ncCH = nA*2 + nG*2 + nC*2 + nU*2
    ncyO = nA*1 + nG*1 + nC*1 + nU*1
    ncCHOH = nA*1 + nG*1 + nC*1 + nU*1
    nCH2 = nA*1 + nG*1 + nC*1 + nU*1
    nPO4 = nA*1 + nG*1 + nC*1 + nU*1
    naN = nG*1 + nC*1
    ncNH = nU*1 + nG*1
    ncCO = nC*1 + nG*1 + nU*2
    naC = nU*1
    nCH3 = nU*1

    names = []
    _solvents = []

    for i in 1:length(solvents)
        if solvents[i] == "water"
            append!(_solvents, [("water",["H2O"=>1])])
            append!(names, ["water"])
        elseif solvents[i] == "PEG6k"
            append!(_solvents, [("PEG6k",["CH2OH"=>2,"cO"=>135,"CH2OE"=>270])])
            append!(names, ["PEG6k"])
        elseif solvents[i] == "ethanol"
            append!(_solvents, [("ethanol",["CH3"=>1,"CH2OH"=>1])])
            append!(names, ["ethanol"])
        elseif solvents[i] == "isopropanol"
            append!(_solvents, [("isopropanol",["CH3"=>2,"CHOH"=>1])])
            append!(names, ["isopropanol"])
        elseif solvents[i] == "acetone"
            append!(_solvents, [("acetone",["CH3COCH3"=>1])])
            append!(names, ["acetone"])
        else
            append!(_solvents, solvents[i]) # Default to water if unknown
            append!(names, [solvents[i][1]])
        end
    end

    _ions = []
    append!(_ions, [("mRNA", [
                "aCH"=>naCH,
                "aN2"=>naN2,
                "aCNH2"=>naCNH2,
                "afCp"=>nafC,
                "aCH_5"=>naCH_5,
                "aN_5"=>naN_5,
                "cN"=>ncN,
                "cCH"=>ncCH,
                "cyO"=>ncyO,
                "cCHOHn"=>ncCHOH,
                "CH2"=>nCH2,
                "PO4-"=>nPO4,
                "aN"=>naN,
                "cNH"=>ncNH,
                "cC=O"=>ncCO,
                "aCCH"=>naC,
                "CH3"=>nCH3])])
    append!(_ions, [("potassium",["K+"=>1])])
    append!(names, ["mRNA", "potassium"])
    for i in 1:length(salts)
        if salts[i] == "NaCl"
            append!(_ions, [("sodium",["Na+"=>1])])
            append!(_ions, [("chloride",["Cl-"=>1])])
            append!(names, ["sodium", "chloride"])
        elseif salts[i] == "NaAc"
            append!(_ions, [("sodium",["Na+"=>1])])
            append!(_ions, [("acetate",["CH3"=>1,"COO-"=>1])])
            append!(names, ["sodium", "acetate"])
        elseif salts[i] == "NH4Ac"
            append!(_ions, [("ammonium",["NH4+"=>1])])
            append!(_ions, [("acetate",["CH3"=>1,"COO-"=>1])])
            append!(names, ["ammonium", "acetate"])
        else
            error("Salt $salts[i] not recognized. Supported salts are: NaCl, NaAc, NH4Ac.")
        end
    end

    paramdir = joinpath(@__DIR__, "..", "parameters")
    paramfiles = joinpath.(paramdir, ["params_like.csv","params_charges.csv","params_unlike.csv","params_assoc.csv"])
    charge = vcat(zeros(Int,length(_solvents)), [-nPO4, 1, 1, -1])
    fluid = SAFTgammaEMie(_solvents,_ions;neutralmodel_userlocations=paramfiles,
                                        ionmodel_userlocations=paramfiles,
                                        RSPmodel_userlocations = [joinpath(paramdir,"params_charges.csv")],
                                        charge = charge,)

    if !isnothing(solid)
        Gform = solid[1]*1e3*nPO4
        Hform = solid[2]*1e3*nPO4
        nx = Int(round(solid[3]*nPO4))
        nw = Int(round(solid[4]*nPO4))
    elseif "Na" == salts[1][1:end-2]
        Gform = 17.7124*1e3*nPO4
        Hform = 33.8432*1e3*nPO4
        nx = Int(round(0.50716*nPO4))
        nw = Int(round(4.8986*nPO4))
    else
        Gform = 17.5102*1e3*nPO4
        Hform = 23.3154*1e3*nPO4
        nx = Int(round(0.7334*nPO4))
        nw = Int(round(14.9451*nPO4))
    end
    solid = SolidKs(["water",salts[1][1:end-2]*".mRNA"]; userlocations = (;
                # Gform = [0.0,
                # 38181000.0],
                Gform = [0.0,
                Gform],
                # Hform = [0.0,
                # 82000000.0],
                Hform = [0.0,
                Hform],
                Tref =  [298.15,
                        298.15]))

    mapping = [(("water",1),)=>(("water",1)),
            (("water", nw),("mRNA",1),("potassium",nPO4-nx),(salts[1][1:end-2],nx))=>((salts[1][1:end-2]*".mRNA",1))]
    

    model = CompositeModel(names; mapping = mapping, fluid = fluid,
                                            solid = solid)

    id = (model.fluid.neutralmodel.groups.flattenedgroups .!= _ions[3][2][end][1]) .&
         (model.fluid.neutralmodel.groups.flattenedgroups .!= _ions[4][2][end][1]) .&
         (model.fluid.neutralmodel.groups.flattenedgroups .!= "H2O") .&
         (model.fluid.neutralmodel.groups.flattenedgroups .!= "PO4-") .&
         (model.fluid.neutralmodel.groups.flattenedgroups .!= "CH2OH") .&
         (model.fluid.neutralmodel.groups.flattenedgroups .!= "CH2OE") .&
         (model.fluid.neutralmodel.groups.flattenedgroups .!= "cO")

    model.fluid.neutralmodel.params.epsilon.values[end-1,id] .*= 0.
    
    model.fluid.neutralmodel.params.epsilon.values[id,end-1] .*= 0.
    

    if _ions[4][1] == "chloride"
        model.fluid.neutralmodel.params.epsilon.values[end,id] .*= 0.
        
        model.fluid.neutralmodel.params.epsilon.values[id,end] .*= 0.
        
    end
    return model
end


function sle_solubility_mrna(model::CompositeModel{F,S},p,T,z;solute=nothing,x0=nothing) where F <: EoSModel where S <: Clapeyron.SolidKsModel
    mapping = model.mapping
    if isnothing(solute)
        solute = model.solid.components
    end
    p = p*one(eltype(model))
    T = T*one(eltype(model))
    sol = zeros(length(solute),length(model.components))
    idxs = convert(Vector{Int},indexin(solute,model.solid.components))
    idx_sol = zeros(Bool,length(model.solid.components))
    idx_sol[idxs] .= true
    for i in 1:length(solute)
        idx_sol_s = zeros(Bool,length(model.solid.components))
        idx_sol_s[model.solid.components .==solute[i]] .= true

        #TODO: express this in terms of melting_temperature

        idx_sol_l = zeros(Bool,length(model.fluid.components))
        solute_l = mapping[idx_sol_s][1]
        ν_l = [solute_l[1][i][2] for i in 1:length(solute_l[1])]
        solute_l = [solute_l[1][i][1] for i in 1:length(solute_l[1])]

        solute_l = ["mRNA","FLuc","COVID","seqA","seqB","seqC","poly(A)","poly(G)","poly(C)","poly(U)","potassium"]
        for i in solute_l
            idx_sol_l[model.fluid.components .== i] .= true
        end

        solute_l_eq = ["mRNA","FLuc","COVID","seqA","seqB","seqC","poly(A)","poly(G)","poly(C)","poly(U)","potassium","sodium","ammonium","water"]
        idx_sol_l_eq = zeros(Bool,length(model.fluid.components))
        for i in solute_l_eq
            idx_sol_l_eq[model.fluid.components .== i] .= true
        end
        idx_solv = zeros(Bool,length(model.fluid.components))
        # if length(solute_l) == length(model)
        #     idx_solv[findfirst(idx_sol_l)] = true
        # else
            idx_solv[.!(idx_sol_l)] .= true
        # end

        solid_r,idx_sol_r = index_reduction(model.solid,idx_sol_s)
        μsol = chemical_potential(solid_r,p,T,[1.])
        # println(μsol)

        if isnothing(x0)
            x0 = Clapeyron.x0_sle_solubility(model,p,T,z,idx_solv,idx_sol_l,ν_l,μsol)
        end
        # idx_solv = model.components 
        # idx_solv =Bool.([1,0,0,1,1])
        # idx_sol_l =Bool.([0,1,1,0,0])

        f!(F,x) = obj_sle_solubility_mrna(F,model,p,T,z,exp10(x[1]),idx_sol_l,idx_sol_s,idx_solv,idx_sol_l_eq,ν_l)
        results = Clapeyron.Solvers.nlsolve(f!,x0,LineSearch(Newton()),NEqOptions(f_abstol=1e-6,f_reltol=1e-8),ForwardDiff.Chunk{1}())
        # println(results)
        sol[i,:] .=  z
        sol[i,.!(idx_solv)] .= exp10(Solvers.x_sol(results)[1]).*[1,sum(ν_l[2:end])]
        sol[i,:] ./= sum(sol[i,:])
    end
    if length(solute) == 1
        return sol[1,:]
    else
        return sol
    end
end

function obj_sle_solubility_mrna(F,model::CompositeModel{L,S},p,T,zsolv,solu,idx_sol_l,idx_sol_s,idx_solv,idx_sol_l_eq,ν_l) where L <: EoSModel where S <: Clapeyron.SolidKsModel
    # println(Clapeyron.Solvers.primalval(solu))
    z = zeros(typeof(solu),length(model.fluid))
    z .= deepcopy(zsolv)
    z[.!(idx_solv)] .= solu.*[1,sum(ν_l[2:end])]
        z[idx_sol_l] .+= solu*[1,sum(ν_l[2:end])]
    # z[idx_solv] .= zsolv[idx_solv]
    z ./= sum(z)
    # println(Clapeyron.Solvers.primalval.(z))
     

    if typeof(model.fluid) <: ESElectrolyteModel
        zref = ones(length(model.fluid))*1e-30
        idx_water = find_water_indx(model.fluid)
        zref[idx_water] = 1.0
        # zref = zeros(length(model.fluid))

        # ineutral = model.fluid.charge .== 0
        # zref[.!(ineutral)] .= 1e-30
        # zref[ineutral] .= zsolv[ineutral]
        zref ./= sum(zref)
        vref = volume(model.fluid,p,T,zref)
        μref = Clapeyron.VT_chemical_potential_res(model.fluid,vref,T,zref) .- Rgas()*T*log(vref*p/(Rgas()*T*sum(zref)))


        v = volume(model.fluid,p,T,z; vol0 = vref)
        μ = Clapeyron.VT_chemical_potential_res(model.fluid,v,T,z) .- Rgas()*T*log(v*p/(Rgas()*T*sum(z))) + Rgas()  * T * log.(z)

        μliq = (μ - μref)[idx_sol_l_eq]
    else
        pure   = split_pure_model(model.fluid)
        μ_mixt = chemical_potential(model.fluid, p, T, z)
        μ_ref = gibbs_free_energy.(pure, p, T)
        μliq = (μ_mixt - μ_ref)[idx_sol_l]
    end

    solid_r,idx_sol_r = index_reduction(model.solid,idx_sol_s)
    μsol = chemical_potential(solid_r,p,T,[1.])
    μliq = sum(μliq.*ν_l)
    F[1] = (μliq - μsol[1])/(Rgas()*T)

    # println(Clapeyron.Solvers.primalval(F[1]))
    return F
end