using Combinatorics, LinearAlgebra, SparseArrays, Arpack, DelimitedFiles, Distributions, Random
include("generate_basis.jl")
include("all_functions.jl")

#--------------fixed parameters----------------
Morb = 10
Npara = Int(Morb ÷ 2)
Nparb = 1# impurity
Vaa = 0.0
J = 1.0
pbc = false
hardcore = true
if hardcore; Vaa = 0.0; end
opp_hop = false
delta_b = 0.0

tlist = 10 .^ range(-3, 6, length=800)

#outdir = "stark_data/M10/entropy"
#if !isdir(outdir) mkpath(outdir) end

#--------------basis and initial state----------------
b_tc, tfb, torbs, D_tc = generate_basis_two_component(Npara, Nparb, Morb, hardcore)
psi0, i0 = create_density_wave_initial_state(Morb, b_tc, D_tc)

#--------------entropy decomposition----------------
function entanglement_entropy_half_chain_decomposed(psi_t, M, basis, D, tlist)
    cut = M ÷ 2
    left_orbitals= 1:(2*cut)
    right_orbitals = (2*cut + 1):(2*M)

    left_basis= [basis[i, left_orbitals]  for i in 1:D]
    right_basis = [basis[i, right_orbitals] for i in 1:D]

    left_tags  = Dict{Vector{Int}, Int}()
    right_tags = Dict{Vector{Int}, Int}()

    for i in 1:D
        if !haskey(left_tags, left_basis[i])
            left_tags[left_basis[i]] = length(left_tags) + 1
        end
        if !haskey(right_tags, right_basis[i])
            right_tags[right_basis[i]] = length(right_tags) + 1
        end
    end

    dim_left= length(left_tags)
    dim_right = length(right_tags)

    left_indx= [left_tags[left_basis[i]]  for i in 1:D]
    right_indx = [right_tags[right_basis[i]] for i in 1:D]

    Nleft_of_index = zeros(Int, dim_left)
    for (config, idx) in left_tags
        Nleft_of_index[idx] = sum(config)
    end

    S_tot = zeros(length(tlist))
    S_N   = zeros(length(tlist))
    S_C   = zeros(length(tlist))

    for t in 1:length(tlist)
        psi = psi_t[:, t]
        psi_matrix = zeros(ComplexF64, dim_left, dim_right)

        for i in 1:D
            i_left  = left_indx[i]
            i_right = right_indx[i]
            psi_matrix[i_left, i_right] += psi[i]
        end

        rho_left = psi_matrix * psi_matrix'

        e = real.(eigen(rho_left).values)
        S = 0.0
        for p in e
            if p > 1e-10
                S -= p * log(p)
            end
        end
        S_tot[t] = S

        diag_probs = real.(diag(rho_left))
        p_n = Dict{Int, Float64}()
        for idx in 1:dim_left
            n = Nleft_of_index[idx]
            p_n[n] = get(p_n, n, 0.0) + diag_probs[idx]
        end

        Sn = 0.0
        for p in values(p_n)
            if p > 1e-10
                Sn -= p * log(p)
            end
        end
        S_N[t] = Sn
        S_C[t] = S_tot[t] - S_N[t]
    end

    return S_tot, S_N, S_C
end
function build_impurity_tags(basis, D, Morb)
    bath_cols = 1:2:(2*Morb)
    imp_cols= 2:2:(2*Morb)

    bath_config = [basis[i, bath_cols] for i in 1:D]
    imp_config= [basis[i, imp_cols]  for i in 1:D]

    bath_tags = Dict{Vector{Int}, Int}()
    imp_tags  = Dict{Vector{Int}, Int}()
    for i in 1:D
        if !haskey(bath_tags, bath_config[i]); bath_tags[bath_config[i]] = length(bath_tags) + 1; end
        if !haskey(imp_tags,  imp_config[i]);  imp_tags[imp_config[i]]  = length(imp_tags)  + 1; end
    end

    bath_indx = [bath_tags[bath_config[i]] for i in 1:D]
    imp_indx= [imp_tags[imp_config[i]]   for i in 1:D]
    dim_bath = length(bath_tags)
    dim_imp= length(imp_tags)

    return bath_indx, imp_indx, dim_bath, dim_imp
end

function impurity_rho(psi, bath_indx, imp_indx, dim_bath, dim_imp, D)
    C = zeros(ComplexF64, dim_imp, dim_bath)   # C[imp_id, bath_id] = c_{bath,imp}
    for i in 1:D
        C[imp_indx[i], bath_indx[i]] += psi[i]
    end
    return C * C'   
end

function compute_S_imp_vs_time(Morb, basis, D, tlist, psi_t)
    bath_indx, imp_indx, dim_bath, dim_imp = build_impurity_tags(basis, D, Morb)

    Nt = length(tlist)
    S_vs_t = zeros(Float64, Nt)
    S_gap_vs_t = zeros(Float64, Nt)

    for t in 1:Nt
        rho_imp = impurity_rho(psi_t[:, t], bath_indx, imp_indx, dim_bath, dim_imp, D)
        rho_imp ./= tr(rho_imp)
        Energies = eigvals(Hermitian(rho_imp))
        Energies_sorted = sort(Energies, rev=true)
        S_gap_vs_t[t] = Energies_sorted[1] - Energies_sorted[2]

        S = 0.0
        for p in Energies
            if p > 1e-14
                S -= p * log(p)
            end
        end
        S_vs_t[t] = S
    end
    return S_vs_t, S_gap_vs_t
end
