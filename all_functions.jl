##<<<<<<<<<<<<<<<<<<<<<,kinetic term>>>>>>>>>>>>>>>>>>>>>>>>>>
function kinetic_component(M, basis, tag_fock_basis, tag_orbital, D, J, g, pbc::Bool, hardcore::Bool,species::Symbol, opp_hop::Bool)
    tag_to_index = Dict(round(Int, tag_fock_basis[i]) => i for i in eachindex(tag_fock_basis))
    Hkin = spzeros(Float64, D, D)

    if species == :a
        offset = 0  # a at positions 1,3,5,...
    else
        offset = 1  # b at positions 2,4,6,...
    end

    for u in 1:D
        A = basis[u, :]
        for i in 1:M
            i_site = 2*(i - 1) + 1 + offset
            # ---------- LEFT HOP ----------
            left_site = i_site - 2
            if pbc && left_site < 1
                left_site = 2*M - (1 - offset)  # a=2M-1, b=2M
            end
            if left_site >= 1
                if A[i_site] >= 1
                    if !hardcore || A[left_site] == 0
                        A_new = copy(A)
                        A_new[i_site] -= 1
                        A_new[left_site] += 1
                        tag = round(Int, dot(A_new, tag_orbital))

                        if haskey(tag_to_index, tag)
                            v = tag_to_index[tag]
                            if hardcore
                                coeff = 1.0
                            else
                                coeff = sqrt(A[i_site] * (A[left_site] + 1))
                            end

                            if species == :a
                                Hkin[u, v] += -J * exp(+g) * coeff
                            else  # species == :b
                                if opp_hop == true #species b oppsotie hop condition
                                    Hkin[u, v] += -J * coeff#hermitian
                                else
                                    Hkin[u, v] += -J * coeff #hermitian
                                end
                            end
                        end
                    end
                end
            end

            # ---------- RIGHT HOP ----------
            right_site = i_site + 2
            if pbc && right_site > 2*M
                right_site = offset + 1  # a=1, b=2
            end
            if right_site <= 2*M
                if A[i_site] >= 1
                    if !hardcore || A[right_site] == 0
                        A_new = copy(A)
                        A_new[i_site] -= 1
                        A_new[right_site] += 1
                        tag = round(Int, dot(A_new, tag_orbital))

                        if haskey(tag_to_index, tag)
                            v = tag_to_index[tag]
                            if hardcore
                                coeff = 1.0
                            else
                                coeff = sqrt(A[i_site] * (A[right_site] + 1))
                            end

                            if species == :a
                                Hkin[u, v] += -J * exp(-g) * coeff
                            else  # species == :b
                                if opp_hop == true #species b oppsotie hop condition
                                    Hkin[u, v] += -J* coeff
                                else
                                    Hkin[u, v] += -J * coeff
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return Hkin
end
#-----------------interspecies energy term----------
function interspecies_hamiltonian(Morb, basis, D, Uab)
    # H_int = spzeros(Float64, D, D)
    H_int = spzeros(Complex{Float64}, D, D)
    for u in 1:D
        A = basis[u, :]
        interspecies_energy = 0.0
        for i in 1:Morb
            n_a = A[2i - 1]  # a-particles at sites 2i-1
            n_b = A[2i]# b-particles at sites 2i
            # interspecies_energy += (-im)*Uab * n_a * n_b
            interspecies_energy += Uab.*n_a.*n_b
        end
        H_int[u, u] = interspecies_energy
    end
    return H_int
end
#<<<<<<<<<<<<<intra-species interaction term<<<<<<<<<<
function onsite_hamiltonian(Morb, basis, D, Vaa)
    # H_onsite = spzeros(Float64, D, D)
    H_onsite = spzeros(Complex{Float64}, D, D)
    for u in 1:D
        A = basis[u, :]
        onsite_energy = 0.0
        for i in 1:Morb
            n_a = A[2i - 1]  # a-particles at sites 2i-1
            n_b = A[2i]  # b-particles at sites 2i
            # onsite_energy += (-im)*(Vaa) * ((n_a *(n_a - 1)) + (n_b*(n_b -1)))
            onsite_energy += (Vaa) * ((n_a *(n_a - 1)) + (n_b*(n_b -1)))

        end
        H_onsite[u, u] = onsite_energy
    end
    return H_onsite
end
function tilt_hamiltonian_two_component(Morb, basis, D, delta_a, delta_b)
    H_tilt = spzeros(Float64, D, D)
    for u in 1:D
        A = basis[u, :]
        tilt_energy = 0.0
        for i in 1:Morb
            n_a = A[2*i - 1]
            n_b = A[2*i]
            tilt_energy += i * (delta_a * n_a + delta_b * n_b)
            # tilt_energy += 0.5 * (i/Morb)^2 * (n_a) #only 'a' subjected to this potential
        end
        H_tilt[u, u] = tilt_energy
    end
    return H_tilt
end
##>>>>>>>>>------full Hamiltonian-------------
function build_full_hamiltonian(Morb, basis, tag_fock_basis, tag_orbital, D_tc, Uab, Vaa, J, g, pbc, hardcore,
    delta_a, delta_b, opp_hop)
    Hkin_a = kinetic_component(Morb, basis, tag_fock_basis,tag_orbital, D_tc, J, g, pbc, hardcore, :a, opp_hop)
    Hkin_b = kinetic_component(Morb, basis, tag_fock_basis,tag_orbital, D_tc, J, g, pbc, hardcore, :b, opp_hop)
    H_interspecies = interspecies_hamiltonian(Morb, basis, D_tc, Uab)  #Onsite interaction term
    H_tilt_two_comp = tilt_hamiltonian_two_component(Morb, basis, D_tc, delta_a, delta_b)
    H_onsite = onsite_hamiltonian(Morb, basis, D_tc, Vaa)
    H_kinetic = Hkin_a + Hkin_b;
    if hardcore
        H = H_kinetic .+ H_interspecies .+ H_tilt_two_comp;  #Full Hamiltonian 
    else
        H = H_kinetic .+ H_interspecies .+ H_tilt_two_comp .+ H_onsite;  #Full Hamiltonian 
    end
    return H;
end
# ##>>>>>>>>>>>>>>>>>>>>>density calculation<<<<<<<<<<<<<<<<<<<<<<<<<<<
function compute_site_densities(Morb, basis, tlist, psi_t)
    Nt = length(tlist)
    n_a = zeros(Float64, Nt, Morb)
    n_b = zeros(Float64, Nt, Morb)
    for t in 1:Nt
        psi= psi_t[:, t]
        norm_sq = real(dot(psi, psi))
        probs   = abs2.(psi) ./ norm_sq
        for i in 1:Morb
            n_a[t, i] = sum(probs .* basis[:, 2i-1])
            n_b[t, i] = sum(probs .* basis[:, 2i])
        end
    end
    return n_a, n_b
end

function calculate_densities(Morb, basis, D, psi_gs)
    density_a = zeros(Float64, Morb)
    density_b = zeros(Float64, Morb)
    psi = vec(psi_gs)
    norm_sq = real(dot(psi, psi))
    for u in 1:D
        prob = abs2(psi[u]) / norm_sq
        for i in 1:Morb
            
            density_a[i] += prob * basis[u, 2i-1]
            density_b[i] += prob * basis[u, 2i]
        end
    end
    return density_a, density_b
end
#>>>>>>>>>>>>>>>>>>>>>>>n_tot and sz<<<<<<<<<<<<<<<<<<<
function plot_n_tot_sz_bar!(Morb,dens_a,dens_b)
    sites = 1:Morb
    n_tot = dens_a .+ dens_b;
    sz = (dens_a .- dens_b)./2; 
    p = plot(layout=(2,2), size=(500,500), sharex=true) 
    bar!(p[1], sites, n_tot, label=L"N_{tot}",color= :cyan, alpha=1, 
    linewidth=2,ylims=(0,(maximum(n_tot) +1))) #ylabel=L"\langle n_{\sigma,i} \rangle|_{\sigma=a,b}")
    plot!(p[1], ylabel=L"N_{tot}= N_a + N_b", legend=:best)
    bar!(p[2], sites, sz, alpha=1, linewidth=2,color=:orange,label=L"S_z")
    plot!(p[2], ylabel=L"2~S_z = N_a - N_b"
    ,legend=true, ylims=(minimum(sz) -0.2,maximum(sz) +0.2))
    bar!(p[3], sites, dens_a, ylabel=L"N_a",alpha=1, linewidth=2,color=:red,label=L"N_a",
    ylims=(0,(maximum(n_tot) +1)),xlabel=L"Sites")
    bar!(p[4], sites, dens_b, alpha=1, linewidth=2,color=:green,label=L"N_b")
    plot!(p[4], ylabel=L"N_b", ylims=(0,(maximum(n_tot) +1)),
    legend=:best,xlabel=L"Sites")
    title_str = latexstring("\\Delta_a = $delta_a, U_{ab} = $Uab, g = $g")
    plot!(p, plot_title=title_str, xticks=1:2:Morb)
end

##>>>>>>>>>>>>>>>>>>entropy <<<<<<<<<<<<<<<<<<<<
function entanglement_entropy_half_chain(psi_t, M, basis, D, tlist) 
    cut_position = M÷2
    left_orbitals = 1:(2*cut_position)
    right_orbitals = (2*cut_position + 1):(2*M)
    
    left_basis = [basis[i, left_orbitals] for i in 1:D]
    right_basis = [basis[i, right_orbitals] for i in 1:D]
    
    left_tags = Dict{Vector{Int}, Int}()
    right_tags = Dict{Vector{Int}, Int}()
    
    for i in 1:D
        if !haskey(left_tags, left_basis[i])
            left_tags[left_basis[i]] = length(left_tags) + 1
        end
        if !haskey(right_tags, right_basis[i])
            right_tags[right_basis[i]] = length(right_tags) + 1
        end
    end
    
    dim_left = length(left_tags)
    dim_right = length(right_tags)
    S_half_chain = zeros(length(tlist))
    
    left_indx = [left_tags[left_basis[i]] for i in 1:D]
    right_indx = [right_tags[right_basis[i]] for i in 1:D]
    
    for t in 1:length(tlist)
        psi = psi_t[:, t]
        psi_matrix = zeros(ComplexF64, dim_left, dim_right)
        
        for i in 1:D
            i_left = left_indx[i]
            i_right = right_indx[i]
            psi_matrix[i_left, i_right] += psi[i]
        end
        
        rho_left = psi_matrix * psi_matrix'
        
        e = eigen(rho_left).values
        S = 0.0
        for p in real(e)
            if p > 1e-10
                S -= p * log(p)
            end
        end
        S_half_chain[t] = S
    end
    return S_half_chain
end
###>>>>>>>>>>>>>>>>>>>>>>DYNAMICS<<<<<<<<<<<<<<<<<<<

function calculate_expectation_values_from_matrix(Morb, b_tc, D_tc, psi_t, tlist)
    num_time_points = length(tlist)
    n_tot_time = zeros(Morb, num_time_points)
    sz_time = zeros(Morb, num_time_points)
    density_a_time = zeros(Morb, num_time_points)
    density_b_time = zeros(Morb, num_time_points)
    
    for j in 1:num_time_points
        psi_tj = psi_t[:, j]
        norm_sq = real(dot(psi_tj, psi_tj))
        
        density_a_t = zeros(Morb)
        density_b_t = zeros(Morb)
        
        for u in 1:D_tc
            A = b_tc[u, :]
            prob = abs2(psi_tj[u]) / norm_sq

            for i in 1:Morb
                density_a_t[i] += prob * A[2i-1]
                density_b_t[i] += prob * A[2i]
            end
        end
        
        n_tot_time[:, j] = density_a_t .+ density_b_t
        sz_time[:, j] = (density_a_t .- density_b_t) ./ 2
        density_a_time[:, j] = density_a_t
        density_b_time[:, j] = density_b_t
    end
    return n_tot_time, sz_time, density_a_time, density_b_time
end

#>>>>>>>>>>>>>>>>>>>Biorthogonal eigenvectors>>><<<<<<<<<<<<<<<<

function calculate_biorthogonal_states(H)
    D_tc = size(H,2)
    H_L = H'
    E_L, V_L = eigen(Matrix(H_L))
    E_R, V_R = eigen(Matrix(H))

    idx_R = sortperm(real.(E_R))
    E_R_sorted = E_R[idx_R]
    V_R_sorted = V_R[:,idx_R]

    idx_L = sortperm(real.(E_L))
    E_L_sorted = E_L[idx_L]
    V_L_sorted = V_L[:, idx_L]
    for n in 1:D_tc
        overlap = dot(V_L_sorted[:, n], V_R_sorted[:, n])
        V_R_sorted[:, n] = V_R_sorted[:, n] / sqrt(overlap)
        V_L_sorted[:, n] = V_L_sorted[:, n] / conj(sqrt(overlap)) #important to normalise
        # V_L_sorted[:, n] ./= overlap
    end
    return E_L_sorted, E_R_sorted, V_L_sorted, V_R_sorted
end

function time_evolve(D, E, V_R, V_L, tlist, psi0)
    Nt = length(tlist)
    psi_t_list = Vector{Vector{ComplexF64}}(undef, Nt)
    
    norms = [dot(V_L[:, n], V_R[:, n]) for n in 1:D]
    overlaps = [dot(V_L[:, n], psi0) / norms[n] for n in 1:D]
    
    for ii in 1:Nt
        t = tlist[ii]
        psi = zeros(ComplexF64, D)
        for n in 1:D
            psi .+= overlaps[n] * exp(-im * E[n] * t) * V_R[:, n]
        end
        norm_t = sqrt(real(dot(psi, psi)))
        psi_t_list[ii] = psi./ norm_t
        # psi_t_list[ii] = psi
    end
    return hcat(psi_t_list...)

end

#>>>>>>>>>>>imbalance calculation<<<<<<<<<<<<<<<
function create_density_wave_initial_state(Morb, basis, D)
    state = zeros(Int, 2*Morb)
    psi0 = zeros(ComplexF64, D)
    # center = Morb ÷ 2 
    # center = Morb ÷ 2 + 1
    center_candidates = [Morb÷2, Morb÷2+1, Morb÷2-1]
    center = 0
    for c in center_candidates
        if 1 <= c <= Morb && mod(c, 2) == 0  # even site = bath empty
            center = c
            break
        end
    end  
    state[2*center] = 1
    flag = false
    for i in 1:Morb
        if mod(i, 2) == 1  
            state[2i-1] = 1
        end
    end
    for u in 1:D
        if basis[u, :] == state
            psi0[u] = 1.0
            flag = true
            println("Initial state found at basis index u = ", u)
            break
        end
    end
    if !flag
        println("state not found, check M, Na and Nb, returning zeros")
    end
    return psi0, center
end
##
function calculate_imbalance_a(Morb, basis, D, psi_t, tlist)
    # odd sites = initially occupied, even sites = initially empty
    Nt = length(tlist)
    I_t = zeros(Nt)
    
    for j in 1:Nt
        psi = psi_t[:, j]
        norm_sq = real(dot(psi, psi))
        No = 0.0  # occupation on odd sites (initially occupied)
        Ne = 0.0  # occupation on even sites (initially empty)
        
        for u in 1:D
            A = basis[u, :]
            prob = abs2(psi[u]) / norm_sq
            for i in 1:Morb
                na = A[2i-1]  # a-species only
                if mod(i, 2) == 1  # odd site 
                    No += prob * na
                else               # even site 
                    Ne += prob * na
                end
            end
        end
        I_t[j] = (No - Ne) / (No + Ne)
    end
    return I_t
end
function plot_density(i,den_bos, den_imp)
    n_tot = den_bos + den_imp
    Plots.plot(i, den_bos, marker=:circle, markercolor=:green, lw=3.0, ls=:solid, 
    lc =:green, label=L"a", legend=true, ylims=(0,(maximum(n_tot) +1)))
    Plots.plot!(i, den_imp, marker=:square, markercolor=:red, lw=3.0, ls=:solid, lc =:red,
    ylabel=L"n", xlabel=L"i", size=(450,400),ylims=(0,(maximum(n_tot) +1)),
    label=L"b", legend=true)
end
##>>>>>>>>>>>>>OLD stuff>>>>>>>>>>>>>>>


