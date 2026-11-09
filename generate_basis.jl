###<<<<<<<,generate orbitals unique numbers>>>>>>>>>>>>>>
function get_tag_orbital_two_component(Morb)
    tag_precision = 1e12
    tag_orbital = round.(sqrt.((1:2*Morb) .* 100 .+ 3) .* tag_precision)
    return tag_orbital
end
##>>>>>>>>>>>>>>>>>>generate one comp basis>>>>>>>>>>>>>
function generate_basis(Npar, Morb, hardcore::Bool)
   if hardcore
        D = binomial(Morb, Npar)
    else
        D = binomial(Npar + Morb - 1, Npar)
    end
    basis = zeros(Int, D, Morb)
    basis_temp = zeros(Int, Morb)
    if hardcore
        # Hardcore bosons: occupation is 0 or 1
        basis_temp[1:Npar] .= 1
        for u in 1:D
            basis[u, :] .= basis_temp
            k = findlast(i -> basis_temp[i] == 1 && i < Morb - (Npar - sum(basis_temp[1:i])), 1:Morb)
            if isnothing(k)
                break
            end
            basis_temp[k] = 0
            right = Npar - sum(basis_temp[1:k])
            basis_temp[(k + 1):(k + right)] .= 1
            basis_temp[(k + right + 1):end] .= 0
        end
    else
        # Softcore bosons: can be upto N at one site
        basis_temp[1] = Npar
        for u in 1:D
            basis[u, :] = basis_temp
            k = findlast(x -> x > 0, basis_temp[1:(Morb - 1)])
            if !isnothing(k)
                basis_new = zeros(Int, Morb)
                basis_new[1:(k - 1)] .= basis_temp[1:(k - 1)]
                basis_new[k] = basis_temp[k] - 1
                basis_new[k + 1] = Npar - sum(basis_new[1:k])
                basis_temp .= basis_new
            end
        end
    end
    return basis, D
end
###<<<<<<<<,generate two comp basis>>>>>>>>>>>>>>>>
function generate_basis_two_component(Na, Nb, Morb, hardcore::Bool)
    basis_a, Da = generate_basis(Na, Morb, hardcore)
    basis_b, Db = generate_basis(Nb, Morb, hardcore)
    D_tc = Da * Db
    basis_tc = zeros(Int, D_tc, 2*Morb)
    tag_fock_basis = zeros(Int, D_tc)

    tag_orbital = get_tag_orbital_two_component(Morb)

    count = 1
    for ua in 1:Da
        for ub in 1:Db
            for i in 1:Morb
                basis_tc[count, 2i - 1] = basis_a[ua, i]
                basis_tc[count, 2i]= basis_b[ub, i]
            end
            tag_fock_basis[count] = dot(basis_tc[count, :], tag_orbital)
            count += 1
        end
    end

    return basis_tc, tag_fock_basis, tag_orbital, D_tc
end
