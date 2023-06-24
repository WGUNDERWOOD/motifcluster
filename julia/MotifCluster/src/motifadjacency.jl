"""
Build a motif adjacency matrix from an adjacency matrix.
Entry (`i, j`) of a motif adjacency matrix is the
sum of the weights of all motifs containing both
nodes `i` and `j`.
"""
function build_motif_adjacency_matrix(adj_mat::AbstractArray{<:Real}, motif_name::String;
        motif_type::String = "struc", mam_weight_type::String = "unweighted")
    adj_mat_sparse = sparse(adj_mat)
    if motif_name == "Ms"
        return mam_Ms(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "Md"
        return mam_Md(adj_mat_sparse, mam_weight_type)
    elseif motif_name == "M1"
        return mam_M1(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M2"
        return mam_M2(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M3"
        return mam_M3(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M4"
        return mam_M4(adj_mat_sparse, mam_weight_type)
    elseif motif_name == "M5"
        return mam_M5(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M6"
        return mam_M6(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M7"
        return mam_M7(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M8"
        return mam_M8(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M9"
        return mam_M9(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M10"
        return mam_M10(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M11"
        return mam_M11(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M12"
        return mam_M12(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "M13"
        return mam_M13(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "Mcoll"
        return mam_Mcoll(adj_mat_sparse, motif_type, mam_weight_type)
    elseif motif_name == "Mexpa"
        return mam_Mexpa(adj_mat_sparse, motif_type, mam_weight_type)
    end
end

function mam_Ms(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            motif_adj_mat = J .+ J'
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            motif_adj_mat = Js .+ Js'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            G = build_G(adj_mat)
            motif_adj_mat = G .+ G'
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            motif_adj_mat = Gs .+ Gs'
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            motif_adj_mat = G .+ G'
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            motif_adj_mat = Gs .+ Gs'
        end
    end
    return motif_adj_mat
end

function mam_Md(adj_mat, mam_weight_type)
    if mam_weight_type == "unweighted"
        Jd = build_Jd(adj_mat)
        motif_adj_mat = Jd
    elseif mam_weight_type == "mean"
        Gd = build_Gd(adj_mat)
        motif_adj_mat = Gd ./ 2
    elseif mam_weight_type == "product"
        Gp = build_Gp(adj_mat)
        motif_adj_mat = Gp
    end
    return motif_adj_mat
end

function mam_M1(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            C = J' .* (J * J)
            motif_adj_mat = C .+ C'
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            C = Js' .* (Js * Js)
            motif_adj_mat = C .+ C'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            C = J' .* (J * G) .+ J' .* (G * J)
            C += G' .* (J * J)
            motif_adj_mat = (C + C') / 3
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            C = Js' .* (Js * Gs) + Js' .* (Gs * Js)
            C += Gs' .* (Js * Js)
            motif_adj_mat = (C + C') / 3
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            C = G' .* (G * G)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            C = Gs' .* (Gs * Gs)
            motif_adj_mat = C + C'
        end
    end
    return motif_adj_mat
end

function mam_M2(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            C = J' .* (Jd * J) + J' .* (J * Jd)
            C += Jd .* (J * J)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Jd = build_Jd(adj_mat)
            C = Js' .* (Jd * Js) + Js' .* (Js * Jd)
            C += Jd .* (Js * Js)
            motif_adj_mat = C + C'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            G = build_G(adj_mat)
            C = J' .* (Jd * G) + J' .* (Gd * J)
            C += G' .* (Jd * J)
            C = C + J' .* (J * Gd) + J' .* (G * Jd)
            C += G' .* (J * Jd)
            C = C + Jd .* (J * G) + Jd .* (G * J) + Gd .* (J * J)
            motif_adj_mat = (C + C') / 4
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Jd = build_Jd(adj_mat)
            Gs = build_Gs(adj_mat)
            Gd = build_Gd(adj_mat)
            C = Js' .* (Jd * Gs) + Js' .* (Gd * Js)
            C += Gs' .* (Jd * Js)
            C = C + Js' .* (Js * Gd)
            C += Js' .* (Gs * Jd) + Gs' .* (Js * Jd)
            C = C + Jd .* (Js * Gs) + Jd .* (Gs * Js) + Gd .* (Js * Js)
            motif_adj_mat = (C + C') / 4
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Gp = build_Gp(adj_mat)
            C = G' .* (Gp * G) + G' .* (G * Gp)
            C += Gp .* (G * G)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Gp = build_Gp(adj_mat)
            C = Gs' .* (Gp * Gs) + Gs' .* (Gs * Gp)
            C += Gp .* (Gs * Gs)
            motif_adj_mat = C + C'
        end
    end
    return motif_adj_mat
end

function mam_M3(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            C = J .* (Jd * Jd) + Jd .* (Jd * J) + Jd .* (J * Jd)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Jd = build_Jd(adj_mat)
            C = Js .* (Jd * Jd) + Jd .* (Jd * Js) + Jd .* (Js * Jd)
            motif_adj_mat = C + C'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            G = build_G(adj_mat)
            C = J .* (Jd * Gd) + J .* (Gd * Jd) + G .* (Jd * Jd)
            C = C + Jd .* (Jd * G) + Jd .* (Gd * J) + Gd .* (Jd * J)
            C = C + Jd .* (J * Gd) + Jd .* (G * Jd) + Gd .* (J * Jd)
            motif_adj_mat = (C + C') / 5
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Jd = build_Jd(adj_mat)
            Gs = build_Gs(adj_mat)
            Gd = build_Gd(adj_mat)
            C = Js .* (Jd * Gd) + Js .* (Gd * Jd) + Gs .* (Jd * Jd)
            C = C + Jd .* (Jd * Gs) + Jd .* (Gd * Js) + Gd .* (Jd * Js)
            C = C + Jd .* (Js * Gd) + Jd .* (Gs * Jd) + Gd .* (Js * Jd)
            motif_adj_mat = (C + C') / 5
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Gp = build_Gp(adj_mat)
            C = G .* (Gp * Gp) + Gp .* (Gp * G) + Gp .* (G * Gp)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Gp = build_Gp(adj_mat)
            C = Gs .* (Gp * Gp) + Gp .* (Gp * Gs) + Gp .* (Gs * Gp)
            motif_adj_mat = C + C'
        end
    end
    return motif_adj_mat
end

function mam_M4(adj_mat, mam_weight_type)
    if mam_weight_type == "unweighted"
        Jd = build_Jd(adj_mat)
        motif_adj_mat = Jd .* (Jd * Jd)
    elseif mam_weight_type == "mean"
        Jd = build_Jd(adj_mat)
        Gd = build_Gd(adj_mat)
        motif_adj_mat = (Jd .* (Jd * Gd) + Jd .* (Gd * Jd) + Gd .* (Jd * Jd)) / 6
    elseif mam_weight_type == "product"
        Gp = build_Gp(adj_mat)
        motif_adj_mat = Gp .* (Gp * Gp)
    end
    return motif_adj_mat
end

function mam_M5(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            C = J .* (J * J) + J .* (J * J') + J .* (J' * J)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            C = Js .* (Js * Js) + Js .* (Js * Js')
            C += Js .* (Js' * Js)
            motif_adj_mat = C + C'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            C = J .* (J * G) + J .* (G * J) + G .* (J * J)
            C = C + J .* (J * G') + J .* (G * J')
            C += G .* (J * J')
            C = C + J .* (J' * G) + J .* (G' * J)
            C += G .* (J' * J)
            motif_adj_mat = (C + C') / 3
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            C = Js .* (Js * Gs) + Js .* (Gs * Js) + Gs .* (Js * Js)
            C = C + Js .* (Js * Gs')
            C += Js .* (Gs * Js') + Gs .* (Js * Js')
            C = C + Js .* (Js' * Gs)
            C += Js .* (Gs' * Js) + Gs .* (Js' * Js)
            motif_adj_mat = (C + C') / 3
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            C = G .* (G * G) + G .* (G * G')
            C += G .* (G' * G)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            C = Gs .* (Gs * Gs) + Gs .* (Gs * Gs')
            C += Gs .* (Gs' * Gs)
            motif_adj_mat = C + C'
        end
    end
    return motif_adj_mat
end

function mam_M6(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            C = J .* (J * Jd)
            Cprime = Jd .* (J' * J)
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Jd = build_Jd(adj_mat)
            C = Js .* (Js * Jd)
            Cprime = Jd .* (Js' * Js)
            motif_adj_mat = C + C' + Cprime
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            G = build_G(adj_mat)
            C = J .* (J * Gd) + J .* (G * Jd) + G .* (J * Jd)
            Cprime = Jd .* (J' * G) + Jd .* (G' * J)
            Cprime += Gd .* (J' * J)
            motif_adj_mat = (C + C' + Cprime) / 4
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            C = Js .* (Js * Gd) + Js .* (Gs * Jd) + Gs .* (Js * Jd)
            Cprime = Jd .* (Js' * Gs)
            Cprime += Jd .* (Gs' * Js) + Gd .* (Js' * Js)
            motif_adj_mat = (C + C' + Cprime) / 4
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Gp = build_Gp(adj_mat)
            C = G .* (G * Gp)
            Cprime = Gp .* (G' * G)
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Gp = build_Gp(adj_mat)
            C = Gs .* (Gs * Gp)
            Cprime = Gp .* (Gs' * Gs)
            motif_adj_mat = C + C' + Cprime
        end
    end
    return motif_adj_mat
end

function mam_M7(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            C = J .* (Jd * J)
            Cprime = Jd .* (J * J')
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Jd = build_Jd(adj_mat)
            C = Js .* (Jd * Js)
            Cprime = Jd .* (Js * Js')
            motif_adj_mat = C + C' + Cprime
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            G = build_G(adj_mat)
            C = J .* (Jd * G) + J .* (Gd * J) + G .* (Jd * J)
            Cprime = Jd .* (J * G') + Jd .* (G * J')
            Cprime += Gd .* (J * J')
            motif_adj_mat = (C + C' + Cprime) / 4
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            C = Js .* (Jd * Gs) + Js .* (Gd * Js) + Gs .* (Jd * Js)
            Cprime = Jd .* (Js * Gs') + Jd .* (Gs * Js')
            Cprime += Gd .* (Js * Js')
            motif_adj_mat = (C + C' + Cprime) / 4
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Gp = build_Gp(adj_mat)
            C = G .* (Gp * G)
            Cprime = Gp .* (G * G')
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Gp = build_Gp(adj_mat)
            C = Gs .* (Gp * Gs)
            Cprime = Gp .* (Gs * Gs')
            motif_adj_mat = C + C' + Cprime
        end
    end
    return motif_adj_mat
end

function mam_M8(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Id = build_Id(adj_mat)
            C = a_b_one(J, J) - J .* (J)
            Cprime = J' * J - Id .* (J' * J)
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Je = build_Je(adj_mat)
            C = a_b_one(Js, Js) - Js .* (Js * Je)
            Cprime = Js' * Js - Je .* (Js' * Js)
            motif_adj_mat = C + C' + Cprime
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = a_b_one(J, G) - J .* (G) + a_b_one(G, J) - G .* (J)
            Cprime = J' * G - Id .* (J' * G)
            Cprime = Cprime + G' * J - Id .* (G' * J)
            motif_adj_mat = (C + C' + Cprime) / 2
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = a_b_one(Js, Gs) - Js .* (Gs * Je)
            C = C + a_b_one(Gs, Js) - Gs .* (Js * Je)
            Cprime = Js' * Gs - Je .* (Js' * Gs)
            Cprime = Cprime + Gs' * Js - Je .* (Gs' * Js)
            motif_adj_mat = (C + C' + Cprime) / 2
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = a_b_one(G, G) - G .* (G)
            Cprime = G' * G - Id .* (G' * G)
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = a_b_one(Gs, Gs) - Gs .* (Gs * Je)
            Cprime = Gs' * Gs - Je .* (Gs' * Gs)
            motif_adj_mat = C + C' + Cprime
        end
    end
    return motif_adj_mat
end

function mam_M9(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Id = build_Id(adj_mat)
            C = a_one_b(J, J') - 2 * J .* (J') + J * J
            C = C - Id .* (J * J) + a_b_one(J, J')
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Je = build_Je(adj_mat)
            C = a_one_b(Js, Js') - Js .* (Je * Js')
            C = C + Js * Js - Je .* (Js * Js)
            C = C + a_b_one(Js, Js') - Js .* (Js' * Je)
            motif_adj_mat = C + C'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = a_one_b(J, G') - 2 * J .* (G') + J * G
            C = C + a_one_b(G, J') - 2 * G .* (J') + G * J
            C = C - Id .* (J * G) + a_b_one(J, G')
            C = C - Id .* (G * J) + a_b_one(G, J')
            motif_adj_mat = (C + C') / 2
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = a_one_b(Js, Gs') - Js .* (Je * Gs')
            C = C + a_one_b(Gs, Js') - Gs .* (Je * Js')
            C = C + Js * Gs - Je .* (Js * Gs)
            C = C + a_b_one(Js, Gs') - Js .* (Gs' * Je)
            C = C + Gs * Js - Je .* (Gs * Js)
            C = C + a_b_one(Gs, Js') - Gs .* (Js' * Je)
            motif_adj_mat = (C + C') / 2
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = a_one_b(G, G') - 2 * G .* (G') + G * G
            C = C - Id .* (G * G) + a_b_one(G, G')
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = a_one_b(Gs, Gs') - Gs .* (Je * Gs')
            C = C + Gs * Gs - Je .* (Gs * Gs)
            C = C + a_b_one(Gs, Gs') - Gs .* (Gs' * Je)
            motif_adj_mat = C + C'
        end
    end
    return motif_adj_mat
end

function mam_M10(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Id = build_Id(adj_mat)
            C = a_one_b(J, J) - J .* (J)
            Cprime = J * J' - Id .* (J * J')
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Je = build_Je(adj_mat)
            C = a_one_b(Js, Js) - Js .* (Je * Js)
            Cprime = Js * Js' - Je .* (Js * Js')
            motif_adj_mat = C + C' + Cprime
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = a_one_b(J, G) - J .* (G) + a_one_b(G, J) - G .* (J)
            Cprime = J * G' - Id .* (J * G')
            Cprime = Cprime + G * J' - Id .* (G * J')
            motif_adj_mat = (C + C' + Cprime) / 2
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = a_one_b(Js, Gs) - Js .* (Je * Gs)
            C = C + a_one_b(Gs, Js) - Gs .* (Je * Js)
            Cprime = Js * Gs' - Je .* (Js * Gs')
            Cprime = Cprime + Gs * Js' - Je .* (Gs * Js')
            motif_adj_mat = (C + C' + Cprime) / 2
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = a_one_b(G, G) - G .* (G)
            Cprime = G * G' - Id .* (G * G')
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = a_one_b(Gs, Gs) - Gs .* (Je * Gs)
            Cprime = Gs * Gs' - Je .* (Gs * Gs')
            motif_adj_mat = C + C' + Cprime
        end
    end
    return motif_adj_mat
end

function mam_M11(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            Jd = build_Jd(adj_mat)
            Id = build_Id(adj_mat)
            J = build_J(adj_mat)
            C = a_b_one(Jd, J) - Jd .* (J)
            C = C + Jd * J - Id .* (Jd * J)
            C = C + a_b_one(J, Jd) - J .* (Jd)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Jd = build_Jd(adj_mat)
            Je = build_Je(adj_mat)
            Js = build_Js(adj_mat)
            C = a_b_one(Jd, Js) - Jd .* (Js * Je)
            C = C + Jd * Js - Je .* (Jd * Js)
            C = C + a_b_one(Js, Jd) - Js .* (Jd * Je)
            motif_adj_mat = C + C'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            Id = build_Id(adj_mat)
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            C = a_b_one(Jd, G) - Jd .* (G) + a_b_one(Gd, J) - Gd .* (J)
            C = C + Jd * G - Id .* (Jd * G) + Gd * J - Id .* (Gd * J)
            C = C + a_b_one(J, Gd) - J .* (Gd) + a_b_one(G, Jd) - G .* (Jd)
            motif_adj_mat = (C + C') / 3
        elseif motif_type == "struc"
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            Je = build_Je(adj_mat)
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            C = a_b_one(Jd, Gs) - Jd .* (Gs * Je)
            C = C + a_b_one(Gd, Js) - Gd .* (Js * Je)
            C = C + Jd * Gs - Je .* (Jd * Gs) + Gd * Js - Je .* (Gd * Js)
            C = C + a_b_one(Js, Gd) - Js .* (Gd * Je)
            C = C + a_b_one(Gs, Jd) - Gs .* (Jd * Je)
            motif_adj_mat = (C + C') / 3
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            Gp = build_Gp(adj_mat)
            Id = build_Id(adj_mat)
            G = build_G(adj_mat)
            C = a_b_one(Gp, G) - Gp .* (G)
            C = C + Gp * G - Id .* (Gp * G)
            C = C + a_b_one(G, Gp) - G .* (Gp)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Gp = build_Gp(adj_mat)
            Je = build_Je(adj_mat)
            Gs = build_Gs(adj_mat)
            C = a_b_one(Gp, Gs) - Gp .* (Gs * Je)
            C = C + Gp * Gs - Je .* (Gp * Gs)
            C = C + a_b_one(Gs, Gp) - Gs .* (Gp * Je)
            motif_adj_mat = C + C'
        end
    end
    return motif_adj_mat
end

function mam_M12(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            Jd = build_Jd(adj_mat)
            Id = build_Id(adj_mat)
            J = build_J(adj_mat)
            C = a_one_b(Jd, J) - Jd .* (J)
            C = C + J * Jd - Id .* (J * Jd)
            C = C + a_one_b(J, Jd) - J .* (Jd)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Jd = build_Jd(adj_mat)
            Je = build_Je(adj_mat)
            Js = build_Js(adj_mat)
            C = a_one_b(Jd, Js) - Jd .* (Je * Js)
            C = C + Js * Jd - Je .* (Js * Jd)
            C = C + a_one_b(Js, Jd) - Js .* (Je * Jd)
            motif_adj_mat = C + C'
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            Id = build_Id(adj_mat)
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            C = a_one_b(Jd, G) - Jd .* (G) + a_one_b(Gd, J) - Gd .* (J)
            C = C + J * Gd - Id .* (J * Gd) + G * Jd - Id .* (G * Jd)
            C = C + a_one_b(J, Gd) - J .* (Gd) + a_one_b(G, Jd) - G .* (Jd)
            motif_adj_mat = (C + C') / 3
        elseif motif_type == "struc"
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            Je = build_Je(adj_mat)
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            C = a_one_b(Jd, Gs) - Jd .* (Je * Gs)
            C = C + a_one_b(Gd, Js) - Gd .* (Je * Js)
            C = C + Js * Gd - Je .* (Js * Gd) + Gs * Jd - Je .* (Gs * Jd)
            C = C + a_one_b(Js, Gd) - Js .* (Je * Gd)
            C = C + a_one_b(Gs, Jd) - Gs .* (Je * Jd)
            motif_adj_mat = (C + C') / 3
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            Gp = build_Gp(adj_mat)
            Id = build_Id(adj_mat)
            G = build_G(adj_mat)
            C = a_one_b(Gp, G) - Gp .* (G)
            C = C + G * Gp - Id .* (G * Gp)
            C = C + a_one_b(G, Gp) - G .* (Gp)
            motif_adj_mat = C + C'
        elseif motif_type == "struc"
            Gp = build_Gp(adj_mat)
            Je = build_Je(adj_mat)
            Gs = build_Gs(adj_mat)
            C = a_one_b(Gp, Gs) - Gp .* (Je * Gs)
            C = C + Gs * Gp - Je .* (Gs * Gp)
            C = C + a_one_b(Gs, Gp) - Gs .* (Je * Gp)
            motif_adj_mat = C + C'
        end
    end
    return motif_adj_mat
end

function mam_M13(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            Jd = build_Jd(adj_mat)
            Id = build_Id(adj_mat)
            C = a_b_one(Jd, Jd) - Jd .* (Jd)
            Cprime = Jd * Jd - Id .* (Jd * Jd)
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Jd = build_Jd(adj_mat)
            Je = build_Je(adj_mat)
            C = a_b_one(Jd, Jd) - Jd .* (Jd * Je)
            Cprime = Jd * Jd - Je .* (Jd * Jd)
            motif_adj_mat = C + C' + Cprime
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            Id = build_Id(adj_mat)
            C = a_b_one(Jd, Gd) - Jd .* (Gd) + a_b_one(Gd, Jd) - Gd .* (Jd)
            C = C + Jd * Gd - Id .* (Jd * Gd)
            motif_adj_mat = (C + C') / 4
        elseif motif_type == "struc"
            Jd = build_Jd(adj_mat)
            Gd = build_Gd(adj_mat)
            Je = build_Je(adj_mat)
            C = a_b_one(Jd, Gd) - Jd .* (Gd * Je)
            C = C + a_b_one(Gd, Jd) - Gd .* (Jd * Je)
            C = C + Jd * Gd - Je .* (Jd * Gd)
            motif_adj_mat = (C + C') / 4
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            Gp = build_Gp(adj_mat)
            Id = build_Id(adj_mat)
            C = a_b_one(Gp, Gp) - Gp .* (Gp)
            Cprime = Gp * Gp - Id .* (Gp * Gp)
            motif_adj_mat = C + C' + Cprime
        elseif motif_type == "struc"
            Gp = build_Gp(adj_mat)
            Je = build_Je(adj_mat)
            C = a_b_one(Gp, Gp) - Gp .* (Gp * Je)
            Cprime = Gp * Gp - Je .* (Gp * Gp)
            motif_adj_mat = C + C' + Cprime
        end
    end
    return motif_adj_mat
end

function mam_Mcoll(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Id = build_Id(adj_mat)
            C = J * J' - Id .* (J * J')
            motif_adj_mat = C
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Je = build_Je(adj_mat)
            C = Js * Js' - Je .* (Js * Js')
            motif_adj_mat = C
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = J * G' - Id .* (J * G')
            C = C + G * J' - Id .* (G * J')
            motif_adj_mat = C / 2
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = Js * Gs' - Je .* (Js * Gs')
            C = C + Gs * Js' - Je .* (Gs * Js')
            motif_adj_mat = C / 2
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = G * G' - Id .* (G * G')
            motif_adj_mat = C
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = Gs * Gs' - Je .* (Gs * Gs')
            motif_adj_mat = C
        end
    end
    return motif_adj_mat
end

function mam_Mexpa(adj_mat, motif_type, mam_weight_type)
    if mam_weight_type == "unweighted"
        if motif_type == "func"
            J = build_J(adj_mat)
            Id = build_Id(adj_mat)
            C = J' * J - Id .* (J' * J)
            motif_adj_mat = C
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Je = build_Je(adj_mat)
            C = Js' * Js - Je .* (Js' * Js)
            motif_adj_mat = C
        end
    elseif mam_weight_type == "mean"
        if motif_type == "func"
            J = build_J(adj_mat)
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = J' * G - Id .* (J' * G)
            C = C + G' * J - Id .* (G' * J)
            motif_adj_mat = C / 2
        elseif motif_type == "struc"
            Js = build_Js(adj_mat)
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = Js' * Gs - Je .* (Js' * Gs)
            C = C + Gs' * Js - Je .* (Gs' * Js)
            motif_adj_mat = C / 2
        end
    elseif mam_weight_type == "product"
        if motif_type == "func"
            G = build_G(adj_mat)
            Id = build_Id(adj_mat)
            C = G' * G - Id .* (G' * G)
            motif_adj_mat = C
        elseif motif_type == "struc"
            Gs = build_Gs(adj_mat)
            Je = build_Je(adj_mat)
            C = Gs' * Gs - Je .* (Gs' * Gs)
            motif_adj_mat = C
        end
    end
    return motif_adj_mat
end
