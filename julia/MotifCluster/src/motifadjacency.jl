"""
Build a motif adjacency matrix from an adjacency matrix.
Entry (`i, j`) of a motif adjacency matrix is the
sum of the weights of all motifs containing both
nodes `i` and `j`.
"""
function build_motif_adjacency_matrix(adj_mat::AbstractArray{<:Real}, motif_name::String,
        motif_type::String, mam_weight_type::String, mam_method::String)
    if motif_name == "Ms"
        return mam_Ms(adj_mat, motif_type, mam_weight_type)
    elseif motif_name == "Md"
        return mam_Md(adj_mat, mam_weight_type)
    elseif motif_name == "M1"
        return mam_M1(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M2"
        return mam_M2(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M3"
        return mam_M3(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M4"
        return mam_M4(adj_mat, mam_weight_type, mam_method)
    elseif motif_name == "M5"
        return mam_M5(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M6"
        return mam_M6(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M7"
        return mam_M7(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M8"
        return mam_M8(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M9"
        return mam_M9(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M10"
        return mam_M10(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M11"
        return mam_M11(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M12"
        return mam_M12(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "M13"
        return mam_M13(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "Mcoll"
        return mam_Mcoll(adj_mat, motif_type, mam_weight_type, mam_method)
    elseif motif_name == "Mexpa"
        return mam_Mexpa(adj_mat, motif_type, mam_weight_type, mam_method)
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

function mam_M1(adj_mat, motif_type, mam_weight_type, mam_method)
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

#=

function mam_M2(adj_mat, motif_type, mam_weight_type, mam_method): # pylint: disable=too-many-statements


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J' * (Jd @ J) + J' * (J @ Jd)
C += Jd * (J @ J)
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js' * (Jd @ Js) + Js' * (Js @ Jd)
C += Jd * (Js @ Js)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J' * (Jd @ G) + J' * (Gd @ J)
C += G' * (Jd @ J)
C = C + J' * (J @ Gd) + J' * (G @ Jd)
C += G' * (J @ Jd)
C = C + Jd * (J @ G) + Jd * (G @ J) + Gd * (J @ J)
motif_adj_mat = (C + C') / 4

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
Gs = build_Gs(adj_mat)
Gd = build_Gd(adj_mat)
C = Js' * (Jd @ Gs) + Js' * (Gd @ Js)
C += Gs' * (Jd @ Js)
C = C + Js' * (Js @ Gd)
C += Js' * (Gs @ Jd) + Gs' * (Js @ Jd)
C = C + Jd * (Js @ Gs) + Jd * (Gs @ Js) + Gd * (Js @ Js)
motif_adj_mat = (C + C') / 4

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G' * (Gp @ G) + G' * (G @ Gp)
C += Gp * (G @ G)
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs' * (Gp @ Gs) + Gs' * (Gs @ Gp)
C += Gp * (Gs @ Gs)
motif_adj_mat = C + C'

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J'.multiply(Jd * J) + J'.multiply(J * Jd)
C += Jd.multiply(J * J)
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js'.multiply(Jd * Js) + Js'.multiply(Js * Jd)
C += Jd.multiply(Js * Js)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J'.multiply(Jd * G) + J'.multiply(Gd * J)
C += G'.multiply(Jd * J)
C = C + J'.multiply(J * Gd) + J'.multiply(G * Jd)
C += G'.multiply(J * Jd)
C = C + Jd.multiply(J * G) + Jd.multiply(G * J) + Gd.multiply(J * J)
motif_adj_mat = (C + C') / 4

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
Gs = build_Gs(adj_mat)
Gd = build_Gd(adj_mat)
C = Js'.multiply(Jd * Gs) + Js'.multiply(Gd * Js)
C += Gs'.multiply(Jd * Js)
C = C + Js'.multiply(Js * Gd)
C += Js'.multiply(Gs * Jd) + Gs'.multiply(Js * Jd)
C = C + Jd.multiply(Js * Gs) + Jd.multiply(Gs * Js) + Gd.multiply(Js * Js)
motif_adj_mat = (C + C') / 4

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G'.multiply(Gp * G) + G'.multiply(G * Gp)
C += Gp.multiply(G * G)
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs'.multiply(Gp * Gs) + Gs'.multiply(Gs * Gp)
C += Gp.multiply(Gs * Gs)
motif_adj_mat = C + C'

return motif_adj_mat


function mam_M3(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J * (Jd @ Jd) + Jd * (Jd @ J) + Jd * (J @ Jd)
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js * (Jd @ Jd) + Jd * (Jd @ Js) + Jd * (Js @ Jd)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J * (Jd @ Gd) + J * (Gd @ Jd) + G * (Jd @ Jd)
C = C + Jd * (Jd @ G) + Jd * (Gd @ J) + Gd * (Jd @ J)
C = C + Jd * (J @ Gd) + Jd * (G @ Jd) + Gd * (J @ Jd)
motif_adj_mat = (C + C') / 5

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
Gs = build_Gs(adj_mat)
Gd = build_Gd(adj_mat)
C = Js * (Jd @ Gd) + Js * (Gd @ Jd) + Gs * (Jd @ Jd)
C = C + Jd * (Jd @ Gs) + Jd * (Gd @ Js) + Gd * (Jd @ Js)
C = C + Jd * (Js @ Gd) + Jd * (Gs @ Jd) + Gd * (Js @ Jd)
motif_adj_mat = (C + C') / 5

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G * (Gp @ Gp) + Gp * (Gp @ G) + Gp * (G @ Gp)
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs * (Gp @ Gp) + Gp * (Gp @ Gs) + Gp * (Gs @ Gp)
motif_adj_mat = C + C'

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J.multiply(Jd * Jd) + Jd.multiply(Jd * J) + Jd.multiply(J * Jd)
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js.multiply(Jd * Jd) + Jd.multiply(Jd * Js) + Jd.multiply(Js * Jd)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J.multiply(Jd * Gd) + J.multiply(Gd * Jd) + G.multiply(Jd * Jd)
C = C + Jd.multiply(Jd * G) + Jd.multiply(Gd * J) + Gd.multiply(Jd * J)
C = C + Jd.multiply(J * Gd) + Jd.multiply(G * Jd) + Gd.multiply(J * Jd)
motif_adj_mat = (C + C') / 5

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
Gs = build_Gs(adj_mat)
Gd = build_Gd(adj_mat)
C = Js.multiply(Jd * Gd) + Js.multiply(Gd * Jd) + Gs.multiply(Jd * Jd)
C = C + Jd.multiply(Jd * Gs) + Jd.multiply(Gd * Js) + Gd.multiply(Jd * Js)
C = C + Jd.multiply(Js * Gd) + Jd.multiply(Gs * Jd) + Gd.multiply(Js * Jd)
motif_adj_mat = (C + C') / 5

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G.multiply(Gp * Gp) + Gp.multiply(Gp * G) + Gp.multiply(G * Gp)
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs.multiply(Gp * Gp) + Gp.multiply(Gp * Gs) + Gp.multiply(Gs * Gp)
motif_adj_mat = C + C'

return motif_adj_mat


function mam_M4(adj_mat, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
Jd = build_Jd(adj_mat)
motif_adj_mat = Jd * (Jd @ Jd)

if mam_weight_type == "mean"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
motif_adj_mat = (Jd * (Jd @ Gd) + Jd * (Gd @ Jd) + Gd * (Jd @ Jd)) / 6

if mam_weight_type == "product"
Gp = build_Gp(adj_mat)
motif_adj_mat = Gp * (Gp @ Gp)

if mam_method == "sparse"
if mam_weight_type == "unweighted"
Jd = build_Jd(adj_mat)
motif_adj_mat = Jd.multiply(Jd * Jd)

if mam_weight_type == "mean"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
motif_adj_mat = (Jd.multiply(Jd * Gd) + Jd.multiply(Gd * Jd) + Gd.multiply(Jd * Jd)) / 6

if mam_weight_type == "product"
Gp = build_Gp(adj_mat)
motif_adj_mat = Gp.multiply(Gp * Gp)

return motif_adj_mat


function mam_M5(adj_mat, motif_type, mam_weight_type, mam_method)

if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
C = J * (J @ J) + J * (J @ J') + J * (J' @ J)
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
C = Js * (Js @ Js) + Js * (Js @ Js')
C += Js * (Js' @ Js)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
G = build_G(adj_mat)
C = J * (J @ G) + J * (G @ J) + G * (J @ J)
C = C + J * (J @ G') + J * (G @ J')
C += G * (J @ J')
C = C + J * (J' @ G) + J * (G' @ J)
C += G * (J' @ J)
motif_adj_mat = (C + C') / 3

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
C = Js * (Js @ Gs) + Js * (Gs @ Js) + Gs * (Js @ Js)
C = C + Js * (Js @ Gs')
C += Js * (Gs @ Js') + Gs * (Js @ Js')
C = C + Js * (Js' @ Gs)
C += Js * (Gs' @ Js) + Gs * (Js' @ Js)
motif_adj_mat = (C + C') / 3

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
C = G * (G @ G) + G * (G @ G')
C += G * (G' @ G)
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
C = Gs * (Gs @ Gs) + Gs * (Gs @ Gs')
C += Gs * (Gs' @ Gs)
motif_adj_mat = C + C'

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
C = J.multiply(J * J) + J.multiply(J * J') + J.multiply(J' * J)
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
C = Js.multiply(Js * Js) + Js.multiply(Js * Js')
C += Js.multiply(Js' * Js)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
G = build_G(adj_mat)
C = J.multiply(J * G) + J.multiply(G * J) + G.multiply(J * J)
C = C + J.multiply(J * G') + J.multiply(G * J')
C += G.multiply(J * J')
C = C + J.multiply(J' * G) + J.multiply(G' * J)
C += G.multiply(J' * J)
motif_adj_mat = (C + C') / 3

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
C = Js.multiply(Js * Gs) + Js.multiply(Gs * Js) + Gs.multiply(Js * Js)
C = C + Js.multiply(Js * Gs')
C += Js.multiply(Gs * Js') + Gs.multiply(Js * Js')
C = C + Js.multiply(Js' * Gs)
C += Js.multiply(Gs' * Js) + Gs.multiply(Js' * Js)
motif_adj_mat = (C + C') / 3

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
C = G.multiply(G * G) + G.multiply(G * G')
C += G.multiply(G' * G)
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
C = Gs.multiply(Gs * Gs) + Gs.multiply(Gs * Gs')
C += Gs.multiply(Gs' * Gs)
motif_adj_mat = C + C'

return motif_adj_mat


function mam_M6(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J * (J @ Jd)
Cprime = Jd * (J' @ J)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js * (Js @ Jd)
Cprime = Jd * (Js' @ Js)
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J * (J @ Gd) + J * (G @ Jd) + G * (J @ Jd)
Cprime = Jd * (J' @ G) + Jd * (G' @ J)
Cprime += Gd * (J' @ J)
motif_adj_mat = (C + C' + Cprime) / 4

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
C = Js * (Js @ Gd) + Js * (Gs @ Jd) + Gs * (Js @ Jd)
Cprime = Jd * (Js' @ Gs)
Cprime += Jd * (Gs' @ Js) + Gd * (Js' @ Js)
motif_adj_mat = (C + C' + Cprime) / 4

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G * (G @ Gp)
Cprime = Gp * (G' @ G)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs * (Gs @ Gp)
Cprime = Gp * (Gs' @ Gs)
motif_adj_mat = C + C' + Cprime

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J.multiply(J * Jd)
Cprime = Jd.multiply(J' * J)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js.multiply(Js * Jd)
Cprime = Jd.multiply(Js' * Js)
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J.multiply(J * Gd) + J.multiply(G * Jd) + G.multiply(J * Jd)
Cprime = Jd.multiply(J' * G) + Jd.multiply(G' * J)
Cprime += Gd.multiply(J' * J)
motif_adj_mat = (C + C' + Cprime) / 4

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
C = Js.multiply(Js * Gd) + Js.multiply(Gs * Jd) + Gs.multiply(Js * Jd)
Cprime = Jd.multiply(Js' * Gs)
Cprime += Jd.multiply(Gs' * Js) + Gd.multiply(Js' * Js)
motif_adj_mat = (C + C' + Cprime) / 4

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G.multiply(G * Gp)
Cprime = Gp.multiply(G' * G)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs.multiply(Gs * Gp)
Cprime = Gp.multiply(Gs' * Gs)
motif_adj_mat = C + C' + Cprime

return motif_adj_mat


function mam_M7(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J * (Jd @ J)
Cprime = Jd * (J @ J')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js * (Jd @ Js)
Cprime = Jd * (Js @ Js')
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J * (Jd @ G) + J * (Gd @ J) + G * (Jd @ J)
Cprime = Jd * (J @ G') + Jd * (G @ J')
Cprime += Gd * (J @ J')
motif_adj_mat = (C + C' + Cprime) / 4

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
C = Js * (Jd @ Gs) + Js * (Gd @ Js) + Gs * (Jd @ Js)
Cprime = Jd * (Js @ Gs') + Jd * (Gs @ Js')
Cprime += Gd * (Js @ Js')
motif_adj_mat = (C + C' + Cprime) / 4

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G * (Gp @ G)
Cprime = Gp * (G @ G')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs * (Gp @ Gs)
Cprime = Gp * (Gs @ Gs')
motif_adj_mat = C + C' + Cprime

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
C = J.multiply(Jd * J)
Cprime = Jd.multiply(J * J')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
Jd = build_Jd(adj_mat)
C = Js.multiply(Jd * Js)
Cprime = Jd.multiply(Js * Js')
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
G = build_G(adj_mat)
C = J.multiply(Jd * G) + J.multiply(Gd * J) + G.multiply(Jd * J)
Cprime = Jd.multiply(J * G') + Jd.multiply(G * J')
Cprime += Gd.multiply(J * J')
motif_adj_mat = (C + C' + Cprime) / 4

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
C = Js.multiply(Jd * Gs) + Js.multiply(Gd * Js) + Gs.multiply(Jd * Js)
Cprime = Jd.multiply(Js * Gs') + Jd.multiply(Gs * Js')
Cprime += Gd.multiply(Js * Js')
motif_adj_mat = (C + C' + Cprime) / 4

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Gp = build_Gp(adj_mat)
C = G.multiply(Gp * G)
Cprime = Gp.multiply(G * G')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Gp = build_Gp(adj_mat)
C = Gs.multiply(Gp * Gs)
Cprime = Gp.multiply(Gs * Gs')
motif_adj_mat = C + C' + Cprime

return motif_adj_mat


function mam_M8(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
C = J * (J @ Jn)
Cprime = Jn * (J' @ J)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
J0 = build_J0(adj_mat)
C = Js * (Js @ J0)
Cprime = J0 * (Js' @ Js)
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
G = build_G(adj_mat)
C = J * (G @ Jn) + G * (J @ Jn)
Cprime = Jn * (J' @ G) + Jn * (G' @ J)
motif_adj_mat = (C + C' + Cprime) / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = Js * (Gs @ J0) + Gs * (Js @ J0)
Cprime = J0 * (Js' @ Gs) + J0 * (Gs' @ Js)
motif_adj_mat = (C + C' + Cprime) / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Jn = build_Jn(adj_mat)
C = G * (G @ Jn)
Cprime = Jn * (G' @ G)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = Gs * (Gs @ J0)
Cprime = J0 * (Gs' @ Gs)
motif_adj_mat = C + C' + Cprime

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_b_one(J, J) - J.multiply(J)
Cprime = J' * J - Id.multiply(J' * J)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_b_one(Js, Js) - Js.multiply(Js * Je)
Cprime = Js' * Js - Je.multiply(Js' * Js)
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_b_one(J, G) - J.multiply(G) + mcut._a_b_one(G, J) - G.multiply(J)
Cprime = J' * G - Id.multiply(J' * G)
Cprime = Cprime + G' * J - Id.multiply(G' * J)
motif_adj_mat = (C + C' + Cprime) / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_b_one(Js, Gs) - Js.multiply(Gs * Je)
C = C + mcut._a_b_one(Gs, Js) - Gs.multiply(Js * Je)
Cprime = Js' * Gs - Je.multiply(Js' * Gs)
Cprime = Cprime + Gs' * Js - Je.multiply(Gs' * Js)
motif_adj_mat = (C + C' + Cprime) / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_b_one(G, G) - G.multiply(G)
Cprime = G' * G - Id.multiply(G' * G)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_b_one(Gs, Gs) - Gs.multiply(Gs * Je)
Cprime = Gs' * Gs - Je.multiply(Gs' * Gs)
motif_adj_mat = C + C' + Cprime

return motif_adj_mat


function mam_M9(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
C = J * (Jn @ J') + Jn * (J @ J)
C += J * (J' @ Jn)
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
J0 = build_J0(adj_mat)
C = Js * (J0 @ Js') + J0 * (Js @ Js)
C += Js * (Js' @ J0)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
G = build_G(adj_mat)
C = J * (Jn @ G') + G * (Jn @ J')
C = C + Jn * (J @ G) + Jn * (G @ J)
C = C + J * (G' @ Jn) + G * (J' @ Jn)
motif_adj_mat = (C + C') / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = Js * (J0 @ Gs') + Gs * (J0 @ Js')
C = C + J0 * (Js @ Gs) + J0 * (Gs @ Js)
C = C + Js * (Gs' @ J0) + Gs * (Js' @ J0)
motif_adj_mat = (C + C') / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Jn = build_Jn(adj_mat)
C = G * (Jn @ G') + Jn * (G @ G) + G * (G' @ Jn)
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = Gs * (J0 @ Gs') + J0 * (Gs @ Gs)
C += Gs * (Gs' @ J0)
motif_adj_mat = C + C'

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_one_b(J, J') - 2 * J.multiply(J') + J * J
C = C - Id.multiply(J * J) + mcut._a_b_one(J, J')
motif_adj_mat = C + C'

if motif_type == "struc"
Js = build_Js(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_one_b(Js, Js') - Js.multiply(Je * Js')
C = C + Js * Js - Je.multiply(Js * Js)
C = C + mcut._a_b_one(Js, Js') - Js.multiply(Js' * Je)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_one_b(J, G') - 2 * J.multiply(G') + J * G
C = C + mcut._a_one_b(G, J') - 2 * G.multiply(J') + G * J
C = C - Id.multiply(J * G) + mcut._a_b_one(J, G')
C = C - Id.multiply(G * J) + mcut._a_b_one(G, J')
motif_adj_mat = (C + C') / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_one_b(Js, Gs') - Js.multiply(Je * Gs')
C = C + mcut._a_one_b(Gs, Js') - Gs.multiply(Je * Js')
C = C + Js * Gs - Je.multiply(Js * Gs)
C = C + mcut._a_b_one(Js, Gs') - Js.multiply(Gs' * Je)
C = C + Gs * Js - Je.multiply(Gs * Js)
C = C + mcut._a_b_one(Gs, Js') - Gs.multiply(Js' * Je)
motif_adj_mat = (C + C') / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_one_b(G, G') - 2 * G.multiply(G') + G * G
C = C - Id.multiply(G * G) + mcut._a_b_one(G, G')
motif_adj_mat = C + C'

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_one_b(Gs, Gs') - Gs.multiply(Je * Gs')
C = C + Gs * Gs - Je.multiply(Gs * Gs)
C = C + mcut._a_b_one(Gs, Gs') - Gs.multiply(Gs' * Je)
motif_adj_mat = C + C'

return motif_adj_mat


function mam_M10(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
C = J * (Jn @ J)
Cprime = Jn * (J @ J')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
J0 = build_J0(adj_mat)
C = Js * (J0 @ Js)
Cprime = J0 * (Js @ Js')
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
G = build_G(adj_mat)
C = J * (Jn @ G) + G * (Jn @ J)
Cprime = Jn * (J @ G') + Jn * (G @ J')
motif_adj_mat = (C + C' + Cprime) / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = Js * (J0 @ Gs) + Gs * (J0 @ Js)
Cprime = J0 * (Js @ Gs') + J0 * (Gs @ Js')
motif_adj_mat = (C + C' + Cprime) / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Jn = build_Jn(adj_mat)
C = G * (Jn @ G)
Cprime = Jn * (G @ G')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = Gs * (J0 @ Gs)
Cprime = J0 * (Gs @ Gs')
motif_adj_mat = C + C' + Cprime

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_one_b(J, J) - J.multiply(J)
Cprime = J * J' - Id.multiply(J * J')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Js = build_Js(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_one_b(Js, Js) - Js.multiply(Je * Js)
Cprime = Js * Js' - Je.multiply(Js * Js')
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_one_b(J, G) - J.multiply(G) + mcut._a_one_b(G, J) - G.multiply(J)
Cprime = J * G' - Id.multiply(J * G')
Cprime = Cprime + G * J' - Id.multiply(G * J')
motif_adj_mat = (C + C' + Cprime) / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_one_b(Js, Gs) - Js.multiply(Je * Gs)
C = C + mcut._a_one_b(Gs, Js) - Gs.multiply(Je * Js)
Cprime = Js * Gs' - Je.multiply(Js * Gs')
Cprime = Cprime + Gs * Js' - Je.multiply(Gs * Js')
motif_adj_mat = (C + C' + Cprime) / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_one_b(G, G) - G.multiply(G)
Cprime = G * G' - Id.multiply(G * G')
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_one_b(Gs, Gs) - Gs.multiply(Je * Gs)
Cprime = Gs * Gs' - Je.multiply(Gs * Gs')
motif_adj_mat = C + C' + Cprime

return motif_adj_mat


function mam_M11(adj_mat, motif_type, mam_weight_type, mam_method): # pylint: disable=too-many-statements


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Jn = build_Jn(adj_mat)
J = build_J(adj_mat)
C = Jd * (J @ Jn) + Jn * (Jd @ J) + J * (Jd @ Jn)
motif_adj_mat = C + C'

if motif_type == "struc"
Jd = build_Jd(adj_mat)
J0 = build_J0(adj_mat)
Js = build_Js(adj_mat)
C = Jd * (Js @ J0) + J0 * (Jd @ Js) + Js * (Jd @ J0)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Jn = build_Jn(adj_mat)
J = build_J(adj_mat)
G = build_G(adj_mat)
C = Jd * (G @ Jn) + Gd * (J @ Jn)
C = C + Jn * (Jd @ G) + Jn * (Gd @ J)
C = C + J * (Gd @ Jn) + G * (Jd @ Jn)
motif_adj_mat = (C + C') / 3

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
J0 = build_J0(adj_mat)
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
C = Jd * (Gs @ J0) + Gd * (Js @ J0)
C = C + J0 * (Jd @ Gs) + J0 * (Gd @ Js)
C = C + Js * (Gd @ J0) + Gs * (Jd @ J0)
motif_adj_mat = (C + C') / 3

if mam_weight_type == "product"
if motif_type == "func"
Gp = build_Gp(adj_mat)
Jn = build_Jn(adj_mat)
G = build_G(adj_mat)
C = Gp * (G @ Jn) + Jn * (Gp @ G) + G * (Gp @ Jn)
motif_adj_mat = C + C'

if motif_type == "struc"
Gp = build_Gp(adj_mat)
J0 = build_J0(adj_mat)
Gs = build_Gs(adj_mat)
C = Gp * (Gs @ J0) + J0 * (Gp @ Gs) + Gs * (Gp @ J0)
motif_adj_mat = C + C'

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Id = build_Id(adj_mat)
J = build_J(adj_mat)
C = mcut._a_b_one(Jd, J) - Jd.multiply(J)
C = C + Jd * J - Id.multiply(Jd * J)
C = C + mcut._a_b_one(J, Jd) - J.multiply(Jd)
motif_adj_mat = C + C'

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Je = build_Je(adj_mat)
Js = build_Js(adj_mat)
C = mcut._a_b_one(Jd, Js) - Jd.multiply(Js * Je)
C = C + Jd * Js - Je.multiply(Jd * Js)
C = C + mcut._a_b_one(Js, Jd) - Js.multiply(Jd * Je)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Id = build_Id(adj_mat)
J = build_J(adj_mat)
G = build_G(adj_mat)
C = mcut._a_b_one(Jd, G) - Jd.multiply(G) + mcut._a_b_one(Gd, J) - Gd.multiply(J)
C = C + Jd * G - Id.multiply(Jd * G) + Gd * J - Id.multiply(Gd * J)
C = C + mcut._a_b_one(J, Gd) - J.multiply(Gd) + mcut._a_b_one(G, Jd) - G.multiply(Jd)
motif_adj_mat = (C + C') / 3

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Je = build_Je(adj_mat)
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
C = mcut._a_b_one(Jd, Gs) - Jd.multiply(Gs * Je)
C = C + mcut._a_b_one(Gd, Js) - Gd.multiply(Js * Je)
C = C + Jd * Gs - Je.multiply(Jd * Gs) + Gd * Js - Je.multiply(Gd * Js)
C = C + mcut._a_b_one(Js, Gd) - Js.multiply(Gd * Je)
C = C + mcut._a_b_one(Gs, Jd) - Gs.multiply(Jd * Je)
motif_adj_mat = (C + C') / 3

if mam_weight_type == "product"
if motif_type == "func"
Gp = build_Gp(adj_mat)
Id = build_Id(adj_mat)
G = build_G(adj_mat)
C = mcut._a_b_one(Gp, G) - Gp.multiply(G)
C = C + Gp * G - Id.multiply(Gp * G)
C = C + mcut._a_b_one(G, Gp) - G.multiply(Gp)
motif_adj_mat = C + C'

if motif_type == "struc"
Gp = build_Gp(adj_mat)
Je = build_Je(adj_mat)
Gs = build_Gs(adj_mat)
C = mcut._a_b_one(Gp, Gs) - Gp.multiply(Gs * Je)
C = C + Gp * Gs - Je.multiply(Gp * Gs)
C = C + mcut._a_b_one(Gs, Gp) - Gs.multiply(Gp * Je)
motif_adj_mat = C + C'

return motif_adj_mat


function mam_M12(adj_mat, motif_type, mam_weight_type, mam_method): # pylint: disable=too-many-statements


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Jn = build_Jn(adj_mat)
J = build_J(adj_mat)
C = Jd * (Jn @ J) + Jn * (J @ Jd) + J * (Jn @ Jd)
motif_adj_mat = C + C'

if motif_type == "struc"
Jd = build_Jd(adj_mat)
J0 = build_J0(adj_mat)
Js = build_Js(adj_mat)
C = Jd * (J0 @ Js) + J0 * (Js @ Jd) + Js * (J0 @ Jd)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Jn = build_Jn(adj_mat)
G = build_G(adj_mat)
J = build_J(adj_mat)
C = Jd * (Jn @ G) + Gd * (Jn @ J)
C = C + Jn * (J @ Gd) + Jn * (G @ Jd)
C = C + J * (Jn @ Gd) + G * (Jn @ Jd)
motif_adj_mat = (C + C') / 3

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
J0 = build_J0(adj_mat)
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
C = Jd * (J0 @ Gs) + Gd * (J0 @ Js)
C = C + J0 * (Js @ Gd) + J0 * (Gs @ Jd)
C = C + Js * (J0 @ Gd) + Gs * (J0 @ Jd)
motif_adj_mat = (C + C') / 3

if mam_weight_type == "product"
if motif_type == "func"
Gp = build_Gp(adj_mat)
Jn = build_Jn(adj_mat)
G = build_G(adj_mat)
C = Gp * (Jn @ G) + Jn * (G @ Gp) + G * (Jn @ Gp)
motif_adj_mat = C + C'

if motif_type == "struc"
Gp = build_Gp(adj_mat)
J0 = build_J0(adj_mat)
Gs = build_Gs(adj_mat)
C = Gp * (J0 @ Gs) + J0 * (Gs @ Gp) + Gs * (J0 @ Gp)
motif_adj_mat = C + C'

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Id = build_Id(adj_mat)
J = build_J(adj_mat)
C = mcut._a_one_b(Jd, J) - Jd.multiply(J)
C = C + J * Jd - Id.multiply(J * Jd)
C = C + mcut._a_one_b(J, Jd) - J.multiply(Jd)
motif_adj_mat = C + C'

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Je = build_Je(adj_mat)
Js = build_Js(adj_mat)
C = mcut._a_one_b(Jd, Js) - Jd.multiply(Je * Js)
C = C + Js * Jd - Je.multiply(Js * Jd)
C = C + mcut._a_one_b(Js, Jd) - Js.multiply(Je * Jd)
motif_adj_mat = C + C'

if mam_weight_type == "mean"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Id = build_Id(adj_mat)
J = build_J(adj_mat)
G = build_G(adj_mat)
C = mcut._a_one_b(Jd, G) - Jd.multiply(G) + mcut._a_one_b(Gd, J) - Gd.multiply(J)
C = C + J * Gd - Id.multiply(J * Gd) + G * Jd - Id.multiply(G * Jd)
C = C + mcut._a_one_b(J, Gd) - J.multiply(Gd) + mcut._a_one_b(G, Jd) - G.multiply(Jd)
motif_adj_mat = (C + C') / 3

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Je = build_Je(adj_mat)
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
C = mcut._a_one_b(Jd, Gs) - Jd.multiply(Je * Gs)
C = C + mcut._a_one_b(Gd, Js) - Gd.multiply(Je * Js)
C = C + Js * Gd - Je.multiply(Js * Gd) + Gs * Jd - Je.multiply(Gs * Jd)
C = C + mcut._a_one_b(Js, Gd) - Js.multiply(Je * Gd)
C = C + mcut._a_one_b(Gs, Jd) - Gs.multiply(Je * Jd)
motif_adj_mat = (C + C') / 3

if mam_weight_type == "product"
if motif_type == "func"
Gp = build_Gp(adj_mat)
Id = build_Id(adj_mat)
G = build_G(adj_mat)
C = mcut._a_one_b(Gp, G) - Gp.multiply(G)
C = C + G * Gp - Id.multiply(G * Gp)
C = C + mcut._a_one_b(G, Gp) - G.multiply(Gp)
motif_adj_mat = C + C'

if motif_type == "struc"
Gp = build_Gp(adj_mat)
Je = build_Je(adj_mat)
Gs = build_Gs(adj_mat)
C = mcut._a_one_b(Gp, Gs) - Gp.multiply(Je * Gs)
C = C + Gs * Gp - Je.multiply(Gs * Gp)
C = C + mcut._a_one_b(Gs, Gp) - Gs.multiply(Je * Gp)
motif_adj_mat = C + C'

return motif_adj_mat


function mam_M13(adj_mat, motif_type, mam_weight_type, mam_method)

if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Jn = build_Jn(adj_mat)
C = Jd * (Jd @ Jn)
Cprime = Jn * (Jd @ Jd)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Jd = build_Jd(adj_mat)
J0 = build_J0(adj_mat)
C = Jd * (Jd @ J0)
Cprime = J0 * (Jd @ Jd)
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Jn = build_Jn(adj_mat)
Gd = build_Gd(adj_mat)
C = Jd * (Gd @ Jn) + Gd * (Jd @ Jn) + Jn * (Jd @ Gd)
motif_adj_mat = (C + C') / 4

if motif_type == "struc"
Jd = build_Jd(adj_mat)
J0 = build_J0(adj_mat)
Gd = build_Gd(adj_mat)
C = Jd * (Gd @ J0) + Gd * (Jd @ J0) + J0 * (Jd @ Gd)
motif_adj_mat = (C + C') / 4

if mam_weight_type == "product"
if motif_type == "func"
Gp = build_Gp(adj_mat)
Jn = build_Jn(adj_mat)
C = Gp * (Gp @ Jn)
Cprime = Jn * (Gp @ Gp)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gp = build_Gp(adj_mat)
J0 = build_J0(adj_mat)
C = Gp * (Gp @ J0)
Cprime = J0 * (Gp @ Gp)
motif_adj_mat = C + C' + Cprime

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_b_one(Jd, Jd) - Jd.multiply(Jd)
Cprime = Jd * Jd - Id.multiply(Jd * Jd)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_b_one(Jd, Jd) - Jd.multiply(Jd * Je)
Cprime = Jd * Jd - Je.multiply(Jd * Jd)
motif_adj_mat = C + C' + Cprime

if mam_weight_type == "mean"
if motif_type == "func"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_b_one(Jd, Gd) - Jd.multiply(Gd) + mcut._a_b_one(Gd, Jd) - Gd.multiply(Jd)
C = C + Jd * Gd - Id.multiply(Jd * Gd)
motif_adj_mat = (C + C') / 4

if motif_type == "struc"
Jd = build_Jd(adj_mat)
Gd = build_Gd(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_b_one(Jd, Gd) - Jd.multiply(Gd * Je)
C = C + mcut._a_b_one(Gd, Jd) - Gd.multiply(Jd * Je)
C = C + Jd * Gd - Je.multiply(Jd * Gd)
motif_adj_mat = (C + C') / 4

if mam_weight_type == "product"
if motif_type == "func"
Gp = build_Gp(adj_mat)
Id = build_Id(adj_mat)
C = mcut._a_b_one(Gp, Gp) - Gp.multiply(Gp)
Cprime = Gp * Gp - Id.multiply(Gp * Gp)
motif_adj_mat = C + C' + Cprime

if motif_type == "struc"
Gp = build_Gp(adj_mat)
Je = build_Je(adj_mat)
C = mcut._a_b_one(Gp, Gp) - Gp.multiply(Gp * Je)
Cprime = Gp * Gp - Je.multiply(Gp * Gp)
motif_adj_mat = C + C' + Cprime

return motif_adj_mat


function mam_Mcoll(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
C = Jn * (J @ J')
motif_adj_mat = C

if motif_type == "struc"
Js = build_Js(adj_mat)
J0 = build_J0(adj_mat)
C = J0 * (Js @ Js')
motif_adj_mat = C

if mam_weight_type == "mean"
if motif_type == "func"
G = build_G(adj_mat)
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
C = Jn * (J @ G') + Jn * (G @ J')
motif_adj_mat = C / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = J0 * (Js @ Gs') + J0 * (Gs @ Js')
motif_adj_mat = C / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Jn = build_Jn(adj_mat)
C = Jn * (G @ G')
motif_adj_mat = C

if motif_type == "struc"
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = J0 * (Gs @ Gs')
motif_adj_mat = C

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Id = build_Id(adj_mat)
C = J * J' - Id.multiply(J * J')
motif_adj_mat = C

if motif_type == "struc"
Js = build_Js(adj_mat)
Je = build_Je(adj_mat)
C = Js * Js' - Je.multiply(Js * Js')
motif_adj_mat = C

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = J * G' - Id.multiply(J * G')
C = C + G * J' - Id.multiply(G * J')
motif_adj_mat = C / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = Js * Gs' - Je.multiply(Js * Gs')
C = C + Gs * Js' - Je.multiply(Gs * Js')
motif_adj_mat = C / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = G * G' - Id.multiply(G * G')
motif_adj_mat = C

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = Gs * Gs' - Je.multiply(Gs * Gs')
motif_adj_mat = C

return motif_adj_mat


function mam_Mexpa(adj_mat, motif_type, mam_weight_type, mam_method)


if mam_method == "dense"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
C = Jn * (J' @ J)
motif_adj_mat = C

if motif_type == "struc"
Js = build_Js(adj_mat)
J0 = build_J0(adj_mat)
C = J0 * (Js' @ Js)
motif_adj_mat = C

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
Jn = build_Jn(adj_mat)
G = build_G(adj_mat)
C = Jn * (J' @ G) + Jn * (G' @ J)
motif_adj_mat = C / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = J0 * (Js' @ Gs) + J0 * (Gs' @ Js)
motif_adj_mat = C / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Jn = build_Jn(adj_mat)
C = Jn * (G' @ G)
motif_adj_mat = C

if motif_type == "struc"
Gs = build_Gs(adj_mat)
J0 = build_J0(adj_mat)
C = J0 * (Gs' @ Gs)
motif_adj_mat = C

if mam_method == "sparse"
if mam_weight_type == "unweighted"
if motif_type == "func"
J = build_J(adj_mat)
Id = build_Id(adj_mat)
C = J' * J - Id.multiply(J' * J)
motif_adj_mat = C

if motif_type == "struc"
Js = build_Js(adj_mat)
Je = build_Je(adj_mat)
C = Js' * Js - Je.multiply(Js' * Js)
motif_adj_mat = C

if mam_weight_type == "mean"
if motif_type == "func"
J = build_J(adj_mat)
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = J' * G - Id.multiply(J' * G)
C = C + G' * J - Id.multiply(G' * J)
motif_adj_mat = C / 2

if motif_type == "struc"
Js = build_Js(adj_mat)
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = Js' * Gs - Je.multiply(Js' * Gs)
C = C + Gs' * Js - Je.multiply(Gs' * Js)
motif_adj_mat = C / 2

if mam_weight_type == "product"
if motif_type == "func"
G = build_G(adj_mat)
Id = build_Id(adj_mat)
C = G' * G - Id.multiply(G' * G)
motif_adj_mat = C

if motif_type == "struc"
Gs = build_Gs(adj_mat)
Je = build_Je(adj_mat)
C = Gs' * Gs - Je.multiply(Gs' * Gs)
motif_adj_mat = C

return motif_adj_mat
=#
