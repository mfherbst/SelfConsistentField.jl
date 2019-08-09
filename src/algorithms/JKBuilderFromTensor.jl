using TensorOperations

"""
Two electron term builder, which builds J + K combined
from an electron repulsion eri tensor, which is stored completely in memory.
"""
struct JKBuilderFromTensor <: SelfConsistentField.TwoElectronBuilder
    eri::AbstractArray{Float64, 4}  # make generic in T
end

"""
Build a two electron term from a passed density.

If the density object is a single matrix, then it is interpreted as the alpha density
and 2*density is used to build J and density is used to build K.
A single AbstractArray will be returned

If the density object is a stack of two matrices, then the first density is interpreted
as the alpha density and the second as the beta density. The sum is used for J
and the individual terms for K. A tuple of AbstractArrays will be returned.
"""
function add_2e_term!(out::AbstractArray,
                      builder::JKBuilderFromTensor,
                      density::AbstractArray,
                      total_density::AbstractArray)
    n_bas, _, n_spin = size(density)
    Dtot = total_density
    eri = builder.eri

    @assert n_spin == 1 || n_spin == 2
    @assert size(Dtot) == (n_bas, n_bas)
    @assert size(density) == (n_bas, n_bas, n_spin)

    if n_spin == 1
        Da = view(density, :, :, 1)
        @tensor begin
            JKa[μ,ν] := eri[α,β,μ,ν] * Dtot[α,β] - eri[μ,β,α,ν] * Da[α,β]
        end
        out[:, :, 1] += JKa
    elseif n_spin == 2
        Da = view(density, :, :, 1)
        Db = view(density, :, :, 2)
        @tensor begin
            J[μ,ν] := eri[α,β,μ,ν] * Dtot[α,β]
        end
        @tensor begin
            JKa[μ,ν] := J[μ,ν] - eri[μ,β,α,ν] * Da[α,β]
            JKb[μ,ν] := J[μ,ν] - eri[μ,β,α,ν] * Db[α,β]
        end
        out[:, :, 1] += JKa
        out[:, :, 2] += JKb
    end
end
