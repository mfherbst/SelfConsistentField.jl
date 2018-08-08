function compute_pulay_error(fock::AbstractArray{T,3},
                             density::AbstractArray{T,3},
                             overlap::AbstractArray{T,2}) where T
    n_bas, _, n_spin = size(density)
    @assert n_spin == 1 || n_spin == 2

    # TODO We should avoid the memory allocation
    error = Array{eltype(density)}(undef, n_bas, n_bas, n_spin)
    for σ in 1:n_spin
        error[:, :, σ] = compute_pulay_error(view(fock, :, :, σ),
                                             view(density, :, :, σ),
                                             overlap)
    end
    return error
end

function compute_pulay_error(fock::AbstractArray{T,2},
                             density::AbstractArray{T,2},
                             overlap::AbstractArray{T,2}) where T
    n_bas, _ = size(density)
    @assert size(fock) == (n_bas, n_bas)
    @assert size(density) == (n_bas, n_bas)
    @assert size(overlap) == (n_bas, n_bas)
    return overlap * density * fock - fock * density * overlap
end
