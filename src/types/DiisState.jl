"""
    Container type for the history of one spin type for cDIIS
"""
mutable struct DiisHistory
    iterate::CircularBuffer
    error::CircularBuffer
    iterationhistory::CircularBuffer
    density::CircularBuffer
    energies::CircularBuffer
    n_diis_size::Int

    function DiisHistory(n_diis_size::Int)
        new(CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{Dict}(n_diis_size),
            n_diis_size
        )
    end
end

"""
    Pushes current iterate and error matrices to historys of both spin types
"""
function push_iterate!(algorithm::Algorithm, history::DiisHistory, iterate::AbstractArray, error::Union{AbstractArray,Nothing} = nothing, density::Union{AbstractArray,Nothing} = nothing, energies::Union{Dict,Nothing} = nothing)
    pushfirst!(history.iterate,  iterate)

    # Push difference to previous iterate if no error given
    if needs_error(algorithm)
        pushfirst!(history.error,
                   error == nothing ? iterate - history.iterate[1] : error)
    end
    needs_density(algorithm) ? pushfirst!(history.density, density) : 0
    energies != nothing ? pushfirst!(history.energies, energies) : 0
end

"""
    Removes 'count' Items from DiisHistory.
    Warning: Assumes all buffers are of the same length.
"""
function purge_from_history(algorithm::Algorithm, history::DiisHistory, count::Int)
    if count == 0
        return
    end
    @assert length(history.iterate) == length(history.iterationhistory)
    for i in 1:count
        pop!(history.iterate)
        pop!(history.iterationhistory)
        needs_error(algorithm) && pop!(history.error)
        needs_density(algorithm) && pop!(history.density)
        pop!(history.energies)
    end
end

