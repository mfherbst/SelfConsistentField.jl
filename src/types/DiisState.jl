"""
    Container type for the state of one spin type for cDIIS
"""
mutable struct DiisState
    iterate::CircularBuffer
    error::CircularBuffer
    iterationstate::CircularBuffer
    density::CircularBuffer
    energies::CircularBuffer
    n_diis_size::Int

    function DiisState(n_diis_size::Int)
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
    Pushes current iterate and error matrices to states of both spin types
"""
function push_iterate!(algorithm::Algorithm, state::DiisState, iterate::AbstractArray, error::Union{AbstractArray,Nothing} = nothing, density::Union{AbstractArray,Nothing} = nothing, energies::Union{Dict,Nothing} = nothing)
    pushfirst!(state.iterate,  iterate)

    # Push difference to previous iterate if no error given
    if needs_error(algorithm)
        pushfirst!(state.error,
                   error == nothing ? iterate - state.iterate[1] : error)
    end
    needs_density(algorithm) ? pushfirst!(state.density, density) : 0
    energies != nothing ? pushfirst!(state.energies, energies) : 0
end

"""
    Removes 'count' Items from DiisState.
    Warning: Assumes all buffers are of the same length.
"""
function purge_from_state(algorithm::Algorithm, state::DiisState, count::Int)
    if count == 0
        return
    end
    @assert length(state.iterate) == length(state.iterationstate)
    for i in 1:count
        pop!(state.iterate)
        pop!(state.iterationstate)
        needs_error(algorithm) && pop!(state.error)
        needs_density(algorithm) && pop!(state.density)
        pop!(state.energies)
    end
end

