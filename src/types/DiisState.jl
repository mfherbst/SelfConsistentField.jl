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
function push_iterate!(acc::Accelerator, state::DiisState, iterate::AbstractArray, error::Union{AbstractArray,Nothing} = nothing, density::Union{AbstractArray,Nothing} = nothing, energies::Union{Dict,Nothing} = nothing)
    pushfirst!(state.iterate,  iterate)

    # Push difference to previous iterate if no error given
    pushfirst!(state.error,
               error != nothing ? error : iterate - state.iterate[1])
    needs_density(acc) ? pushfirst!(state.density, density) : 0
    energies != nothing ? pushfirst!(state.energies, energies) : 0
end

"""
    Removes 'count' Items from DiisState.
    Warning: Assumes you have already pushed a new iterate to the state.
"""
function purge_from_state(acc::Accelerator, state::DiisState, count::Int)
    if count == 0
        return
    end
    for i in 1:count + 1
        pop!(state.iterate)
        pop!(state.iterationstate)
        needs_error(acc) & (state.error != nothing) ? pop!(state.error) : 0
        needs_density(acc) & (state.density != nothing) ? pop!(state.density) : 0
    end
end

