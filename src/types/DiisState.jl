"""
    Container type for the state of one spin type for cDIIS
"""
mutable struct DiisState
    iterate::CircularBuffer
    error::CircularBuffer
    iterationstate::CircularBuffer
    n_diis_size::Int

    function DiisState(n_diis_size::Int)
        new(CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            CircularBuffer{AbstractArray}(n_diis_size),
            n_diis_size
        )
    end
end

"""
    Pushes current iterate and error matrices to states of both spin types
"""
function push_iterate!(state::DiisState, iterate::AbstractArray, error::Union{AbstractArray,Nothing} = nothing)
    pushfirst!(state.iterate,  iterate)

    # Push difference to previous iterate if no error given
    pushfirst!(state.error,
               error != nothing ? error : iterate - state.iterate[1])
end

"""
    Removes 'count' Items from DiisState.
    Warning: Assumes you have already pushed a new iterate to the state.
"""
function purge_from_state(state::DiisState, count::Int)
    for i in 1:count + 1
        pop!(state.iterate)
        pop!(state.error)
        pop!(state.iterationstate)
    end
end

