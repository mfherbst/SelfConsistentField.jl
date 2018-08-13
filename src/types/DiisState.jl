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
    Obtaining views for specific spins
"""
function spin(obj::AbstractArray, dim::Int)
    view(obj, ntuple(x -> Colon(), ndims(obj) - 1)..., dim)
end

"""
    Get number of spins in object
"""
function spincount(obj::AbstractArray)
    size(obj, ndims(obj))
end

"""
    Iterator for indices of spins
"""
function spinloop(obj::Union{AbstractArray, Accelerator})
    typeof(obj) == Accelerator ?
        (1:spincount(obj.state.iterate)) :
        (1:spincount(obj))
end

"""
    remove 'count' Items from DiisState
"""
function purge_from_state(state::DiisState, count::Int)
    for i in 1:2*count
        pop!(state.iterate)
        pop!(state.error)
        pop!(state.iterationstate)
    end
end

