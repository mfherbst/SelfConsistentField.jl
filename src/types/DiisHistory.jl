"""
    Container type for the history of one spin type for cDIIS
"""
struct DiisHistory
    iterate::CircularBuffer
    error::Union{Nothing,CircularBuffer}
    iterationhistory::CircularBuffer
    density::Union{Nothing,CircularBuffer}
    energies::Union{Nothing,CircularBuffer}
    n_diis_size::Int
end

function DiisHistory(n_diis_size::Int; need_error = false, need_density = false, need_energies = false, params...)
    DiisHistory(CircularBuffer{AbstractArray}(n_diis_size),
        need_error ? CircularBuffer{AbstractArray}(n_diis_size) : nothing,
        CircularBuffer{AbstractArray}(n_diis_size),
        need_density ? CircularBuffer{AbstractArray}(n_diis_size) : nothing,
        need_energies ? CircularBuffer{Dict}(n_diis_size) : nothing,
        n_diis_size
    )
end

function copy(cb::CircularBuffer{T}) where {T}
    append!(CircularBuffer{T}(cb.capacity), cb)
end

copy(::Nothing) = nothing

function copy(dh::DiisHistory)
    DiisHistory(copy(dh.iterate), copy(dh.error), copy(dh.iterationhistory), copy(dh.density), copy(dh.energies), dh.n_diis_size)
end

"""
    Pushes current iterate and error matrices to histories of both spin types
"""
function push_iterate!(history::DiisHistory, iterate::AbstractArray, error::Union{AbstractArray,Nothing} = nothing, density::Union{AbstractArray,Nothing} = nothing, energies::Union{Dict,Nothing} = nothing)
    pushfirst!(history.iterate,  iterate)

    # Push difference to previous iterate if no error given
    if history.error != nothing
        pushfirst!(history.error,
                   error == nothing ? iterate - history.iterate[1] : error)
    end
    history.density != nothing && pushfirst!(history.density, density)
    history.energies != nothing && pushfirst!(history.energies, energies)
    return history
end

"""
    Loop over Spin submatrices and push them to the relevant DiisHistory struct
"""
function push_iterate(history::Tuple{DiisHistory,DiisHistory}, state::ScfIterState)
    push_spinblock!(storage, σ) = push_iterate!(storage,
                           spin(get_iterate_matrix(state), σ),
                           spin(state.error_pulay, σ),
                           spin(state.density, σ),
                           state.energies)

    if spincount(state) == 1
        new_historyα = push_spinblock!(copy(history[1]), 1)
        new_history = new_historyα, new_historyα
    else
        new_history = push_spinblock!(copy(history[1]), 1), push_spinblock!(copy(history[2]), 2)
    end
    return new_history
end

"""
    Removes 'count' Items from DiisHistory.
    Warning: Assumes all buffers are of the same length.
"""
function purge_from_history!(history::DiisHistory, count::Int)
    for i in 1:count
        length(history.iterate) > 0 && pop!(history.iterate)
        length(history.iterationhistory) > 0 && pop!(history.iterationhistory)
        history.error != nothing && length(history.error) > 0 && pop!(history.error)
        history.density != nothing && length(history.density) > 0 && pop!(history.density)
        history.energies != nothing && length(history.energies) > 0 && pop!(history.energies)
    end
end

