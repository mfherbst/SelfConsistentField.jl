function is_closed_shell(obj)
    return obj.n_elec[1] == obj.n_elec[2]
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
