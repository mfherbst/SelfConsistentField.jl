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
function spincount(iterate::Iterate)
    spincount(iterate.fock)
end

"""
    Iterator for indices of spins
"""
function spinloop(obj::AbstractArray)
    1:spincount(obj)
end
function spinloop(iterate::Iterate)
    1:spincount(iterate.fock)
end
