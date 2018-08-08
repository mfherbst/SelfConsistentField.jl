import Base.+

# TODO For the future
"""
Struct representing a density
"""
struct Density
    density::Union{Nothing,AbstractArray}
    orbcoeff::Union{Nothing,Array{Float64}}
    occupation::Union{Nothing,Array{Float64}}
end

function +(d1::Density, d2::Density)
    if d1.density != nothing and d2.density != nothing
        sum = d1.density + d2.density
    else
        sum = nothing
    end

    # TODO I hope the next statement is not linear, but constant in time
    if d1.orbcoeff == d2.orbcoeff
        return Density(sum, d2.orbcoeff,
                       d1.occupation + d2.occupation)
    else
        return Density(sum, nothing, nothing)
    end
end


function *(v::Number, d::Density)
    if d.density != nothing
        newdens = d.density * v
    else
        newdens = nothing
    end
    if d.occupation != nothing
        newocc = v * d.occupation
    else
        newocc = nothing
    end
    return Density(newdens, d.orbcoeff, newocc)
end
