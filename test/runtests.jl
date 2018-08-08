@static if VERSION < v"0.7-"
    const Test = Base.Test
end

using Test
using SelfConsistentField

include("functionality_hartee_fock.jl")
