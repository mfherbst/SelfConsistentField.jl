@static if VERSION < v"0.7-"
    const Test = Base.Test
end

using Test
using SelfConsistentField

@testset "SelfConsistentField" begin
    #include("functionality_hartee_fock.jl")
    include("types.jl")
end
