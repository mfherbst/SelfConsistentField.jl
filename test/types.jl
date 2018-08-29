using SelfConsistentField: FockIterState, DiisState,
          ScfIterState, get_iterate_matrix, update_iterate_matrix,
          LogLevel, Parameters, ReportMessage, push_iterate!,
          Algorithm, purge_from_state, ScfConvergence

@testset "FockIterState structure" begin
    state = FockIterState([1,2,3], [4,5,6], Dict("a" => 3.0, "b" => 4.0), [7, 8, 9], [10, 11, 12], [13, 14, 15])
    @test state isa ScfIterState
    @test state.fock == [1,2,3] && state.error_pulay == [4,5,6] &&
      state.energies == Dict("a" => 3.0, "b" => 4.0) && state.orbcoeff == [7,8,9] &&
      state.orben == [10, 11, 12] && state.density == [13, 14, 15]
    @test get_iterate_matrix(state) == [1,2,3]
    updated_iterate = update_iterate_matrix(state, [42,63,78])
    @test updated_iterate.fock == [42,63,78] && updated_iterate.error_pulay == [4,5,6] &&
      updated_iterate.energies == Dict("a" => 3.0, "b" => 4.0) && updated_iterate.orbcoeff == [7,8,9] &&
      updated_iterate.orben == [10, 11, 12] && updated_iterate.density == [13, 14, 15]
end

@testset "Report" begin
    @test LogLevel(:memory => Set([:debug, :warn])) isa Dict{Symbol, Set}
    @test Parameters(:max_diis_size => 5) isa Dict{Symbol, Any}
    rmess =  ReportMessage("testmsg")
    @test rmess.msg == "testmsg"
    @test rmess.data == nothing
    # TODO not done yet
end

struct TestAlgorithm <: Algorithm end
@testset "push/purge_state" begin
    need_error = true
    need_density = true

    SelfConsistentField.needs_error(::TestAlgorithm) = need_error
    SelfConsistentField.needs_density(::TestAlgorithm) = need_density

    function testentry()
        testenergies = Dict( "e0" => rand(), "e1" => rand(), "e2" => rand(), "total" => rand())
        (rand(6,6), rand(6,6), rand(6,6), testenergies)
    end
    function check_entries(state::DiisState, entrylist::Vector, range::Vector{Int})
        for i in eachindex(range)
            @test state.iterate[i] == entrylist[range[i]][1]
            @test SelfConsistentField.needs_error(TestAlgorithm()) ? state.error[i] == entrylist[range[i]][2] : true
            @test SelfConsistentField.needs_density(TestAlgorithm()) ? state.density[i] == entrylist[range[i]][3] : true
            @test state.energies[i] == entrylist[range[i]][4]
        end
    end
    function check_diisstate_length(dstate::DiisState, expected_length::Int)
        @test expected_length == length(dstate.iterate) == length(dstate.energies)
        if SelfConsistentField.needs_error(TestAlgorithm())
            @test expected_length == length(dstate.error)
        else
            @test dstate.error == []
        end
        if SelfConsistentField.needs_density(TestAlgorithm())
            @test expected_length == length(dstate.density)
        else
            @test dstate.density == []
        end
    end
    function run_pushiterate_tests()
        dstate = DiisState(5)
        check_diisstate_length(dstate, 0)

        entrylist = map(x -> testentry(), 1:7)

        @testset "push 1:3" begin
            map(i -> push_iterate!(TestAlgorithm(), dstate, entrylist[i]...), 1:3)
            map(i -> push!(dstate.iterationstate, rand(5)), 1:3)
            check_entries(dstate, entrylist, [3,2,1])
            check_diisstate_length(dstate, 3)
        end

        @testset "purge 1" begin
            purge_from_state(TestAlgorithm(), dstate, 1)
            check_diisstate_length(dstate, 2)
        end

        @testset "push 4:6" begin
            map(i -> push_iterate!(TestAlgorithm(), dstate, entrylist[i]...), 4:6)
            map(i -> push!(dstate.iterationstate, rand(5)), 4:6)
            check_entries(dstate, entrylist, [6,5,4,3,2])
            check_diisstate_length(dstate, 5)
        end

        @testset "push 7" begin
            map(i -> push_iterate!(TestAlgorithm(), dstate, entrylist[i]...), 7:7)
            map(i -> push!(dstate.iterationstate, rand(5)), 7:7)
            check_entries(dstate, entrylist, [7,6,5,4,3])
            check_diisstate_length(dstate, 5)
        end
    end

    @testset "needs_error == needs_density == true" begin
        need_error = true
        need_density = true
        run_pushiterate_tests()
    end
    @testset "needs_error == false, needs_density == true" begin
        need_error = false
        need_density = true
        run_pushiterate_tests()
    end
    @testset "needs_error == true, needs_density == false" begin
        need_error = true
        need_density = false
        run_pushiterate_tests()
    end
    @testset "needs_error == false, needs_density == false" begin
        need_error = false
        need_density = false
        run_pushiterate_tests()
    end
end

@testset "ScfConvergence structure" begin
    conv = ScfConvergence(3.0, Dict("total" => 4.4, "e1" => 5.7), false)
    @test conv.error_norm == 3.0
    @test conv.energy_change == Dict("total" => 4.4, "e1" => 5.7)
    @test conv.is_converged == false
end
