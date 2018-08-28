@testset "FockIterState structure" begin
    state = FockIterState([1,2,3], [4,5,6], Dict("a" => 3.0, "b" => 4.0), [7, 8, 9], [10, 11, 12])
    @test state <: ScfIterState
    @test state.fock == [1,2,3] && state.error_pulay == [4,5,6] &&
      state.energies == Dict("a" => 3.0, "b" => 4.0) && fock.orben == [7,8,9] &&
      density == [10,11,12]
    @test get_iterate_matrix(state) == [1,2,3]
    @test update_iterate_matrix(state, [42,63,78]) == FockIterState([42,63,78], [4,5,6], Dict("a" => 3.0, "b" => 4.0), [7, 8, 9], [10, 11, 12])
end

@testset "Report" begin
    @test Dict(:memory => Set(:debug, :warn)) isa LogLevel
    @test Dict(:max_diis_size => 5) isa Parameters
    rmess =  ReportMessage("testmsg")
    @test rmess.msg == "testmsg"
    @test rmess.data == nothing
    # TODO not done yet
end

@testset "push/purge_state" begin
    need_error = true
    need_density = true

    struct TestAlgorithm end
    needs_error(::TestAlgorithm) = need_error

    function testentry()
        testenergies = Dict( "e0" => rand(), "e1" => rand(), "e2" => rand(), "total" => rand())
        (rand(6,6), rand(6,6), rand(6,6), testenergies)
    end
    function check_entries(dstate::DiisState, entrylist::Vector, range::Vector{Int})
        for i in eachindex(range)
            @test state.iterate[i] == entrylist[range[i]][1]
            @test state.error[i] == needs_error(TestAlgorithm()) ? entrylist[range[i]][2] : nothing
            @test state.density[i] == needs_density(TestAlgorithm()) ? entrylist[range[i]][3] : nothing
            @test state.energies[i] == entrylist[range[i]][4]
        end
    end
    function check_diisstate_length(dstate::DiisState, expected_length::Int)
        @test expected_length == length(dstate.iterate) == length(dstate.energies)
        if needs_error(TestAlgorithm())
            @test expected_length == length(dstate.error)
        else
            @test dstate.error == nothing
        end
        if needs_density(TestAlgorithm())
            @test expected_length == length(dstate.density)
        else
            @test dstate.density == nothing
        end
    end
    function run_pushiterate_tests()
        dstate = DIISstate(5)
        check_diisstate_length(dstate.iterate, 0)

        entrylist = map!(x -> testentry(), 1:8)
        map(i -> push_iterate!(TestAlgorithm(), dstate, entrylist[i]...), 1:3)
        check_entries(dstate, entrylist, [1,2,3])
        check_diisstate_length(dstate.iterate, 3)
        purge_from_state(TestAlgorithm(), dstate, 1)
        check_diisstate_length(dstate.iterate, 2)
        map(i -> push_iterate!(TestAlgorithm(), dstate, entrylist[i]...), 4:6)
        check_entries(dstate, entrylist, [6,5,4,2,1])
        map(i -> push_iterate!(TestAlgorithm(), dstate, entrylist[i]...), 7:7)
        check_entries(dstate, entrylist, [7,6,5,4,2])
    end

    need_error = true
    need_density = true
    run_pushiterate_tests()
    need_error = false
    need_density = true
    run_pushiterate_tests()
    need_error = true
    need_density = false
    run_pushiterate_tests()
    need_error = false
    need_density = false
    run_pushiterate_tests()
end
