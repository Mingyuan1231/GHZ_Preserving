@testset "GHZState Construction" begin
    s = GHZState(3, 2)
    @test s.qubit_num == 3
    @test s.ghz_num == 2
    @test length(s.phases) == 6
    @test all(x -> x == false, s.phases)
end
