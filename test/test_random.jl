@testset "Random GHZ Objects" begin
    s = rand(GHZState, 3, 2)
    @test isa(s, GHZState)
    g = rand(GHZGate, 3)
    @test isa(g, GHZGate)
    m = rand(GHZMeasure, 3)
    @test isa(m, GHZMeasure)
end
