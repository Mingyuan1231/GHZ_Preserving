@testset "GHZGate Application" begin
    s = GHZState(3, 2)
    g = GHZGate(3, 1, fill(1,2), fill(1,3), 1, 2)
    s2 =apply!(s,g)
    @test isa(s2, GHZState)
end
