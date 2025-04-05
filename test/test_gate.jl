@testset "GHZGate Application" begin
    s = GHZState(3, 2)
    g = GHZPreserving.GHZGate(3, 1, fill(1,2), fill(1,3), 1, 2)
    s2 =GHZPreserving.apply!(s,g)
    @test isa(s2, GHZState)
end


s=GHZState(3, 2)
g=GHZPreserving.GHZGate(3, 1, [1,2], [1,2,3], 1, 2)
GHZPreserving.apply!(s,g) # test apply with array indices, should work without error
