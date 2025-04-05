@testset "Pauli Noise Application" begin
    s = GHZState(3, 1)
    noisy_op = PauliNoiseOp(3, 1, 0.3, 0.3, 0.3)
    s2 = apply!(s, noisy_op)
    @test isa(s2, GHZState)
end
