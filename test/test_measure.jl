@testset "GHZ Measurement" begin
    s = GHZState(3, 1, BitVector([true, false, true]))
    op = GHZMeasure(3, 1, 1)
    new_s, result = Measure!(s, op)
    @test isa(new_s, GHZState)
    @test result == [true] || result == [false]
end
