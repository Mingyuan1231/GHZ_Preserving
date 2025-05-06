using Test
using GHZPreserving

@testset "Perfect gate translation via toQCcircuit (Hgroup only)" begin
        #generate series of random gates in GHZPreserving format(only H group is considered here for simplicity)
    gates=[]
    for _ in 1:10
        i = rand(1:3)
        j = rand(setdiff(1:3, [i]))  # j ≠ i
        gate = Hgroup{3}(rand(1:6), i, j)
        push!(gates, gate)
    end

    #transfer to QC format
    qcgates=toQCcircuit(gates)

    #generate random GHZ states, 3 qubits, 3 raw states, fidelity 0.7
    state=rand(GHZState,3,3,0.7)

    #transfer to QC format
    qcstate=Stabilizer(state)

    #apply gates to the state in both formats
    for g in gates
        apply!(state, g)
    end

    for g in qcgates
        apply!(qcstate, g)
    end

    # Check if both representations are equivalent
    @test canonicalize!(qcstate)==canonicalize!(Stabilizer(state))
end