using Test
using GHZPreserving
using StatsBase
using HypothesisTests

@testset "NoisyGate GHZPreserving vs QuantumClifford consistency" begin
    #set a simple series of gates
    gates=[Hgroup{3}(2, 1, 2)]

    qcgates=toQCcircuit(gates)
    noise = PauliNoise(1/5)
    qcgates_noisy = [NoisyGate(g, noise) for g in qcgates] 

    state_cache=[]
    qcstate_cache=[]
    
    for i in 1:100000
        #generate a GHZ state, 3 qubits, 3 raw states, fidelity 1.0
        state=GHZState(3,2)
        qcstate=Stabilizer(state)

        gn=NoisyGate(gates[1], noise)

        for g in gates
            apply!(state, gn)
        end
  
        idx=GHZPreserving.ghzstateindex(state,3)
        push!(state_cache, idx)

        for g in qcgates_noisy
            apply!(qcstate, g)
        end
 
        qcidx=GHZPreserving.ghzstateindex(qcstate,3)
        push!(qcstate_cache, qcidx)
    end

    ks_noise= ApproximateTwoSampleKSTest(Float64.(state_cache), Float64.(qcstate_cache))

    @test pvalue(ks_noise) > 0.05
end
