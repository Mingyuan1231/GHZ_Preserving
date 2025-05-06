using Test
using GHZPreserving
using StatsBase, Distributions

@testset "Noisy GHZ Measurement vs QuantumClifford Equivalence" begin
    #start with simple GHZ state, 3 qubits, 2 raw states
    state=GHZState(3,2)
    qcstate=MixedDestabilizer(state)
    p=0.8
    #this time we do noisy measurement with 0.8 success rate, we want to measure 1st state in X basis
    ghznoisymeasure=NoisyGHZMeasureNoisyReset(GHZMeasure(3, 1, 1),p)
    #since we need to deal with noisy, toQCcircuit will reture series of projectors.
    qcnoisymeaure=toQCcircuit(ghznoisymeasure)

    sample_number=100000
    #ghz part
    ghz_result=[]
    for i in 1:sample_number
        state,ghz_res=apply!(state, ghznoisymeasure)
        push!(ghz_result, ghz_res)
    end

    mean(ghz_result)    #this return 1-p, due to ghz preserving convention.


    #qc part
    qc_result=[]

    for j in 1:sample_number
        qc_res=[]

        for i in 1:length(qcnoisymeaure)
            qcstate, res=projectrand!(qcstate, qcnoisymeaure[i])
            push!(qc_res, res)
        end

        #check results
        if qcnoisymeaure[1] isa sMX
            #0x02 is -1, 0x00 is +1
            eigen_results = map(re -> re == 0x02 ? -1 : 1, qc_res)
            if rand()>p
                eigen_results*=-1   #flip result
            end
            qc = Int(prod(eigen_results) == 1)
        else
            # Convert 0x00 → 0, 0x02 → 1
            bitvals = map(r -> r == 0x02 ? 1 : 0, qc_res)

            # Add correlated measurement errors between adjacent qubits
            rber = Bernoulli(1 - p)
            error_phase = [rand(rber) for _ in 1:length(bitvals)]
            error_correlation = [error_phase[i] ⊻ error_phase[i+1] for i in 1:length(bitvals)-1]

            # Inject the errors
            bitvals_noisy = [bitvals[i] ⊻ error_correlation[i] for i in 1:length(error_correlation)+1]

            # Check if pairwise identical
            qc = Int(all(i -> bitvals_noisy[i] == bitvals_noisy[i+1], 1:length(bitvals_noisy)-1))
        end
        push!(qc_result, qc)
    end
    mean(qc_result)   #this return p, due to qc convention, and both result agree with each other.

    #test 5 sigma
    @test abs(mean(ghz_result)[1]+mean(qc_result)-1)<5*1/sqrt(sample_number)  
end
