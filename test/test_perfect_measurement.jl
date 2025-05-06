using Test
using GHZPreserving
using StatsBase

@testset "Perfect GHZ Measurement vs QC Measurement" begin
    #start with simple GHZ state, 3 qubits, 2 raw states
    state=GHZState(3,2)
    qcstate=MixedDestabilizer(state)
    #perfect measurement, we want to measure 1st state in X basis
    ghzmeasure=GHZMeasure(3, 1, 1)
    qcmeaure=toQCcircuit(ghzmeasure)

    measured_state, ghz_result=measure!(state, ghzmeasure)
    #we get ghz_result 0 for + phases and 1 for - phases, 

    #qc part need to keep the measurement result for each qubit
    qc_result=[]

    for i in 1:length(qcmeaure)
        qcstate, qc_res=projectrand!(qcstate, qcmeaure[i])
        push!(qc_result, qc_res)
    end

    #check results
    if qcmeaure[1] isa sMX
        eigen_results = map(re -> re == 0x02 ? -1 : 1, qc_result)
        qc = Int(prod(eigen_results) == 1)
    else
        qc=Int(all(x -> x == qc_result[1], qc_result))
    end

    #we will get qc=1 for measurement result as +1 and qc=0 for measurement result as -1.
    @test Int(ghz_result[1])==1-qc   
    #the result agree with each other, although qc is 1 and ghz is 0, 
    #but this is due to different convention, the meaning is the same.
end 