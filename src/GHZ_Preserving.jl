module GHZPreserving

using QuantumClifford
using LinearAlgebra
using Random
using Distributions

export GHZState, GHZMeasure, GHZGate, Measure!, rand, tensor!,
    Hgroup, Fgroup, PauliGroup,
    CNOT, CZ, F1, depolarize!, 
    PauliNoiseOp, NoisyMeasure, NoisyMeasureNoisyReset

"""
convert bit(phase) to int
"""
function bit_to_int(bits)
    reduce(⊻,(bit<<(index-1) for (index,bit) in enumerate(reverse(bits)))) + 1 # +1 so that we use julia indexing convenctions
end

"""
convert int to bit(phase)
"""
function int_to_bit(int,digits)     #digits should be 2n, as dealing with 2 input state, n qubit each
    int = int - 1 # -1 so that we use julia indexing conventions
    Bool[int>>shift&0x1 for shift in digits-1:-1:0]
end

################
#Generate permutation based on number of qubits
################

"""
For input n, generate all possible states for n qubit state 
(2^n states in total)
This funtion return all states in the vector form, the order is ++..++, ++..+- , ... to --..--, with binary carrying from last digit.
"""
function generate_state(n) #generate all possible states for n qubit state
    state=canonicalize!(ghz(n))
    states=[]
    for i in 0:2^n-1
        newstate = copy(state)
        phases(newstate) .= 2*reverse(digits(i;base=2, pad=n)) #reverse for common convention.
        push!(states, newstate)
    end
    return states   #return all states in the vector form, the order is ++..++, ++..+- , ... to --..--, with binary carrying from last digit.
end

"""
Restore the original order of phase, beacuse canonicalize will change the order of phase
"""
function adjust_order(phase,n)
    restored=Vector(undef, 2n)
    restored[1]=phase[1]
    restored[2:n]=phase[3:n+1]
    restored[n+1]=phase[2]
    restored[n+2:2n]=phase[n+2:2n]
    return restored
end

"""
Generate all gates for H gorup
There are six gates in H group, which are SWAP, CNOT12, INV_CNOT21, INV_CNOT12, Identity, CNOT21
The idx of these gate already be found in the twoqubitgates, which are [2,3,62,63,121,151]
This function return all these gates in the vector form in the order of [SWAP,CNOT12,INV_CNOT21,INV_CNOT12,Identity,CNOT21]
"""
function generate_Hgroup() #generate H group 
    sameinx=[2,3,62,63,121,151] #in the order of [SWAP,CNOT12,INV_CNOT21,INV_CNOT12,Identity,CNOT21]
    samegate=[]
    for i in 1:6
        push!(samegate,twoqubitgates[sameinx[i]])
    end
    return samegate
end

"""
Generate permutation tuple for each gate in H group
This will return a 2D array, each row is the permutation tuple for each gate in H group
1st idx is from 1 to 6, indicate which gate in H group
2nd idx is from 1 to 2^n, indicate the input state
the element in these two idx is the output state.
"""
function generate_Hperm(n,gate,states) #arg n is the number of qubits, gate is the gate in H group
    idx_homo_gate=[Vector{Int}() for _ in 1:6]
    total_state=2^n
    for k in 1:6    #samegate loop
        for i in 1:total_state
            for j in 1:total_state
                temp_state=states[i] ⊗ states[j]
                for t in 1:n
                    apply!(temp_state, gate[k], [t,t+n])
                end
                canonicalize!(temp_state)
                #restore the original order of phase
                phase=adjust_order(temp_state.tab.phases,n)
                #find the matching index in phase_perm
                idx=(bit_to_int(phase)+1)÷2 #convention check
                push!(idx_homo_gate[k],idx)
            end
        end
    end
    return idx_homo_gate    #permutaion index in the order of [SWAP,CNOT12,INV_CNOT21,INV_CNOT12,IDENTITY,CNOT21]
end

"""
Generate all gates for F group
There are eight gates in F group which are [Phase1,CZ+Phase1, Identity, CZ, Phase2, CZ+Phase2, BothPhase, CZ+BothPhase]
The idx of these gate already be found in the twoqubitgates, which are [20,88,121,189,241,309,620,688]
This function return all these gates in the vector form in the order of [Phase1,CZ+Phase1, Identity, CZ, Phase2, CZ+Phase2, BothPhase, CZ+BothPhase]
"""
function generate_Fgroup() #generate F group
    bi_inx=[20,88,121,189,241,309,620,688] # in the order of [Phase1,CZ+Phase1, Identity, CZ, Phase2, CZ+Phase2, BothPhase, CZ+BothPhase]
    bi_gate=[]
    for i in 1:8
        push!(bi_gate,twoqubitgates[bi_inx[i]])
    end
    return bi_gate
end

"""
generate permutation tuple each gate in F group and every nodes combinations
This will turn a 3D array, need 3 idx to specifies
1st idx is from 1 to n-1, indicates which nodes it acting on(Alice, Bob...)
2nd idx is from 1 to 8, indicates which gate in F group
3rd idx is from 1 to 2^n, indicates the input state
the element in these three idx is the output state.
"""
function generate_Fperm(n,bigate,states) # generate F group permutation for n qubit state
    idx_phase_gate = [] #store all permutaion
    total_state=2^n

    for t in 1:n-1      #for n qubit, there is n-1 F group application
        current_perm = [Int[] for _ in 1:8]

        for k in 1:8    #bigate loop
            for i in 1:total_state
                for j in 1:total_state
                    temp_state=states[i] ⊗ states[j]
                    apply!(temp_state, bigate[k], [t,t+n])
                    apply!(temp_state, bigate[k], [t+1,t+n+1])
                    canonicalize!(temp_state)
                    #restore the original order of phase
                    phase=adjust_order(temp_state.tab.phases,n)
                    #find the matching index in phase_perm
                    idx=(bit_to_int(phase)+1)÷2 #convention check
                    push!(current_perm[k],idx)
                end
            end
        end
        
        push!(idx_phase_gate,current_perm) # in the order of [Phase1,CZ+Phase1, Identity, CZ, Phase2, CZ+Phase2, BothPhase, CZ+BothPhase]
    end
    return idx_phase_gate # in the order of AB, BC, CD...
end

"""
generate permutation tuple for each gate in Pauli group
This will return a 2D array, each row is the permutation tuple for each gate in Pauli group
1st idx is from 1 to 4, indicate which gate in Pauli group, in the order of[X,Y,Z,I]
2nd idx is from 1 to 2^n, indicate the input state
the element in these two idx is the output state.
"""
function generate_Pperm(n,states)
    total_state=2^n
    pauligate=[P"X",P"Y",P"Z",P"I"]
    idx_pauli_gate = [[[] for _ in 1:total_state] for _ in 1:4]
    for k in 1:4
        for i in 1:total_state
            for j in 1:n
                temp_state=copy(states[i])
                apply!(temp_state, pauligate[k], [j])
                canonicalize!(temp_state)
                #find the matching index
                for (idx, state) in enumerate(states)
                    if temp_state.tab.phases == state.tab.phases
                        push!(idx_pauli_gate[k][i], idx)
                        break
                    end
                end
            end
        end
    end
    return idx_pauli_gate
end

#######################
#Initialization
#######################
"""
Initialization
"""
# get all two-qubit gates without phases
const twoqubitgates = collect(enumerate_cliffords(2))
# 6 gate in H group
const hgroup_gate = generate_Hgroup()
# 8 gate in F group
const fgroup_gate = generate_Fgroup()
#Dict for state cache
const State_Cache = Dict{Int, Vector}()
#Dict for H group permutation cache
const HPerm_Cache = Dict{Int, Vector}()
#Dict for F group permutation cache
const FPerm_Cache = Dict{Int, Vector}()
#Dict for Pauli group permutation cache
const PPerm_Cache = Dict{Int, Vector}()


#########################
#Main Struct
#########################

abstract type GHZOp <: QuantumClifford.AbstractCliffordOperator end

"""
For n qubit GHZ state, there are 2^n states, and it is represented by a BitVector of length n.
If the phase is "+" then it is 0, if the phase is "-" then it is 1.
"""
struct GHZState <: QuantumClifford.AbstractStabilizer
    qubit_num::Int  #number of qubits in one ghz state
    ghz_num::Int    #number of ghz states
    phases::BitVector
end
GHZState(n::Integer, m::Integer)=GHZState(n,m,BitVector(falses(n*m)))
GHZState(t::Tuple)=GHZState(BitVector(t))

Base.copy(state::GHZState) = GHZState(state.qubit_num, state.ghz_num, copy(state.phases))
Base.:(==)(l::GHZState, r::GHZState) = l.phases == r.phases

"""
H group
GHZ preserving gate contain 6 gates which can apply homogeneously to all qubits
1st argument, gate index: 1-6 in the order of [Swap,Cnot12,Inv_Cnot21,Inv_Cnot12,Identity,Cnot21]
2nd argument, ghz index: 1-n
3rd argument, ghz index: 1-n
"""
struct Hgroup{N} <: GHZOp
    gate_idx::Int
    ghz_idx1::Int
    ghz_idx2::Int
    perm::Any
end

function Hgroup{N}(g,q1,q2) where N
    if !haskey(State_Cache, N)
        State_Cache[N] = generate_state(N)
    end
    state = State_Cache[N]

    if !haskey(HPerm_Cache, N)
        HPerm_Cache[N]=generate_Hperm(N,hgroup_gate,state)
    end
    perm = HPerm_Cache[N]

    1 <= g <= 6 || throw(ArgumentError("Invalid gate index"))
    #1 <= q1 <= N || throw(ArgumentError("Invalid qubit index"))
    #1 <= q2 <= N || throw(ArgumentError("Invalid qubit index"))
    return Hgroup{N}(g,q1,q2,perm)
end

"""
F group
GHZ preserving gate contain 8 gates which can apply bi-locally to all qubits
1st argument, gate index: 1-8 in the order of [Phase1,CZ+Phase1, Identity, CZ, Phase2, CZ+Phase2, BothPhase, CZ+BothPhase]
2nd argument, ghz index: 1 - n
3rd argument, ghz index: 1 - n
4th argument, node index:1 - (n-1)
"""
struct Fgroup{N} <: GHZOp
    gate_idx::Int
    ghz_idx1::Int
    ghz_idx2::Int
    node_idx::Int
    perm::Any
end

function Fgroup{N}(g,q1,q2,node) where N
    if !haskey(State_Cache, N)
        State_Cache[N] = generate_state(N)
    end
    state = State_Cache[N]

    if !haskey(FPerm_Cache, N)
        FPerm_Cache[N]=generate_Fperm(N,fgroup_gate,state)
    end
    perm = FPerm_Cache[N]

    1 <= g <= 8 || throw(ArgumentError("Invalid gate index"))
    #1 <= q1 <= N || throw(ArgumentError("Invalid qubit index"))
    #1 <= q2 <= N || throw(ArgumentError("Invalid qubit index"))
    1 <= node <= (N-1) || throw(ArgumentError("Invalid node index"))
    return Fgroup{N}(g,q1,q2,node,perm)
end

"""
Pauli group
all 4 Pauli gates can apply as they only change the phases.
1st argument, gate_idx: 1 for X, 2 for Y, 3 for Z, 4 for I
2nd argument, ghz_idx: the state idx acting on
3nd argument, qubit_idx: qubit acting on 1-n
"""
struct PauliGroup{N} <: GHZOp
    gate_idx::Int
    ghz_idx::Int
    qubit_idx::Int
    perm::Any
end

function PauliGroup{N}(p,s,q) where N
    if !haskey(State_Cache, N)
        State_Cache[N] = generate_state(N)
    end
    state = State_Cache[N]

    if !haskey(PPerm_Cache, N)
        PPerm_Cache[N]=generate_Pperm(N,state)
    end
    perm = PPerm_Cache[N]

    1 <= p <= 4 || throw(ArgumentError("The permutation index needs to be between 1 and 4"))
    s>0 || throw(ArgumentError("The Bell pair indices have to be positive integers."))
    return PauliGroup{N}(p,s,q,perm)
end

"""
define "apply!" operation for H group
"""
function QuantumClifford.apply!(s::GHZState, op::Hgroup)
    n=s.qubit_num
    #concatenate the bit arrays from ghz_idx1 and ghz_idx2
    start_idx1=(op.ghz_idx1-1)*n+1
    end_idx1=op.ghz_idx1*n
    start_idx2=(op.ghz_idx2-1)*n+1
    end_idx2=op.ghz_idx2*n
    combined_bits=vcat(s.phases[start_idx1:end_idx1],s.phases[start_idx2:end_idx2])
    phase_idx=bit_to_int(combined_bits)

    #permutation
    new_state_idx=op.perm[op.gate_idx][phase_idx]
    new_phase = int_to_bit(new_state_idx,2n)

    #assign new phase to the state
    s.phases[start_idx1:end_idx1] .= new_phase[1:n]
    s.phases[start_idx2:end_idx2] .= new_phase[(n+1):end]
    return s
end

"""
define "apply!" operation for F group
"""
function QuantumClifford.apply!(s::GHZState, op::Fgroup)
    n=s.qubit_num
    #concatenate the bit arrays from ghz_idx1 and ghz_idx2
    start_idx1=(op.ghz_idx1-1)*n+1
    end_idx1=op.ghz_idx1*n
    start_idx2=(op.ghz_idx2-1)*n+1
    end_idx2=op.ghz_idx2*n
    combined_bits=vcat(s.phases[start_idx1:end_idx1],s.phases[start_idx2:end_idx2])
    phase_idx=bit_to_int(combined_bits)

    #permutation
    new_state_idx=op.perm[op.node_idx][op.gate_idx][phase_idx]
    new_phase = int_to_bit(new_state_idx,2n)

    #assign new phases to the state
    s.phases[start_idx1:end_idx1] .= new_phase[1:n]
    s.phases[start_idx2:end_idx2] .= new_phase[(n+1):end]
    return s
end

"""
define "apply!" operation for Pauli group
"""
function QuantumClifford.apply!(s::GHZState, op::PauliGroup)
    n=s.qubit_num
    start_idx=(op.ghz_idx-1)*n+1
    end_idx=op.ghz_idx*n
    phase_idx=bit_to_int(s.phases[start_idx:end_idx])

    #permutation
    new_state_idx=op.perm[op.gate_idx][phase_idx][op.qubit_idx]
    new_phase = int_to_bit(new_state_idx,n)

    #assign new phases to the state
    s.phases[start_idx:end_idx] .= new_phase
    return s
end

"""
define "*" operation for GHZOp
"""
function Base.:(*)(op::GHZOp, s::GHZState; phases::Bool=true)
    s = copy(s)
    apply!(s,op)
end

"""
define tensor! operation for GHZ state
"""
function tensor!(s1::GHZState, s2::GHZState)
    q1 = s1.qubit_num
    q2 = s2.qubit_num
    
    if q1 != q2
        throw(ArgumentError("Qubit num should be same to preserve"))
    end

    n1 = s1.ghz_num
    n2 = s2.ghz_num

    return GHZState(q1, n1+n2, vcat(s1.phases, s2.phases))
end

################################
#Measurement
################################

"""
Measurement on the GHZ state
1st argument, basis_idx: 1 for X, 2 for Y, 3 for Z, Y is not useful as there is no simple Y measurement in GHZ state
2nd argument, ghz_idx: the index of ghz state to be measured
"""
struct GHZMeasure <: QuantumClifford.AbstractMeasurement
    n::Int
    basis_idx::Int
    ghz_idx::Int
end

function GHZMeasure(n,b,si)
    #b as basis index, si as state index, n as number of qubits
    1 <= b <= 3 || throw(ArgumentError("Invalid basis index"))
    1 <= si <= n || throw(ArgumentError("Invalid state index"))
    return GHZMeasure(n,b,si)
end

"""
define "Measure!" operation for GHZ state
it will return the post measurement state and a probabilistic result
"""
function Measure!(s::GHZState, op::GHZMeasure)
    n=s.qubit_num
    start_idx=(op.ghz_idx-1)*n+1
    end_idx=op.ghz_idx*n

    if op.basis_idx==3  #Z measurement
        measure_result=[s.phases[i] for i in start_idx+1:end_idx] #Z basis measurement result should follow the phase of Z stabilizer
        s.phases[start_idx+1:end_idx] .= 0 #reset the state to 0
    elseif op.basis_idx==1  #X measurement
        #the result of X basis should be the same as Z basis except +,- instead of 0,1
        measure_result=[s.phases[start_idx]] #X basis measurement result should be the phase of 1st stabilizer
        s.phases[start_idx] = 0 #reset the state to 0
    #there is no simple Y measurement in GHZ state.
    end
    return s, measure_result
end

#################################
#Full GHZ preserving gate
#################################
"""
Most general representation of a GHZ preserving gate on n qubits two GHZ states
The general gate consists of:
    number of qubits in the GHZ state,
    one H group gate, 
    a series of F group gates, the length should be equal to n-1,
    a series of Pauli group gates, the length should be equal to n,
    two ghz state idx.
"""
struct GHZGate <: GHZOp
    N::Int                      #number of qubits
    H::Int                      #H group index
    F::Array{Int,1}             #F group index
    Paulis::Array{Int,1}        #Pauli group index
    ghz_idx1::Int               #ghz state index1
    ghz_idx2::Int               #ghz state index2
    function GHZGate(n,h,f,p,i1,i2)
        (1 <= h <= 6) || throw(ArgumentError("Invalid H group index, it should be between 1 and 6"))
        (length(f) == n-1) || throw(ArgumentError("F group series length error, it should be n-1"))
        all(1 <= fi <= 8 for fi in f) || throw(ArgumentError("Invalid F group index, it should be between 1 and 8"))
        (length(p) == n) || throw(ArgumentError("Pauli group series length error, it should be n"))
        all(1 <= pi <= 4 for pi in p) || throw(ArgumentError("Invalid Pauli group index, it should be between 1 and 4"))
        (i1 > 0 && i2 > 0) || throw(ArgumentError("GHZ state index should be positive integers"))
        i1 != i2 || throw(ArgumentError("GHZ state index should be different"))
        new(n,h,f,p,i1,i2)
    end
end

function QuantumClifford.apply!(s::GHZState, g::GHZGate)
    n=s.qubit_num
    m=s.ghz_num
    
    #apply H group
    apply!(s,Hgroup(g.H,g.ghz_idx1,g.ghz_idx2))
    #apply F group
    for i in 1:n-1
        apply!(s,Fgroup(g.F[i],g.ghz_idx1,g.ghz_idx2,i))
    end
    #apply Pauli group
    for i in 1:n
        apply!(s,PauliGroup(g.Paulis[i],g.ghz_idx1,i))
    end
    return s
end

function QuantumClifford._apply!(s::GHZState, g::GHZGate)
    return QuantumClifford.apply!(s, g)
end

#region good operations

#########################
#Typically good operations
#########################

"""
CNOT as a good operation
"""
struct CNOT <: GHZOp
    n::Int
    ghz_idx1::Int
    ghz_idx2::Int
    function CNOT12(n,i1,i2)
        (i1 > 0 && i2 > 0) || throw(ArgumentError("GHZ state index should be positive integers"))
        (1<=i1<=n && 1<=i2<=n) || throw(ArgumentError("Invalid GHZ state index"))
        new(n,i1,i2)
    end
end

function QuantumClifford.apply!(s::GHZState, op::CNOT)
    n=s.qubit_num
    m=s.ghz_num
    #apply H group
    apply!(s,Hgroup(2,op.ghz_idx1,op.ghz_idx2))
    return s
end

"""
CZ as a good operation
"""
struct CZ <: GHZOp
    n::Int
    ghz_idx1::Int
    ghz_idx2::Int
    node::Int
    function CZ12(n,i1,i2,node)
        (i1 > 0 && i2 > 0) || throw(ArgumentError("GHZ state index should be positive integers"))
        (1<=i1<=n && 1<=i2<=n) || throw(ArgumentError("Invalid GHZ state index"))
        (1<=node<=n-1) || throw(ArgumentError("Invalid node index"))
        new(n,i1,i2,node)
    end
end

function QuantumClifford.apply!(s::GHZState, op::CZ)
    n=s.qubit_num
    m=s.ghz_num
    #apply H group
    apply!(s,Fgroup(4,op.ghz_idx1,op.ghz_idx2,op.node))
    return s
end

"""
Phase gate as a good operation
"""
struct F1 <: GHZOp
    n::Int
    ghz_idx1::Int
    ghz_idx2::Int
    node::Int
    function F1(n,i1,i2,node)
        (i1 > 0 && i2 > 0) || throw(ArgumentError("GHZ state index should be positive integers"))
        (1<=i1<=n && 1<=i2<=n) || throw(ArgumentError("Invalid GHZ state index"))
        (1<=node<=n-1) || throw(ArgumentError("Invalid node index"))
        new(n,i1,i2,node)
    end
end

function QuantumClifford.apply!(s::GHZState, op::F1)
    n=s.qubit_num
    m=s.ghz_num
    #apply H group
    apply!(s,Fgroup(1,op.ghz_idx1,op.ghz_idx2,op.node))
    return s
end

#endregion

#region Noisy

############
#Noisy
############

"""
Pauli noise with probability px, py, pz
"""
struct PauliNoiseOp <: GHZOp
    qubit_num::Int
    ghz_idx::Int
    px::Float64
    py::Float64
    pz::Float64
end

function QuantumClifford.apply!(s::GHZState, op::PauliNoiseOp)
    n=op.qubit_num
    i=op.ghz_idx

    for q in 1:n
        r=rand()
        if r<op.px
            apply!(s,PauliGroup{n}(1,i,q))
        elseif r<op.px+op.py
            apply!(s,PauliGroup{n}(2,i,q))
        elseif r<op.px+op.py+op.pz
            apply!(s,PauliGroup{n}(3,i,q))
        end
    end
    return s
end

"""
Depolarize error
"""
function depolarize!(s::GHZState, op::Hgroup, p::Float64,)
    n=s.qubit_num
    i1=op.ghz_idx1
    i2=op.ghz_idx2
    
    for q in 1:n
        if rand()>p
            case1 = rand(1:4)
            apply!(s,PauliGroup{n}(case1,i1,q))
            case2 = rand(1:4)
            apply!(s,PauliGroup{n}(case2,i2,q))
        end
    end
   
end

function depolarize!(s::GHZState, op::Fgroup, p::Float64)
    n=s.qubit_num
    i1=op.ghz_idx1
    i2=op.ghz_idx2
    node=op.node_idx

    if rand()>p
        case1 = rand(1:4)
        apply!(s,PauliGroup{n}(case1,i1,node))
        case2 = rand(1:4)
        apply!(s,PauliGroup{n}(case2,i2,node+1))
    end

end

"""
Noisy measurement
"""
struct NoisyMeasure <: GHZOp
    m::GHZMeasure
    p::Float64
end

"""
Common noisy measurement with probability p giving wrong result
"""
function QuantumClifford.apply!(s::GHZState, op::NoisyMeasure)
    state, result=Measure!(s, op.m)
    n=s.qubit_num

    if op.m.basis_idx==3  #Z measurement
        #error in coincidence measurement.
        r=Bernoulli(op.p)
        measurement_errors = [rand(r) for _ in 1:n]
        error_phase = [i ⊻ j for (i, j) in zip(
            measurement_errors[begin:end-1], measurement_errors[2:end])]

        result .= result .⊻ error_phase
        return state, result
    else #X measurement
        if rand() > op.p
            result .= result .⊻ 1
        end
        return state, result
    end
end


"""
Noisy measurement and Pauli noise after the reset
"""
struct NoisyMeasureNoisyReset
    m::GHZMeasure
    p::Float64
    f_in::Float64
end

#=
function QuantumClifford.applywstatus!(s::GHZState, op::NoisyMeasure)
    state, result=Measure!(s, op.m)
    original_result = copy(result)

    #each qubit have independent filp error with probability p
    result .= result .⊻ (rand(length(result)) .< op.p)
    
    if any(original_result .!= result)
        return state, failure_stat
    else
        return state, continue_stat
    end
end
=#

function QuantumClifford.applywstatus!(s::GHZState, op::NoisyMeasure, f_in::Float64=1.0)
    state, result=Measure!(s, op.m)
    n=s.qubit_num

    if op.m.basis_idx==3  #Z measurement
        #error in coincidence measurement.
        r=Bernoulli(op.p)
        measurement_errors = [rand(r) for _ in 1:n]
        error_phase = [i ⊻ j for (i, j) in zip(
            measurement_errors[begin:end-1], measurement_errors[2:end])]

        result .= result .⊻ error_phase
        
        if result == zeros(Int, n-1)
            #if success, insert a new state
            new_state = rand(GHZState, n, 1, f_in)
            idx=op.m.ghz_idx
            state.phases[(idx-1)*n+1:idx*n] = new_state.phases

            return state, continue_stat
        else
            return state, failure_stat
        end

    else #X measurement
        if rand() > op.p
            result .= result .⊻ 1
        end

        if result == zeros(Int, 1)
            #if success, insert a new state
            new_state = rand(GHZState, n, 1, f_in)
            idx=op.m.ghz_idx
            state.phases[(idx-1)*n+1:idx*n] = new_state.phases

            return state, continue_stat
        else
            return state, failure_stat
        end
    end
end

function QuantumClifford.applywstatus!(s::GHZState, op::NoisyMeasureNoisyReset)
    state, result = Measure!(s, op.m)
    
    #flip error
    flipped =  result .⊻ (rand(length(result)) .< op.p)

    #if flipped apply Pauli noise
    if any(flipped)
        apply!(state, PauliNoiseOp(s.qubit_num,op.m.ghz_idx,op.px,op.py,op.pz))
    end

    state, any(flipped) ? failure_stat : continue_stat 
end

#endregion

#region Random

############
#Random
############

"""
Random GHZ diagonal state, with n qubits and m states
"""
function Random.rand(::Type{GHZState}, n::Int, m::Int)
    return GHZState(n,m,BitVector(rand(Bool,n*m)))
end

"""
Generate random GHZ state with n qubits and m states, with fidelity p
"""
function Random.rand(::Type{GHZState}, n::Int, m::Int, p::Float64)
    phases=BitVector()
    for i in 1:m
        if rand() < p
            #generate target state, all 0
            append!(phases, BitVector(falses(n)))
        else
            #generate random state except all 0
            random_bits = BitVector(rand(Bool,n))
            while all(!x for x in random_bits)  #avoid all 0 state
                random_bits = BitVector(rand(Bool,n))
            end
            append!(phases, random_bits)
        end
    end
    return GHZState(n,m,phases)
end

"""
Random GHZGate
"""
function Random.rand(::Type{GHZGate}, n::Int)
    h=rand(1:6)
    f=rand(1:8,n-1)
    p=rand(1:4,n)
    ghz1=rand(1:n)
    ghz2=rand(1:n)
    while ghz1==ghz2
        ghz2=rand(1:n)
    end
    return GHZGate(n,h,f,p,ghz1,ghz2)
end

"""
Random CNOT on ghz state i and j
"""
function Random.rand(::Type{CNOT}, n::Int)
    i=rand(1:n)
    j=rand(1:n)
    while i==j
        j=rand(1:n)
    end
    return CNOT(n,i,j)
end

"""
Random CZ on ghz state i and j
"""
function Random.rand(::Type{CZ}, n::Int)
    i=rand(1:n)
    j=rand(1:n)
    while i==j
        j=rand(1:n)
    end
    node=rand(1:n-1)
    return CZ(n,i,j,node)
end

"""
Random F1 on ghz state i and j
"""
function Random.rand(::Type{F1}, n::Int)
    i=rand(1:n)
    j=rand(1:n)
    while i==j
        j=rand(1:n)
    end
    node=rand(1:n-1)
    return F1(n,i,j,node)
end

"""
Random measurement on ghz state i
"""
function Random.rand(::Type{GHZMeasure}, n::Int)
    b=rand([1,3])   #no 2 as no Y for now.
    i=rand(1:n)
    return GHZMeasure(n,b,i)
end

#endregion

#region to QC

#################################
#Convension from ghz preserving to QC
#################################

"""
Convert GHZ preserving gate to QC gate
"""
function toQCcircuit end

function toQCcircuit(g::Hgroup)
    return hgroup_gate[g.gate_idx]
end

function toQCcircuit(g::Fgroup)
    return fgroup_gate[g.gate_idx]

end

const pg_qc = (sId1,sX,sZ,sY)

function toQCcircuit(g::PauliGroup)
    return pg_qc[g.gate_idx]
end

function toQCcircuit(g::GHZMeasure)
    m=(sMX, sMY, sMZ)[g.basis_idx]
end

function QuantumClifford.Stabilizer(s::GHZState)
    # make sure the state is in the cache
    if !haskey(State_Cache, s.qubit_num)
        State_Cache[s.qubit_num] = generate_state(s.qubit_num)
    end

    states=State_Cache[s.qubit_num]

    result_states = []
    for i in 1:s.ghz_num
        start_idx=(i-1)*s.qubit_num+1
        end_idx=i*s.qubit_num
        push!(result_states, states[bit_to_int(s.phases[start_idx:end_idx])])
    end

    combined_state=copy(result_states[1])
    for i in 2:length(result_states)
        combined_state=combined_state ⊗ result_states[i]
    end

    return combined_state
end

function QuantumClifford.MixedDestabilizer(s::GHZState)
    MixedDestabilizer(Stabilizer(s))
end
#endregion

end 