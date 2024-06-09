module GHZPreserving

using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using LinearAlgebra
using Random


#convert bit(phase) to int
function bit_to_int(bits)
    reduce(⊻,(bit<<(index-1) for (index,bit) in enumerate(bits))) + 1 # +1 so that we use julia indexing convenctions
end

#convert int to bit(phase)
function int_to_bit(int,digits)     #digits should be 2n, as dealing with 2 input state, n qubit each
    int = int - 1 # -1 so that we use julia indexing conventions
    Bool[int>>shift&0x1 for shift in 0:digits-1]
end

###############
#Initialization
#Generating all the permutations for given n qubit ghz state.
###############

# get all two-qubit gates without phases
twoqubitgates = collect(enumerate_cliffords(2))

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
function generate_Hperm(n,gate) #arg n is the number of qubits, gate is the gate in H group
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
                #find the matching index in phase_perm
                idx=bit_to_int(temp_state.tab.phases)
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
function generate_Fperm(n,bigate) # generate F group permutation for n qubit state
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
                    #find the matching index in phase_perm
                    idx=bit_to_int(temp_state.tab.phases)
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
1st idx is from 1 to 4, indicate which gate in Pauli group, in the order of[I,X,Y,Z]
2nd idx is from 1 to 2^n, indicate the input state
the element in these two idx is the output state.
"""
function generate_Pauli_perm(n)
    total_state=2^n
    pauligate=[P"I",P"X",P"Y",P"Z"]
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

#input number of qubits, user should enter the value here.
n=3
states=generate_state(n)
h_group=generate_Hgroup()
f_group=generate_Fgroup()
h_perm=generate_Hperm(n,h_group)
f_perm=generate_Fperm(n,f_group)
p_perm=generate_Pauli_perm(n)

################
#Struct
################

abstract type GHZOp <: QuantumClifford.AbstractCliffordOperator end

"""
For n qubit GHZ state, there are 2^n states, and it is represented by a BitVector of length n.
If the phase is "+" then it is 0, if the phase is "-" then it is 1.
"""
struct GHZState <: QuantumClifford.AbstractStabilizer
    phases::BitVector
end
GHZState(n::Integer)=GHZState(BitVector(falses(2n)))
GHZState(t::Tuple)=GHZState(BitVector(t))

"""
H group
1st argument, gate index: 1-6 in the order of [SWAP,CNOT12,INV_CNOT21,INV_CNOT12,IDENTITY,CNOT21]
2nd argument, ghz index: 1-n
3rd argument, ghz index: 1-n
"""
struct Hgroup <: GHZOp
    gate_idx::Int
    ghz_idx1::Int
    ghz_idx2::Int
    function Hgroup(g,q1,q2)
        1 <= g <= 6 || throw(ArgumentError("Invalid gate index"))
        1 <= q1 <= n || throw(ArgumentError("Invalid qubit index"))
        1 <= q2 <= n || throw(ArgumentError("Invalid qubit index"))
        new(g,q1,q2)
    end
end

"""
F group
1st argument, gate index: 1-8 in the order of [Phase1,CZ+Phase1, Identity, CZ, Phase2, CZ+Phase2, BothPhase, CZ+BothPhase]
2nd argument, ghz index: 1 - n
3rd argument, ghz index: 1 - n
4th argument, node index:1 - (n-1)
"""
struct Fgroup <: GHZOp
    gate_idx::Int
    ghz_idx1::Int
    ghz_idx2::Int
    node_idx::Int
    function Fgroup(g,q1,q2,node)
        1 <= g <= 8 || throw(ArgumentError("Invalid gate index"))
        1 <= q1 <= n || throw(ArgumentError("Invalid qubit index"))
        1 <= q2 <= n || throw(ArgumentError("Invalid qubit index"))
        1 <= node <= (n-1) || throw(ArgumentError("Invalid node index"))
        new(g,q1,q2)
    end
end

"""
Pauli group
1st argument, gate_idx: 1 for X, 2 for Y, 3 for Z, 4 for I
2nd argument, ghz_idx: the state idx acting on
3nd argument, qubit_idx: qubit acting on 1-n
"""
struct PauliGroup <: GHZOp
    gate_idx::Int
    ghz_idx::Int
    qubit_idx::Int
    function PauliPermutation(p,s)
        1 <= p <= 4 || throw(ArgumentError("The permutation index needs to be between 1 and 4"))
        s>0 || throw(ArgumentError("The Bell pair indices have to be positive integers."))
        new(p,s)
    end
end

"""
define "apply!" operation for H group
"""
function QuantumClifford.apply!(s::GHZState, op::Hgroup)
    #concatenate the bit arrays from ghz_idx1 and ghz_idx2
    start_idx1=(op.ghz_idx1-1)*n+1
    end_idx1=op.ghz_idx1*n
    start_idx2=(op.ghz_idx2-1)*n+1
    end_idx2=op.ghz_idx2*n
    combined_bits=vcat(s.phases[start_idx1:end_idx1],s.phases[start_idx2:end_idx2])
    phase_idx=bit_to_int(combined_bits)

    #permutation
    new_state_idx=h_perm[op.gate_idx][phase_idx]
    new_phase = int_to_bit(new_state_idx)

    #assign new phase to the state
    s.phases[start_idx1:end_idx1] .= new_phase[1:n]
    s.phases[start_idx2:end_idx2] .= new_phase[(n+1):end]
    return state
end

"""
define "apply!" operation for F group
"""
function QuantumClifford.apply!(s::GHZState, op::Fgroup)
    #concatenate the bit arrays from ghz_idx1 and ghz_idx2
    start_idx1=(op.ghz_idx1-1)*n+1
    end_idx1=op.ghz_idx1*n
    start_idx2=(op.ghz_idx2-1)*n-1
    end_idx2=op.ghz_idx2*n
    combined_bits=vcat(s.phases[start_idx1:end_idx1],s.phases[start_idx2:end_idx2])
    phase_idx=bit_to_int(combined_bits)

    #permutation
    new_state_idx=f_perm[op.node_idx][op.gate_idx][phase_idx]
    new_phase = int_to_bit(new_state_idx)

    #assign new phases to the state
    s.phases[start_idx1:end_idx1] .= new_phase[1:n]
    s.phases[start_idx2:end_idx2] .= new_phase[(n+1):end]
    return state
end

"""
define "apply!" operation for Pauli group
"""
function QuantumClifford.apply!(s::GHZState, op::PauliGroup)
    start_idx=(op.ghz_idx-1)*n+1
    end_idx=op.ghz_idx*n
    phase_idx=bit_to_int(s.phases[start_idx:end_idx])

    #permutation
    new_state_idx=p_perm[op.gate_idx][phase_idx]
    new_phase = int_to_bit(new_state_idx)

    #assign new phases to the state
    s.phases[start_idx:end_idx] .= new_phase
    return state
end

"""
define "*" operation for GHZOp
"""
function Base.:(*)(op::GHZOp, s::GHZState; phases::Bool=true)
    s = copy(s)
    apply!(s,op)
end

####################
#Measurement
####################

"""
Measurement on the GHZ state
1st argument, basis_idx: 1 for X, 2 for Y, 3 for Z
2nd argument, ghz_idx: the index of ghz state to be measured
"""
struct GHZMeasure <: QuantumClifford.AbstractMeasurement
    basis_idx::Int
    ghz_idx::Int
    function GHZMeasure(b,q)
        1 <= b <= 3 || throw(ArgumentError("Invalid basis index"))
        1 <= q <= n || throw(ArgumentError("Invalid qubit index"))
        new(b,q)
    end
end

"""
define "Measure!" operation for GHZ state
it will return the post measurement state and a probabilistic result
"""
function Measure!(s::GHZState, op::GHZMeasure)
    start_idx=(op.ghz_idx-1)*n+1
    end_idx=op.ghz_idx*n

    if op.basis_idx==3  #Z measurement
        measure_result=[rand() < 0.5 ? (s.phases[i] ? 1:0) : (s.phases[i] ? 0:1) for i in start_idx:end_idx] #randomly choose the measurement result
        s.phases[start_idx:end_idx] .= 0 #reset the state to 0
    elseif op.basis_idx==1  #X measurement
        #the result of X basis should be the same as Z basis except +,- instead of 0,1
        measure_result=[rand() < 0.5 ? (s.phases[i] ? 1:0) : (s.phases[i] ? 0:1) for i in start_idx:end_idx] #randomly choose the measurement result
        s.phases[start_idx:end_idx] .= 0 #reset the state to 0
    #TODO Y measurement
    end
    return state, measure_result
end

######################
#Full GHZ preserving gate
######################
#I don't think this part in BP gate can be adopted here, for BPGate, the number of qubit is fixed which is 2, while here the number of qubits is not defined before user enter it.

######################
#Noisy
######################

struct PauliNoiseOp <: GHZOp
    ghz_idx::Int
    qubit_idx::Int
    px::Float64
    py::Float64
    pz::Float64
end

function QuantumClifford.apply!(s::GHZState, op::PauliNoiseOp)
    i=op.ghz_idx
    q=op.qubit_idx
    r=rand()
    if r<g.px
        apply!(s,PauliGroup(1,i,q))
    elseif r<g.px+g.py
        apply!(s,PauliGroup(2,i,q))
    elseif r<g.px+g.py+g.pz
        apply!(s,PauliGroup(3,i,q))
    end
    return state
end

"""
Noisy measurement
"""
struct NoisyMeasure <: GHZOp
    m::GHZMeasure
    p::Float64
end

"""
Noisy measurement and Pauli noise after the reset
"""
struct NoisyMeasureNoisyReset
    m::GHZMeasure
    p::Float64
    px::Float64
    py::Float64
    pz::Float64
end

function QuantumClifford.apply!(s::GHZState, op::NoisyMeasure)
    state,result=Measure!(s,op.m)
    state, result⊻(rand()<op.p) ? continue_stat : failure_stat
end

function QuantumClifford.apply!(s::GHZState, op::NoisyMeasureNoisyReset)
    state, result = Measure!(s, op.m)
    cont = result⊻(rand()<op.p)
    cont && apply!(state, PauliNoiseOp(op.m.ghz_idx,op.px,op.py,op.pz))
    state, cont ? continue_stat : failure_stat
end

########################
#Random
########################

"""
Random GHZ diagonal state
"""
function Random.rand(::Type{GHZState}, n::Int)
    return GHZState(BitArray(rand(Bool,2n)))
end

#for rest of random, I don't think it is necessary, as it speficy "good perm", and all my gates are good perm, there is no need to pick it out.

########################
#TO QC
########################

#TODO, this part should be done after verify above content is correct.