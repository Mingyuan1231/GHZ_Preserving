using QuantumClifford
using LinearAlgebra
using Random

# get all two-qubit gates without phases
twoqubitgates = collect(enumerate_cliffords(2))

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
1st idx is from 1 to 4, indicate which gate in Pauli group, in the order of[I,X,Y,Z]
2nd idx is from 1 to 2^n, indicate the input state
the element in these two idx is the output state.
"""
function generate_Pauli_perm(n,states)
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

#######################
#Initialization
#######################

struct Perm{N}
    n::Int
    h::Any
    f::Any
    p::Any

    function Perm{N}() where N
        states = generate_state(N)
        h_group = generate_Hgroup()
        f_group = generate_Fgroup()
        h_perm = generate_Hperm(N, h_group, states)
        f_perm = generate_Fperm(N, f_group, states)
        p_perm = generate_Pauli_perm(N, states)
        new(N, h_perm, f_perm, p_perm)
    end
end



###################
#test
###################

n=3
perm=Perm{n}()
state=generate_state(n)
hgroup=generate_Hgroup()
fgroup=generate_Fgroup()
hperm=generate_Hperm(n,hgroup,state)
fperm=generate_Fperm(n,fgroup,state)
pperm=generate_Pauli_perm(n,state)


check=state[3] ⊗ state[4]
check.tab.phases
adjust_order(check.tab.phases,n)

for i in 1:8
    temp=state[i] ⊗ state[i]
    println(temp.tab.phases)
    println((bit_to_int(temp.tab.phases)+1)÷2)
end

GHZState(b::BitVector)=GHZState(3,2,b)

s=GHZState((0,1,0,1,0,0))

#checking H group
swap=Hgroup(1,1,2)
apply!(s,swap)


cnot12=Hgroup(2,1,2)
apply!(s,cnot12)

inv_cnot21=Hgroup(3,1,2)
int_cnt12=Hgroup(4,1,2)
identity=Hgroup(5,1,2)
cnot21=Hgroup(6,2,1)


n=s.qubit_num
op=identity
#concatenate the bit arrays from ghz_idx1 and ghz_idx2
start_idx1=(op.ghz_idx1-1)*n+1
end_idx1=op.ghz_idx1*n
start_idx2=(op.ghz_idx2-1)*n+1
end_idx2=op.ghz_idx2*n
combined_bits=vcat(s.phases[start_idx1:end_idx1],s.phases[start_idx2:end_idx2])
phase_idx=bit_to_int(combined_bits)

#permutation
new_state_idx=perm.h[op.gate_idx][phase_idx]
new_phase = int_to_bit(new_state_idx,2n)

#assign new phase to the state
s.phases[start_idx1:end_idx1] .= new_phase[1:n]
s.phases[start_idx2:end_idx2] .= new_phase[(n+1):end]

total_state=8
states=generate_state(3)
idx_homo_gate=[]
for i in 1:total_state
    for j in 1:total_state
        temp_state=states[i] ⊗ states[j]
        for t in 1:n
            apply!(temp_state, twoqubitgates[121], [t,t+n])
        end
        canonicalize!(temp_state)
        #find the matching index in phase_perm
        idx=(bit_to_int(temp_state.tab.phases)+1)÷2 #convention check
        push!(idx_homo_gate,idx)
    end
end

st=state[1] ⊗ states[5]

apply!(st,twoqubitgates[121],[1,4])
apply!(st,twoqubitgates[121],[2,5])
apply!(st,twoqubitgates[121],[3,6])
canonicalize!(st)
idx=(bit_to_int(st.tab.phases)+1)÷2

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
        1 <= q1 <= perm.n || throw(ArgumentError("Invalid qubit index"))
        1 <= q2 <= perm.n || throw(ArgumentError("Invalid qubit index"))
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
        1 <= q1 <= perm.n || throw(ArgumentError("Invalid qubit index"))
        1 <= q2 <= perm.n || throw(ArgumentError("Invalid qubit index"))
        1 <= node <= (n-1) || throw(ArgumentError("Invalid node index"))
        new(g,q1,q2,node)
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
    n=s.qubit_num
    #concatenate the bit arrays from ghz_idx1 and ghz_idx2
    start_idx1=(op.ghz_idx1-1)*n+1
    end_idx1=op.ghz_idx1*n
    start_idx2=(op.ghz_idx2-1)*n+1
    end_idx2=op.ghz_idx2*n
    combined_bits=vcat(s.phases[start_idx1:end_idx1],s.phases[start_idx2:end_idx2])
    phase_idx=bit_to_int(combined_bits)

    #permutation
    new_state_idx=perm.h[op.gate_idx][phase_idx]
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
    new_state_idx=perm.f[op.node_idx][op.gate_idx][phase_idx]
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
    new_state_idx=perm.p[op.gate_idx][phase_idx]
    new_phase = int_to_bit(new_state_idx,2n)

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

################################
#Measurement
################################

"""
Measurement on the GHZ state
1st argument, basis_idx: 1 for X, 2 for Y, 3 for Z
2nd argument, ghz_idx: the index of ghz state to be measured
"""
struct GHZMeasure <: QuantumClifford.AbstractMeasurement
    basis_idx::Int
    ghz_idx::Int
    function GHZMeasure(b,si)  #b as basis index, si as state index, n as number of qubits
        1 <= b <= 3 || throw(ArgumentError("Invalid basis index"))
        1 <= si <= perm.n || throw(ArgumentError("Invalid state index"))
        new(b,si)
    end
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
        s.phases[start_idx] .= 0 #reset the state to 0
    #TODO Y measurement
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


############
#Noisy
############

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

#=
function QuantumClifford.apply!(s::GHZState, op::NoisyMeasure)
    state, result=Measure!(s, op.m)
    state, result⊻(rand()<op.p) ? continue_stat : failure_stat
end
=#

function QuantumClifford.apply!(s::GHZState, op::NoisyMeasure)
    state, result=Measure!(s, op.m)
    if rand()<op.p
        result .= result .⊻ true        #flip the result
    end
    return state, result
end

#=
function QuantumClifford.apply!(s::GHZState, op::NoisyMeasureNoisyReset)
    state, result = Measure!(s, op.m)
    cont = result⊻(rand()<op.p)
    cont && apply!(state, PauliNoiseOp(op.m.ghz_idx,op.px,op.py,op.pz))
    state, cont ? continue_stat : failure_stat
end
=#

#not sure here is necessary
#in the BP case, `PauliNoiseOp(idx,px,py,pz)` causes qubit-pair `idx` to flip to one of the other 3 Bell states with probabilities `px`, `py`, `pz` respectively.
#but here we will have 2^n-1 states which it can be filpped into, it is not very smart to assign the probability to each state.
#also, the probabilities of each state is different as some state is more likely to be flipped than others. 1 qubit error is more likely than 2 qubit error and so on.

function QuantumClifford.apply!(s::GHZState, op::NoisyMeasureNoisyReset)
    state, result = Measure!(s, op.m)
    if rand()<op.p
        result .= result .⊻ true        #flip the result
    end
    
    #apply Pauli noise
    if any(result)
        apply!(state, PauliNoiseOp(op.m.ghz_idx, op.px, op.py, op.pz))
    end

    #choose the next state
    continue_stat = result
    failure_stat = !result

    return state, any(result) ? continue_stat : failure_stat
end

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
Generate random GHZ state with n qubits and m states, with fidelity p, rho=p|0⟩⟨0|+(1-p)I/2^n, |0⟩⟨0| stands for all 0 state
"""
function Random.rand(::Type{GHZState}, n::Int, m::Int, p::Float64)
    phases=BitVector()
    for i in 1:m
        if rand() < p
            #generate target state, all 0
            append!(phases, BitVector(falses(n)))
        else
            #generate random state
            append!(phases, BitVector(rand(Bool,n)))
        end
    end
    return GHZState(n,m,phases)
end

############
#test
############

NM=NoisyMeasure(GHZMeasure(3,2),0.5)
apply!(teststate,NM)
teststate
NM.m
state=Measure!(teststate,NM.m)
state,result=apply!(teststate,NM)
state
result

teststate=GHZState(3,2)
testrandstate=rand(GHZState,3,2)

check=GHZState((0,0,0,1,1,0))
check.phases

Measure!(check,GHZMeasure(3,2))
Measure!(testrandstate,Z_measure)


for i in 1:100
    testrandstate=rand(GHZState,3,2,0.25)
    be=copy(testrandstate.phases)
    apply!(testrandstate,cnot12)
    println(be,testrandstate.phases)
end

teststate=generate_state(3)
for i in 1:8
    println(teststate[i].tab.phases)
    println(bit_to_int(teststate[i].tab.phases))
end

bit_to_int([0,0,0,0,0,0,0])
rand(GHZState,3,2,0.5)
Z_measure=GHZMeasure(3,2)

succ_dict=Dict{Tuple{Float64,Float64},Float64}()
fide_dict=Dict{Tuple{Float64,Float64},Float64}()
statecheck=[]
t=10000

for p2 in [0.9,0.95,0.99]
    NM=NoisyMeasure(Z_measure,1-p2)
    for fin in 0:0.05:1
        count=0
        fout=0
        for i in 1:t
            testrandstate=rand(GHZState,3,2,fin)
            apply!(testrandstate,cnot12)
            s,result=apply!(testrandstate,NM)
            if result ==[0,0]
                count+=1
                push!(statecheck,s)
                #count fidelity
                fout+=Int(s.phases[1:3] == [0,0,0])
            end
        end
        succ_dict[(p2,fin)]=count/t
        fide_dict[(p2,fin)]=fout/count
    end
end

using Plots

p2_values = unique([key[1] for key in keys(succ_dict)])

# Initialize the plots
p_success = plot(title="Success Rate vs. fin", xlabel="fin", ylabel="Success Rate", legend=:topleft,xlims=(0,1),ylims=(0,1),aspect_ratio=:equal)
p_fidelity = plot(title="Fidelity vs. fin", xlabel="fin", ylabel="Fidelity", legend=:topleft,xlims=(0,1),ylims=(0,1),aspect_ratio=:equal)

# Plotting for success_rate_dict
for p2 in p2_values
    fin_values = [key[2] for key in keys(succ_dict) if key[1] == p2]
    success_rates = [succ_dict[(p2, fin)] for fin in fin_values]
    
    sorted_indices = sortperm(fin_values)
    plot!(p_success, fin_values[sorted_indices], success_rates[sorted_indices], fillalpha=0.3, label="p2 = $p2")
end

# Plotting for fidelity_dict
for p2 in p2_values
    fin_values = [key[2] for key in keys(fide_dict) if key[1] == p2]
    fidelities = [fide_dict[(p2, fin)] for fin in fin_values]
    
    sorted_indices = sortperm(fin_values)
    plot!(p_fidelity, fin_values[sorted_indices], fidelities[sorted_indices], fillalpha=0.3, label="p2 = $p2")
end

p=plot(p_fidelity, p_success, layout=(1,2))
display(p)


using Makie
using CairoMakie

function makefig_succ()
    f_in=0:0.05:1
    fig=Figure()
    ax=Axis(fig[1,1],aspect=AxisAspect(1),xlabel="F_in",ylabel="succ rate")
    limits!(ax,0,1,0,1)
    ax.xticks=0:0.1:1
    ax.yticks=0:0.1:1
    lines!(ax,f_in,succ_out)
    return fig
end

function makefig_fide()
    f_in=0:0.05:1
    fig=Figure()
    ax=Axis(fig[1,1],aspect=AxisAspect(1),xlabel="F_in",ylabel="succ rate")
    limits!(ax,0,1,0,1)
    ax.xticks=0:0.1:1
    ax.yticks=0:0.1:1
    lines!(ax,f_in,fide_out)
    return fig
end

makefig_succ()
makefig_fide()

