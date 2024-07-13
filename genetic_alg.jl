using .GHZPreserving
using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using Random
using Statistics

mutable struct Individual
    history::String
    qubit_num::Int
    ghz_num::Int
    f_in::Float64
    ops::Vector{Union{Hgroup,Fgroup}} #modify?
    success::Float64
    f_out::Float64
end 

function drop_op(indiv::Individual) 
    new_indiv = deepcopy(indiv)
    deleteat!(new_indiv.ops, rand(1:length(new_indiv.ops)))
    new_indiv.history = "drop_m"
    return new_indiv
end

#TODO adding noise
function gain_op(indiv::Individual)
    new_indiv = deepcopy(indiv)
    n=indiv.qubit_num
    idx=indiv.ghz_num
    if rand() < 0.0 #no measurement adding for now
        #chance to measure
        #rand[1,3] as X or Z measure, no Y for now, measure random ghz state
        measure=GHZMeasure(n,3,rand(1:idx))  #3 tobe edit later as measure basis, Z for now
        gate=NoisyMeasre(measure,1)    #1 tobe edit later as noise level, no noise for now
    elseif rand() < 0.5
        #chance to add H gate to two random ghz state
        perm=randperm(idx)[1:2]
        gate = Hgroup{n}(rand(1:6),perm[1],perm[2])
    else
        #chance to add F gate to two random ghz state at random node.
        perm=randperm(idx)[1:2]
        gate = Fgroup{n}(rand(1:8),perm[1],perm[2],rand(1:idx-1))
    end

    if length(new_indiv.ops) == 0
        push!(new_indiv.ops, gate)
    else
        insert!(new_indiv.ops, rand(1:length(new_indiv.ops)), gate)
    end
    new_indiv.history = "gain_m"
    return new_indiv
end

function swap_op(indiv::Individual)
    new_indiv = deepcopy(indiv)
    ind1, ind2 = rand(1:length(new_indiv.ops)), rand(1:length(new_indiv.ops))
    op1, op2 = new_indiv.ops[ind1], new_indiv.ops[ind2]
    new_indiv.ops[ind1] = op2
    new_indiv.ops[ind2] = op1
    new_indiv.history = "swap_m"
    return new_indiv
end

#TODO adding noise
function mutate(gate::Hgroup{N}) where N
    new_gate=Hgroup{N}(rand(1:6),gate.ghz_idx1,gate.ghz_idx2)
    return new_gate
end

function mutate(gate::Fgroup{N}) where N
    new_gate=Fgroup{N}(rand(1:8),gate.ghz_idx1,gate.ghz_idx2,gate.node_idx)
    return new_gate
end

#not sure this need to be added.
function mutate(gate::PauliGroup{N}) where N
    new_gate=PauliGroup{N}(rand(1:3),gate.ghz_idx,gate.qubit_idx)
    return new_gate
end

function mutate(indiv::Individual)
    new_indiv = deepcopy(indiv)
    new_indiv.ops = [mutate(gate) for gate in new_indiv.ops]
    new_indiv.history = "ops_m"
    return new_indiv
end

function new_child(indiv::Individual, indiv2::Individual, max_ops::Int)
    new_indiv = deepcopy(indiv)
    ops1, ops2 = indiv.ops, indiv2.ops
    if rand() < 0.5
        ops1 = ops1[end:-1:1]
    end
    if length(ops1) == 0 || length(ops2) == 0
        return new_indiv
    end

    if rand() < 0.5
        ops2 = ops2[end:-1:1]
    end
    new_indiv.ops = vcat(ops1[1:rand(1:min(length(ops1), max_ops))], ops2[1:rand(1:length(ops2))])[1:min(end, max_ops)]
    new_indiv.history = "child"
    return new_indiv
end

mutable struct Population
    n::Int
    ghz_num::Int
    qubit_num::Int
    f_in::Float64
    p2::Float64
    η::Float64
    population_size::Int
    starting_pop_multiplier::Int
    max_gen::Int
    max_ops::Int
    starting_ops::Int
    pairs::Int
    children_per_pair::Int
    mutants_per_individual_per_type::Int
    p_lose_operation::Float64
    p_add_operation::Float64
    p_swap_operations::Float64
    p_mutate_operations::Float64
    individuals::Vector{Individual}
    selection_history::Dict{String, Vector{Int64}}
end

function ini_pop!(population::Population)
    population.individuals = [Individual("random", population.qubit_num, population.ghz_num, population.f_in, [], 0.0, 0.0) for i=1:population.population_size*population.starting_pop_multiplier]
    n=population.qubit_num
    idx=population.ghz_num
    for indiv in population.individuals

        for i in rand(1:population.starting_ops)
            if rand() < 0.5
                #chance to add H gate to two random ghz state
                perm=randperm(idx)[1:2]
                gate = Hgroup{n}(rand(1:6),perm[1],perm[2])
            else
                #chance to add F gate to two random ghz state at random node.
                perm=randperm(idx)[1:2]
                gate = Fgroup{n}(rand(1:8),perm[1],perm[2],rand(1:idx-1))
            end
            push!(indiv.ops, gate)
        end
    end
end

function cull!(population::Population)
    population.individuals = population.individuals[1:population.population_size]
end

function sort!(population::Population) 
    Threads.@threads for indiv in population.individuals
        calculate_performance!(indiv) 
    end
    population.individuals = sort(population.individuals, by = x -> x.f_out, rev=true)
end

function step!(population::Population)
    for indiv in population.individuals
        indiv.history = "survivor"
    end

    parents = [(rand(population.individuals), rand(population.individuals)) for i=1:population.pairs]
    for (p1, p2) in parents
        population.individuals = vcat(population.individuals, [new_child(p1, p2, population.max_ops) for j=1:population.children_per_pair])
    end

    for indiv in population.individuals[1:population.population_size]
        population.individuals = vcat(population.individuals, [drop_op(indiv) for i=1:population.mutants_per_individual_per_type if rand() < population.p_lose_operation && length(indiv.ops) > 0])
        population.individuals = vcat(population.individuals, [gain_op(indiv) for i=1:population.mutants_per_individual_per_type if rand() < population.p_add_operation && length(indiv.ops) < population.max_ops])
        population.individuals = vcat(population.individuals, [swap_op(indiv) for i=1:population.mutants_per_individual_per_type if rand() < population.p_swap_operations && length(indiv.ops) > 0])
        population.individuals = vcat(population.individuals, [mutate(indiv) for i=1:population.mutants_per_individual_per_type if rand() < population.p_mutate_operations && length(indiv.ops) > 0])
    end

    sort!(population)
    cull!(population)
end

function run!(population::Population)
    println(Threads.nthreads())
    for hist in ["manual", "survivor", "random", "child", "drop_m", "gain_m", "swap_m", "ops_m"]
        population.selection_history[hist] = Vector{Int64}()
    end
    ini_pop!(population)
    sort!(population)
    cull!(population)
    for i=1:population.max_gen
        step!(population)
        for hist in ["manual", "survivor", "random", "child", "drop_m", "gain_m", "swap_m", "ops_m"]
            push!(population.selection_history[hist], reduce(+, [1 for indiv in population.individuals if indiv.history==hist], init=0))
        end
    end
end

function fidelities(population::Population)
    return [i.fitness for i in population.individuals]
end

function succ_probs(population::Population)
    return [i.performance.success_probability for i in population.individuals]
end

function calculate_performance!(indiv::Individual) 
    n=indiv.qubit_num
    idx=indiv.ghz_num
    count=0 #counting success trials

    #=
    for i=1:t
        state = rand(GHZState,n,idx,indiv.f_in)
        for g in indiv.ops
            if typeof(g) == NoisyMeasure
                s, result=apply!(state,g)
                if result == zeros(Int,n-1)
                    count+=1

                end
            else
                apply!(state,g)
            end
        end

    end
    =#
    t=10000
    fout=0
    state_tobe_measured=rand(1:idx)
    meas=GHZMeasure(n,3,state_tobe_measured)
    NM=NoisyMeasure(meas,1)
    for i in 1:t
        state = rand(GHZState,n,idx,indiv.f_in)
        for g in indiv.ops
            apply!(state,g)
        end

        s, result=apply!(state,NM)

        if result == zeros(Int,n-1)
            count+=1
            fout+=Int(s.phases[1:n] == zeros(Int,n))
        end
    end
    indiv.success=count/t
    indiv.f_out=fout/count
end

               #n,ghz,q,fin,p2,η,size,multi,mgen,mops,start_ops,pairs,children_per_pair,mutants_per_individual_per_type,p_single_operation_mutates,p_lose_operation,p_add_operation,p_swap_operations,p_mutate_operations,individuals,selection_history
POP=Population(5, 3, 3, 0.8, 1, 1, 100, 10, 10, 10, 5, 50, 2, 5, 0.2, 0.2, 0.2, 0.2, [], Dict())
run!(POP)
POP.individuals
sort!(POP)
op=POP.individuals[1].ops
length(op)
POP.individuals[1].success
POP.individuals[1].f_out

#test
test=POP.individuals[1]
calculate_performance!(test)
tn=3 #qubit
ti=2 #states
t=10000

meas=GHZMeasure(3,3,2)  #Z
NM=NoisyMeasure(meas,1)

c=0 #counting success trials
fout=0
for i in 1:t
    state = GHZPreserving.rand(GHZState,tn,ti,0.8)
    for g in test.ops
        apply!(state,g)
    end

    s, result=apply!(state,NM)

    if result == zeros(Int,tn-1)
        c+=1
        fout+=Int(s.phases[1:tn] == zeros(Int,tn))
    end
end
s=c/t
fout=fout/c

