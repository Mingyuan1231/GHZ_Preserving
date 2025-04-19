##
using GHZPreserving
using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using Random
using Statistics

mutable struct Individual
    history::String
    n::Int          #number of raw state
    q::Int          #number of qubit
    k::Int          #number of state to keep
    r::Int          #number of register
    f_in::Float64   #input fidelity
    p2::Float64     #depolarizing noise level
    η::Float64      #measurement noise level
    ops::Vector{Union{Hgroup,Bgroup,NoisyGHZMeasure}}
    success::Float64
    f_out::Vector{Float64}
end 

function drop_op(indiv::Individual) 
    new_indiv = deepcopy(indiv)
    deleteat!(new_indiv.ops, rand(1:length(new_indiv.ops)))
    new_indiv.history = "drop_m"
    return new_indiv
end

function gain_op_H_only(indiv::Individual)
    new_indiv = deepcopy(indiv)
    n=indiv.n
    q=indiv.q
    k=indiv.k
    r=indiv.r

    measure_count = count(op -> isa(op, NoisyGHZMeasure), indiv.ops)

    if rand() < 0.3 && measure_count < (n - k)
        measure = GHZMeasure(n, rand([1, 3]), rand(k+1:r))
        gate = NoisyGHZMeasure(measure, indiv.η)
    else
        perm = randperm(r)[1:2]
        gate = Hgroup{q}(2, perm[1], perm[2])
    end

    if isempty(new_indiv.ops)
        push!(new_indiv.ops, gate)
    else
        insert!(new_indiv.ops, rand(1:length(new_indiv.ops)), gate)
    end

    new_indiv.history = "gain_H"
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

function mutate_H_only(gate)
    if isa(gate, NoisyGHZMeasure)
        current_measure = gate.m
        new_basis = current_measure.basis_idx == 1 ? 3 : 1
        new_measure = GHZMeasure(current_measure.n, new_basis, current_measure.ghz_idx)
        return NoisyGHZMeasure(new_measure, gate.p)
    elseif isa(gate, Hgroup)
        N = typeof(gate).parameters[1]
        idx1, idx2 = gate.ghz_idx1, gate.ghz_idx2
        return Hgroup{N}(2, idx1, idx2)
    else
        return gate
    end
end

function mutate_H_only(indiv::Individual)
    new_indiv = deepcopy(indiv)
    new_indiv.ops = [mutate_H_only(g) for g in new_indiv.ops]
    new_indiv.history = "mutate_H"
    return new_indiv
end


function new_child(indiv::Individual, indiv2::Individual, max_ops::Int)
    n=indiv.n
    q=indiv.q
    k=indiv.k
    r=indiv.r
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

    #random select operation
    combined_ops = vcat(ops1, ops2)
    select_ops = combined_ops[randperm(length(combined_ops))[1:min(length(combined_ops), max_ops)]]
    
    #delete extra measures, the number limit for middle step measure is n-r
    measure_count = count(op -> isa(op,NoisyGHZMeasure), select_ops)
    while measure_count > (n-r)
        for i in eachindex(select_ops)
            if isa(select_ops[i],NoisyGHZMeasure)
                deleteat!(select_ops, i)
                measure_count -= 1
                break
            end
        end
    end

    #add measure to the end of the circuit, always keep first k state, measure the rest
    for i in k+1:r
        basis=rand([1,3]) #1 for X, 3 for Z
        measure=GHZMeasure(q,basis,i)
        gate=NoisyGHZMeasure(measure,indiv.η)
        push!(select_ops, gate)
    end

    new_indiv.ops = select_ops
    new_indiv.history = "child"
    return new_indiv
end


mutable struct Population
    n::Int
    q::Int
    k::Int
    r::Int
    f_in::Float64
    p2::Float64
    η::Float64
    population_size::Int
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

function ini_pop_H_only!(population::Population)
    population.individuals = [Individual("random", population.n, population.q, population.k, population.r, population.f_in, population.p2, population.η, [], 0.0, []) for _ =1:population.population_size]
    n, q, k, r = population.n, population.q, population.k, population.r

    for indiv in population.individuals
        num_ops = rand(1:population.starting_ops - 1)
        for _ in 1:num_ops
            perm = randperm(r)[1:2]
            gate = Hgroup{q}(2, perm[1], perm[2])
            push!(indiv.ops, gate)
        end
        for i in k+1:r
            basis = rand([1,3])
            gate = NoisyGHZMeasure(GHZMeasure(q, basis, i), population.η)
            push!(indiv.ops, gate)
        end
    end
end


function cull!(population::Population)
    population.individuals = population.individuals[1:population.population_size]
end

function sort!(population::Population, current_gen::Int) 
    Threads.@threads for indiv in population.individuals
        calculate_performance!(indiv,current_gen) 
    end
    population.individuals = sort(population.individuals, by = x -> x.f_out, rev=true)

    new_indiv=[]
    prob=rand()
    for i in 1:length(population.individuals)
        if prob < 0.8
            push!(new_indiv,population.individuals[i])
        end
    end
end


function step_H_only!(population::Population, current_gen::Int)
    for indiv in population.individuals
        indiv.history = "survivor"
    end

    parents = [(rand(population.individuals), rand(population.individuals)) for _ in 1:population.pairs]
    for (p1, p2) in parents
        population.individuals = vcat(population.individuals, [new_child(p1, p2, population.max_ops) for _ in 1:population.children_per_pair])
    end

    for indiv in population.individuals[1:population.population_size]
        if rand() < population.p_lose_operation && !isempty(indiv.ops)
            population.individuals = vcat(population.individuals, drop_op(indiv))
        end
        if rand() < population.p_add_operation && length(indiv.ops) < population.max_ops
            population.individuals = vcat(population.individuals, gain_op_H_only(indiv))
        end
        if rand() < population.p_swap_operations && !isempty(indiv.ops)
            population.individuals = vcat(population.individuals, swap_op(indiv))
        end
        if rand() < population.p_mutate_operations && !isempty(indiv.ops)
            population.individuals = vcat(population.individuals, mutate_H_only(indiv))
        end
    end

    sort!(population, current_gen)
    cull!(population)
end


function run_H_only!(population::Population)
    for hist in ["manual", "survivor", "random", "child", "drop_m", "gain_H", "swap_m", "mutate_H"]
        population.selection_history[hist] = Vector{Int64}()
    end
    ini_pop_H_only!(population)
    sort!(population, 1)
    cull!(population)
    for i in 1:population.max_gen
        step_H_only!(population, i)
        for hist in keys(population.selection_history)
            push!(population.selection_history[hist], count(indiv -> indiv.history == hist, population.individuals))
        end
    end
end


#TODO hashing?

function calculate_performance!(indiv::Individual, current_gen::Int) 
    n=indiv.n
    q=indiv.q
    k=indiv.k
    r=indiv.r
    
    t=1000*current_gen #total number of trials
    count = 0 #counting success trials
    fout = zeros(Int, k) #fidelity

    for i in 1:t
        state = rand(GHZState,q,r,indiv.f_in) # only need to consider r raw states, as constraint by register number.
        status = continue_stat
 
        for g in indiv.ops
            if isa(g,NoisyGHZMeasure)
                state, status=applywstatus!(state, g, indiv.f_in)
                #if fail, break the loop
                if status == failure_stat
                    break
                end
                
            else
                apply!(state,g)
                depolarize!(state, g, indiv.p2)
            end
        end

        if status == continue_stat
            count+=1
            for j in 1:k
                start_idx=(j-1)*q+1
                end_idx=j*q
                fout[j]+=Int(state.phases[start_idx:end_idx] == zeros(Int,q))
            end
            #fout+=Int(state.phases[1:q] == zeros(Int,q))
        end
    end

    indiv.success=count/t
    indiv.f_out=fout/count
end

import Base.show
function show(io::IO, hg::Hgroup{N}) where N
    print(io, "Hgroup{$N}($(hg.gate_idx), $(hg.ghz_idx1), $(hg.ghz_idx2))")
end    

##


                #n, q, k, r, fin,   p2,   η, size,mgen,mops,start_ops,pairs,children_per_pair,mutants_per_individual_per_type,p_single_operation_mutates,p_lose_operation,p_add_operation,p_swap_operations,p_mutate_operations,individuals,selection_history
POP_H=Population(4, 3, 1, 3, 0.7, 0.99, 0.99, 30, 20, 20, 5, 50, 2, 5, 0.2, 0.2, 0.2, 0.2, [], Dict())

run_H_only!(POP_H)

idv=deepcopy(POP_H.individuals[1])
length(idv.ops)
for op in idv.ops
    println(op)
end


##
using Plots

fide_data=Dict{Tuple{Int,Int,Float64},Vector}()
succ_data=Dict{Tuple{Int,Int,Float64},Float64}()
circuit_data=[]

for i in 4:8
    for regi in 3:i #register number should equal or less than raw state number, more register is redundant
        for j in 0.7:0.025:1.0
        POP=Population(i, 3, 1, regi, j, 0.99, 0.99, 30, 20, 20, 5, 50, 2, 5, 0.2, 0.2, 0.2, 0.2, [], Dict())
        run_H_only!(POP)
 
        fide_data[(i,regi,j)]=POP.individuals[1].f_out
        succ_data[(i,regi,j)]=POP.individuals[1].success
        push!(circuit_data,POP.individuals[1]) # save the best circuit for plotting later
        end
    end
end

##
open("circuit_data.txt", "w") do io
    for line in circuit_data
        println(io, line)
    end
end
##
