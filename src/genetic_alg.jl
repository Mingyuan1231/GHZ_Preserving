using .GHZPreserving
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
    ops::Vector{Union{Hgroup,Fgroup,NoisyMeasure}}
    success::Float64
    f_out::Vector{Float64}
end 

function drop_op(indiv::Individual) 
    new_indiv = deepcopy(indiv)
    deleteat!(new_indiv.ops, rand(1:length(new_indiv.ops)))
    new_indiv.history = "drop_m"
    return new_indiv
end


function gain_op(indiv::Individual)
    new_indiv = deepcopy(indiv)
    n=indiv.n
    q=indiv.q
    k=indiv.k
    r=indiv.r

    #measure number
    measure_count = count(op -> isa(op,NoisyMeasure), indiv.ops)

    if rand() < 0.2 && measure_count < (n-r) #adding measurement
        #chance to measure
        #rand[1,3] as X or Z measure, no Y for now, measure random ghz state
        measure=GHZMeasure(n,rand([1,3]),rand(k+1:r))  #as keep k result state, so only measure (k+1:r) state. 
        gate=NoisyMeasure(measure, indiv.η)
    elseif rand() < 0.5
        #chance to add H gate to two random ghz state
        perm=randperm(n)[1:2]
        gate = Hgroup{q}(rand(1:6),perm[1],perm[2])
    else
        #chance to add F gate to two random ghz state at random node.
        perm=randperm(n)[1:2]
        gate = Fgroup{q}(rand(1:8),perm[1],perm[2],rand(1:q-1))
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

function mutate(gate)
    #if it is a measure, then mutate to different basis
    if isa(gate, NoisyMeasure)
        current_measure=gate.m
        n=current_measure.n
        idx=current_measure.ghz_idx
        new_basis=current_measure.basis_idx == 1 ? 3 : 1

        new_measure=GHZMeasure(n,new_basis,idx)
        p=gate.p

        return NoisyMeasure(new_measure,p)
    end

    #if it is Hgroup or Fgroup, then mutate to random different gate
    N = typeof(gate).parameters[1]
    idx1=gate.ghz_idx1
    idx2=gate.ghz_idx2
    
    if rand() < 0.5
        new_gate=Hgroup{N}(rand(1:6),idx1,idx2)
    else
        node_idx = typeof(gate) <: Fgroup ? gate.node_idx : rand(1:N-1)  #if Fgroup, keep the node_idx, else randomize, this can be modify to fully randomize
        new_gate=Fgroup{N}(rand(1:8),idx1,idx2,node_idx)
    end

    return new_gate
end

function mutate(indiv::Individual)
    new_indiv = deepcopy(indiv)
    new_indiv.ops = [mutate(gate) for gate in new_indiv.ops]
    new_indiv.history = "ops_m"
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
    measure_count = count(op -> isa(op,NoisyMeasure), select_ops)
    while measure_count > (n-r)
        for i in eachindex(select_ops)
            if isa(select_ops[i],NoisyMeasure)
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
        gate=NoisyMeasure(measure,indiv.η)
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

function ini_pop!(population::Population)
    population.individuals = [Individual("random", population.n, population.q, population.k, population.r, population.f_in, population.p2, population.η, [], 0.0, []) for _ =1:population.population_size]
    n=population.n
    q=population.q
    k=population.k
    r=population.r

    for indiv in population.individuals
        #adding random gates
        num_ops=rand(1:population.starting_ops-1)
        for i in 1:num_ops
            if rand() < 0.5
                #chance to add H gate to two random ghz state
                perm=randperm(n)[1:2]
                gate = Hgroup{q}(rand(1:6),perm[1],perm[2])
            else
                #chance to add F gate to two random ghz state at random node.
                perm=randperm(n)[1:2]
                gate = Fgroup{q}(rand(1:8),perm[1],perm[2],rand(1:q-1))
            end
            push!(indiv.ops, gate)
        end

        #add measure to the end of the circuit, always keep first k state, measure the rest
        for i in k+1:r
            basis=rand([1,3]) #1 for X, 3 for Z
            measure=GHZMeasure(q,basis,i)
            gate=NoisyMeasure(measure,population.η)
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

    # 使用自定义比较逻辑排序，并添加了数组长度的检查
    population.individuals = sort(population.individuals, lt=(x, y) -> begin
        # 检查每个个体的 f_out 是否足够长
        if length(x.f_out) < 2 || length(y.f_out) < 2
            # 如果任一个体的 f_out 不包含至少两个元素，则退回到只比较第一个元素
            return x.f_out[1] > y.f_out[1]
        else
            # 比较第一个fidelity
            if abs(x.f_out[1] - y.f_out[1]) / max(x.f_out[1], y.f_out[1]) < 0.05
                # 如果第一个fidelity差别小于5%，比较第二个
                return x.f_out[2] > y.f_out[2]
            else
                # 否则只比较第一个
                return x.f_out[1] > y.f_out[1]
            end
        end
    end)
end
#=
function sort!(population::Population) 
    Threads.@threads for indiv in population.individuals
        calculate_performance!(indiv) 
    end
    population.individuals = sort(population.individuals, by = x -> x.f_out, rev=true)
end
=#

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

#TODO hashing?

function calculate_performance!(indiv::Individual) 
    n=indiv.n
    q=indiv.q
    k=indiv.k
    r=indiv.r
    
    t=10000 #total number of trials
    count = 0 #counting success trials
    fout = zeros(Int, k) #fidelity

    for i in 1:t
        state = rand(GHZState,q,n,indiv.f_in)
        status = continue_stat
 
        for g in indiv.ops
            if isa(g,NoisyMeasure)
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

function show(io::IO, fg::Fgroup{N}) where N
    print(io, "Fgroup{$N}($(fg.gate_idx), $(fg.ghz_idx1), $(fg.ghz_idx2), $(fg.node_idx))")
end

              #n, q, k, r, fin,p2,η,size,mgen,mops,start_ops,pairs,children_per_pair,mutants_per_individual_per_type,p_single_operation_mutates,p_lose_operation,p_add_operation,p_swap_operations,p_mutate_operations,individuals,selection_history
POP=Population(5, 3, 2, 4, 0.8, 0.99, 1, 20, 10, 20, 5, 50, 2, 5, 0.2, 0.2, 0.2, 0.2, [], Dict())
ini_pop!(POP)
sort!(POP)
cull!(POP)

step!(POP)

length(POP.individuals)
POP.individuals
for i in 1:20
    for j in 1:length(POP.individuals[i].ops)
        println(typeof(POP.individuals[i].ops[j]))
    end
    println(" ")
end

run!(POP)


idv=POP.individuals[1]
length(idv.ops)
for op in idv.ops
    println(op)
end

idv2=POP.individuals[2]
length(idv2.ops)
for op in idv2.ops
    println(op)
end

child=new_child(idv,idv2,10)
child.ops

gain_op(child)
length(child.ops)


for i in 1:20
    calculate_performance!(POP.individuals[i])
    println(POP.individuals[i].success)
    println(POP.individuals[i].f_out)
end


testsucc_dict=Dict{Tuple{Float64,Float64},Float64}()
testfide_dict=Dict{Tuple{Float64,Float64},Float64}()
for fin in 0:0.1:1
    POP=Population(5, 3, 2, 4, fin, 0.99, 1, 20, 10, 20, 5, 50, 2, 5, 0.2, 0.2, 0.2, 0.2, [], Dict())
    run!(POP)
    testsucc_dict[(0.99,fin)]=mean([indiv.success for indiv in POP.individuals])
    testfide_dict[(0.99,fin)]=mean([mean(indiv.f_out) for indiv in POP.individuals])
end


p2_values = [0.99]
for p2 in p2_values
    fin_values = [key[2] for key in keys(testsucc_dict) if key[1] == p2]
    success_rates = [testsucc_dict[(p2, fin)] for fin in fin_values]
    
    sorted_indices = sortperm(fin_values)
    #plot!(p_success, fin_values[sorted_indices], success_rates[sorted_indices], fillalpha=0.3, label="p2 = $p2")
    plot!(ps, fin_values[sorted_indices], success_rates[sorted_indices], fillalpha=0.3, label="p2 = $p2")
end

# Plotting for fidelity_dict
for p2 in p2_values
    fin_values = [key[2] for key in keys(testfide_dict) if key[1] == p2]
    fidelities = [testfide_dict[(p2, fin)] for fin in fin_values]
    
    sorted_indices = sortperm(fin_values)
    #plot!(p_fidelity, fin_values[sorted_indices], fidelities[sorted_indices], fillalpha=0.3, label="p2 = $p2")
    plot!(p, fin_values[sorted_indices], fidelities[sorted_indices], fillalpha=0.3, label="p2 = $p2")
end