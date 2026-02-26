# basis.jl
# 优化点：使用递归(DFS)生成基底，直接剪枝不满足守恒律的分支，避免生成后筛选
# 费米子算符操作
@inline function get_phase_and_state(state::Int, a::Int, b::Int, c::Int, d::Int)
    sign = 1
    # Annihilate d
    if isodd(count_ones(state & ((1<<(d-1))-1))) sign = -sign end
    temp = state ⊻ (1<<(d-1))
    # Annihilate c
    if isodd(count_ones(temp & ((1<<(c-1))-1))) sign = -sign end
    temp = temp ⊻ (1<<(c-1))
    # Create b
    if isodd(count_ones(temp & ((1<<(b-1))-1))) sign = -sign end
    temp = temp | (1<<(b-1))
    # Create a
    if isodd(count_ones(temp & ((1<<(a-1))-1))) sign = -sign end
    temp = temp | (1<<(a-1))
    return temp, sign
end




function generate_basis(p_val::Int, n::Int, target_m1, target_m2)
    # 1. 生成所有单体轨道标签
    orb_labels = Tuple{Float64,Float64,Float64,Float64}[]
    for two_j in 0:p_val
        j = two_j / 2.0; k = p_val/2.0 - j
        for two_m1 in -two_j:2:two_j
            m1 = two_m1 / 2.0
            two_k = Int(p_val - two_j)
            for two_m2 in -two_k:2:two_k
                m2 = two_m2 / 2.0
                push!(orb_labels, (j, m1, k, m2))
            end
        end
    end
    n_orb = length(orb_labels)
    valid_states = Int[]

    # 2. 递归搜索 (DFS)
    # idx: 当前考虑的轨道索引; n_left: 还需要选几个粒子
    function dfs(idx, n_left, current_m1, current_m2, current_state)
        # 成功条件
        if n_left == 0
            # 浮点数比较使用容差
            if abs(current_m1 - target_m1) < 1e-6 && abs(current_m2 - target_m2) < 1e-6
                push!(valid_states, current_state)
            end
            return
        end
        
        # 剪枝：如果剩下的轨道全选也不够，或者已经没有轨道可选
        if idx > n_orb || n_left < 0 return end

        # 路径 A: 选取当前轨道 idx
        label = orb_labels[idx]
        dfs(idx + 1, n_left - 1, 
            current_m1 + label[2], 
            current_m2 + label[4], 
            current_state | (1 << (idx-1)))
        
        # 路径 B: 不选取当前轨道
        dfs(idx + 1, n_left, current_m1, current_m2, current_state)
    end

    # 启动递归
    dfs(1, n, 0.0, 0.0, 0)
    
    # 3. 构建索引 Map
    state_map = Dict{Int, Int}()
    sizehint!(state_map, length(valid_states))
    for (i, s) in enumerate(valid_states)
        state_map[s] = i
    end

    return valid_states, state_map, orb_labels, n_orb
end




# [重要] 此函数返回 5 个值，增加了 basis_occ,用于计算correlation 


function get_basis(p::Int, n::Int, target_m1, target_m2)
    orb_labels = SO_5label(p)
    n_orb = length(orb_labels)
    
    valid_states = Int[]
    basis_occ = Vector{Vector{Int}}() # [新增] 缓存占据轨道列表
    
    for pos in combinations(1:n_orb, n)
        M1 = sum(orb_labels[i][2] for i in pos)
        M2 = sum(orb_labels[i][4] for i in pos)
        
        if abs(M1 - target_m1) < 1e-6 && abs(M2 - target_m2) < 1e-6
            state_int = 0
            for i in pos; state_int |= (1 << (i-1)); end
            push!(valid_states, state_int)
            push!(basis_occ, pos) # [新增] 存储位置列表
        end
    end
    
    state_map = Dict{Int, Int}(s => i for (i, s) in enumerate(valid_states))
    
    # 返回 5 个值
    return valid_states, basis_occ, state_map, orb_labels, n_orb
end


# 用于计算correlation function


function SO_5label(p::Int)
    states = Tuple{Float64,Float64,Float64,Float64}[]
    for two_j in 0:p
        j = two_j / 2.0
        k = p/2.0 - j
        for two_m1 in -two_j:2:two_j
            m1 = two_m1 / 2.0
            for two_m2 in -Int(2*k):2:Int(2*k)
                m2 = two_m2 / 2.0
                push!(states, (j, m1, k, m2))
            end
        end
    end
    return states
end