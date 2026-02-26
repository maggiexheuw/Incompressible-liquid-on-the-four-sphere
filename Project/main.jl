# main_lazy.jl
# 优化策略：
# 1. 移除全量预计算 (对于小系统，预计算开销 > 收益)
# 2. 保留 fast_cg 和 InteractionParams 带来的数学加速
# 3. 保留多线程并行构建
# 4. 这种模式称为 "On-the-fly" 计算，适合 N_orb 较大但填充较少的情况

using LinearAlgebra
using SparseArrays
using Arpack
using Base.Threads
using Printf

# 引入之前的三个优化模块 (确保这三个文件在同一目录)
include("reduce.jl")
include("cg.jl")
include("basis.jl")



function main_optimized(p_val::Int, n::Int, target_m1, target_m2, alpha, beta)
    println("--- 开始计算 (p=$p_val, n=$n) [Lazy Mode] ---")
    
    # [Step 1] 准备参数
    params = get_interaction_params(alpha, beta)
    
    # [Step 2] 生成基矢
    t_basis = time()
    L_basis, state_map, orb_labels, n_orb = generate_basis(p_val, n, target_m1, target_m2)
    dim = length(L_basis)
    println("   Basis Dimension: $dim (耗时 $(round(time()-t_basis, digits=4))s)")
    
    if dim == 0 return [], [], spzeros(0,0) end

    # [Step 3] 预计算配对索引 (这个开销很小且必须做)
    # 我们只建立 (M1, M2) -> [(a,b)...] 的映射，不计算具体矩阵元
    pair_map = Dict{Tuple{Float64,Float64}, Vector{Tuple{Int,Int}}}()
    for i in 1:n_orb, j in i+1:n_orb
        m1 = orb_labels[i][2] + orb_labels[j][2]
        m2 = orb_labels[i][4] + orb_labels[j][4]
        key = (m1, m2)
        if !haskey(pair_map, key) pair_map[key] = Tuple{Int,Int}[] end
        push!(pair_map[key], (i, j))
    end

    # [Step 4] 多线程构建哈密顿量 (On-the-fly 计算)
    println("Building Hamiltonian (Threads: $(Threads.nthreads()))...")
    t_build = time()
    
    n_th = Threads.nthreads()
    rows_th = [Int[] for _ in 1:n_th]
    cols_th = [Int[] for _ in 1:n_th]
    vals_th = [Float64[] for _ in 1:n_th]
    
    # 预估容量
    est_nnz = dim * n * n * 5
    for k in 1:n_th
        sizehint!(rows_th[k], est_nnz ÷ n_th)
        sizehint!(cols_th[k], est_nnz ÷ n_th)
        sizehint!(vals_th[k], est_nnz ÷ n_th)
    end

    Threads.@threads for i in 1:dim
        tid = Threads.threadid()
        state_i = L_basis[i]
        
        # 快速提取占据轨道
        occ_list = Int[]
        ptr = 1
        curr = state_i
        while curr > 0
            tz = trailing_zeros(curr)
            push!(occ_list, ptr + tz)
            curr >>= (tz + 1)
            ptr += (tz + 1)
        end
        n_occ = length(occ_list)

        # 遍历 c, d
        for idx1 in 1:n_occ
            c = occ_list[idx1]
            for idx2 in idx1+1:n_occ
                d = occ_list[idx2]
                
                s_c, s_d = orb_labels[c], orb_labels[d]
                m_key = (s_c[2] + s_d[2], s_c[4] + s_d[4])
                
                if haskey(pair_map, m_key)
                    possible_pairs = pair_map[m_key]
                    mask_rem = state_i ⊻ (1<<(c-1)) ⊻ (1<<(d-1))
                    
                    for (a, b) in possible_pairs
                        # Pauli 检查
                        if (mask_rem & (1<<(a-1))) == 0 && (mask_rem & (1<<(b-1))) == 0
                            
                            s_a, s_b = orb_labels[a], orb_labels[b]
                            
                            # === 关键修改：实时计算，不查表 ===
                            # 因为 fast_cg 极快，这里直接算的代价很低
                            # 且避免了 O(N_orb^4) 的预热时间
                            val_dir = calculate_matrix_element(params, s_a, s_b, s_c, s_d)
                            val_exc = calculate_matrix_element(params, s_a, s_b, s_d, s_c)
                            total_val = -(val_dir - val_exc)
                            
                            if abs(total_val) > 1e-10
                                new_state, sign = get_phase_and_state(state_i, a, b, c, d)
                                if haskey(state_map, new_state)
                                    j = state_map[new_state]
                                    push!(rows_th[tid], j)
                                    push!(cols_th[tid], i)
                                    push!(vals_th[tid], total_val * sign)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # [Step 5] 合并与对角化
    rows = vcat(rows_th...)
    cols = vcat(cols_th...)
    vals = vcat(vals_th...)
    
    # 构造矩阵
    H = sparse(rows, cols, vals, dim, dim)
    H = (H + H') / 2
    
    println("   Hamiltonian Built (耗时 $(round(time()-t_build, digits=4))s)")
    println("Diagonalizing...")
    
    t_diag = time()
    if dim < 2000
        eigen_vals, eigen_vecs = eigen(Matrix(H))
    else
        try
            eigen_vals, eigen_vecs = eigs(H, nev=min(20, dim), which=:SR)
            eigen_vals = real(eigen_vals)
        catch e
            println("Warning: Arpack failed, using dense diagonalization.")
            eigen_vals, eigen_vecs = eigen(Matrix(H))
        end
    end
    println("   Diagonalization done (耗时 $(round(time()-t_diag, digits=4))s)")
    
    return eigen_vals, eigen_vecs,H
end

# === 测试运行 ===
println("Warmup run (compilation)...")
# 第一次运行会包含编译时间，请忽略这行的时间
main_optimized(2, 2, 0.0, 0.0, 1.0, 1.0) 

println("\n=== Real Run (p=3, n=4) ===")
t_start = time()
vals, vecs = main_optimized(3, 4, 0.0, 0.0, 1.0, 4.3)
t_end = time()

println("Total Execution Time: $(round(t_end - t_start, digits=4)) s")
for (i, e) in enumerate(vals[1:min(20, end)])
    @printf("E_%d = %.8f\n", i, e)
end