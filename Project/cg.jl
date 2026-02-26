# cg.jl
# 优化点：修复 factorial 报错，内联 CG 计算，接收 InteractionParams 参数

# 预计算阶乘表：先用 big(i) 计算防止溢出，再转为 Float64 存储
const FACT = Float64[factorial(big(i)) for i in 0:100]

# 阶乘辅助函数：直接查表
@inline function fast_fact(n::Float64)
    return FACT[Int(round(n)) + 1]
end

# 优化后的 CG 系数计算
@inline function fast_cg(j1, m1, j2, m2, j3, m3)
    # 快速守恒律检查
    if abs(m1 + m2 - m3) > 1e-9 return 0.0 end
    if j1 + j2 < j3 || j3 < abs(j1 - j2) return 0.0 end

    a = (2*j3+1) * fast_fact(j1+j2-j3) * fast_fact(j1-j2+j3) * fast_fact(-j1+j2+j3)
    b = fast_fact(j1+j2+j3+1)
    c = fast_fact(j1+m1) * fast_fact(j1-m1) * fast_fact(j2+m2) * fast_fact(j2-m2) * fast_fact(j3+m3) * fast_fact(j3-m3)
    
    min1 = Int(max(0, j2-j3-m1, j1-j3+m2))
    max1 = Int(min(j1+j2-j3, j1-m1, j2+m2))
    
    sum0 = 0.0
    for k in min1:max1
        denom = fast_fact(Float64(k)) * fast_fact(j1+j2-j3-k) * fast_fact(j1-m1-k) * fast_fact(j2+m2-k) * fast_fact(j3-j2+m1+k) * fast_fact(j3-j1-m2+k)
        
        term = 1.0 / denom
        if isodd(k)
            sum0 -= term
        else
            sum0 += term
        end
    end
    return sum0 * sqrt(a * c / b)
end

# 矩阵元计算：参数化 reduced 字典
function calculate_matrix_element(p::InteractionParams, t1, t2, t3, t4)
    alpha, beta = p.alpha, p.beta
    j1, m1, k1, n1 = t1; j2, m2, k2, n2 = t2
    j3, m3, k3, n3 = t3; j4, m4, k4, n4 = t4
    
    element = 0.0
    tol = 1e-9

    # Block 1
    if abs((j1+j2) - (j3+j4)) < tol
        term = beta * get(p.reduced2, (j1,j2), 0.0) * get(p.reduced2, (j3,j4), 0.0)
        if abs(term) > 1e-12
            element += term * fast_cg(j1,m1,j2,m2,j1+j2,m1+m2) * fast_cg(k1,n1,k2,n2,k1+k2,n1+n2) *
                              fast_cg(j3,m3,j4,m4,j3+j4,m3+m4) * fast_cg(k3,n3,k4,n4,k3+k4,n3+n4)
        end
    end
    
    # Block 2
    if abs((k1+k2) - (k3+k4)) < tol && abs(j1-j2) <= j1+j2-1 && abs(m1+m2) <= j1+j2-1
        term = alpha * get(p.reduced, (j1,j2,1), 0.0) * get(p.reduced, (j3,j4,1), 0.0) + 
               beta  * get(p.reduced1, (j1,j2,1), 0.0) * get(p.reduced1, (j3,j4,1), 0.0)
        if abs(term) > 1e-12
            element += term * fast_cg(j1,m1,j2,m2,j1+j2-1,m1+m2) * fast_cg(k1,n1,k2,n2,k1+k2,n1+n2) *
                              fast_cg(j3,m3,j4,m4,j3+j4-1,m3+m4) * fast_cg(k3,n3,k4,n4,k3+k4,n3+n4)
        end
    end

    # Block 3
    if abs((j1+j2) - (j3+j4)) < tol && abs(k1-k2) <= k1+k2-1 && abs(n1+n2) <= k1+k2-1
        term = alpha * get(p.reduced, (j1,j2,2), 0.0) * get(p.reduced, (j3,j4,2), 0.0) + 
               beta  * get(p.reduced1, (j1,j2,2), 0.0) * get(p.reduced1, (j3,j4,2), 0.0)
        if abs(term) > 1e-12
            element += term * fast_cg(j1,m1,j2,m2,j1+j2,m1+m2) * fast_cg(k1,n1,k2,n2,k1+k2-1,n1+n2) *
                              fast_cg(j3,m3,j4,m4,j3+j4,m3+m4) * fast_cg(k3,n3,k4,n4,k3+k4-1,n3+n4)
        end
    end

    # Block 4
    if abs((j1+j2) - (j3+j4 - 1)) < tol && abs(k1-k2) <= k1+k2-1 && abs(n1+n2) <= k1+k2-1 && abs(j3-j4) <= j3+j4-1 && abs(m3+m4) <= j3+j4-1
        term = alpha * get(p.reduced, (j1,j2,2), 0.0) * get(p.reduced, (j3,j4,1), 0.0) + 
               beta  * get(p.reduced1, (j1,j2,2), 0.0) * get(p.reduced1, (j3,j4,1), 0.0)
        if abs(term) > 1e-12
            element += term * fast_cg(j1,m1,j2,m2,j1+j2,m1+m2) * fast_cg(k1,n1,k2,n2,k1+k2-1,n1+n2) *
                              fast_cg(j3,m3,j4,m4,j3+j4-1,m3+m4) * fast_cg(k3,n3,k4,n4,k3+k4,n3+n4)
        end
    end

    # Block 5
    if abs((j1+j2 - 1) - (j3+j4)) < tol && abs(j1-j2) <= j1+j2-1 && abs(m1+m2) <= j1+j2-1 && abs(k3-k4) <= k3+k4-1 && abs(n3+n4) <= k3+k4-1
        term = alpha * get(p.reduced, (j1,j2,1), 0.0) * get(p.reduced, (j3,j4,2), 0.0) + 
               beta  * get(p.reduced1, (j1,j2,1), 0.0) * get(p.reduced1, (j3,j4,2), 0.0)
        if abs(term) > 1e-12
            element += term * fast_cg(j1,m1,j2,m2,j1+j2-1,m1+m2) * fast_cg(k1,n1,k2,n2,k1+k2,n1+n2) *
                              fast_cg(j3,m3,j4,m4,j3+j4,m3+m4) * fast_cg(k3,n3,k4,n4,k3+k4-1,n3+n4)
        end
    end

    return element
end