# reduce.jl
# 优化点：引入 InteractionParams 结构体，解决全局变量导致的类型不稳定问题

struct InteractionParams
    alpha::Float64
    beta::Float64
    reduced::Dict{Tuple{Float64,Float64,Int}, Float64}
    reduced1::Dict{Tuple{Float64,Float64,Int}, Float64}
    reduced2::Dict{Tuple{Float64,Float64}, Float64}
end

function get_interaction_params(alpha_val::Float64, beta_val::Float64)
    # 初始化 reduced 字典
    r = Dict{Tuple{Float64,Float64,Int}, Float64}()
    r[(1.5,1.5,1)] = sqrt(3)/2; r[(1.5,1,1)] = sqrt(5)/4
    r[(1.5,0.5,1)] = sqrt(3)/6; r[(1.5,0.5,2)] = 0.0
    r[(1.5,0,1)] = 0.0;         r[(1.5,0,2)] = 0.0
    r[(1,1.5,1)] = sqrt(5)/4;   r[(1,1,1)] = sqrt(3)/3; r[(1,1,2)] = 0.5
    r[(1,0.5,1)] = sqrt(3)/4;   r[(1,0.5,2)] = sqrt(3)/4
    r[(1,0,2)] = sqrt(3)/6;     r[(1,0,1)] = 0.0
    r[(0.5,1.5,1)] = sqrt(3)/6; r[(0.5,1.5,2)] = 0.0
    r[(0.5,1,1)] = sqrt(3)/4;   r[(0.5,1,2)] = sqrt(3)/4
    r[(0.5,0.5,1)] = 0.5;       r[(0.5,0.5,2)] = sqrt(3)/3
    r[(0.5,0,2)] = sqrt(5)/4
    r[(0,1.5,1)] = 0.0;         r[(0,1.5,2)] = 0.0
    r[(0,1,2)] = sqrt(3)/6;     r[(0,1,1)] = 0.0
    r[(0,0.5,2)] = sqrt(5)/4;   r[(0,0,2)] = sqrt(3)/2

    # 初始化 reduced1 字典
    r1 = Dict{Tuple{Float64,Float64,Int}, Float64}()
    r1[(1.5,1.5,1)] = 0.5;       r1[(1.5,1,1)] = sqrt(3)/4
    r1[(1.5,0.5,1)] = sqrt(3)/6; r1[(1.5,0.5,2)] = 0.0
    r1[(1.5,0,1)] = 0.0;         r1[(1.5,0,2)] = 0.0
    r1[(1,1.5,1)] = sqrt(3)/4;   r1[(1,1,1)] = sqrt(3)/3; r1[(1,1,2)] = -sqrt(3)/2
    r1[(1,0.5,1)] = sqrt(5)/4;   r1[(1,0.5,2)] = -sqrt(5)/4
    r1[(1,0,2)] = -sqrt(3)/6;    r1[(1,0,1)] = 0.0
    r1[(0.5,1.5,1)] = sqrt(3)/6; r1[(0.5,1.5,2)] = 0.0
    r1[(0.5,1,1)] = sqrt(5)/4;   r1[(0.5,1,2)] = -sqrt(5)/4
    r1[(0.5,0.5,1)] = sqrt(3)/2; r1[(0.5,0.5,2)] = -sqrt(3)/3
    r1[(0.5,0,2)] = -sqrt(3)/4
    r1[(0,1.5,1)] = 0.0;         r1[(0,1.5,2)] = 0.0
    r1[(0,1,2)] = -sqrt(3)/6;    r1[(0,1,1)] = 0.0
    r1[(0,0.5,2)] = -sqrt(3)/4;  r1[(0,0,2)] = -0.5

    # 初始化 reduced2 字典
    r2 = Dict{Tuple{Float64,Float64}, Float64}()
    r2[(1.5,1.5)] = 0.0;        r2[(1.5,1)] = 1/sqrt(2)
    r2[(1.5,0.5)] = 1/sqrt(2);  r2[(1.5,0)] = 0.5
    r2[(1,1.5)] = -1/sqrt(2);   r2[(1,1)] = 0.0
    r2[(1,0.5)] = 0.5;          r2[(1,0)] = 1/sqrt(2)
    r2[(0.5,1.5)] = -1/sqrt(2); r2[(0.5,1)] = -0.5
    r2[(0.5,0.5)] = 0.0;        r2[(0.5,0)] = 1/sqrt(2)
    r2[(0,1.5)] = -0.5;         r2[(0,1)] = -1/sqrt(2)
    r2[(0,0.5)] = -1/sqrt(2);   r2[(0,0)] = 0.0

    return InteractionParams(alpha_val, beta_val, r, r1, r2)
end