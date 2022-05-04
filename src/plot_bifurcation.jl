include("StommelModel.jl")
include("NewtonMethod.jl")
include("ContMethod.jl")
include("detect_ranges.jl")

using LinearAlgebra
using LaTeXStrings
using Formatting
using .StommelModel
using .NewtonMethod
using .ContMethod1


function calBifur(;
    coe :: NamedTuple,
    use_Θ :: Bool,
    p_rng :: AbstractArray{Float64, 1},
)
    local m = Model(c=coe.c, μ=coe.μ, ν=coe.ν, p=p_rng[1], ξ=coe.ξ, use_Θ=use_Θ)
    m.x .= [ 0.0,  0.0 ]

    function _F_dFdx_dFdp(x, p)
        return (
            StommelModel.F(m; x=x, p=p[1]),
            StommelModel.dFdx(m; x=x, p=p[1]),
            StommelModel.dFdp(m; x=x, p=p[1]),
        )
    end

    scale = 1e-2
    cmi = ContMethod1.CMInfo(
        Nx = length(m.x),
        F_dFdx_dFdp = _F_dFdx_dFdp,
        mx   = [scale, scale],
        mp   = scale,
        nwt_min_iter = 5,
        nwt_max_iter = 20,
        res  = 1e-10,
        ds_exp_lb = -20.0,
        newton_callback = function(s) end,
        skip_unconverge = true, 
    )

    m.p  = p_rng[1]
    m.x .= 0

    _F_dFdx = function(x)
        return (
            StommelModel.F(m; x=x),
            StommelModel.dFdx(m; x=x),
        )
    end
    

    function record!(arr)
        eg_vals = eigvals(StommelModel.dFdx(m))
        stable = ( all(eg_vals .< 0) ) ? 1.0 : 0.0
        push!(arr, [m.p, m.x[1], m.x[2], StommelModel.cal_ψ(m), stable])
    end 
   
    local ppp = []
    NewtonMethod.doNewtonMethod!(
        m.x,
        _F_dFdx,
        1e-10,
        5,
        20;
    )
    
    # Now do continuition
    setS!(cmi, s=cmi.s, x=m.x, p=[ m.p ])
    setS!(cmi, s=cmi.ṡ, x=[0.0 ; 0.0], p=[ 0.01 ])

    record!(ppp)

    local p_container = [ m.p ] 
    
    for t = 1:10000

        ContMethod1.doContinuition!(cmi)
        setXP!(cmi; s=cmi.s, x=m.x, p=p_container)
        m.p = p_container[1]

        record!(ppp)

        if !( p_rng[1] <= m.p <= p_rng[2] ) || isnan(m.p)
            break
        end

    end
  
    local data = Dict()
    local varnames = [:p, :x, :y, :ψ, :stable]
    for (i, k) in enumerate(varnames)
        data[k] = zeros(Float64, length(ppp))
        for j=1:length(ppp)
            data[k][j] = ppp[j][i]
        end
    end

    return data
end

ξs = [-0.5, 0.5]
p_rng = [0.0, 3.0]


println("Loading PyPlot")
using PyPlot
plt = PyPlot
plt.ion()
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)


ax.set_xlabel(L"p")
ax.set_ylabel(L"\psi")
ax.grid()

cnt = 1
for ξ in ξs, use_Θ in [true, false]

    coe = (
        c = 3.6e3,
        μ = 5.0,
        ν = 1.0,
        ξ = ξ,
    )
        
    data = calBifur(;
        coe = coe,
        use_Θ = use_Θ,
        p_rng = p_rng,
    )

 
    vals, rngs = detectRanges(data[:stable])
    
    for (k, rng) in enumerate(rngs)

        args = Dict()
        if k==1
            args[:label] = format("\$\\xi\$ = {:.1f}, use_\$\\Theta\$ = {}", ξ, use_Θ)
        else
            args[:label] = nothing
        end

        args[:linestyle] = (vals[k] == 0) ? "dashed" : "solid"
        args[:color] = ["black", "red", "blue", "green"][cnt]
        ax.plot(data[:p][rng], data[:ψ][rng]; args...)

    end

    global cnt += 1

end

ax.legend()

sleep(500)
