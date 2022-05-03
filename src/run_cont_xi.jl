include("StommelModel.jl")
include("NewtonMethod.jl")
include("ContMethod.jl")

using LinearAlgebra
using LaTeXStrings
using Formatting
using .StommelModel
using .NewtonMethod
using .ContMethod1

c = 3.6e3
μ = 7.0
ν = 1.0
ξ_vec = collect(range(-1, 10, length=21))
p_rng = [0.0, 3.0]

m = Model(c=c, μ=μ, ν=ν, p=p_rng[1], ξ=0.0)

m.x .= [ 0.0,  0.0 ]


println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true, subplot_kw=Dict("projection" => "3d"))

function _F_dFdx_dFdp(x, p)
    return (
        StommelModel.F(m; x=x, p=p[1]),
        StommelModel.dFdx(m; x=x, p=p[1]),
        StommelModel.dFdp(m; x=x, p=p[1]),
    )
end


scale = 5e-2

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

ax.set_xlabel(L"p")
ax.set_ylabel(L"\xi")
ax.set_zlabel(L"\psi")
ax.grid()

        

# First, find a steady state
#_F_dFdx_dFdp(x, m.p)[1:2]
p_container = [ m.p ] 
for ξ in ξ_vec

    println("Now do ξ = ", ξ)

    m.ξ  = ξ
    m.p  = p_rng[1]
    m.x .= 0

    _F_dFdx = function(x)
        return (
            StommelModel.F(m; x=x),
            StommelModel.dFdx(m; x=x),
        )
    end
    
    local ppp = []

    NewtonMethod.doNewtonMethod!(
        m.x,
        _F_dFdx,
        1e-10,
        5,
        20;
    )
    
    println("Done newton. Solution = ", m.x)

    #ax.scatter(m.p, m.x[2], marker="*", s=50, color="red")

    # Now do continuition
    setS!(cmi, s=cmi.s, x=m.x, p=[ m.p ])
    setS!(cmi, s=cmi.ṡ, x=[0.0 ; 0.0], p=[ 0.01 ])
        
    push!(ppp, [m.p, ξ, m.x[1], m.x[2], StommelModel.cal_ψ(m)])
    #ax.scatter3D(m.p, μ, ppp[end][5], marker="*", s=10, color="red")

    for t = 1:10000
        #print(t, "/10000")

        ContMethod1.doContinuition!(cmi)
        setXP!(cmi; s=cmi.s, x=m.x, p=p_container)
        m.p    = p_container[1]

        push!(ppp, [m.p, ξ, m.x[1], m.x[2], StommelModel.cal_ψ(m)])

        if !( p_rng[1] <= m.p <= p_rng[2] ) || isnan(m.p)
            break
        end

    end
    println()

    println("Plot")
    
    global xxx = zeros(Float64, length(ppp))
    global yyy = copy(xxx)
    global zzz = copy(xxx)

    for i=1:length(ppp)
        xxx[i] = ppp[i][1]
        yyy[i] = ppp[i][2]
        zzz[i] = ppp[i][5]
    end
     

    ax.plot(xxx, yyy, zzz, color="black")
    println("End Plot")


end

plt.show()
