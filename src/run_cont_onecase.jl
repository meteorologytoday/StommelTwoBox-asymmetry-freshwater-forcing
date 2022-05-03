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
ξ = 0.0
p_rng = [0.0, 5.0]

m = Model(c=c, μ=μ, ν=ν, p=p_rng[1], ξ=ξ)

m.x .= [ 1.0,  0.0 ]


println("Loading PyPlot")
using PyPlot
plt = PyPlot
plt.ion()
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

function _F_dFdx_dFdp(x, p)
    return (
        StommelModel.F(m; x=x, p=p[1]),
        StommelModel.dFdx(m; x=x, p=p[1]),
        StommelModel.dFdp(m; x=x, p=p[1]),
    )
end

scale = 5e-3
cmi = ContMethod1.CMInfo(
    Nx = length(m.x),
    F_dFdx_dFdp = _F_dFdx_dFdp,
    mx   = [scale, scale],
    mp   = scale,
    nwt_min_iter = 5,
    nwt_max_iter = 20,
    res  = 1e-10,
    ds_exp_lb = -20.0,
    newton_callback = function(s)
    end,
    skip_unconverge = false, 
)

ax.set_xlabel(L"p")
ax.set_ylabel(L"ψ")
ax.grid()

        
# First, find a steady state
#_F_dFdx_dFdp(x, m.p)[1:2]
p_container = [ m.p ] 
m.p  = p_rng[1]
m.x .= 0

ds = Dict(
    :dx => copy(m.x),
    :dp => [ 0.0 ],
)


_F_dFdx = function(x)
    return (
        StommelModel.F(m; x=x),
        StommelModel.dFdx(m; x=x),
    )
end

ppp = []

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
    
push!(ppp, [m.p, μ, m.x[1], m.x[2], StommelModel.cal_ψ(m), [ ds[:dx][1], ds[:dx][2], ds[:dp][1]] ])
#ax.scatter3D(m.p, μ, ppp[end][5], marker="*", s=10, color="red")

for t = 1:10000
    println(t, "/10000")

    ContMethod1.doContinuition!(cmi)
    setXP!(cmi; s=cmi.s, x=m.x, p=p_container)
    m.p    = p_container[1]
    
    setXP!(cmi; s=cmi.ṡ, x=ds[:dx], p=ds[:dp])

    push!(ppp, [m.p, μ, m.x[1], m.x[2], StommelModel.cal_ψ(m), [ ds[:dx][1], ds[:dx][2], ds[:dp][1]] ])

    if !( p_rng[1] <= m.p <= p_rng[2] ) || isnan(m.p)
        break
    end

end
println()

println("Plot")

global ps = zeros(Float64, length(ppp))
global ψs = copy(ps)
global ys = copy(ps)
global dys = copy(ps)
global dps = copy(ps)


for i=1:length(ppp)
    ps[i] = ppp[i][1]
    ys[i] = ppp[i][4]
    ψs[i] = ppp[i][5]
    dys[i] = ppp[i][6][2]
    dps[i] = ppp[i][6][3]
end
 

ax.plot(ps, ψs, color="black")

#plt.show(block=false)

#=
for i=2:length(ppp)

    if mod(i, 5) != 0
        continue
    end

    println("Plot vector of i = $i")

    l = (dps[i]^2 + dys[i]^2)^0.5
    dp = dps[i] / l
    dy = dys[i] / l
    ax.plot([xxx[i], xxx[i]+dp], [yyy[i], yyy[i]+dy])
    sleep(0.1)
end
sleep(500)
=# 

readline()
println("End Plot")


