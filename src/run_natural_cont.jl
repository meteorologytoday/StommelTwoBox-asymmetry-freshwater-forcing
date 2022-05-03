include("StommelModel.jl")
include("NewtonMethod.jl")
include("ContMethod.jl")

using LinearAlgebra
using LaTeXStrings
using Formatting
using .StommelModel
using .NewtonMethod
using .ContMethod1

c = 3.6e5
μ = 7.0
ν = 1.0
ξ = 0.0

p_vec = collect(range(0.0, 2.0, length=201))

m = Model(c=c, μ=μ, ν=1.0, p=p_vec[1], ξ=0.0)

m.x .= [ 1.0,  0.0 ]
ppp = []

for p in p_vec

    println("Now do p = ", p)
    
    m.p  = p

    _F_dFdx = function(x)
        return (
            StommelModel.F(m; x=x, p=m.p),
            StommelModel.dFdx(m; x=x, p=m.p),
        )
    end



    NewtonMethod.doNewtonMethod!(
        m.x,
        _F_dFdx,
        1e-10,
        5,
        20;
    )
    
    push!(ppp, [m.p, m.μ, m.x[1], m.x[2], StommelModel.cal_ψ(m; x=m.x, p=m.p)])

end

for p in reverse(p_vec)

    println("Now do p = ", p)

    m.p  = p

    _F_dFdx = function(x)
        return (
            StommelModel.F(m; x=x, p=m.p),
            StommelModel.dFdx(m; x=x, p=m.p),
        )
    end



    NewtonMethod.doNewtonMethod!(
        m.x,
        _F_dFdx,
        1e-10,
        5,
        20;
        verbose = true
    )
    
    push!(ppp, [m.p, m.μ, m.x[1], m.x[2], StommelModel.cal_ψ(m; x=m.x, p=m.p)])

end


println("F : ", StommelModel.F(m; x=m.x, p=m.p))

global ps = zeros(Float64, length(ppp))
global ys = copy(ps)
global ψs = copy(ps)

for i=1:length(ppp)
    ps[i] = ppp[i][1]
    ys[i] = ppp[i][4]
    ψs[i] = ppp[i][5]
end
 

println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("Done")

fig, ax = plt.subplots(1, 1, constrained_layout=true)

ax.set_xlabel(L"$ p $")
ax.set_ylabel(L"$ y $")
ax.grid()




ax.plot(ps, ys, color="black")

plt.show()
