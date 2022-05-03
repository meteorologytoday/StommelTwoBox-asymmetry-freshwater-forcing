include("StommelModel.jl")
include("NewtonMethod.jl")

using .StommelModel
using .NewtonMethod

α = 3.6e3
μ = 7.5^0.5
p_vec = range(0, 2, length=101)

m = Model(α=α, μ=μ, p=p_vec[1])

using PyPlot

fig, ax = subplots(1, 1, constrained_layout=true)
    
PyPlot.show()

ax.set_xlim([p_vec[1], p_vec[end]])
ax.set_ylim([0, 1.5])
ax.grid()
for i=1:2

    _p_vec = (i == 1) ? p_vec : reverse(p_vec)
    c = (i == 1) ? "red" : "blue"
    for p in _p_vec

        m.p = p


        _F_dFdx(x) = ( 
            StommelModel.F(x, m.α, m.μ, m.p),
            StommelModel.dFdx(x, m.α, m.μ, m.p),
        )


        NewtonMethod.doNewtonMethod!(
            m.x,
            _F_dFdx,
            1e-10,
            5,
            20;
        )
        ax.scatter(p, m.x[2], marker="o", s=10, color=c)
        sleep(0.001)
    end

end
