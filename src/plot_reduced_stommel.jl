using Formatting
println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("done")

using Roots

fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=true)

p = 1.2
ξ = 0.0
μ = 5.0
ν = 1.0

example_coe = (p, ξ, μ, ν)


ys = collect(range(-1, 2, length=301))
_dydt(y, p, ξ, μ, ν) = p - ( 1 + abs( μ * (1 - y) - ν * p * ξ ) ) * y

_dydt_upward(y, p, ξ, μ, ν) = p - ( 1 + ( μ * (1 - y) - ν * p * ξ ) ) * y

dydt(p, ξ, μ, ν) = [ _dydt(y, p, ξ, μ, ν) for y in ys ]


for ξ in [1.0, 2.0, -1.0, 0.0]

    if ξ != 0.0
        ax.plot(ys, dydt(p, ξ, μ, ν), "-", color="#cccccc", zorder=1)
    end
    y_kink = 1 - ν / μ * p * ξ 
    ax.text(y_kink + ((ξ==0) ? 0.1 : 0), _dydt(y_kink, p, ξ, μ, ν) + 0.1, "\$ \\xi = $(format("{:d}", ξ)) \$", color= ( (ξ == 0) ? "black" : "#aaaaaa" ), ha="center", va="bottom")
end

ax.plot([-5, 5], [0, 0], color="black", linewidth=1)
ax.plot([0, 0], [-5, 5], color="black", linewidth=1)
ax.plot(ys, dydt(example_coe...), "k-", zorder=5)

x0 = 0
x1 = 1
y0 = _dydt(x0, example_coe...)
y1 = _dydt(x1, example_coe...)
ax.scatter([x0, x1], [y0, y1], s=50, marker="o", color="red", zorder=10)
ax.annotate("\$ A_1 = \\left(0, p\\right) \$", xy=(x0, y0),  xycoords="data",
            xytext=(0.3, 1.25), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=15,
)

ax.annotate("\$ A_2 = \\left(1 - \\frac{\\nu p \\xi}{\\mu}, \\left(1 + \\frac{\\nu}{\\mu} \\xi \\right) p - 1\\right) \$", xy=(x1, y1),  xycoords="data",
            xytext=(1, 1), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=12
)

#ax.text(x0 - 0.05, y0, "\$ \\left(0, p\\right) \$", va="center", ha="center", fontsize=15)
#ax.text(x1, y1 + 0.25, "\$ \\left(1 + \\frac{\\nu}{\\mu} \\xi \\right) p - 1 \$", va="center", ha="center", fontsize=15)

slope = ( _dydt(0, example_coe...) - _dydt(1, example_coe...)) / (0 - 1)
ax.plot([0, 2], [_dydt(0, example_coe...), _dydt(0, example_coe...) + 2 * slope], color="red", linestyle="dashed")

# Dot steady states

f(y) = _dydt(y, example_coe...)
y_ss    = [ find_zero(f,  init_pt) for init_pt in [0.2, 0.8, 1.1] ]
dydt_ss = f.(y_ss)
ax.scatter(y_ss, dydt_ss, s=50, marker="x", color="red", zorder=20)

ax.plot([y_ss[1], y_ss[1]], [ 0, -0.5 ], color="#888888", linestyle="dashed")
ax.plot([y_ss[3], y_ss[3]], [ 0, -0.5 ], color="#888888", linestyle="dashed")
ax.text(y_ss[1], -0.5, "thermo", color="red", ha="center", va="top")
ax.text(y_ss[3], -0.5, "halo", color="red", ha="center", va="top")

# Dot the lowest point on the parabola
y_low = (1 + μ - ν * p * ξ) / (2 * μ)
dydt_low = f(y_low)
ax.scatter([y_low], [dydt_low], s=50, marker="o", color="red", zorder=20)

ax.annotate("\$ A_3 = \\left( \\frac{1 + \\mu - \\nu p \\xi}{2 \\mu} , p - \\frac{\\left(1 + \\mu - \\nu p \\xi\\right)^2}{4 \\mu}\\right) \$", xy=(y_low, dydt_low),  xycoords="data",
            xytext=(y_low, -1.3), textcoords="data",
            arrowprops=Dict("facecolor" => "black", "shrink" => 0.15 , "headwidth" => 5.0, "width" => 0.5),
            horizontalalignment="center", verticalalignment="center", fontsize=15
)




ax.set_ylabel("\$\\mathrm{d} y / \\mathrm{d}\\tau\$")
ax.set_xlabel("\$y\$")

ax.grid()
ax.set_ylim([-2, 1.5])
ax.set_xlim([-0.1, 1.4])

fig.savefig("figure-reduced_stommel.png", dpi=300)

plt.show()

