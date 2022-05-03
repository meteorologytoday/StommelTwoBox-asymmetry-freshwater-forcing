include("detect_ranges.jl")

using Formatting

println("Loading JLD2")
using JLD2

println("Loading data...")
datum = [ JLD2.load("output_$i.jld2") for i=0:2]


select_ξ = [
    [0,],
    [0, -0.8, 0.8],
    [0, -0.8, 0.8],
]

function find_ξ_idx(ξ_vec, ξ)
    
    dist = abs.(ξ_vec .- ξ)
    _, idx = findmin(dist)
    
    return idx
 
end




println("Loading PyPlot")
using PyPlot
plt = PyPlot
println("done")

fig, ax = plt.subplots(1, length(datum), figsize=(3 * length(datum), 3), constrained_layout=true)


for (i, data) in enumerate(datum)

    println("Data i=$i")

    ξs = select_ξ[i]

    println("selected: ", ξs)

    for (j, ξ) in enumerate(ξs)

        ξ_vec = data["ξ_vec"]
        ξ_idx = find_ξ_idx(ξ_vec, ξ)
        println("Selected ξ = $ξ => idx = $ξ_idx")
        
        d = data["scans"][ξ_idx]

        vals, rngs = detectRanges(d["stability"])
        
        for (k, rng) in enumerate(rngs)

            args = Dict()
            if k==1
                args[:label] = format("\$ \\xi = {:.2f} \$", ξ)
            else
                args[:label] = nothing
            end

            println("rng: ", rng)
            args[:linestyle] = (vals[k] == 0) ? "dashed" : "solid"
            args[:color] = ["black", "red", "blue", "green"][j]
            ax[i].plot(d["p"][rng], d["q"][rng]; args...)

        end
    end

    ax[i].set_title("\$ \\Theta_{$i} \$")
    ax[i].legend()

    ax[i].set_xlabel("p") 
    ax[i].set_ylabel("q") 
end




plt.show()



