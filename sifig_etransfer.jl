using ElectrochemicalKinetics
# using DelimitedFiles
# using Serialization
# using DataFrames
using CairoMakie
using Colors

I_max = 300
colors = sequential_palette(357, 12)[3:end]

α_vals = 0.05:0.1:0.95

theme = Theme(fontsize = 22,
            linewidth = 4,
            font = "Noto")
set_theme!(theme)
tls = 20
lw = 3
f = Figure()
ax = Axis(f[1,1]; xlabel="x", xgridvisible=false, ygridvisible=false, ylabel=L"\textrm{I }[I_0]")
xlims!(ax, (0, 1))
ylims!(ax, (0, (I_max-5)/900))

for i in 1:length(α_vals)
    α = α_vals[i]
    filename = "./data/sifig_etransfer/$(α).txt"
    # println(α)
    # bv = ButlerVolmer(100, α)
    # pbs, I = phase_diagram(bv, I_max=I_max, I_step=0.25)
    # open(filename,"w") do io
    #     writedlm(io, [pbs I])
    # end
    data = readdlm(filename)
    pbs = data[:,1]
    I = data[:,2]
    lines!(ax, pbs, I./900, label="α=$α", color=colors[i])
end

axislegend()

save("sifig_etransfer.png", f)
f