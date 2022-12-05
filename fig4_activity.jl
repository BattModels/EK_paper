using ElectrochemicalKinetics
using CairoMakie
using DelimitedFiles

bv = ButlerVolmer(100, 0.5)

I_step = 1
I_max = 330
T = 330

a1(x) = 1 .- x
a2(x) = (1 .- x).^2
a3(x) = (1 .- x).^(2/3)
activity_fs = [a3, a1, a2] # so they'll plot in a visually easier order
a_fn_strs = Dict(a1=>"a1", a2=>"a2", a3=>"a3")
a_labels = Dict(a1=>L"a(x) = 1 - x", a2=>L"a(x) = (1-x)^2", a3=>L"a(x) = (1-x)^{2/3}")

# generate phase maps
pbs0 = find_phase_boundaries(0, bv, T=T, guess=[0.01, 0.99])

for a in activity_fs
    println(a_fn_strs[a])
    pbs, I = phase_diagram(bv, start_guess=pbs0, T=T, I_step=I_step, I_max=I_max, activity_function_o=a, activity_function_r=a, warn=false)
    open("data/fig4/$(a_fn_strs[a]).txt", "w") do io
        writedlm(io, [pbs I])
    end
end

# plotting...
theme = Theme(fontsize = 22,
            linewidth = 4,
            font = "Noto")
set_theme!(theme)

f = Figure(resolution=(500,450))
ax = Axis(f[1,1]; xlabel="x", xgridvisible=false, ygridvisible=false, ylabel=L"\textrm{I }[I_0]")
ylims!(ax, (0, (I_max-1)/900))

for ac in activity_fs
    fname = "data/fig4/$(a_fn_strs[ac]).txt"
    data = readdlm(fname)
    pbs = data[:,1]
    I = data[:,2]
    lines!(ax, pbs, I./900, label=a_labels[ac])
end

axislegend(ax)

save("fig4.png", f)
f