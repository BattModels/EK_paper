"""
Generate the TOC figure. This is mostly copy/pasted, with light modifications, from the scripts for figure 1 and the SI figure.
"""

using ElectrochemicalKinetics
using CairoMakie
using DelimitedFiles
using DataFrames
using Interpolations

# model/data options
x = 0:0.002:0.999
I₀ = 900
λ = 0.3
m = Marcus(I₀, λ) # inverts at a current of I₀ and voltage of λ
bv = ButlerVolmer(I₀/9, 0.5) # aligns with that Marcus model at low overpotential and 300K
model = m
T = 375
Ω = 0.1
I_force = I₀/40
pbs₀ = find_phase_boundaries(I_force, model, Ω=Ω, T=T, warn=false, lin_thresh=2)
x₀ = 0.84
x0_ind = searchsorted(x, x₀).start
I_calc = I_force/(1-x₀)

# compute things
µth_vals = µ_thermo(x, Ω=Ω, T=T)
µ = µ_kinetic(I_force, model, Ω=Ω, T=T)
µ_vals = µ.(x)
Δµ = µ_vals[x0_ind]-µth_vals[x0_ind]
g = g_kinetic(I_force, model, Ω=Ω, T=T)
Δg = g(x₀) - g_thermo(x₀, Ω=Ω, T=T)

# make figure
theme = Theme(fontsize = 24,
            linewidth = 6,
            font = "Noto")
set_theme!(theme)

ts = 24 # font size for when overriding

f = Figure(resolution=100 .* (7, 7))

µ_ax = Axis(f[1:2,1:6],
    xgridvisible = false,
    xtickalign=1,
    xticks = ([],[]),
    xticklabelsvisible = false,
    yticks = ([], []),
    yticklabelsvisible = false,
    ygridvisible = false,
    ylabel = "µ",
    xlabel="x"
    )

# shaded region to show integration of µ
band!(µ_ax, x[1:x0_ind], µth_vals[1:x0_ind], µ_vals[1:x0_ind], color=(:green, 0.3))

# µ values (plotted second so they're on top of shading)
lines!(µ_ax, x, µ_vals, label=L"\mu_{\mathrm{kinetic}}(x; I)", color=:orange)
lines!(µ_ax, x, µth_vals, label=L"\mu_{\mathrm{thermo}}(x)", linestyle=:dot, color=:steelblue3)

# tweaks/annotations
limits!(µ_ax, (0,1), (-0.4Ω, 2.3Ω))
axislegend(µ_ax, position=:rt, labelsize=32)
arrows!(µ_ax, [x₀, x₀], [µth_vals[x0_ind], µ_vals[x0_ind]-1e-3], [0, 0], [Δµ-2e-3, -Δµ+2e-3], color=:mediumpurple4, linewidth=3, arrowsize=10)
text!(µ_ax, x₀+0.01, Ω/4, text=L"\Delta \mu(x_0)", color=:mediumpurple4, textsize=ts)

# add Tafel plot inset
inset_ax = Axis(f[1:2,1:6],
    width=Relative(0.55),
    height=Relative(0.53),
    halign=0.22,
    valign=0.94,
    backgroundcolor=:gray95,
    yscale=log10,
    xticks = ([], []),
    # xticklabelsize=tls-2,
    yticks=([], []),
    ylabel=L"\log(I)",
    ylabelsize=ts,
)

Vmax = 0.2
V = 1e-3:1e-3:Vmax
lines!(inset_ax, V, m(V, T=T), color=:mediumorchid4, label="Marcus")
lines!(inset_ax, V, bv(V, T=T), color=(:teal, 0.3), label="Butler-Volmer")
limits!(inset_ax, (0,Vmax), (0.01*model.A,3*model.A))

# tweaks/annotations for Tafel inset...
inset_ax.xlabel="V"
inset_ax.xgridvisible=false
inset_ax.ygridvisible=false

axislegend(inset_ax, position=:rb)#, labelsize=16)

arrows!(inset_ax, [0, 0.92*Δµ], [I_calc, I_calc], [0.92*Δµ-1e-3, -0.92*Δµ+1e-3], [0,0], color=:mediumpurple4, linewidth=3, arrowsize=10)

text!(inset_ax, 0.003, 1.04*I_calc, text=L"\Delta \mu(x_0)", color=:mediumpurple4, textsize=ts)

textx = 0.34
texty = -0.01
text!(µ_ax, textx, texty, text=L"\Delta g(x_0)", color=:green, textsize=ts)
arrows!(µ_ax, [textx+0.12], [texty+0.01], [0.13], [0.012], color=:green4, linewidth=2)

# 3D phase diagram on the other side
T_vals = 270:10:590
pb_ax_1 = Axis3(f[3,4:6], azimuth=3pi/4, elevation=pi/12, aspect=(10,8,6), xlabel="x", ylabel="T [K]", zlabel="I", xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, zticksvisible=false, zticklabelsvisible=false, xgridvisible=false, ygridvisible=false, zgridvisible=false, xlabeloffset=10, ylabeloffset=10, zlabeloffset=5, title="Marcus")
pb_ax_2 = Axis3(f[3,1:3], azimuth=3pi/4, elevation=pi/12, aspect=(10,8,6), xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, zticksvisible=false, zticklabelsvisible=false, xgridvisible=false, ygridvisible=false, zgridvisible=false, xlabelvisible=false, ylabelvisible=false, zlabelvisible=false, title="Butler-Volmer", titlecolor=(:black,0.5), xspinecolor_1=(:black,0.5), xspinecolor_2=(:black,0.5), xspinecolor_3=(:black,0.5), xlabelcolor=(:black,0.5), yspinecolor_1=(:black,0.5), yspinecolor_2=(:black,0.5), yspinecolor_3=(:black,0.5), ylabelcolor=(:black,0.5), zspinecolor_1=(:black,0.5), zspinecolor_2=(:black,0.5), zspinecolor_3=(:black,0.5), zlabelcolor=(:black,0.5),)
zlims!.([pb_ax_1, pb_ax_2], Ref((0, 1000)))

data = readdlm("./data/sifig_3d/Marcus_clean.txt")
df = DataFrame(data, [:x, :T, :I])
regular_x = 0.0:0.03:0.99
regular_T = T_vals
grid_I = Matrix{Union{Float64,Missing}}(undef, length(regular_x), length(regular_T))
for i in 1:length(T_vals)
    T = regular_T[i]
    # println(T)
    subset = df[df.T .== T, :]
    # interpolate on an irregular grid, give zero outside provided domain
    interp = extrapolate(interpolate((subset.x,), subset.I, Gridded(Linear())), 0)
    grid_I[:, i] .= interp(regular_x)
end
wireframe!(pb_ax_1, regular_x, regular_T, grid_I, linewidth=2, color=:mediumorchid4)

data = readdlm("./data/sifig_3d/ButlerVolmer_clean.txt")
df = DataFrame(data, [:x, :T, :I])
regular_x = 0.0:0.03:0.99
regular_T = T_vals
grid_I = Matrix{Union{Float64,Missing}}(undef, length(regular_x), length(regular_T))
for i in 1:length(T_vals)
    T = regular_T[i]
    # println(T)
    subset = df[df.T .== T, :]
    # interpolate on an irregular grid, give zero outside provided domain
    interp = extrapolate(interpolate((subset.x,), subset.I, Gridded(Linear())), 0)
    grid_I[:, i] .= interp(regular_x)
end
wireframe!(pb_ax_2, regular_x, regular_T, grid_I, linewidth=2, color=(:teal, 0.15))

save("tocfig.png", f)
f