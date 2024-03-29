"""
Code to generate Figure 1 from the manuscript, which serves as the "explainer" figure for how the phase diagrams are constructed. It plots the thermodynamic and kinetic versions of both the chemical potential and Gibbs free energy, as well as the phase diagram for both the intercalation and deintercalation reaction directions. The Tafel plot is inset into the chemical potential plot. Voltages are normalized by kT, energies by the interaction parameter Ω, and currents by I₀, the prefactor of the Marcus model (and maximum current it can attain).
"""

using ElectrochemicalKinetics
using CairoMakie
using DelimitedFiles

# plotting options
theme = Theme(fontsize = 22,
            linewidth = 4,
            font = "Noto")
set_theme!(theme)
tls = 20

# model/data options
x = 0:0.002:0.999
I₀ = 900
λ = 0.3
model = Marcus(I₀, λ) # inverts at a current of I₀ and voltage of λ
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
f = Figure(resolution=(1250,900))
grid = f[1:6,1] = GridLayout()

# first, chemical potential plot
µ_ax = Axis(grid[1:3,1], # 1:3 is so relative sizes get set well automatically
    xgridvisible = false,
    xtickalign=1,
    xticks = ([0,0.5,x₀,1.0],["0","0.5","x₀","1.0"]),
    xticklabelsvisible = false,
    yticks = ([0, Ω, 2*Ω], ["0", "Ω", "2Ω"]),
    yticklabelsize = tls,
    ygridvisible = false,
    ylabel = "chemical potential"
    )

# shaded region to show integration of µ
band!(µ_ax, x[1:x0_ind], µth_vals[1:x0_ind], µ_vals[1:x0_ind], color=(:green, 0.3))

# µ values (plotted second so they're on top of shading)
lines!(µ_ax, x, µ_vals, label=L"\mu_{\mathrm{kin}}(x;\text{ }I=A^{\text{M}}/40)", color=:orange)
lines!(µ_ax, x, µth_vals, label=L"\mu_{\mathrm{thermo}}(x)", linestyle=:dot, color=:steelblue3)

# tweaks/annotations
limits!(µ_ax, (0,1), (-0.4Ω, 2.3Ω))
axislegend(µ_ax, position=:rt)
arrows!(µ_ax, [x₀, x₀], [µth_vals[x0_ind], µ_vals[x0_ind]-1e-3], [0, 0], [Δµ-2e-3, -Δµ+2e-3], color=:mediumpurple4, linewidth=3, arrowsize=10)
text!(µ_ax, x₀+0.01, Ω/4, text=L"\Delta \mu(x_0)", color=:mediumpurple4)

# add Tafel plot inset
inset_ax = Axis(grid[1:3,1],
    width=Relative(0.44),
    height=Relative(0.5),
    halign=0.24,
    valign=0.94,
    backgroundcolor=:gray95,
    yscale=log10,
    xticks = ([0, kB*T, 2kB*T, 3kB*T, 4kB*T], ["0", L"k_{\text{B}}T", L"2k_{\text{B}}T", L"3k_{\text{B}}T", L"4k_{\text{B}}T"]),
    xticklabelsize=tls-2,
    yticks=([I_calc], [L"\frac{A^{\text{M}}/40}{1-x_0}"]),
    yticklabelsize=tls
)

V = 1e-3:1e-3:0.15
lines!(inset_ax, V, model(V, T=T), color=:lightcoral, label="Rate Relationship")
limits!(inset_ax, (0,4kB*T), (0.01*model.A,model.A))

# tweaks/annotations for Tafel inset...
inset_ax.xlabel="V"
inset_ax.xgridvisible=false
inset_ax.ygridvisible=false

axislegend(inset_ax, position=:rb, labelsize=16)

arrows!(inset_ax, [0, Δµ], [I_calc, I_calc], [Δµ-1e-3, -Δµ+1e-3], [0,0], color=:mediumpurple4, linewidth=3, arrowsize=10)

text!(inset_ax, 0.005, 1.04*I_calc, text=L"\Delta \mu(x_0)", color=:mediumpurple4)

textx = 0.33
texty = -0.01
text!(µ_ax, textx, texty, text=L"\Delta g(x_0)", color=:green)
arrows!(µ_ax, [textx+0.13], [texty+0.01], [0.13], [0.012], color=:green4, linewidth=2)

# next, Gibbs free energy subplot
g_ax = Axis(grid[4:6,1], 
    xticks = ([0,0.5,x₀,1.0],["0","0.5","x₀","1.0"]),
    xticklabelsize = tls,
    yticks = ([0, Ω/5, 2Ω/5], ["0", "Ω/5", "2Ω/5"]),
    ylabel = "molar Gibbs free energy",
    yticklabelsize = tls,
    xgridvisible = false,
    ygridvisible = false
    )

limits!(g_ax, (0,1), (0.17Ω, 0.47Ω))

lines!(g_ax, x[1:end-1], g.(x[1:end-1]), label=L"g_{\mathrm{kin}}(x;\text{ }I=A^{\text{M}}/40)", color=:orange)
lines!(g_ax, x, g_thermo(x, Ω=Ω, T=T), label=L"g_{\mathrm{thermo}}(x)", linestyle=:dot, color=:steelblue3)

# again, some tweaks and annotations
axislegend(g_ax, position=:lt)

arrows!(g_ax, [x₀, x₀], [g_thermo(x₀, Ω=Ω, T=T), g(x₀)], [0,0], [Δg-1e-4, -Δg+1e-4], color=:green, linewidth=3, arrowsize=10)
text!(g_ax, x₀+0.015, 0.33Ω, text=L"\Delta g(x_0)", color=:green)


lines!(g_ax, pbs₀, g(pbs₀), linestyle=:dash, color=:gray25, linewidth=2)
text!(g_ax, 0.44, 0.273Ω, text="common tangent", rotation=0.53, color=:gray25, textsize=20)

# now the phase map
pb_grid = f[1:6,2] = GridLayout()

pb_ax_1 = Axis(pb_grid[3:5,1], 
    xticksvisible = false,
    xticklabelsvisible = false,
    yticks = ([0, 0.05I₀, 0.1I₀], [L"0",L"A^{\text{M}}/20",L"A^{\text{M}}/10"]),
    yticklabelsize = tls,
    xgridvisible = false,
    ygridvisible = false,
    ylabel = "Current"
    )

limits!(pb_ax_1, 0,1, 0, I₀/10)

# this is the code that constructed the phase diagram and saved the data. I precomputed and saved it so iterating on the figure design would be faster, but it should all run just fine
# pbs_top, I_top = phase_diagram(model; I_step=0.0005I₀, I_max=0.08I₀, Ω=Ω, warn=false, T=T, lin_thresh=2)
# pbs_bottom, I_bottom = phase_diagram(model; I_start=-0.002I₀, I_step=0.0001I₀, I_max=-0.0001I₀, Ω=Ω, warn=false, T=T, start_guess=[0.058, 0.978], lin_thresh=2)

# pbs = vcat(pbs_bottom[1:Int(length(pbs_bottom)/2)], pbs_top, pbs_bottom[Int(length(pbs_bottom)/2)+1:end])
# I = vcat(I_bottom[1:Int(length(I_bottom)/2)], I_top, I_bottom[Int(length(I_bottom)/2)+1:end])

# open("data/fig1.txt", "w") do io
#     writedlm(io, [pbs I])
# end

# read in the saved data
data = readdlm("./data/fig1.txt")
pbs = data[:,1]
I = data[:,2]

# plot and shade in phase map
lines!(pb_ax_1, pbs, I, color=:darkgray)
band!(pb_ax_1, pbs, zeros(length(pbs)).-I_force, I, color=(:darkgrey, 0.3))

# now for the stripping direction
pb_ax_2 = Axis(pb_grid[6,1], 
    xlabel = "x",
    xticklabelsize = tls,
    yticks = ([-0.002I₀, 0],[L"-0.002A^{\text{M}}",""]),
    yticklabelsize = tls,
    xgridvisible = false,
    ygridvisible = false,
    )

limits!(pb_ax_2, 0,1, -0.002I₀, 0)
lines!(pb_ax_2, pbs, I, color=:darkgray)
band!(pb_ax_2, pbs, zeros(length(pbs)).-I_force, I, color=(:darkgrey, 0.3))

# some overall layout tweaks
rowsize!(pb_grid, 6, Relative(0.08))
rowsize!(grid, 1, Relative(0.2))
rowsize!(grid, 2, Relative(0.2))
rowsize!(pb_grid, 1, 0.15)
colsize!(f.layout, 1, Relative(0.55))
rowgap!(grid, 3, -40)
rowgap!(pb_grid, 5, 10)

# label subparts of figure
Label(grid[3,1], "a)", 
    font = "TeX Gyre Heros Bold",
    textsize=26, 
    tellwidth=false,
    tellheight=false,
    justification=:left,
    padding = (0,670,0,70)
    )

Label(grid[6,1], "b)", 
    font = "TeX Gyre Heros Bold",
    textsize=26, 
    tellwidth=false,
    tellheight=false,
    justification=:left,
    padding = (0,670,0,70)
    )

Label(pb_grid[6,1], "c)", 
    font = "TeX Gyre Heros Bold",
    textsize=26, 
    tellwidth=false,
    tellheight=false,
    justification=:left,
    padding = (0,610,50,0)
    )

# save and display
save("explainer_figure.png", f)
f