"""
Code to generate Figure 3 from the manuscript, which plots phase diagrams for the two integral-based models (with the asymptotic MHC approximation included for comparison) with three different reorganization energies λ, given in units of the thermal energy.
"""

using ElectrochemicalKinetics
using DelimitedFiles
using CairoMakie

kT = 298 * kB
λ₁ = 2 * kT
λ₂ = 5 * kT
λ₃ = 10 * kT

mhckv1 = MarcusHushChidseyDOS(15000, λ₁, "./data/fig3/Li_100_dos.txt")
mhc1 = MarcusHushChidsey(150000, λ₁)
amhc1 = AsymptoticMarcusHushChidsey(150000, λ₁)

mhckv2 = MarcusHushChidseyDOS(15000, λ₂, "./data/fig3/Li_100_dos.txt")
mhc2 = MarcusHushChidsey(150000, λ₂)
amhc2 = AsymptoticMarcusHushChidsey(150000, λ₂)

mhckv3 = MarcusHushChidseyDOS(15000, λ₃, "./data/fig3/Li_100_dos.txt")
mhc3 = MarcusHushChidsey(150000, λ₃)
amhc3 = AsymptoticMarcusHushChidsey(150000, λ₃)

pbs0 = find_phase_boundaries(0, amhc1, guess=[0.01, 0.99])

# code to build phase diagrams takes a long time (~hours in total) to run for these integral models, data is saved in files (with some slight refinement at higher current resolution to smooth out the tops) for plotting

# pbs_mhckv1, I_mhckv1 = phase_diagram(mhckv1, start_guess=pbs0, I_step=250, verbose=true)
# open("data/fig3/pbs_I_mhcd_lambda_2kT.txt", "w") do io
#     writedlm(io, [pbs_mhckv1 I_mhckv1])
# end
data = readdlm("data/fig3/pbs_I_mhcd_lambda_2kT.txt")
pbs_mhckv1 = data[:,1]
I_mhckv1 = data[:,2]

# pbs_mhc1, I_mhc1 = phase_diagram(mhc1, start_guess=pbs0, I_step=250, verbose=true)
# open("data/fig3/pbs_I_mhc_lambda_2kT.txt", "w") do io
#     writedlm(io, [pbs_mhc1 I_mhc1])
# end
data = readdlm("data/fig3/pbs_I_mhc_lambda_2kT.txt")
pbs_mhc1 = data[:,1]
I_mhc1 = data[:,2]

# pbs_mhckv2, I_mhckv2 = phase_diagram(mhckv2, start_guess=pbs0, I_step=100, verbose=true)
# open("data/fig3/pbs_I_mhcd_lambda_5kT.txt", "w") do io
#     writedlm(io, [pbs_mhckv2 I_mhckv2])
# end
data = readdlm("data/fig3/pbs_I_mhcd_lambda_5kT.txt")
pbs_mhckv2 = data[:,1]
I_mhckv2 = data[:,2]

# pbs_mhc2, I_mhc2 = phase_diagram(mhc2, start_guess=pbs0, I_step=200, verbose=true)
# open("data/fig3/pbs_I_mhc_lambda_5kT.txt", "w") do io
#     writedlm(io, [pbs_mhc2 I_mhc2])
# end
data = readdlm("data/fig3/pbs_I_mhc_lambda_5kT.txt")
pbs_mhc2 = data[:,1]
I_mhc2 = data[:,2]

# pbs_mhckv3, I_mhckv3 = phase_diagram(mhckv3, start_guess=pbs0, I_step=50, verbose=true)
# open("data/fig3/pbs_I_mhcd_lambda_10kT.txt", "w") do io
#     writedlm(io, [pbs_mhckv3 I_mhckv3])
# end
data = readdlm("data/fig3/pbs_I_mhcd_lambda_10kT.txt")
pbs_mhckv3 = data[:,1]
I_mhckv3 = data[:,2]

# pbs_mhc3, I_mhc3 = phase_diagram(mhc3, start_guess=pbs0, I_step=50, verbose=true)
# open("data/fig3/pbs_I_mhc_lambda_10kT.txt", "w") do io
#     writedlm(io, [pbs_mhc3 I_mhc3])
# end
data = readdlm("data/fig3/pbs_I_mhc_lambda_10kT.txt")
pbs_mhc3 = data[:,1]
I_mhc3 = data[:,2]

# pbs_amhc1, I_amhc1 = phase_diagram(amhc1, I_step=10, start_guess=pbs0, tol=1e-2)
# open("data/fig3/pbs_I_amhc_lambda_2kT.txt", "w") do io
#     writedlm(io, [pbs_amhc1 I_amhc1])
# end
data = readdlm("data/fig3/pbs_I_amhc_lambda_2kT.txt")
pbs_amhc1 = data[:,1]
I_amhc1 = data[:,2]

# pbs_amhc2, I_amhc2 = phase_diagram(amhc2, I_step=10, start_guess=pbs0, tol=1e-2)
# open("data/fig3/pbs_I_amhc_lambda_5kT.txt", "w") do io
#     writedlm(io, [pbs_amhc2 I_amhc2])
# end
data = readdlm("data/fig3/pbs_I_amhc_lambda_5kT.txt")
pbs_amhc2 = data[:,1]
I_amhc2 = data[:,2]

# pbs_amhc3, I_amhc3 = phase_diagram(amhc3, I_step=10, start_guess=pbs0, tol=1e-2)
# open("data/fig3/pbs_I_amhc_lambda_10kT.txt", "w") do io
#     writedlm(io, [pbs_amhc3 I_amhc3])
# end
data = readdlm("data/fig3/pbs_I_amhc_lambda_10kT.txt")
pbs_amhc3 = data[:,1]
I_amhc3 = data[:,2]

theme = Theme(fontsize = 22,
            linewidth = 4,
            font = "Noto")
set_theme!(theme)

model_colors = Dict(AsymptoticMarcusHushChidsey=>:goldenrod, MarcusHushChidsey=>:tomato2, MarcusHushChidseyDOS=>:steelblue4)

fig = Figure(resolution=(1300,500))
grid = fig[1, 1:2] = GridLayout()

ax_args = Dict(:xlabel=>"x", :xgridvisible=>false, :ygridvisible=>false)
ax1 = Axis(grid[1,1]; title=L"λ=2k_{\text{B}}T", ylabel=L"\textrm{I }[\text{units of }A^{\text{M}}]", ax_args...)
ax2 = Axis(grid[1,2]; title=L"λ=5k_{\text{B}}T", yticklabelsvisible=false, ax_args...)
ax3 = Axis(grid[1,3]; title=L"λ=10k_{\text{B}}T", yticklabelsvisible=false, ax_args...)

# nondimensionalize by the equivalent Marcus model I₀
lines!(ax1, pbs_amhc1, I_amhc1./9000, label="asymptotic MHC", color=model_colors[AsymptoticMarcusHushChidsey])
lines!(ax1, pbs_mhckv1, I_mhckv1./9000, label="MHC+DOS (Li metal)", color=model_colors[MarcusHushChidseyDOS])
lines!(ax1, pbs_mhc1, I_mhc1./9000, label="MHC", color=model_colors[MarcusHushChidsey])

lines!(ax2, pbs_amhc2, I_amhc2./9000, label="asymptotic MHC", color=model_colors[AsymptoticMarcusHushChidsey])
lines!(ax2, pbs_mhckv2, I_mhckv2./9000, label="MHC+DOS (Li metal)", color=model_colors[MarcusHushChidseyDOS])
lines!(ax2, pbs_mhc2, I_mhc2./9000, label="MHC", color=model_colors[MarcusHushChidsey])

lines!(ax3, pbs_amhc3, I_amhc3./9000, label="asymptotic MHC", color=model_colors[AsymptoticMarcusHushChidsey])
lines!(ax3, pbs_mhckv3, I_mhckv3./9000, label="MHC+DOS (Li metal)", color=model_colors[MarcusHushChidseyDOS])
lines!(ax3, pbs_mhc3, I_mhc3./9000, label="MHC", color=model_colors[MarcusHushChidsey])

axislegend(ax3)
ylims!.([ax1, ax2, ax3], Ref((0,7000/9000)))

save("fig3.png", fig)
fig