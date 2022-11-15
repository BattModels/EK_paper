using ElectrochemicalKinetics
using DelimitedFiles
using CairoMakie

λ₁ = 2 * 298 * kB
λ₂ = 5 * 298 * kB
λ₃ = 10 * 298 * kB

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

# pbs_amhc1, I_amhc1 = phase_diagram(amhc1, I_step=20, start_guess=pbs0, tol=1e-2)
# pbs_amhc2, I_amhc2 = phase_diagram(amhc2, I_step=20, start_guess=pbs0, tol=1e-2)
# pbs_amhc3, I_amhc3 = phase_diagram(amhc3, I_step=20, start_guess=pbs0, tol=1e-2)

theme = Theme(fontsize = 22,
            linewidth = 4,
            font = "Noto")
set_theme!(theme)

fig = Figure(resolution=(1300,500))
grid = fig[1, 1:2] = GridLayout()

ax_args = Dict(:xlabel=>"x", :xgridvisible=>false, :ygridvisible=>false)
ax1 = Axis(grid[1,1]; title=L"λ=2k_{\text{B}}T", ylabel="I", ax_args...)
ax2 = Axis(grid[1,2]; title=L"λ=5k_{\text{B}}T", yticklabelsvisible=false, ax_args...)
ax3 = Axis(grid[1,3]; title=L"λ=10k_{\text{B}}T", yticklabelsvisible=false, ax_args...)

lines!(ax1, pbs_amhc1, I_amhc1, label="asymptotic MHC")
lines!(ax1, pbs_mhckv1, I_mhckv1, label="MHC+DOS (Li metal)")
lines!(ax1, pbs_mhc1, I_mhc1, label="MHC")

lines!(ax2, pbs_amhc2, I_amhc2, label="asymptotic MHC")
lines!(ax2, pbs_mhckv2, I_mhckv2, label="MHC+DOS (Li metal)")
lines!(ax2, pbs_mhc2, I_mhc2, label="MHC")

lines!(ax3, pbs_amhc3, I_amhc3, label="asymptotic MHC")
lines!(ax3, pbs_mhckv3, I_mhckv3, label="MHC+DOS (Li metal)")
lines!(ax3, pbs_mhc3, I_mhc3, label="MHC")

axislegend(ax3)
ylims!.([ax1, ax2, ax3], Ref((0,7000)))

fig