using ElectrochemicalKinetics
using DelimitedFiles
using Serialization
using DataFrames
using CairoMakie

λ = .224

# obviously end member that's "all vacancies" has µ=0
# we can choose the other end to be 0 (effectively setting a reference potential)
muoA = 0
muoB = 0

# prefactors chosen such that they all coincide at low overpotential
bv = ButlerVolmer(100, 0.5)
m = Marcus(900, λ) # max current ~900
amhc = AsymptoticMarcusHushChidsey(15000, λ)
# mhc = MarcusHushChidsey(15000, λ)

models = [bv, m, amhc]
# plot_models(models) # if you want to check Tafel plots
T_vals = [250, 300, 450]
Ω_vals = [0.075, 0.1] # "within the range" for interfaces between Li metal and SEI materials

T_list = []
Ω_list = []
pbs_list = []
step_list = []

for Ω in Ω_vals
    for T in T_vals
        push!(T_list, T)
        push!(Ω_list, Ω)
        pbs0 = []
        pbs0 = find_phase_boundaries(0, bv, Ω=Ω, muoA=muoA, muoB=muoB, T=T, guess=[0.01, 0.99])
        push!(pbs_list, pbs0)
    end
end

df_start = DataFrame([T_list, Ω_list, pbs_list], [:T, :Ω, :pb0])

model_type_list = []
T_list = []
Ω_list = []
data_list = []

model_text = Dict(:ButlerVolmer => "Butler-Volmer", :Marcus => "Marcus", :AsymptoticMarcusHushChidsey => "asymptotic MHC")

for row in eachrow(df_start)
    println(row)
    start_guess = row.pb0
    if !isapprox(start_guess[1], start_guess[2]) && !any(isnan.(start_guess))
        for m in models
            model_type = nameof(typeof(m))
            print(model_type)
            # figure out I_step, I_max
            pbs1 = find_phase_boundaries(1, m, T=row.T, Ω=row.Ω, muoA=muoA, muoB=muoB, guess=start_guess)
            steps = abs.(pbs1 .- row.pb0)
            I_step = 1
            I_max = 250
            if any(steps.<0.001)
                I_max = 700
            end

            # build phase diagram
            pbs, I = phase_diagram(m, start_guess=start_guess, T=row.T, Ω=row.Ω, muoA=muoA, muoB=muoB, I_step=I_step, I_max=I_max, warn=false, tol=1e-2)
        
            push!(model_type_list, model_text[model_type])
            push!(T_list, row.T)
            push!(Ω_list, row.Ω)
            push!(data_list, [pbs I])
        end
    end
end

# this is not the most elegant way to store this data but it works...
df_full = DataFrame([model_type_list, T_list, Ω_list, data_list], [:model, :T, :Ω, :pb_data])

convert_data(mat) = mat[:, 1], mat[:,2]

theme = Theme(fontsize = 22,
            linewidth = 4,
            font = "Noto")
set_theme!(theme)
tls = 20
lw = 3

f = Figure(resolution=(900, 900))
g = GridLayout()
T_inds = 1:length(T_vals)
Ω_inds = 1:length(Ω_vals)
ax_args = Dict(:xgridvisible => false,
    :ygridvisible => false,
    :xticklabelsize => tls,
    :yticklabelsize => tls)
axes = [Axis(f; ax_args...) for i in 1:length(Ω_inds)+1, j in T_inds]
g[1:length(Ω_inds)+1, T_inds] = axes
f.layout[1,1] = g

for T_ind in T_inds
    for Ω_ind in Ω_inds
        T = T_vals[T_ind]
        Ω = Ω_vals[Ω_ind]
        df = subset(df_full, :T=>x->x.==T, :Ω=>x->x.==Ω)
        for r in eachrow(df)
            data = convert_data(r.pb_data)
            lines!(g[Ω_ind+1, T_ind], data...; label=string(r.model))
        end
    end
end

axes[2,1].yticks = ([0,50,100,150,200], ["0","50","100","150","200"])
axes[3,1].yticks = ([0,100,200,300,400,500], ["0","100","200","300","400", "500"])
xlims!.(axes[2:3,:], Ref((0,1)))
ylims!.(axes[3,:], Ref((0, 530)))
ylims!.(axes[2,:], Ref((0, 210)))

setproperty!.(axes[1,:], :xticklabelsvisible, false)
setproperty!.(axes[:,2:end], :yticklabelsvisible, false)

Label(g[1, 4], "Rate models", tellheight=false, width=5, rotation=3pi/2, padding=(-25,0,0,0), font = "TeX Gyre Heros Bold",)
Label(g[2, 4], "Ω=$(Ω_vals[1])", tellheight=false, width=5, rotation=3pi/2, padding=(-25,0,0,0), font = "TeX Gyre Heros Bold",)
Label(g[3, 4], "Ω=$(Ω_vals[2])", tellheight=false, width=5, rotation=3pi/2, padding=(-25,0,0,0), font = "TeX Gyre Heros Bold",)

for j in 1:length(T_vals)
    axes[1,j].title = "T=$(T_vals[j])"
end

axes[2,1].ylabel = "Current"
axes[3,2].xlabel = "x"
colgap!(g, 30)
rowgap!(g, 15)

axislegend(axes[3,3], position=:ct)

V = 0.001:0.01:0.35
for T_ind in T_inds
    T = T_vals[T_ind]
    lines!(axes[1,T_ind], V, abs.(bv(V, T=T)))
    lines!(axes[1,T_ind], V, abs.(m(V, T=T)))
    lines!(axes[1,T_ind], V, abs.(amhc(V, T=T)))
end

for i in 1:length(T_vals)
    axes[1,i].yscale = log10
    axes[1,i].xticklabelsvisible = true
    axes[1,i].limits = (0,0.34,6,2e5)
end
axes[1,2].xlabel = "V"

text!(axes[2,3], 0.1, 50, text="(uniform mixing at\n all x and I)", justification=:center)

save("T_omega_sweep.png", f)
f
