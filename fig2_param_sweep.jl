using ElectrochemicalKinetics
using DelimitedFiles
using Serialization
using DataFrames
using CairoMakie

λ = .3

# prefactors chosen such that they all coincide at low overpotential at 300K
bv = ButlerVolmer(100, 0.5)
m = Marcus(900, λ) # max current ~900
amhc = AsymptoticMarcusHushChidsey(15000, λ)

models = [bv, m, amhc]
T_vals = [275, 350, 500]
Ω_vals = [0.075, 0.1]

I_step = 0.2
I_max = 250

# start out by getting thermodynamic phase boundaries, this will speed up beginning of phase diagram construction later
T_list = []
Ω_list = []
pbs_list = []

for Ω in Ω_vals
    for T in T_vals
        push!(T_list, T)
        push!(Ω_list, Ω)
        pbs0 = []
        pbs0 = find_phase_boundaries(0, bv, Ω=Ω, T=T, guess=[0.01, 0.99])
        push!(pbs_list, pbs0)
    end
end

df_start = DataFrame([T_list, Ω_list, pbs_list], [:T, :Ω, :pb0])

# now we can actually build the phase maps

# this is just so the legend in the plot is nicer-looking later
model_text = Dict(:ButlerVolmer => "Butler-Volmer", :Marcus => "Marcus", :AsymptoticMarcusHushChidsey => "asymptotic MHC")

# iterate through the parameters and build phase maps
# this is the expensive bit...I'm saving out the data so we can skip to plotting, but leaving the code that generated it here
# for row in eachrow(df_start)
#     println(row)
#     start_guess = row.pb0
#     if !isapprox(start_guess[1], start_guess[2]) && !any(isnan.(start_guess))
#         for m in models
#             model_type = nameof(typeof(m))
#             print(model_type)

#             # build phase diagram
#             pbs, I = phase_diagram(m, start_guess=start_guess.+[0.005,-0.005], T=row.T, Ω=row.Ω, I_step=I_step, I_max=I_max, warn=false, tol=1e-2)
        
#             open("data/fig2/$(nameof(typeof(m)))_T$(row.T)_Omega$(row.Ω).txt", "w") do io
#                 writedlm(io, [pbs I])
#             end
#         end
#     end
# end

# set up for plotting...
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
        # df = subset(df_full, :T=>x->x.==T, :Ω=>x->x.==Ω)
        # for r in eachrow(df)
        #     lines!(g[Ω_ind+1, T_ind], r.pbs, r.I; label=string(r.model))
        # end
        for m in models
            fname = "data/fig2/$(nameof(typeof(m)))_T$(T)_Omega$(Ω).txt"
            if isfile(fname)
                data = readdlm("data/fig2/$(nameof(typeof(m)))_T$(T)_Omega$(Ω).txt")
                pbs = data[:,1]
                I = data[:,2]
                lines!(g[Ω_ind+1, T_ind], pbs, I, label=model_text[nameof(typeof(m))])
            end
        end
    end
end

# y axes should line up across each row, then we can just label the leftmost one
# axes[2,1].yticks = ([0,50,100,150], ["0","50","100","150"])
# axes[3,1].yticks = ([0,100,200,300], ["0","100","200","300"])
xlims!.(axes[2:3,:], Ref((0,1)))
ylims!.(axes[3,:], Ref((0, 220)))
ylims!.(axes[2,:], Ref((0, 110)))
setproperty!.(axes[:,2:end], :yticklabelsvisible, false)

# add parameter and axis labels...
Label(g[1, 4], "Rate models", tellheight=false, width=5, rotation=3pi/2, padding=(-25,0,0,0), font = "TeX Gyre Heros Bold",)
Label(g[2, 4], "Ω=$(Ω_vals[1])", tellheight=false, width=5, rotation=3pi/2, padding=(-25,0,0,0), font = "TeX Gyre Heros Bold",)
Label(g[3, 4], "Ω=$(Ω_vals[2])", tellheight=false, width=5, rotation=3pi/2, padding=(-25,0,0,0), font = "TeX Gyre Heros Bold",)

for j in 1:length(T_vals)
    axes[1,j].title = "T=$(T_vals[j])"
end

axes[2,1].ylabel = "Current"
axes[3,2].xlabel = "x"

# tweak spacing
colgap!(g, 30)
rowgap!(g, 15)

axislegend(axes[3,3], position=:ct)

# now add the Tafel plots along the top row
V = 0.001:0.01:0.35
for T_ind in T_inds
    T = T_vals[T_ind]
    lines!(axes[1,T_ind], V, abs.(bv(V, T=T)))
    lines!(axes[1,T_ind], V, abs.(m(V, T=T)))
    lines!(axes[1,T_ind], V, abs.(amhc(V, T=T)))
end

for i in 1:length(T_vals)
    axes[1,i].yscale = log10
    axes[1,i].limits = (0,0.34,6,2e5)
end
axes[1,2].xlabel = "V"

text!(axes[2,3], 0.1, 50, text="(uniform mixing at\n all x and I)", justification=:center)

save("fig2.png", f)
f