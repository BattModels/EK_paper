using ElectrochemicalKinetics
# using GLMakie
using CairoMakie
using DelimitedFiles
using DataFrames
using Interpolations

bv = ButlerVolmer(300, 0.5)
m = Marcus(2700, λ) # max current ~900
amhc = AsymptoticMarcusHushChidsey(45000, λ)
models = [bv, m, amhc]
model_colors = Dict(bv=>:teal, m=>:mediumorchid4, amhc=>:goldenrod)

T_vals = 270:10:590

# for m in models
#     println(nameof(typeof(m)))
#     x = Float32[]
#     y = Float32[]
#     z = Float32[]

#     # NB that this data required a bit of manual cleaning afterwards to make sure it was monotonic for the interpolants below, the results of this cleaning are saved in the repo
#     for T in T_vals
#         I_step = round(300000 * T^(-2), digits=1)
#         println(T, " ", I_step)
#         pbs, I = phase_diagram(m, muoA=0, muoB=0, I_max=1000, I_step=I_step, T=T, start_guess=[0.01, 0.99], tol=1e-2)
#         T_vec = T .* ones(length(I))
#         append!(x, pbs)
#         append!(y, T_vec)
#         append!(z, I)
#     end

#     open("./data/sifig_3d/$(nameof(typeof(m))).txt","w") do io
#         writedlm(io, [x y z])
#     end
# end

f = Figure(resolution=(1500, 400), figure_padding= (10,50,15,15))
ax_args = Dict(:azimuth=>3pi/4, :elevation=>pi/12, :aspect=>(10,8,6), :ylabel=>"T [K]", :zlabel=>"I [arb]")
ax1 = Axis3(f[1,1]; title="Butler-Volmer", ylabel=L"\textrm{I }[I_0]", ax_args...)
ax2 = Axis3(f[1,2]; title="Marcus", ax_args...)
ax3 = Axis3(f[1,3]; title="asymptotic MHC", ax_args...)
axes = [ax1, ax2, ax3]
zlims!.(axes, Ref((0,1000)))

for (m, ax) in zip(models, axes)
    println(nameof(typeof(m)))
    data = readdlm("./data/sifig_3d/$(nameof(typeof(m)))_clean.txt")
    df = DataFrame(data, [:x, :T, :I])

    # to get a draggable scatterplot with GLMakie...
    # f = Figure()
    # a = Axis3(f[1,1], ylabel="T", zlabel="I")
    # limits!(a, (0,1), (250, 600), (0, 1000))
    # scatter!(a, df.x, df.T, df.I)
    # f

    # construct 2D interpolants to individual maps and then build a regular grid in x, T and call interpolant to get I
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

    wireframe!(ax, regular_x, regular_T, grid_I, linewidth=2, color=model_colors[m])
end

save("3d.png", f)
f