using ElectrochemicalKinetics
# using GLMakie
using CairoMakie
using DelimitedFiles
using DataFrames
using Interpolations

bv = ButlerVolmer(300, 0.5)
T_vals = 250:10:590

# x = Float32[]
# y = Float32[]
# z = Float32[]

# for T in T_vals
#     I_step = round(300000 * T^(-2), digits=1)
#     println(T, " ", I_step)
#     pbs, I = phase_diagram(bv, muoA=0, muoB=0, I_max=1000, I_step=I_step, T=T, start_guess=[0.01, 0.99], tol=1e-2)
#     T_vec = T .* ones(length(I))
#     append!(x, pbs)
#     append!(y, T_vec)
#     append!(z, I)
# end

# open("./data/bv3d.txt","w") do io
#     writedlm(io, [x y z])
# end

data = readdlm("./data/bv3d.txt")
df = DataFrame(data, [:x, :T, :I])

# to get a draggable scatterplot with GLMakie...
# f = Figure()
# a = Axis3(f[1,1], ylabel="T", zlabel="I")
# limits!(a, (0,1), (250, 600), (0, 1000))
# scatter!(a, df.x, df.T, df.I)
# f

# to get `surface` to work, construct 2D interpolants to individual maps and then build a regular grid in x, T and call interpolant to get I
# this also required a bit of manual cleanup of the bv3d.txt file to get rid of nonmonotonicities, but not too hard to find with the `issorted` function...
regular_x = 0.0:0.03:0.99
regular_T = T_vals
grid_I = Matrix{Union{Float64,Missing}}(undef, 34, 35)

for i in 1:35
    T = regular_T[i]
    subset = df[df.T .== T, :]
    # interpolate on an irregular grid, give zero outside provided domain
    interp = extrapolate(interpolate((subset.x,), subset.I, Gridded(Linear())), 0)
    grid_I[:, i] .= interp(regular_x)
end

f = Figure(figure_padding= (10,40,10,10))
ax = Axis3(f[1,1], azimuth=3pi/4, elevation=pi/12, aspect=(10, 8, 5), ylabel="T [K]", zlabel="I [arb]")
wireframe!(regular_x, regular_T, grid_I, linewidth=2)
zlims!(ax, (0,1000))

save("3d.png", f)
f