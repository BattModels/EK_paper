# first bit cribbed from PkgBenchmark
pkgid = Base.identify_package("ElectrochemicalKinetics")
pkgfile = Base.locate_package(pkgid)
pkgdir = normpath(joinpath(dirname(pkgfile), ".."))
script = joinpath(pkgdir, "benchmark", "bench_suite.jl")

# define suite and pull the tests we need to run
include(script)
op_suite = suite[@tagged "overpotential fitting"]

# run the benchmarks
tune!(op_suite)
results = run(op_suite)