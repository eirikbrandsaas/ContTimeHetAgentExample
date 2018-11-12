using Parameters, LinearAlgebra, SparseArrays, Plots

include("cont_time_model.jl")
res=Solve_model(delta=1000.0, imp=1, r = 0.01) #solve the model. flip imp=0 if you're a masochist
prim = Primitives() #initialize primtives

#graphing
#value functions
Plots.plot(prim.a_grid, res.val_func, legend=:bottomright, label=["Low Income" "High Income"],
xlabel="Assets", lw=2)
Plots.savefig("cont_time_val_functions.png")

#stationary distribution
Plots.plot(prim.a_grid, res.stat_dist, ylims=(0, 2), label=["Low Income" "High Income"],
xlabel="Assets", ylabel="Conditional Distribution", lw=2)
Plots.savefig("cont_time_stat_dist.png")

#asset drift
Plots.plot(prim.a_grid, [res.adot, zeros(prim.Na)],  label=["Low Income" "High Income" "Zero Line"], xlabel="Assets",
ylabel="Change in Assets", lw=2, ls = [:solid :solid :dash])
Plots.savefig("cont_time_a_dot.png")


#####
