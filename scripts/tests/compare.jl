using Plots
using Serialization


# h from numerical scheme
nx  = 301
lx  = 30
x_num = LinRange(-lx/2, lx/2, nx)
h_num = deserialize("h_huppert_2D_49.jls")
mid_y = nx ÷ 2
h_quer = h_num[:, mid_y]



# h from self-similar solution
lx = 10
nx = 1000
x_sim = LinRange(0, lx, nx)
h_sim = deserialize("h_huppert_selfsim.jls")

# Create the initial plot
p = plot(x_sim, h_sim, label="Self-Similar Solution", color=:blue, linewidth=2, size=(1000, 500))

# 2. overlay second graph
p = plot!(x_num[mid_y : end], h_quer[mid_y : end], label="Numerical Solution", color=:red, linewidth=2, linestyle=:dash)

display(p)
savefig("compare_fig.png")

