using DifferentialEquations
import CSV
using DataFrames

function sirModel!(du, u, p, t)
    S, I, R = u
    Beta, Gamma, N = p
    du[1] = -Beta * S * I / N
    du[2] = Beta * S * I / N - Gamma * I
    du[3] = Gamma * I
end

function jlFunc1()
    # Toy initial conditions
    Beta  = 0.3
    Gamma = 0.1
    N     = 1000.0
    S0    = 990.0
    I0    = 10.0
    R0    = 0.0
    tmax  = 160.0

    R_null = Beta / Gamma
    u0    = [S0, I0, R0]
    tspan = (0.0, tmax)
    p     = (Beta, Gamma, N)

    prob = ODEProblem(sirModel!, u0, tspan, p)
    sol  = solve(prob, Tsit5())

    t = sol.t
    S = [u[1] for u in sol.u]
    I = [u[2] for u in sol.u]
    R = [u[3] for u in sol.u]

    out_df = DataFrame(time = t, S = S, I = I, R = R)

    outName = "func1_output.csv"
    CSV.write(outName, out_df)
    faasr_put_file(outName, outName)
end
