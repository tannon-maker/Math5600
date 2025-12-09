using Plots
using FundamentalsNumericalComputation

"""
things to do: 
    1) fix the dimension mismathc eror in the initialization of our IVP problem
        ideally, our ivp problem spits out an array of problems that we can plug into our rk4 function
    2) tweak our rk4 algorithm to handle a system of four ODES
        this can be done with a nested for loop in the rk4 function definition, we just need to figure out how to tweak the output arrays properly
    4) plot this bitch

"""

function rk4(ivp, n)
    # Time discretization.
    a, b = ivp.tspan
    h = (b - a) / n
    t = [a + i * h for i in 0:n]

    # Initialize output.
    #lets figure out how to make the output an array of arrays
    u = fill(float(ivp.u0), n + 1)

    # Time stepping.
    for i in 1:n
        k₁ = h * ivp.f(u[i], ivp.p, t[i])
        k₂ = h * ivp.f(u[i] + k₁ / 2, ivp.p, t[i] + h / 2)
        k₃ = h * ivp.f(u[i] + k₂ / 2, ivp.p, t[i] + h / 2)
        k₄ = h * ivp.f(u[i] + k₃, ivp.p, t[i] + h)
        u[i+1] = u[i] + (k₁ + 2(k₂ + k₃) + k₄) / 6
    end
    return t, u
end

function error_analysis(true_values, approx_values)
    error = []
    for i in 1:approx_values
        difference = abs(true_values[i] - approx_values[i])
        append!(error, difference)
    end
    return error
end

#this should work to iterate our rk4 method over the tuple of IVP problems

function rk4_tuple_call(ivp_tuple, n)
    """
    this inputs a tuple of ivp problems, and an n time step size, and outputs a time tuple and an approximate solution tuple 
    """
    #all the inner arrays here should be the same thing its just easier to get one for each for loop iteration
    time_tuple = [[], [], [], []]
    #sol_tuple[i] will hold the ith ivp solutions
    sol_tuple = [[], [], [], []]
    for problem in eachindex(ivp_tuple)
        t, u = rk4(ivp_tuple[problem], n)
        time_tuple[problem] = t
        sol_tuple[problem] = u
    end
    return time_tuple, sol_tuple
end





#change this if needed
time_span = (0.0, 3.0)
G = 9.8
# we can change the values later later
#each vetor = [first pendulum value, second pendulum value]
length = [2, 2]
mass = [2, 2]


function double_pendulum(du, u, p, t)
    p = (mass[1], mass[2], length[1], length[2], G)
    #omega_1 = theta_1 prime
    du[1] = u[1]

    diff = u[3] - u[1]
    denominator_1 = (mass[1] + mass[2]) * length[1] - mass[2] * length[1] * cos(diff) * cos(diff)
    #omega_1 prime
    du[2] = (
        (
            mass[2] * length[1] * u[2] * u[2] * sin(diff) * cos(diff) +
            mass[2] * G * sin(u[3]) * cos(diff) +
            mass[2] * length[2] * u[4] * u[4] * sin(diff) - (mass[1] + mass[2]) * G * sin(u[1])
        ) / denominator_1
    )
    #omega_2 = theta_2 prime
    du[3] = u[4]

    denominator_2 = (length[2] / length[1]) * denominator_1
    #omega_2 prime
    du[4] = (
        (
            -mass[2] * length[2] * u[4] * u[4] * sin(diff) * cos(diff) +
            (mass[1] + mass[2]) * G * sin(u[1]) * cos(diff) -
            (mass[1] + mass[2]) * length[1] * u[2] * u[2] * sin(diff) - (mass[1] + mass[2]) * G * sin(u[3])
        ) / denominator_2
    )
    #fuck around with the return value?
    return nothing
end

#try messing around with the initial conditions
initial_theta_one = 90.0
initial_theta_two = 90.0
intial_vel_one = 0.0
intial_vel_two = 0.0
u0 = [initial_theta_one, intial_vel_one, initial_theta_two, intial_vel_two]

p = (mass[1], mass[2], length[1], length[2], G)
#here is our IVP
#need to figure out how to make IVP a tuple of length four to plug into our rk4 function
IVP = ODEProblem(double_pendulum, deg2rad.(u0), time_span, p)

