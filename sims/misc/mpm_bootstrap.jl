# Bootstrap - following the slides 

using DataFrames, CSV, Distributions, Random
using StatsPlots, VegaLite
Random.seed!(123)

# Simulate a population
d = Normal(2.0, 3.0)
x = rand(d, 10000)

# Draw an observed sample 
N = 100
s = sample(x, N, replace = true)
density(s, label = "Sample Density")
m = mean(s)
vline!([m], label = "Observed Mean")
T = sqrt(length(s)) * (m - 2)
sd = std(s) 

# Parametric Bootstrap
B = 1000

function normal_para_boot(B, m, sd)
    b_df = DataFrame(draw = 0, value = 0.0)
    for b in 1:B 
        b_s = DataFrame(value = mean(rand(Normal(m, sd), length(s))))
        b_s[!, :draw] .= b
        append!(b_df, b_s)
    end
    delete!(b_df, 1)
    return(b_df)
end

@time b_df = normal_para_boot(B, m, sd)
density(b_df.value)

# Compare to non-parametric version
function boot(s, B)
    b2_df = DataFrame(draw = 0, value = 0.0)
    for b in 1:B
        s_b = sample(s, length(s), replace = true)
        b_s = DataFrame(value = mean(s_b))
        b_s[!, :draw] .= b
        append!(b2_df, b_s)
    end
    delete!(b2_df, 1)
    return(b2_df)
end
@time b2_df = boot(s, B)

density!(b2_df.value)

# Generate some data for regression bootstrap
data = DataFrame(X = s, 
                 Y = s * 3 + rand(Normal(0, 1), length(s)))

function test_loop()
    for i in 1:1000
        