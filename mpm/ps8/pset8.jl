# pset 8

#############################################################################
# Question 1
#############################################################################

using CSV, DataFrames, Optim, Distributions
using Gadfly

dir = "/Users/tombearpark/Documents/princeton/1st_year/term2/" * 
              "ECO518_Metrics2/mpm/excercises/ps8/"

# load and format the data 
df = DataFrame(CSV.File(joinpath(dir, "mroz.csv")))
insertcols!(df, 1, :c => ones(nrow(df)))
Y = convert(Vector, df.part)
X = convert(Matrix, df[!, Not(:part)])

# return the log likelihood for the probit model
function probit_log_lik(beta, X, Y)
	Xhat = X * beta
	-sum(Y .* log.(cdf(Normal(), μ)) .+ (1 .- Y) .* log.(1 .- cdf(Normal(), μ)))
end

# initialise with OLS
beta0 = (X'X)^(-1) * X' * Y

# Run the optimisation
@time MLE = optimize(beta -> probit_log_lik(beta, X, Y), beta0);
Optim.minimizer(MLE)

#############################################################################
# Question  3
#############################################################################

df = DataFrame(CSV.File(joinpath(dir, "fish.csv")))
rename!(df, Symbol.(replace.(string.(names(df)), Ref(r" "=>""))))

Y = convert(Vector, df.logq)
X = convert(Vector, df.logp)
Z = convert(Matrix, df[:,[:mixed, :stormy]])


function estimate_2SLS(X, Y, Z)
    
    # calculate beta
    beta = (X' * Z * (Z' * Z)^(-1) * Z' * X)^(-1) * X' * Z * (Z' * Z)^(-1) * Z' * Y
    
    # estimate variance 
    u = Y - X * beta
    omega = Matrix()
