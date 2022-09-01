############################################################################################
# Constants
"RNG for sampling based solutions"
const _RNG = Random.Xoshiro(123)    # shorthand
Random.seed!(_RNG, 1)

"Tolerance for stochastic solutions"
const _TOL = 1.0e-6

"Number of samples"
N = 10^3
