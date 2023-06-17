# computing the nodes x and weights w for generalized Gauss-Laguerre quadrature
# \int_0^{+\infty} x^\alpha e^{-x} f(x) d x
using FastGaussQuadrature
using Test, LinearAlgebra, Random, SpecialFunctions
using MAT

n = 100; # number of nodes
alpha = 11; # exponent
x, w = gausslaguerre(n,alpha);

# Create a dictionary to hold the vectors
dict = Dict("x" => x, "w" => w)

# Save the dictionary to a MAT file
matwrite("gl_nodes_weights.mat", dict)