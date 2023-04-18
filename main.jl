#=
#########  ##   ##  ######  ######  ######
#########  #### ##  ##  ##    ##    ##
#########  ## ####  ##  ##    ##    ####
#########  ##   ##  ######    ##    ######
=#

cd(@__DIR__)
include("libs.jl")
include("SymmetryAdaptedBasis.jl")
##

n = 2 # number of variables
r = 4 # degree

initializeBasis(n,r)
