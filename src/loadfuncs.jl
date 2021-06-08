using Distributed
using UnicodePlots

@everywhere using Distributions
@everywhere using RCall
# @everywhere using HDF5
@everywhere using JLD2
@everywhere using DataFrames
@everywhere using LinearAlgebra
@everywhere using SharedArrays
# @everywhere using EcologicalNetworks



if homedir() == "/home/z840"
  @everywhere include("$(homedir())/herbforaging/src/trait_and_rate_functions.jl");
  @everywhere include("$(homedir())/herbforaging/src/withindaysim_singleres.jl");
  @everywhere include("$(homedir())/herbforaging/src/acrossdaysim_singleres.jl");
  @everywhere include("$(homedir())/herbforaging/src/smartpath.jl");
else
  @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/trait_and_rate_functions.jl");
  @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/withindaysim_singleres.jl");
  @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/acrossdaysim_singleres.jl");
  @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/smartpath.jl");
end
