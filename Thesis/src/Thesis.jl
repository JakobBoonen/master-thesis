module Thesis
    using Makie, DelimitedFiles, LinearRegression, LsqFit, LinearAlgebra
    using GLMakie, SparseArrays, FFTW
    using SpecialFunctions

    # using CairoMakie
    using IterTools


    export OscillatorParam,SimParam,InitParam,ModelParam
    export plot

    include("params.jl")
    include("1D_wave.jl")
    include("plots.jl")
    # include("2D_wave.jl")
end # module