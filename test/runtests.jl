using Test
using BAT4Epidemiology

Test.@testset "BAT4Epidemiology" begin
    include("test_SEIRModel.jl")
    #include("test_CPSEIRModel.jl")
end