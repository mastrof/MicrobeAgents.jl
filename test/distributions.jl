using MicrobeAgents, Test
using Random

@testset "Distributions" begin
    @test_throws ArgumentError Arccos(0.5,0.1) # must be a < b
    @test_throws ArgumentError Arccos(-5,0.5) # must be -1 ≤ a ≤ 1
    @test_throws ArgumentError Arccos(-1,2) # must be -1 ≤ b ≤ 1
    d = Arccos()
    # if unspecified a=-1 and b=1
    @test d.a == -1 && d.b == 1
    d = Arccos(-0.5, 0.75)
    @test d.a == -0.5 && d.b == 0.75
    Random.seed!(96)
    x = rand(Arccos())
    y = acos(rand(Random.seed!(96))*2-1)
    @test x == y
    # works with custom rng
    x = rand(MersenneTwister(8), Arccos())
    y = acos(rand(MersenneTwister(8))*2-1)
    @test x == y
end
