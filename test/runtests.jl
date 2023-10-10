using DispersionRelations
using Test

@testset "DispersionRelations.jl" begin
    t = 0.0:0.01:5.0

    @testset "purely growing" begin
        E = exp.(0.04 * t)
        line, γ = fit_pure_growth_rate(t, E)
        @test γ ≈ 0.04
        @test line ≈ E[1] * exp.(γ * t)

        E = exp.(1.53 * t)
        line, γ = fit_pure_growth_rate(t, E)
        @test γ ≈ 1.53

        E = exp.(2.31 * t) + 0.4 * cos.(t)
        line, γ = fit_pure_growth_rate(t, E, time_range=[1.0, 5.0])
        @test isapprox(γ, 2.31, rtol=0.01)
    end

    @testset "complex frequency" begin
        # In all these tests, we square the real part of the oscillating value,
        # so the frequency must be multiplied by two in our assertions.
        
        E = abs2.(real.(exp.(-im * (6.34 + im*1.53) * t)))
        line, ω = fit_complex_frequency(t, E)
        @test isapprox(ω, 2 * (6.34 + im*1.53), rtol=0.01)

        @testset "too few peaks" begin
            E = abs2.(real.(exp.(-im * (0.98 + im*1.53) * t)))
            
            @test_throws ErrorException("Cannot automatically fit complex frequency with fewer than 3 peaks to work with") fit_complex_frequency(t, E)
        end

        @testset "Exactly three peaks" begin
            E = abs2.(real.(exp.(-im * (1.7 + im*1.53) * t)))
            line, ω = fit_complex_frequency(t, E)
            @test isapprox(ω, 2 * (1.7 + im*1.53), rtol=0.01)
        end

        @testset "Manually specified two peaks" begin
            E = abs2.(real.(exp.(-im * (0.98 + im*0.87) * t)))
            line, ω = fit_complex_frequency(t, E, use_peaks=1:2)
            @test isapprox(ω, 2 * (0.98 + im*0.87), rtol=0.01)
        end
    end
end
