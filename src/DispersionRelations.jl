module DispersionRelations

using LinearAlgebra

export fit_complex_frequency, fit_pure_growth_rate

function find_peaks(E)
    vcat(false, false, (diff(sign.(diff(log.(E)))) .== -2))
end

"""
fit_complex_frequency(t, E, use_peaks=nothing)

Find the best fit to an exponential model for the given energy trace.
That is, this function finds a fit of the form E ≈ E₀ * exp(iωt).

If `use_peaks` is provided, it serves as a range or vector of which peaks to use
during the fitting process.
If it is omitted, the function uses a heuristic to find the range of consecutive
peaks which provides the most stable fit.

Return a pair `(line, ω)` where
    - `line` is a flat line that passes through the peaks used for the fit
    - `ω` is the complex frequency of the estimated exponential
"""
function fit_complex_frequency(t, E, use_peaks=nothing)
    peaks = (diff(sign.(diff(log.(E)))) .== -2)

    j = sum(peaks)

    fit_slope(range) = begin
        t̂ = t[3:end][peaks][range]
        X = hcat(t̂, ones(length(range)))
        Ê = E[3:end][peaks][range] .|> log
        y = X \ Ê
    end

    fit_freq(range) = begin
        t̂ = t[3:end][peaks][range]
        X = hcat(1:length(range), ones(length(range)))
        y = X \ t̂
    end

    find_best_fit(fitfunc) = begin
        changes = fill(Inf, j, j)
        fits = Array{Vector{Float64}}(undef, j, j)

        for first in 1:j-1
            for last in (first+2):j
                y = fitfunc(first:last)
                fits[first, last] = y

                y_prev = fitfunc((first+1):last)
                y_next = fitfunc(first:(last-1))
                rel_prev = norm(y - y_prev) / max(norm(y), norm(y_prev))
                rel_next = norm(y - y_next) / max(norm(y), norm(y_next))

                changes[first, last] = max(rel_prev, rel_next)
            end
        end
        bestJ = argmin(changes)
        fits[bestJ]
    end

    γ, intercept = if isnothing(use_peaks) 
        find_best_fit(fit_slope)
    else
        fit_slope(use_peaks)
    end

    line = exp.(γ * t .+ intercept)

    invfreq, intercept = if isnothing(use_peaks)
        find_best_fit(fit_freq)
    else
        fit_freq(use_peaks)
    end
    freq = 2π/invfreq

    line, freq + im*γ
end

"""
find_pure_growth_slope(t, E, time_range=nothing)

Estimate an exponential fit of the form E ≈ E₀ * exp(γt).

Return `(line, γ)`, where
    - `γ` is the best fit to the exponential growth rate,
    - `line` is `E₀ * exp.(γ .* t)`.
"""
function fit_pure_growth_rate(t, E, time_range=nothing)
    if !isnothing(time_range)
        range = searchsortedfirst(t, first(time_range)):searchsortedlast(t, last(time_range))
    else
        range = Colon()
    end
    t̂ = t[range]
    Ê = E[range]
    X = hcat(t̂, ones(length(t̂)))
    y = X \ (log.(Ê))
    γ, intercept = y
    line = exp.(γ * t .+ intercept)
    return line, γ
end


end
