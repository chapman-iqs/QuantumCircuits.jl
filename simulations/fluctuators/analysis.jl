using Suppressor: @suppress_err
@suppress_err using Plots, FFTW

function plot_fluctuators(ts, χtotal)
	return plot(ts, χtotal, xlabel = "t (μs)", legend=:none, title=string("total fluctuation"), titlefontsize=13)
end

function PSD(ts, series)
	t0 = first(ts)
	dt = ts[2] - ts[1]
	tf = last(ts)
	fs = 1/dt          # Sampling rate

	t = t0:dt:tf

	F = fftshift(fft(series))
	freqs = fftshift(fftfreq(length(t), fs))

	df = step(freqs)
	psd = map(f -> abs(f)^2 / 2df, F)

	# return only positive frequencies
	i = findfirst(freqs .> 0)

	return (freqs[i:end], psd[i:end])
end # function PSD

function plotPSD(freqs, PSD; scale=20e3, logscale=true)
	xaxis = logscale ? :log10 : :none
	yaxis = logscale ? :log10 : :none
	positive_freqs = range(0, last(freqs), step=step(freqs)) ./ 1000 # to convert from MHz to kHz
	(m, b) = fitPSD(freqs, PSD)
	title = string("PSD, y = ", m, "x + ", b)
	p = plot(freqs, PSD, title=title, xlim=(0, 20), legend=:none, xlabel="f (kHz)",
				xaxis=xaxis, yaxis=yaxis)

	# plot!(positive_freqs, scale ./positive_freqs, xlims=[1,20], xaxis=xaxis, yaxis=yaxis, linewidth=5)
	return p
end


"""
Uses LsqFit package to do a linear fit to log-log data. The PSD for a (1/f)^α spectrum looks
like a line of slope α in log-log space.
"""
function fitPSD(freqs, psd)
	logfreqs = log10.(freqs)
	logPSD = log10.(psd)

	model(f, p) = p[1] .* f .+ p[2]
	p0 = [-2, 8]

	fit = curve_fit(model, logfreqs, logPSD, p0)
	param = fit.param
	return param
end
