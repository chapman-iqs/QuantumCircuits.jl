@userplot BlochWireFrame
@recipe function f(wf::BlochWireFrame; mesh=30, ax=false, viewϕ=0)

	# Wire frame coordinates ---------------------------------------------------

	x(θ, ϕ) = sin(θ) * cos(ϕ + viewϕ)
	y(θ, ϕ) = sin(θ) * sin(ϕ + viewϕ)
	z(θ, ϕ) = cos(θ)

	θs = range(0, 2π, length=mesh)
	ϕs = range(0, π, length=mesh) .+ viewϕ


	# Plot wireframe -----------------------------------------------------------


	legend := false
	linecolor := "steelblue"
	linewidth := 0.5
	linealpha := 1
	seriestype := path3d
	aspect_ratio := 1.0
	size --> (400,400)

	# Longitudes
	for ϕ in ϕs
		@series begin
			[x(θ, ϕ) for θ in θs], [y(θ, ϕ) for θ in θs], [z(θ, ϕ) for θ in θs]
		end
	end

	# Latitudes
	for θ in θs
		@series begin
			[x(θ, ϕ) for ϕ in ϕs], [y(θ, ϕ) for ϕ in ϕs], [z(θ, ϕ) for ϕ in ϕs]
		end
	end


	# Plot reference axes ------------------------------------------------------
	linewidth := 3
	colors = [palette(:tab10)[i] for i in 1:3]

	if ax

		linecolor := colors[1]
		@series begin
			[0, cos(viewϕ)], [0, sin(viewϕ)], [0, 0]
		end

		linecolor := colors[2]
		@series begin
			[0, -sin(viewϕ)], [0, cos(viewϕ)], [0, 0]
		end

		linecolor := colors[3]
		@series begin
			[0, 0], [0, 0], [0, 1]
		end
	end
end # wireframe

@userplot BlochSphere
@recipe function f(bs::BlochSphere; mesh=30, ax=false, viewϕ=0, vec=nothing, blochmark=false, blochmarkcolor="white")

	args = bs.args
	timed = length(args) == 4

	xss, yss, zss = timed ? args[2:4] : args
	i = timed ? bs.args[1] : length(xss)

	xs = xss[1:i] .* cos(viewϕ) .- yss[1:i] .* sin(viewϕ)
	ys = xss[1:i] .* sin(viewϕ) .+ yss[1:i] .* cos(viewϕ)
	zs = zss[1:i]


	# Plot trajectory ----------------------------------------------------------

	# marker := nothing
	linecolor --> "black" # makes blue only if color is not specified
	linewidth --> :2

	@series begin
		xs, ys, zs
	end


	# Wire frame coordinates ---------------------------------------------------

	x(θ, ϕ) = sin(θ) * cos(ϕ + viewϕ)
	y(θ, ϕ) = sin(θ) * sin(ϕ + viewϕ)
	z(θ, ϕ) = cos(θ)

	θs = range(0, 2π, length=mesh)
	ϕs = range(0, π, length=mesh) .+ viewϕ


	# Plot wireframe -----------------------------------------------------------


	legend := false
	linecolor := "steelblue"
	linewidth := 0.5
	linealpha := 1
	seriestype := path3d
	aspect_ratio := 1.0
	size --> (400,400)

	# Longitudes
	for ϕ in ϕs
		@series begin
			[x(θ, ϕ) for θ in θs], [y(θ, ϕ) for θ in θs], [z(θ, ϕ) for θ in θs]
		end
	end

	# Latitudes
	for θ in θs
		@series begin
			[x(θ, ϕ) for ϕ in ϕs], [y(θ, ϕ) for ϕ in ϕs], [z(θ, ϕ) for ϕ in ϕs]
		end
	end


	# Plot reference axes ------------------------------------------------------
	linewidth := 3
	colors = [palette(:tab10)[i] for i in 1:3]

	if ax

		linecolor := colors[1]
		@series begin
			[0, cos(viewϕ)], [0, sin(viewϕ)], [0, 0]
		end

		linecolor := colors[2]
		@series begin
			[0, -sin(viewϕ)], [0, cos(viewϕ)], [0, 0]
		end

		linecolor := colors[3]
		@series begin
			[0, 0], [0, 0], [0, 1]
		end
	end


	# Plot optional Bloch vector input by user --------------------------------

	if vec != nothing

		(xvv, yvv, zvv) = vec

		xv = xvv * cos(viewϕ) - yvv * sin(viewϕ)
		yv = xvv .* sin(viewϕ) .+ yvv * cos(viewϕ)
		zv = zvv

		linewidth := 2
		linecolor := "red"

		@series begin
			[0, xv], [0, yv], [0, zv]
		end

		marker := (:circle, 5)
		markercolor := "red"
		@series begin
			[xv], [yv], [zv]
		end

	end

	# Plot optional marker ------------------------------------------------

	if blochmark

		marker := (:circle, 3.5)
		markercolor := blochmarkcolor
		@series begin
			[last(xs)], [last(ys)], [last(zs)]
		end

	end
end # blochsphere
