### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ ec7ef8e2-f9a8-11ef-2e7c-0f562efc23a1
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()

	# 
	using Plots
	using HighEnergyTools.Optim
	using HighEnergyTools.Parameters
	using HighEnergyTools.FHist
	using UnROOT
	using HadronicLineshapes
end

# ╔═╡ 59dd328c-6fac-4bf1-907d-10f3d96ec1a6
data_location = joinpath(@__DIR__, "..", "data");

# ╔═╡ aaad123c-17a4-4098-a3bb-ce48cb4426fc
theme(:boxed, ylims=(0, :auto))

# ╔═╡ d0f0452f-5835-4bb6-aee8-6dcfa9a2e2a2
md"""
## Read data
hist from file
"""

# ╔═╡ 5b5b5bbb-ba60-493b-b9d7-9e5d1d7d1f8f
unbinned_data = let
	nS = 1000
	nB = 10000
	support = (0.3, 1.3)
	
	data1 = randn(nS) .* 0.05 .+ 1.0;
	data2 = rand(nB) .* (support[2]-support[1]) .+ support[1]
	# 
	filter(vcat(data1, data2)) do x
		support[1] < x < support[2]
	end
end

# ╔═╡ b1b050a0-f283-45cd-a5df-14eb5bdff035
specific_hist = Hist1D(unbinned_data, binedges=range(0.3,1.3, 100));

# ╔═╡ 7cd14ac4-b688-4336-8532-5c52c3f9600c
scatter(bincenters(specific_hist), specific_hist.bincounts, yerr=sqrt.(specific_hist.bincounts), ylim=(0,:auto))

# ╔═╡ 4edee112-1541-4ea6-a6a2-5a738742760f
md"""
## Build model

```math
I = |\mathrm{BW} e^{i\phi} + c|^2 + b^2
```
"""

# ╔═╡ ff42bba7-a2c2-4951-a9ce-ba4058578b39
begin
	@with_kw struct BWCoherentIncorent{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕ::T # phase of BW
		c::T # coherent backgrond
		b::T # incoh
	end
	function intensity(model::BWCoherentIncorent, m)
		@unpack a0, m0, Γ0, ϕ = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack c,b = model
		abs2(a * cis(ϕ) + c) + abs2(b)
	end
	function BWCoherentIncorent(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncorent)}(values)
		BWCoherentIncorent(; pars...)
	end
end

# ╔═╡ 4a20200b-5de2-480f-b21a-6d83cd5721bf
let # test of the model intensity
	_model = BWCoherentIncorent(; a0=1.0, m0=1.0, Γ0=0.2, ϕ = 0.0, c=0.5, b=1.0)
	@assert intensity(_model, _model.m0) ≈ 1^2+ 0.5^2 + 1
	# 
	_model = BWCoherentIncorent(; a0=1.0, m0=1.0, Γ0=0.2, ϕ = π/2, c=0.5, b=1.0)
	@assert intensity(_model, _model.m0) ≈ 0.5^2 + 1
end

# ╔═╡ 1fdfa846-e0e0-4a2f-881b-49cbcf1c1289
md"""
## Define mismatch

the loss function
"""

# ╔═╡ 0f79a81d-b54e-477d-9aee-c0750d142bf5
function compute_loss(model, data)
	@unpack xv, yv, δyv = data
	yv_bar = intensity.(Ref(model), xv)
	chi2v = (yv - yv_bar).^2 ./ δyv.^2
	chi2 = sum(chi2v)
	return chi2
end

# ╔═╡ 7ac3cb4f-b874-4cb9-9da8-c874943adbf2
md"""
## Adjust

for initial parameters
"""

# ╔═╡ 5f4a600a-eca9-4d42-b48a-307838a14b25
initial_pars = (; a0=1.1, m0=1.0, Γ0=0.2, ϕ=1.1, c=0.3, b=1.2)

# ╔═╡ 6349f6d9-6a29-4356-bc72-d180d6b8585f
initial_model = BWCoherentIncorent(; initial_pars...)

# ╔═╡ ec85b957-7845-4d78-90de-def823fc0428
let
	xlims = (0.3, 1.5)
	plot(ylims=(0,:auto))
	plot!(xlims...) do m
		intensity(initial_model, m)
	end
	plot!(xlims..., fill=0, alpha=0.3) do m
		intensity(BWCoherentIncorent(initial_model; c=0, b=0), m)
	end
end

# ╔═╡ 111dbc4b-66b5-4cd1-bab4-f025262fb963
md"""
## Fit
"""

# ╔═╡ fa95a3e2-be56-42ba-a528-91fe73071bff
const specific_data = (; 
	xv = bincenters(specific_hist),
	yv = specific_hist.bincounts,
	δyv = sqrt.(specific_hist.bincounts))

# ╔═╡ 22ea68a4-d8e1-4b66-add4-2d270a614e23
loss(pars) = compute_loss(BWCoherentIncorent(pars), specific_data)

# ╔═╡ b7fcac1a-5e12-4c2e-a380-bab8ce0859e9
@assert loss(collect(initial_pars)) isa Number

# ╔═╡ aef67dfa-6775-4220-9a9b-151891609bdf
fit_result = optimize(loss, collect(initial_pars), BFGS())

# ╔═╡ 46b0f126-f613-48ba-aba6-be965c66f865
fit_result.minimum

# ╔═╡ 88f88f6e-ff7f-441d-96ae-69269db978b7
best_model = BWCoherentIncorent(fit_result.minimizer)

# ╔═╡ f0615806-2d24-44b1-b162-49f8759a5f98
let
	xlims = (specific_hist.binedges[1][1], specific_hist.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist), specific_hist.bincounts, yerr=sqrt.(specific_hist.bincounts), ylim=(0,:auto))
	plot!(m->intensity(best_model, m), xlims...) 
end

# ╔═╡ Cell order:
# ╠═ec7ef8e2-f9a8-11ef-2e7c-0f562efc23a1
# ╠═59dd328c-6fac-4bf1-907d-10f3d96ec1a6
# ╠═aaad123c-17a4-4098-a3bb-ce48cb4426fc
# ╟─d0f0452f-5835-4bb6-aee8-6dcfa9a2e2a2
# ╠═5b5b5bbb-ba60-493b-b9d7-9e5d1d7d1f8f
# ╠═b1b050a0-f283-45cd-a5df-14eb5bdff035
# ╠═7cd14ac4-b688-4336-8532-5c52c3f9600c
# ╟─4edee112-1541-4ea6-a6a2-5a738742760f
# ╠═ff42bba7-a2c2-4951-a9ce-ba4058578b39
# ╠═4a20200b-5de2-480f-b21a-6d83cd5721bf
# ╟─1fdfa846-e0e0-4a2f-881b-49cbcf1c1289
# ╠═0f79a81d-b54e-477d-9aee-c0750d142bf5
# ╟─7ac3cb4f-b874-4cb9-9da8-c874943adbf2
# ╠═5f4a600a-eca9-4d42-b48a-307838a14b25
# ╠═6349f6d9-6a29-4356-bc72-d180d6b8585f
# ╠═ec85b957-7845-4d78-90de-def823fc0428
# ╟─111dbc4b-66b5-4cd1-bab4-f025262fb963
# ╠═fa95a3e2-be56-42ba-a528-91fe73071bff
# ╠═22ea68a4-d8e1-4b66-add4-2d270a614e23
# ╠═b7fcac1a-5e12-4c2e-a380-bab8ce0859e9
# ╠═aef67dfa-6775-4220-9a9b-151891609bdf
# ╠═46b0f126-f613-48ba-aba6-be965c66f865
# ╠═88f88f6e-ff7f-441d-96ae-69269db978b7
# ╠═f0615806-2d24-44b1-b162-49f8759a5f98
