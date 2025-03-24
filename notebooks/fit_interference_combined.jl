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
data_location72 = joinpath(@__DIR__, "..", "data72");

# ╔═╡ 463d3c0b-606d-42d1-9ad1-105495e3b2b1
data_location0 = joinpath(@__DIR__, "..", "data0");

# ╔═╡ 7e14c341-1c17-4145-8068-e874e7be7cdd
data_location180 = joinpath(@__DIR__, "..", "data180");

# ╔═╡ 81b2984d-20a8-4f96-bd5e-c923d9d00c57
file_path72 = joinpath(data_location72, "unbinned_data72.txt");

# ╔═╡ d1088a6f-d962-490c-ad50-45fcc8774ebc
file_path0 = joinpath(data_location0, "unbinned_data0.txt");

# ╔═╡ 5df7c82a-cf79-4de5-9ac8-2e2bf771bee1
file_path180 = joinpath(data_location180, "unbinned_data180.txt");

# ╔═╡ 379b4299-6df7-4216-ab35-3be396ecc6c3
data_points72 = readlines(file_path72) .|> x -> parse(Float64, x);

# ╔═╡ 927d3bc7-5f39-451a-bf92-6680a670ea7c
data_points0 = readlines(file_path0) .|> x -> parse(Float64, x);

# ╔═╡ 680fa5ff-b7b7-4fae-a144-82b94acb6d0e
data_points180 = readlines(file_path180) .|> x -> parse(Float64, x);

# ╔═╡ aaad123c-17a4-4098-a3bb-ce48cb4426fc
theme(:boxed, ylims=(0, :auto))

# ╔═╡ d0f0452f-5835-4bb6-aee8-6dcfa9a2e2a2
md"""
## Read data
hist from file
"""

# ╔═╡ 5b5b5bbb-ba60-493b-b9d7-9e5d1d7d1f8f
#= unbinned_data = let
	nS = 1000
	nB = 10000
	support = (0.3, 1.3)
	
	data1 = randn(nS) .* 0.05 .+ 1.0;
	data2 = rand(nB) .* (support[2]-support[1]) .+ support[1]
	# 
	filter(vcat(data1, data2)) do x
		support[1] < x < support[2]
	end
end =#

# ╔═╡ b1b050a0-f283-45cd-a5df-14eb5bdff035
#specific_hist = Hist1D(unbinned_data, binedges=range(0.3,1.3, 100));

# ╔═╡ a5a8b38c-d510-4341-a57f-4dc219ed2483
specific_hist72 = Hist1D(data_points72, binedges=range(4.44, 4.49, 500));

# ╔═╡ f735c700-13f2-4e5e-8afa-b4cc6f0aae19
specific_hist0 = Hist1D(data_points0, binedges=range(4.44, 4.49, 500));

# ╔═╡ 561dba23-d72d-4bd3-bc7f-2af365c4da9b
specific_hist180 = Hist1D(data_points180, binedges=range(4.44, 4.49, 500));

# ╔═╡ 7cd14ac4-b688-4336-8532-5c52c3f9600c
scatter(bincenters(specific_hist72), specific_hist72.bincounts, yerr=sqrt.(specific_hist72.bincounts), ylim=(0,:auto))

# ╔═╡ 2cc85370-c066-45e3-8a51-ce566ac592a0
scatter(bincenters(specific_hist0), specific_hist0.bincounts, yerr=sqrt.(specific_hist0.bincounts), ylim=(0,:auto))

# ╔═╡ 290b2ea3-3aa2-452e-aa45-22a02f42ccba
scatter(bincenters(specific_hist180), specific_hist180.bincounts, yerr=sqrt.(specific_hist180.bincounts), ylim=(0,:auto))

# ╔═╡ 4edee112-1541-4ea6-a6a2-5a738742760f
md"""
## Build model

```math
I = |\mathrm{BW} e^{i\phi} + c|^2 + b^2
```
"""

# ╔═╡ ff42bba7-a2c2-4951-a9ce-ba4058578b39
begin
	@with_kw struct BWCoherentIncoherent72{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of BW
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param1
		#β::T # bkg param2
		#γ::T # bkg param3
	end
	function intensity(model::BWCoherentIncoherent72, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = 72
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent72(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent72)}(values)
		BWCoherentIncoherent72(; pars...)
	end
end

# ╔═╡ 44cb256e-2eb2-4d8c-b3e8-2e694292e2f5
begin
	@with_kw struct BWCoherentIncoherent0{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of BW
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param1
		#β::T # bkg param2
		#γ::T # bkg param3
	end
	function intensity(model::BWCoherentIncoherent0, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = 0
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent0(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent0)}(values)
		BWCoherentIncoherent0(; pars...)
	end
end

# ╔═╡ f1e5fa90-e711-4d12-8c37-081ab5f5de67
begin
	@with_kw struct BWCoherentIncoherent180{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of BW
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param1
		#β::T # bkg param2
		#γ::T # bkg param3
	end
	function intensity(model::BWCoherentIncoherent180, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = 180
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent180(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent180)}(values)
		BWCoherentIncoherent180(; pars...)
	end
end

# ╔═╡ b1f19e69-330b-41c0-bfc7-d49a2cd4fa74
#=begin
	@with_kw struct BWCoherentIncoherent0{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		#ϕ::T # phase of BW
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param1
		β::T # bkg param2
		#γ::T # bkg param3
	end
	function intensity(model::BWCoherentIncoherent0, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α, β = model
		bkg = α*m + β 
		@unpack c,b = model
		cB = c*bkg
		iB = b*bkg
		abs2(a * cis(0.0-π) + cB) + abs2(iB)
	end
	function BWCoherentIncoherent0(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent0)}(values)
		BWCoherentIncoherent0(; pars...)
	end
end=#

# ╔═╡ ed58cc37-a624-4870-a6ce-d029246a5a31
#=begin
	@with_kw struct BWCoherentIncoherent180{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		#ϕ::T # phase of BW
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param1
		β::T # bkg param2
		#γ::T # bkg param3
	end
	function intensity(model::BWCoherentIncoherent180, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α, β = model
		bkg = α*m + β 
		@unpack c,b = model
		cB = c*bkg
		iB = b*bkg
		abs2(a * cis(3.1416-π) + cB) + abs2(iB)
	end
	function BWCoherentIncoherent180(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent180)}(values)
		BWCoherentIncoherent180(; pars...)
	end
end=#

# ╔═╡ 4a20200b-5de2-480f-b21a-6d83cd5721bf
#let # test of the model intensity
	#_model = BWCoherentIncorent(; a0=1.0, m0=1.0, Γ0=0.2, ϕ = 0.0, c=0.5, b=1.0)
	#@assert intensity(_model, _model.m0) ≈ 1^2+ 0.5^2 + 1
	# 
	#_model = BWCoherentIncorent(; a0=1.0, m0=1.0, Γ0=0.2, ϕ = π/2, c=0.5, b=1.0)
	#@assert intensity(_model, _model.m0) ≈ 0.5^2 + 1
#end

# ╔═╡ 1fdfa846-e0e0-4a2f-881b-49cbcf1c1289
md"""
## Define mismatch

the loss function
"""

# ╔═╡ 0f79a81d-b54e-477d-9aee-c0750d142bf5
#=function compute_loss(model, data)
	@unpack xv, yv, δyv = data
	yv_bar = intensity.(Ref(model), xv)
	chi2v = (yv - yv_bar).^2 ./ δyv.^2
	chi2 = sum(chi2v)
	return chi2
end=#

# ╔═╡ ca2dfa8c-6f4c-4c03-9818-8e373f024cae
#=function compute_loss(model, datasets)
    chi2_total = 0.0
    for data in datasets
        @unpack xv, yv, δyv = data
        yv_bar = intensity.(Ref(model), xv)
        chi2_total += sum((yv - yv_bar).^2 ./ δyv.^2)
    end
    return chi2_total
end=#

# ╔═╡ 6db468ec-724b-4546-9a29-899b8da8706a
function compute_loss(models, datasets)
    chi2_total = 0.0
    for (model, data) in zip(models, datasets)
        @unpack xv, yv, δyv = data
        yv_bar = intensity.(Ref(model), xv)
        chi2_total += sum((yv - yv_bar).^2 ./ δyv.^2)
    end
    return chi2_total
end

# ╔═╡ 7ac3cb4f-b874-4cb9-9da8-c874943adbf2
md"""
## Adjust

for initial parameters
"""

# ╔═╡ 5f4a600a-eca9-4d42-b48a-307838a14b25
initial_pars = (; a0=47.0, m0=4.461, Γ0=0.001, c=5.0, b=10.0, α=0.1, ϕB=3.14)

# ╔═╡ 6349f6d9-6a29-4356-bc72-d180d6b8585f
#initial_model = BWCoherentIncoherent(; initial_pars...)

# ╔═╡ ec85b957-7845-4d78-90de-def823fc0428
#=let
	xlims = (4.44, 4.49)
	plot(ylims=(0,:auto))
	plot!(xlims...) do m
		intensity(initial_model, m)
	end
	plot!(xlims..., fill=0, alpha=0.3) do m
		intensity(BWCoherentIncoherent(initial_model; c=0, b=0), m)
	end
end=#

# ╔═╡ 111dbc4b-66b5-4cd1-bab4-f025262fb963
md"""
## Fit
"""

# ╔═╡ fa95a3e2-be56-42ba-a528-91fe73071bff
const specific_data72 = (; 
	xv = bincenters(specific_hist72),
	yv = specific_hist72.bincounts,
	δyv = sqrt.(specific_hist72.bincounts))

# ╔═╡ 6f7247ae-8c94-4104-ac35-23665410d475
const specific_data0 = (; 
	xv = bincenters(specific_hist0),
	yv = specific_hist0.bincounts,
	δyv = sqrt.(specific_hist0.bincounts))

# ╔═╡ dbfad0c7-d086-4f4e-82ef-5834deeacca0
const specific_data180 = (; 
	xv = bincenters(specific_hist180),
	yv = specific_hist180.bincounts,
	δyv = sqrt.(specific_hist180.bincounts))

# ╔═╡ 02b27549-178f-4b21-8048-5c4c2a48a937
specific_datasets = [specific_data72, specific_data0, specific_data180]

# ╔═╡ 6dfbda41-b976-4138-8834-fad3155e19fa
specific_models(pars) = [BWCoherentIncoherent72(pars), BWCoherentIncoherent0(pars), BWCoherentIncoherent180(pars)]

# ╔═╡ 22ea68a4-d8e1-4b66-add4-2d270a614e23
loss(pars) = compute_loss(specific_models(pars), specific_datasets)

# ╔═╡ b7fcac1a-5e12-4c2e-a380-bab8ce0859e9
@assert loss(collect(initial_pars)) isa Number

# ╔═╡ aef67dfa-6775-4220-9a9b-151891609bdf
fit_result = optimize(loss, collect(initial_pars), BFGS())

# ╔═╡ 46b0f126-f613-48ba-aba6-be965c66f865
fit_result.minimum

# ╔═╡ 88f88f6e-ff7f-441d-96ae-69269db978b7
best_model72 = BWCoherentIncoherent72(fit_result.minimizer)

# ╔═╡ 2306ae6a-281f-45c9-8f67-2f73da145935
best_model0 = BWCoherentIncoherent0(fit_result.minimizer)

# ╔═╡ 3c554d00-f68a-40d5-88ac-f7a9912f3dce
best_model180 = BWCoherentIncoherent180(fit_result.minimizer)

# ╔═╡ 1837cf4c-18a9-4735-b066-eaca988b92d4
ratio = best_model180.c / best_model180.b

# ╔═╡ e9c36bcf-c852-46eb-8ea7-f33208dbddf1
#angle = best_model180.ϕB % 2π

# ╔═╡ f0615806-2d24-44b1-b162-49f8759a5f98
let
	xlims = (specific_hist72.binedges[1][1], specific_hist72.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist72), specific_hist72.bincounts, yerr=sqrt.(specific_hist72.bincounts), ylim=(0,:auto))
	plot!(m->intensity(best_model72, m), xlims...) 
end

# ╔═╡ 7c1a91f4-73e8-4e72-b012-17c02e1aab9b
let
	xlims = (specific_hist0.binedges[1][1], specific_hist0.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist0), specific_hist0.bincounts, yerr=sqrt.(specific_hist0.bincounts), ylim=(0,:auto))
	plot!(m->intensity(best_model0, m), xlims...) 
end

# ╔═╡ 3aa658ec-479e-476a-ad53-ae431d674152
let
	xlims = (specific_hist180.binedges[1][1], specific_hist180.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist180), specific_hist180.bincounts, yerr=sqrt.(specific_hist180.bincounts), ylim=(0,:auto))
	plot!(m->intensity(best_model180, m), xlims...) 
end

# ╔═╡ Cell order:
# ╠═ec7ef8e2-f9a8-11ef-2e7c-0f562efc23a1
# ╠═59dd328c-6fac-4bf1-907d-10f3d96ec1a6
# ╠═463d3c0b-606d-42d1-9ad1-105495e3b2b1
# ╠═7e14c341-1c17-4145-8068-e874e7be7cdd
# ╠═81b2984d-20a8-4f96-bd5e-c923d9d00c57
# ╠═d1088a6f-d962-490c-ad50-45fcc8774ebc
# ╠═5df7c82a-cf79-4de5-9ac8-2e2bf771bee1
# ╠═379b4299-6df7-4216-ab35-3be396ecc6c3
# ╠═927d3bc7-5f39-451a-bf92-6680a670ea7c
# ╠═680fa5ff-b7b7-4fae-a144-82b94acb6d0e
# ╠═aaad123c-17a4-4098-a3bb-ce48cb4426fc
# ╟─d0f0452f-5835-4bb6-aee8-6dcfa9a2e2a2
# ╠═5b5b5bbb-ba60-493b-b9d7-9e5d1d7d1f8f
# ╠═b1b050a0-f283-45cd-a5df-14eb5bdff035
# ╠═a5a8b38c-d510-4341-a57f-4dc219ed2483
# ╠═f735c700-13f2-4e5e-8afa-b4cc6f0aae19
# ╠═561dba23-d72d-4bd3-bc7f-2af365c4da9b
# ╠═7cd14ac4-b688-4336-8532-5c52c3f9600c
# ╠═2cc85370-c066-45e3-8a51-ce566ac592a0
# ╠═290b2ea3-3aa2-452e-aa45-22a02f42ccba
# ╟─4edee112-1541-4ea6-a6a2-5a738742760f
# ╠═ff42bba7-a2c2-4951-a9ce-ba4058578b39
# ╠═44cb256e-2eb2-4d8c-b3e8-2e694292e2f5
# ╠═f1e5fa90-e711-4d12-8c37-081ab5f5de67
# ╠═b1f19e69-330b-41c0-bfc7-d49a2cd4fa74
# ╠═ed58cc37-a624-4870-a6ce-d029246a5a31
# ╠═4a20200b-5de2-480f-b21a-6d83cd5721bf
# ╟─1fdfa846-e0e0-4a2f-881b-49cbcf1c1289
# ╠═0f79a81d-b54e-477d-9aee-c0750d142bf5
# ╠═ca2dfa8c-6f4c-4c03-9818-8e373f024cae
# ╠═6db468ec-724b-4546-9a29-899b8da8706a
# ╟─7ac3cb4f-b874-4cb9-9da8-c874943adbf2
# ╠═5f4a600a-eca9-4d42-b48a-307838a14b25
# ╠═6349f6d9-6a29-4356-bc72-d180d6b8585f
# ╠═ec85b957-7845-4d78-90de-def823fc0428
# ╟─111dbc4b-66b5-4cd1-bab4-f025262fb963
# ╠═fa95a3e2-be56-42ba-a528-91fe73071bff
# ╠═6f7247ae-8c94-4104-ac35-23665410d475
# ╠═dbfad0c7-d086-4f4e-82ef-5834deeacca0
# ╠═02b27549-178f-4b21-8048-5c4c2a48a937
# ╠═6dfbda41-b976-4138-8834-fad3155e19fa
# ╠═22ea68a4-d8e1-4b66-add4-2d270a614e23
# ╠═b7fcac1a-5e12-4c2e-a380-bab8ce0859e9
# ╠═aef67dfa-6775-4220-9a9b-151891609bdf
# ╠═46b0f126-f613-48ba-aba6-be965c66f865
# ╠═88f88f6e-ff7f-441d-96ae-69269db978b7
# ╠═2306ae6a-281f-45c9-8f67-2f73da145935
# ╠═3c554d00-f68a-40d5-88ac-f7a9912f3dce
# ╠═1837cf4c-18a9-4735-b066-eaca988b92d4
# ╠═e9c36bcf-c852-46eb-8ea7-f33208dbddf1
# ╠═f0615806-2d24-44b1-b162-49f8759a5f98
# ╠═7c1a91f4-73e8-4e72-b012-17c02e1aab9b
# ╠═3aa658ec-479e-476a-ad53-ae431d674152
