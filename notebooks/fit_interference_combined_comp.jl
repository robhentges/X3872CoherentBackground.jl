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
data_location108m = joinpath(@__DIR__, "..", "data108m");

# ╔═╡ 463d3c0b-606d-42d1-9ad1-105495e3b2b1
data_location36 = joinpath(@__DIR__, "..", "data36");

# ╔═╡ 7e14c341-1c17-4145-8068-e874e7be7cdd
data_location144 = joinpath(@__DIR__, "..", "data144");

# ╔═╡ 81b2984d-20a8-4f96-bd5e-c923d9d00c57
file_path108m = joinpath(data_location108m, "unbinned_data108m.txt");

# ╔═╡ d1088a6f-d962-490c-ad50-45fcc8774ebc
file_path36 = joinpath(data_location36, "unbinned_data36.txt");

# ╔═╡ 5df7c82a-cf79-4de5-9ac8-2e2bf771bee1
file_path144 = joinpath(data_location144, "unbinned_data144.txt");

# ╔═╡ 379b4299-6df7-4216-ab35-3be396ecc6c3
data_points108m = readlines(file_path108m) .|> x -> parse(Float64, x);

# ╔═╡ 927d3bc7-5f39-451a-bf92-6680a670ea7c
data_points36 = readlines(file_path36) .|> x -> parse(Float64, x);

# ╔═╡ 680fa5ff-b7b7-4fae-a144-82b94acb6d0e
data_points144 = readlines(file_path144) .|> x -> parse(Float64, x);

# ╔═╡ aaad123c-17a4-4098-a3bb-ce48cb4426fc
theme(:boxed, ylims=(0, :auto))

# ╔═╡ d0f0452f-5835-4bb6-aee8-6dcfa9a2e2a2
md"""
## Read data
hist from file
"""

# ╔═╡ a5a8b38c-d510-4341-a57f-4dc219ed2483
specific_hist108m = Hist1D(data_points108m, binedges=range(4.44, 4.49, 500));

# ╔═╡ f735c700-13f2-4e5e-8afa-b4cc6f0aae19
specific_hist36 = Hist1D(data_points36, binedges=range(4.44, 4.49, 500));

# ╔═╡ 561dba23-d72d-4bd3-bc7f-2af365c4da9b
specific_hist144 = Hist1D(data_points144, binedges=range(4.44, 4.49, 500));

# ╔═╡ 7cd14ac4-b688-4336-8532-5c52c3f9600c
#scatter(bincenters(specific_hist72), specific_hist72.bincounts, yerr=sqrt.(specific_hist72.bincounts), ylim=(0,:auto))

# ╔═╡ 2cc85370-c066-45e3-8a51-ce566ac592a0
#scatter(bincenters(specific_hist0), specific_hist0.bincounts, yerr=sqrt.(specific_hist0.bincounts), ylim=(0,:auto))

# ╔═╡ 290b2ea3-3aa2-452e-aa45-22a02f42ccba
#scatter(bincenters(specific_hist180), specific_hist180.bincounts, yerr=sqrt.(specific_hist180.bincounts), ylim=(0,:auto))

# ╔═╡ 4edee112-1541-4ea6-a6a2-5a738742760f
md"""
## Build model

```math
I = |\mathrm{BW} e^{i\phi} + c|^2 + b^2
```
"""

# ╔═╡ ff42bba7-a2c2-4951-a9ce-ba4058578b39
begin
	@with_kw struct BWCoherentIncoherent108m{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent108m, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = -108
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent108m(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent108m)}(values)
		BWCoherentIncoherent108m(; pars...)
	end
end

# ╔═╡ 44cb256e-2eb2-4d8c-b3e8-2e694292e2f5
begin
	@with_kw struct BWCoherentIncoherent36{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent36, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = 36
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent36(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent36)}(values)
		BWCoherentIncoherent36(; pars...)
	end
end

# ╔═╡ f1e5fa90-e711-4d12-8c37-081ab5f5de67
begin
	@with_kw struct BWCoherentIncoherent144{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent144, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = 144
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent144(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent144)}(values)
		BWCoherentIncoherent144(; pars...)
	end
end

# ╔═╡ 1fdfa846-e0e0-4a2f-881b-49cbcf1c1289
md"""
## Define mismatch

the loss function
"""

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
initial_pars = (; a0=47.0, m0=4.461, Γ0=0.001, c=5.0, b=10.0, α=0.1, ϕB=-0.1)

# ╔═╡ 6349f6d9-6a29-4356-bc72-d180d6b8585f
#initial_model = BWCoherentIncoherent(; initial_pars...)

# ╔═╡ 111dbc4b-66b5-4cd1-bab4-f025262fb963
md"""
## Fit
"""

# ╔═╡ fa95a3e2-be56-42ba-a528-91fe73071bff
const specific_data108m = (; 
	xv = bincenters(specific_hist108m),
	yv = specific_hist108m.bincounts,
	δyv = sqrt.(specific_hist108m.bincounts))

# ╔═╡ 6f7247ae-8c94-4104-ac35-23665410d475
const specific_data36 = (; 
	xv = bincenters(specific_hist36),
	yv = specific_hist36.bincounts,
	δyv = sqrt.(specific_hist36.bincounts))

# ╔═╡ dbfad0c7-d086-4f4e-82ef-5834deeacca0
const specific_data144 = (; 
	xv = bincenters(specific_hist144),
	yv = specific_hist144.bincounts,
	δyv = sqrt.(specific_hist144.bincounts))

# ╔═╡ 02b27549-178f-4b21-8048-5c4c2a48a937
specific_datasets = [specific_data108m, specific_data36, specific_data144]

# ╔═╡ 6dfbda41-b976-4138-8834-fad3155e19fa
specific_models(pars) = [BWCoherentIncoherent108m(pars), BWCoherentIncoherent36(pars), BWCoherentIncoherent144(pars)]

# ╔═╡ 22ea68a4-d8e1-4b66-add4-2d270a614e23
loss(pars) = compute_loss(specific_models(pars), specific_datasets)

# ╔═╡ b7fcac1a-5e12-4c2e-a380-bab8ce0859e9
@assert loss(collect(initial_pars)) isa Number

# ╔═╡ aef67dfa-6775-4220-9a9b-151891609bdf
fit_result = optimize(loss, collect(initial_pars), BFGS())

# ╔═╡ 46b0f126-f613-48ba-aba6-be965c66f865
fit_result.minimum

# ╔═╡ 88f88f6e-ff7f-441d-96ae-69269db978b7
best_model108m = BWCoherentIncoherent108m(fit_result.minimizer)

# ╔═╡ 2306ae6a-281f-45c9-8f67-2f73da145935
best_model36 = BWCoherentIncoherent36(fit_result.minimizer)

# ╔═╡ 3c554d00-f68a-40d5-88ac-f7a9912f3dce
best_model144 = BWCoherentIncoherent144(fit_result.minimizer)

# ╔═╡ 1837cf4c-18a9-4735-b066-eaca988b92d4
ratio = best_model144.c / best_model144.b

# ╔═╡ e9c36bcf-c852-46eb-8ea7-f33208dbddf1
angle = best_model144.ϕB % 2π

# ╔═╡ f0615806-2d24-44b1-b162-49f8759a5f98
let
	xlims = (specific_hist108m.binedges[1][1], specific_hist108m.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist108m), specific_hist108m.bincounts, yerr=sqrt.(specific_hist108m.bincounts), ylim=(0,:auto))
	plot!(m->intensity(best_model108m, m), xlims...) 
end

# ╔═╡ 7c1a91f4-73e8-4e72-b012-17c02e1aab9b
let
	xlims = (specific_hist36.binedges[1][1], specific_hist36.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist36), specific_hist36.bincounts, yerr=sqrt.(specific_hist36.bincounts), ylim=(0,:auto))
	plot!(m->intensity(best_model36, m), xlims...) 
end

# ╔═╡ 3aa658ec-479e-476a-ad53-ae431d674152
let
	xlims = (specific_hist144.binedges[1][1], specific_hist144.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist144), specific_hist144.bincounts, yerr=sqrt.(specific_hist144.bincounts), ylim=(0,:auto))
	plot!(m->intensity(best_model144, m), xlims...) 
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
# ╟─1fdfa846-e0e0-4a2f-881b-49cbcf1c1289
# ╠═6db468ec-724b-4546-9a29-899b8da8706a
# ╟─7ac3cb4f-b874-4cb9-9da8-c874943adbf2
# ╠═5f4a600a-eca9-4d42-b48a-307838a14b25
# ╠═6349f6d9-6a29-4356-bc72-d180d6b8585f
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
