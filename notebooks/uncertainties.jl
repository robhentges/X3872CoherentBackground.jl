### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ ec7ef8e2-f9a8-11ef-2e7c-0f562efc23a1
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
	Pkg.add("ForwardDiff")
	Pkg.add("FiniteDiff")

	# 
	using Plots
	using HighEnergyTools.Optim
	using HighEnergyTools.Parameters
	using HighEnergyTools.FHist
	using UnROOT
	using HadronicLineshapes
	using ForwardDiff
	using FiniteDiff
	using LinearAlgebra
end

# ╔═╡ 54943562-e7ae-4529-bd54-31840fb96a36
data_location144m = joinpath(@__DIR__, "..", "data144m");

# ╔═╡ 59dd328c-6fac-4bf1-907d-10f3d96ec1a6
data_location108m = joinpath(@__DIR__, "..", "data108m");

# ╔═╡ a3645377-b459-43c8-922d-3c55fc6b3eb3
data_location72m = joinpath(@__DIR__, "..", "data72m");

# ╔═╡ b3256657-1252-442d-a5d6-c22bcc417509
data_location36m = joinpath(@__DIR__, "..", "data36m");

# ╔═╡ cb7cf756-0197-4a5d-a02f-65f9fa27ff40
data_location0 = joinpath(@__DIR__, "..", "data0");

# ╔═╡ 463d3c0b-606d-42d1-9ad1-105495e3b2b1
data_location36 = joinpath(@__DIR__, "..", "data36");

# ╔═╡ efb7c8db-023d-43b8-95cf-a3ca899b434f
data_location72 = joinpath(@__DIR__, "..", "data72");

# ╔═╡ a0797ad5-2108-4714-b51f-49805d762a26
data_location108 = joinpath(@__DIR__, "..", "data108");

# ╔═╡ 7e14c341-1c17-4145-8068-e874e7be7cdd
data_location144 = joinpath(@__DIR__, "..", "data144");

# ╔═╡ 9aad20b9-54cf-4584-b180-0ab8feb63f7d
data_location180 = joinpath(@__DIR__, "..", "data180");

# ╔═╡ bda2fdc3-fae6-449c-ab34-2eda148d826f
file_path144m = joinpath(data_location144m, "unbinned_data144m.txt");

# ╔═╡ 81b2984d-20a8-4f96-bd5e-c923d9d00c57
file_path108m = joinpath(data_location108m, "unbinned_data108m.txt");

# ╔═╡ fab8d06d-9a0a-450f-8e9f-b434c3655989
file_path72m = joinpath(data_location72m, "unbinned_data72m.txt");

# ╔═╡ 8a2ba304-b4f5-4ad5-a12e-c17d2e56bc39
file_path36m = joinpath(data_location36m, "unbinned_data36m.txt");

# ╔═╡ aad83f72-0873-4d46-bcc7-10b1961e9844
file_path0 = joinpath(data_location0, "unbinned_data0.txt");

# ╔═╡ d1088a6f-d962-490c-ad50-45fcc8774ebc
file_path36 = joinpath(data_location36, "unbinned_data36.txt");

# ╔═╡ e19f626c-546b-44ab-9953-d46bc635d315
file_path72 = joinpath(data_location72, "unbinned_data72.txt");

# ╔═╡ 2b1c7c62-284a-44e4-b242-1a069c03f178
file_path108 = joinpath(data_location108, "unbinned_data108.txt");

# ╔═╡ 5df7c82a-cf79-4de5-9ac8-2e2bf771bee1
file_path144 = joinpath(data_location144, "unbinned_data144.txt");

# ╔═╡ 75c2cdfc-2d9e-4c30-b28e-75843fecacb0
file_path180 = joinpath(data_location180, "unbinned_data180.txt");

# ╔═╡ fca21510-369d-48fd-8905-c68c79561602
data_points144m = readlines(file_path144m) .|> x -> parse(Float64, x);

# ╔═╡ 379b4299-6df7-4216-ab35-3be396ecc6c3
data_points108m = readlines(file_path108m) .|> x -> parse(Float64, x);

# ╔═╡ b03d25f4-96e9-49ad-869a-e22dcb2303f2
data_points72m = readlines(file_path72m) .|> x -> parse(Float64, x);

# ╔═╡ 88a03e86-5470-4ba6-bc26-e774603978cb
data_points36m = readlines(file_path36m) .|> x -> parse(Float64, x);

# ╔═╡ 0bd617f3-942b-4555-9b57-9a7f9a44c51e
data_points0 = readlines(file_path0) .|> x -> parse(Float64, x);

# ╔═╡ 927d3bc7-5f39-451a-bf92-6680a670ea7c
data_points36 = readlines(file_path36) .|> x -> parse(Float64, x);

# ╔═╡ 9a6e973f-d6c8-4aae-8711-70d7fe5e12d8
data_points72 = readlines(file_path72) .|> x -> parse(Float64, x);

# ╔═╡ 09807f07-1677-40dd-80e1-db5b52b22781
data_points108 = readlines(file_path108) .|> x -> parse(Float64, x);

# ╔═╡ 680fa5ff-b7b7-4fae-a144-82b94acb6d0e
data_points144 = readlines(file_path144) .|> x -> parse(Float64, x);

# ╔═╡ 73dec207-1f13-4138-b551-eeb6a93e63e7
data_points180 = readlines(file_path180) .|> x -> parse(Float64, x);

# ╔═╡ aaad123c-17a4-4098-a3bb-ce48cb4426fc
theme(:boxed, ylims=(0, :auto))

# ╔═╡ d0f0452f-5835-4bb6-aee8-6dcfa9a2e2a2
md"""
## Read data
hist from file
"""

# ╔═╡ e873383c-163b-423e-a9ab-2588180e325a
specific_hist144m = Hist1D(data_points144m, binedges=range(4.44, 4.49, 500));

# ╔═╡ a5a8b38c-d510-4341-a57f-4dc219ed2483
specific_hist108m = Hist1D(data_points108m, binedges=range(4.44, 4.49, 500));

# ╔═╡ 65b20902-c77c-4d15-8ad4-a5596c049d57
specific_hist72m = Hist1D(data_points72m, binedges=range(4.44, 4.49, 500));

# ╔═╡ 3bc4bbf8-90fa-4486-810f-3b39c35756a8
specific_hist36m = Hist1D(data_points36m, binedges=range(4.44, 4.49, 500));

# ╔═╡ 665ff6c9-d847-42a8-96f7-77a36f55352d
specific_hist0 = Hist1D(data_points0, binedges=range(4.44, 4.49, 500));

# ╔═╡ f735c700-13f2-4e5e-8afa-b4cc6f0aae19
specific_hist36 = Hist1D(data_points36, binedges=range(4.44, 4.49, 500));

# ╔═╡ 912d1ddd-88d5-4671-bcf1-9dbf631c28bd
specific_hist72 = Hist1D(data_points72, binedges=range(4.44, 4.49, 500));

# ╔═╡ 74df8a15-21f7-4487-ba46-80a3bd2d35cc
specific_hist108 = Hist1D(data_points108, binedges=range(4.44, 4.49, 500));

# ╔═╡ 561dba23-d72d-4bd3-bc7f-2af365c4da9b
specific_hist144 = Hist1D(data_points144, binedges=range(4.44, 4.49, 500));

# ╔═╡ 8b0dbd51-6813-4f89-8610-86c1cfa8bd48
specific_hist180 = Hist1D(data_points180, binedges=range(4.44, 4.49, 500));

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

# ╔═╡ cbbb18a8-730b-472c-91bc-75a3fb4bff55
begin
	@with_kw struct BWCoherentIncoherent144m{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent144m, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = -144
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent144m(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent144m)}(values)
		BWCoherentIncoherent144m(; pars...)
	end
end

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

# ╔═╡ db382f7b-656d-4346-9c72-9b00ef814faa
begin
	@with_kw struct BWCoherentIncoherent72m{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent72m, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = -72
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent72m(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent72m)}(values)
		BWCoherentIncoherent72m(; pars...)
	end
end

# ╔═╡ 1cc514de-e15d-47e3-a5d5-9a73df7a4609
begin
	@with_kw struct BWCoherentIncoherent36m{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent36m, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = -36
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent36m(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent36m)}(values)
		BWCoherentIncoherent36m(; pars...)
	end
end

# ╔═╡ d10eed86-6ec4-48cd-8514-717c02a6f694
begin
	@with_kw struct BWCoherentIncoherent0{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent0, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
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

# ╔═╡ 1b2f3ce9-948d-4915-9870-c893e8df2bf0
begin
	@with_kw struct BWCoherentIncoherent72{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent72, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
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

# ╔═╡ d32a091e-212a-40e0-821d-30e24d69dfa4
begin
	@with_kw struct BWCoherentIncoherent108{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent108, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
		@unpack c,b,ϕB = model
		cB = c*bkg
		iB = b*bkg
		Δϕ_deg = 108
		Δϕ = Δϕ_deg/180*π
		abs2(a * cis(Δϕ) + cB * cis(ϕB)) + abs2(iB)
	end
	function BWCoherentIncoherent108(values::Vector)
		pars = NamedTuple{fieldnames(BWCoherentIncoherent108)}(values)
		BWCoherentIncoherent108(; pars...)
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

# ╔═╡ 7f061105-29de-4ffd-bcd9-c64bba01e55c
begin
	@with_kw struct BWCoherentIncoherent180{T}
		a0::T # norm of BW
		m0::T # mass of BW
		Γ0::T # width of BW
		ϕB::T # phase of background
		c::T # coherent backgrond
		b::T # incoh
		α::T # bkg param
	end
	function intensity(model::BWCoherentIncoherent180, m)
		@unpack a0, m0, Γ0 = model
		BW = BreitWigner(m0, Γ0)
		a = BW(m^2)*m0*Γ0 * a0
		@unpack α = model
		bkg = α*m^2 + 1 
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

# ╔═╡ b7c7b6dd-467c-43d4-957e-86119c239530
const specific_data144m = (; 
	xv = bincenters(specific_hist144m),
	yv = specific_hist144m.bincounts,
	δyv = sqrt.(specific_hist144m.bincounts))

# ╔═╡ fa95a3e2-be56-42ba-a528-91fe73071bff
const specific_data108m = (; 
	xv = bincenters(specific_hist108m),
	yv = specific_hist108m.bincounts,
	δyv = sqrt.(specific_hist108m.bincounts))

# ╔═╡ 39868266-88df-4748-a8c9-013473303f81
const specific_data72m = (; 
	xv = bincenters(specific_hist72m),
	yv = specific_hist72m.bincounts,
	δyv = sqrt.(specific_hist72m.bincounts))

# ╔═╡ 5db8a51a-42cf-46f1-96cf-75ff7ade52fc
const specific_data36m = (; 
	xv = bincenters(specific_hist36m),
	yv = specific_hist36m.bincounts,
	δyv = sqrt.(specific_hist36m.bincounts))

# ╔═╡ 47fec08a-889e-40bf-8f28-2d1f2cd9da86
const specific_data0 = (; 
	xv = bincenters(specific_hist0),
	yv = specific_hist0.bincounts,
	δyv = sqrt.(specific_hist0.bincounts))

# ╔═╡ 6f7247ae-8c94-4104-ac35-23665410d475
const specific_data36 = (; 
	xv = bincenters(specific_hist36),
	yv = specific_hist36.bincounts,
	δyv = sqrt.(specific_hist36.bincounts))

# ╔═╡ 7ad95e1d-6bb7-400e-b245-c2d00673b778
const specific_data72 = (; 
	xv = bincenters(specific_hist72),
	yv = specific_hist72.bincounts,
	δyv = sqrt.(specific_hist72.bincounts))

# ╔═╡ 15b7bc0e-1b15-44d1-af87-12d9113f28c0
const specific_data108 = (; 
	xv = bincenters(specific_hist108),
	yv = specific_hist108.bincounts,
	δyv = sqrt.(specific_hist108.bincounts))

# ╔═╡ dbfad0c7-d086-4f4e-82ef-5834deeacca0
const specific_data144 = (; 
	xv = bincenters(specific_hist144),
	yv = specific_hist144.bincounts,
	δyv = sqrt.(specific_hist144.bincounts))

# ╔═╡ 81cf986a-7a84-4866-ad34-f30db9eee54b
const specific_data180 = (; 
	xv = bincenters(specific_hist180),
	yv = specific_hist180.bincounts,
	δyv = sqrt.(specific_hist180.bincounts))

# ╔═╡ 02b27549-178f-4b21-8048-5c4c2a48a937
specific_datasets = [specific_data144m, specific_data108m, specific_data72m, specific_data36m, specific_data0, specific_data36, specific_data72, specific_data108, specific_data144, specific_data180]

# ╔═╡ 6dfbda41-b976-4138-8834-fad3155e19fa
specific_models(pars) = [BWCoherentIncoherent144m(pars), BWCoherentIncoherent108m(pars), BWCoherentIncoherent72m(pars), BWCoherentIncoherent36m(pars), BWCoherentIncoherent0(pars), BWCoherentIncoherent36(pars), BWCoherentIncoherent72(pars), BWCoherentIncoherent108(pars), BWCoherentIncoherent144(pars), BWCoherentIncoherent180(pars)]

# ╔═╡ 22ea68a4-d8e1-4b66-add4-2d270a614e23
loss(pars) = compute_loss(specific_models(pars), specific_datasets)

# ╔═╡ b7fcac1a-5e12-4c2e-a380-bab8ce0859e9
@assert loss(collect(initial_pars)) isa Number

# ╔═╡ aef67dfa-6775-4220-9a9b-151891609bdf
fit_result = optimize(loss, collect(initial_pars), BFGS())

# ╔═╡ 46b0f126-f613-48ba-aba6-be965c66f865
fit_result.minimum

# ╔═╡ 3a71d178-dedf-4808-95dc-121bb57489ba
best_model144m = BWCoherentIncoherent144m(fit_result.minimizer)

# ╔═╡ 88f88f6e-ff7f-441d-96ae-69269db978b7
best_model108m = BWCoherentIncoherent108m(fit_result.minimizer)

# ╔═╡ 5717c7c9-3c7a-4cc7-b1ab-ed96d6517e94
best_model72m = BWCoherentIncoherent72m(fit_result.minimizer)

# ╔═╡ 139fae33-e3fb-4116-abaa-b188790ba6aa
best_model36m = BWCoherentIncoherent36m(fit_result.minimizer)

# ╔═╡ c26ed24e-24a0-4646-8863-7c015ab21273
best_model0 = BWCoherentIncoherent0(fit_result.minimizer)

# ╔═╡ 2306ae6a-281f-45c9-8f67-2f73da145935
best_model36 = BWCoherentIncoherent36(fit_result.minimizer)

# ╔═╡ c9f9972d-a2b5-4b76-b3d0-edbd6639b6f0
best_model72 = BWCoherentIncoherent72(fit_result.minimizer)

# ╔═╡ 9b65a7da-bf38-430b-a6d0-ab37cea9284a
best_model108 = BWCoherentIncoherent108(fit_result.minimizer)

# ╔═╡ 3c554d00-f68a-40d5-88ac-f7a9912f3dce
best_model144 = BWCoherentIncoherent144(fit_result.minimizer)

# ╔═╡ b7d844b6-9f30-4292-977c-3454398fd63a
best_model180 = BWCoherentIncoherent180(fit_result.minimizer)

# ╔═╡ 1837cf4c-18a9-4735-b066-eaca988b92d4
ratio = best_model144.c / best_model144.b

# ╔═╡ 9fdc8553-14b4-4922-ad20-871d16f03ddd
coh = best_model144.c / (best_model144.b + best_model144.c)

# ╔═╡ 520d9885-c1ed-47a3-b853-dd657ed45315
incoh = best_model144.b / (best_model144.b + best_model144.c)

# ╔═╡ e9c36bcf-c852-46eb-8ea7-f33208dbddf1
bkg_angle = best_model144.ϕB % π

# ╔═╡ 9f28f268-c44a-4992-a9e4-a0472d329989
bkg_angle_deg = (best_model144.ϕB % π)*180/π

# ╔═╡ be6cdff9-1a47-4b38-813a-a8f16eb9e03e
hessian_matrix = FiniteDiff.finite_difference_hessian(loss, fit_result.minimizer)

# ╔═╡ e511b7e8-ac45-4c8a-a0d5-d8bfd73a43fd
cov_matrix = inv(hessian_matrix)

# ╔═╡ de900b00-4666-4f8c-a3e0-6998e64c5a10
param_uncertainties = sqrt.(diag(cov_matrix))

# ╔═╡ df58764e-6ee5-48dd-86ec-cd4496cd5f57
let
	xlims = (specific_hist144m.binedges[1][1], specific_hist144m.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist144m), specific_hist144m.bincounts, yerr=sqrt.(specific_hist144m.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model144m, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent144m(best_model144m; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model144m, m) - intensity(BWCoherentIncoherent144m(best_model144m; c=0), m) + intensity(BWCoherentIncoherent144m(best_model144m; b=0, c=0), best_model144m.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model144m, m) - intensity(BWCoherentIncoherent144m(best_model144m; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ f0615806-2d24-44b1-b162-49f8759a5f98
let
	xlims = (specific_hist108m.binedges[1][1], specific_hist108m.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist108m), specific_hist108m.bincounts, yerr=sqrt.(specific_hist108m.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model108m, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent108m(best_model108m; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model108m, m) - intensity(BWCoherentIncoherent108m(best_model108m; c=0), m) + intensity(BWCoherentIncoherent108m(best_model108m; b=0, c=0), best_model108m.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model108m, m) - intensity(BWCoherentIncoherent108m(best_model108m; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ d2637790-ac20-4800-ba6a-16bd3f22d8e0
let
	xlims = (specific_hist72m.binedges[1][1], specific_hist72m.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist72m), specific_hist72m.bincounts, yerr=sqrt.(specific_hist72m.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model72m, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent72m(best_model72m; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model72m, m) - intensity(BWCoherentIncoherent72m(best_model72m; c=0), m) + intensity(BWCoherentIncoherent72m(best_model72m; b=0, c=0), best_model72m.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model72m, m) - intensity(BWCoherentIncoherent72m(best_model72m; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ ac8705dc-0d20-4701-ab31-ff1c9a6095b2
let
	xlims = (specific_hist36m.binedges[1][1], specific_hist36m.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist36m), specific_hist36m.bincounts, yerr=sqrt.(specific_hist36m.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model36m, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent36m(best_model36m; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model36m, m) - intensity(BWCoherentIncoherent36m(best_model36m; c=0), m) + intensity(BWCoherentIncoherent36m(best_model36m; b=0, c=0), best_model36m.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model36m, m) - intensity(BWCoherentIncoherent36m(best_model36m; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ b3cef07e-087d-404e-aba7-71d16cdf530e
let
	xlims = (specific_hist0.binedges[1][1], specific_hist0.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist0), specific_hist0.bincounts, yerr=sqrt.(specific_hist0.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model0, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent0(best_model0; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model0, m) - intensity(BWCoherentIncoherent0(best_model0; c=0), m) + intensity(BWCoherentIncoherent0(best_model0; b=0, c=0), best_model0.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model0, m) - intensity(BWCoherentIncoherent0(best_model0; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ 7c1a91f4-73e8-4e72-b012-17c02e1aab9b
let
	xlims = (specific_hist36.binedges[1][1], specific_hist36.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist36), specific_hist36.bincounts, yerr=sqrt.(specific_hist36.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model36, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent36(best_model36; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model36, m) - intensity(BWCoherentIncoherent36(best_model36; c=0), m) + intensity(BWCoherentIncoherent36(best_model36; b=0, c=0), best_model36.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model36, m) - intensity(BWCoherentIncoherent36(best_model36; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ 33afd217-abcf-4326-982b-3c606ba618b0
let
	xlims = (specific_hist72.binedges[1][1], specific_hist72.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist72), specific_hist72.bincounts, yerr=sqrt.(specific_hist72.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model72, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent72(best_model72; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model72, m) - intensity(BWCoherentIncoherent72(best_model72; c=0), m) + intensity(BWCoherentIncoherent72(best_model72; b=0, c=0), best_model72.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model72, m) - intensity(BWCoherentIncoherent72(best_model72; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ 8896fdc0-222c-4a05-97b4-d6094bb89dd8
let
	xlims = (specific_hist108.binedges[1][1], specific_hist108.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist108), specific_hist108.bincounts, yerr=sqrt.(specific_hist108.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model108, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent108(best_model108; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model108, m) - intensity(BWCoherentIncoherent108(best_model108; c=0), m) + intensity(BWCoherentIncoherent108(best_model108; b=0, c=0), best_model108.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model108, m) - intensity(BWCoherentIncoherent108(best_model108; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ 3aa658ec-479e-476a-ad53-ae431d674152
let
	xlims = (specific_hist144.binedges[1][1], specific_hist144.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist144), specific_hist144.bincounts, yerr=sqrt.(specific_hist144.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model144, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent144(best_model144; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model144, m) - intensity(BWCoherentIncoherent144(best_model144; c=0), m) + intensity(BWCoherentIncoherent144(best_model144; b=0, c=0), best_model144.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model144, m) - intensity(BWCoherentIncoherent144(best_model144; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ 143549a1-4fa3-4f01-877c-f9b913f568d9
let
	xlims = (specific_hist180.binedges[1][1], specific_hist180.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist180), specific_hist180.bincounts, yerr=sqrt.(specific_hist180.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model180, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent180(best_model180; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model180, m) - intensity(BWCoherentIncoherent180(best_model180; c=0), m) + intensity(BWCoherentIncoherent180(best_model180; b=0, c=0), best_model180.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model180, m) - intensity(BWCoherentIncoherent180(best_model180; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461 GeV")
	plot!(legend=:right)
end

# ╔═╡ a79425db-298e-445c-8b04-594e7b6fbdbf
block_of_text = """
let
	xlims = (specific_hist72m.binedges[1][1], specific_hist72m.binedges[1][end])
	plot()
	scatter!(bincenters(specific_hist72m), specific_hist72m.bincounts, yerr=sqrt.(specific_hist72m.bincounts), ylim=(0,3000))
	plot!(m->intensity(best_model72m, m), xlims..., lw=2)
	plot!(xlims..., fill=0, alpha=0.3, label="BW") do m
		intensity(BWCoherentIncoherent72m(best_model72m; c=0, b=0), m)
	end
	plot!(xlims..., lw=2, label="coh bkg") do m
        intensity(best_model72m, m) - intensity(BWCoherentIncoherent72m(best_model72m; c=0), m) + intensity(BWCoherentIncoherent72m(best_model72m; b=0, c=0), best_model72m.m0)
    end
	plot!(xlims..., lw=2, color=:red, label="inc bkg") do m
        intensity(best_model72m, m) - intensity(BWCoherentIncoherent72m(best_model72m; b=0), m)
    end
	vline!([4.461], color=:black, linestyle=:dash, lw=0.5, label="m = 4.461")
end
"""

# List containing ten different replacement elements
#=replacement_list = ["101m", "102m", "103m", "104m", "105m", "106m", "107m", "109m", "110m", "111m"]

# Loop through the replacement list and generate modified versions of the text
output_texts = []
for replacement in replacement_list
    modified_text = replace(block_of_text, "108m" => replacement)
    push!(output_texts, modified_text)
end

# Join all modified versions into a single output with separation
final_output = join(output_texts, "\n\n")

println(final_output)=#

# ╔═╡ 595ee368-7d36-40ca-9e97-a1be6984e58f
replacement_list = ["144m", "108m", "72m", "36m", "0", "36", "72", "108", "144", "180"]

# ╔═╡ 481bbab7-9e77-48bb-a2dc-76f8d7a01eac
output_texts = []

# ╔═╡ 758d2de9-af63-4da1-8a55-752c83c60eb3
for replacement in replacement_list
    modified_text = replace(block_of_text, "72m" => replacement)
    push!(output_texts, modified_text)
end

# ╔═╡ 8c94905d-afa9-439c-9a27-eadb756c4838
final_output = join(output_texts, "\n\n")

# ╔═╡ e040db72-52fa-4898-a36c-7efa021e4fad
println(final_output)

# ╔═╡ Cell order:
# ╠═ec7ef8e2-f9a8-11ef-2e7c-0f562efc23a1
# ╠═54943562-e7ae-4529-bd54-31840fb96a36
# ╠═59dd328c-6fac-4bf1-907d-10f3d96ec1a6
# ╠═a3645377-b459-43c8-922d-3c55fc6b3eb3
# ╠═b3256657-1252-442d-a5d6-c22bcc417509
# ╠═cb7cf756-0197-4a5d-a02f-65f9fa27ff40
# ╠═463d3c0b-606d-42d1-9ad1-105495e3b2b1
# ╠═efb7c8db-023d-43b8-95cf-a3ca899b434f
# ╠═a0797ad5-2108-4714-b51f-49805d762a26
# ╠═7e14c341-1c17-4145-8068-e874e7be7cdd
# ╠═9aad20b9-54cf-4584-b180-0ab8feb63f7d
# ╠═bda2fdc3-fae6-449c-ab34-2eda148d826f
# ╠═81b2984d-20a8-4f96-bd5e-c923d9d00c57
# ╠═fab8d06d-9a0a-450f-8e9f-b434c3655989
# ╠═8a2ba304-b4f5-4ad5-a12e-c17d2e56bc39
# ╠═aad83f72-0873-4d46-bcc7-10b1961e9844
# ╠═d1088a6f-d962-490c-ad50-45fcc8774ebc
# ╠═e19f626c-546b-44ab-9953-d46bc635d315
# ╠═2b1c7c62-284a-44e4-b242-1a069c03f178
# ╠═5df7c82a-cf79-4de5-9ac8-2e2bf771bee1
# ╠═75c2cdfc-2d9e-4c30-b28e-75843fecacb0
# ╠═fca21510-369d-48fd-8905-c68c79561602
# ╠═379b4299-6df7-4216-ab35-3be396ecc6c3
# ╠═b03d25f4-96e9-49ad-869a-e22dcb2303f2
# ╠═88a03e86-5470-4ba6-bc26-e774603978cb
# ╠═0bd617f3-942b-4555-9b57-9a7f9a44c51e
# ╠═927d3bc7-5f39-451a-bf92-6680a670ea7c
# ╠═9a6e973f-d6c8-4aae-8711-70d7fe5e12d8
# ╠═09807f07-1677-40dd-80e1-db5b52b22781
# ╠═680fa5ff-b7b7-4fae-a144-82b94acb6d0e
# ╠═73dec207-1f13-4138-b551-eeb6a93e63e7
# ╠═aaad123c-17a4-4098-a3bb-ce48cb4426fc
# ╟─d0f0452f-5835-4bb6-aee8-6dcfa9a2e2a2
# ╠═e873383c-163b-423e-a9ab-2588180e325a
# ╠═a5a8b38c-d510-4341-a57f-4dc219ed2483
# ╠═65b20902-c77c-4d15-8ad4-a5596c049d57
# ╠═3bc4bbf8-90fa-4486-810f-3b39c35756a8
# ╠═665ff6c9-d847-42a8-96f7-77a36f55352d
# ╠═f735c700-13f2-4e5e-8afa-b4cc6f0aae19
# ╠═912d1ddd-88d5-4671-bcf1-9dbf631c28bd
# ╠═74df8a15-21f7-4487-ba46-80a3bd2d35cc
# ╠═561dba23-d72d-4bd3-bc7f-2af365c4da9b
# ╠═8b0dbd51-6813-4f89-8610-86c1cfa8bd48
# ╠═7cd14ac4-b688-4336-8532-5c52c3f9600c
# ╠═2cc85370-c066-45e3-8a51-ce566ac592a0
# ╠═290b2ea3-3aa2-452e-aa45-22a02f42ccba
# ╟─4edee112-1541-4ea6-a6a2-5a738742760f
# ╠═cbbb18a8-730b-472c-91bc-75a3fb4bff55
# ╠═ff42bba7-a2c2-4951-a9ce-ba4058578b39
# ╠═db382f7b-656d-4346-9c72-9b00ef814faa
# ╠═1cc514de-e15d-47e3-a5d5-9a73df7a4609
# ╠═d10eed86-6ec4-48cd-8514-717c02a6f694
# ╠═44cb256e-2eb2-4d8c-b3e8-2e694292e2f5
# ╠═1b2f3ce9-948d-4915-9870-c893e8df2bf0
# ╠═d32a091e-212a-40e0-821d-30e24d69dfa4
# ╠═f1e5fa90-e711-4d12-8c37-081ab5f5de67
# ╠═7f061105-29de-4ffd-bcd9-c64bba01e55c
# ╟─1fdfa846-e0e0-4a2f-881b-49cbcf1c1289
# ╠═6db468ec-724b-4546-9a29-899b8da8706a
# ╟─7ac3cb4f-b874-4cb9-9da8-c874943adbf2
# ╠═5f4a600a-eca9-4d42-b48a-307838a14b25
# ╠═6349f6d9-6a29-4356-bc72-d180d6b8585f
# ╟─111dbc4b-66b5-4cd1-bab4-f025262fb963
# ╠═b7c7b6dd-467c-43d4-957e-86119c239530
# ╠═fa95a3e2-be56-42ba-a528-91fe73071bff
# ╠═39868266-88df-4748-a8c9-013473303f81
# ╠═5db8a51a-42cf-46f1-96cf-75ff7ade52fc
# ╠═47fec08a-889e-40bf-8f28-2d1f2cd9da86
# ╠═6f7247ae-8c94-4104-ac35-23665410d475
# ╠═7ad95e1d-6bb7-400e-b245-c2d00673b778
# ╠═15b7bc0e-1b15-44d1-af87-12d9113f28c0
# ╠═dbfad0c7-d086-4f4e-82ef-5834deeacca0
# ╠═81cf986a-7a84-4866-ad34-f30db9eee54b
# ╠═02b27549-178f-4b21-8048-5c4c2a48a937
# ╠═6dfbda41-b976-4138-8834-fad3155e19fa
# ╠═22ea68a4-d8e1-4b66-add4-2d270a614e23
# ╠═b7fcac1a-5e12-4c2e-a380-bab8ce0859e9
# ╠═aef67dfa-6775-4220-9a9b-151891609bdf
# ╠═46b0f126-f613-48ba-aba6-be965c66f865
# ╠═3a71d178-dedf-4808-95dc-121bb57489ba
# ╠═88f88f6e-ff7f-441d-96ae-69269db978b7
# ╠═5717c7c9-3c7a-4cc7-b1ab-ed96d6517e94
# ╠═139fae33-e3fb-4116-abaa-b188790ba6aa
# ╠═c26ed24e-24a0-4646-8863-7c015ab21273
# ╠═2306ae6a-281f-45c9-8f67-2f73da145935
# ╠═c9f9972d-a2b5-4b76-b3d0-edbd6639b6f0
# ╠═9b65a7da-bf38-430b-a6d0-ab37cea9284a
# ╠═3c554d00-f68a-40d5-88ac-f7a9912f3dce
# ╠═b7d844b6-9f30-4292-977c-3454398fd63a
# ╠═1837cf4c-18a9-4735-b066-eaca988b92d4
# ╠═9fdc8553-14b4-4922-ad20-871d16f03ddd
# ╠═520d9885-c1ed-47a3-b853-dd657ed45315
# ╠═e9c36bcf-c852-46eb-8ea7-f33208dbddf1
# ╠═9f28f268-c44a-4992-a9e4-a0472d329989
# ╠═be6cdff9-1a47-4b38-813a-a8f16eb9e03e
# ╠═e511b7e8-ac45-4c8a-a0d5-d8bfd73a43fd
# ╠═de900b00-4666-4f8c-a3e0-6998e64c5a10
# ╠═df58764e-6ee5-48dd-86ec-cd4496cd5f57
# ╠═f0615806-2d24-44b1-b162-49f8759a5f98
# ╠═d2637790-ac20-4800-ba6a-16bd3f22d8e0
# ╠═ac8705dc-0d20-4701-ab31-ff1c9a6095b2
# ╠═b3cef07e-087d-404e-aba7-71d16cdf530e
# ╠═7c1a91f4-73e8-4e72-b012-17c02e1aab9b
# ╠═33afd217-abcf-4326-982b-3c606ba618b0
# ╠═8896fdc0-222c-4a05-97b4-d6094bb89dd8
# ╠═3aa658ec-479e-476a-ad53-ae431d674152
# ╠═143549a1-4fa3-4f01-877c-f9b913f568d9
# ╠═a79425db-298e-445c-8b04-594e7b6fbdbf
# ╠═595ee368-7d36-40ca-9e97-a1be6984e58f
# ╠═481bbab7-9e77-48bb-a2dc-76f8d7a01eac
# ╠═758d2de9-af63-4da1-8a55-752c83c60eb3
# ╠═8c94905d-afa9-439c-9a27-eadb756c4838
# ╠═e040db72-52fa-4898-a36c-7efa021e4fad
