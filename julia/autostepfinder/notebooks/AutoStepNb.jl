### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ac949aea-be3f-404d-8a8b-22c25b26ee38
using Pkg; Pkg.activate("/Users/pabloherrero/sabat/AutoStepFinder/julia/autostepfinder")

# ╔═╡ a03bd74c-dc79-4f88-a142-815a6e94954c
begin
	Pkg.add("DelimitedFiles")
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	Pkg.add("Images")
	Pkg.add("ImageView")
	Pkg.add("Colors")
	Pkg.add("FixedPointNumbers")
	Pkg.add("Statistics")
	Pkg.add("Revise")
	#Pkg.add("GR")
end

# ╔═╡ 980dcf83-1e5c-4559-a42b-b02b0db094b3
# ╠═╡ show_logs = false
begin
	using DelimitedFiles
	using Plots
	using PlutoUI
	using Images
	using ImageView
	using Colors
	using FixedPointNumbers
	using Statistics
	using Revise
	#using GR
end

# ╔═╡ ee654cd8-b63a-4e5e-96f8-3104bea4fa60
Pkg.instantiate()

# ╔═╡ 16cbcb03-5dc5-4219-8337-c3024d805338
# ╠═╡ show_logs = false
Pkg.resolve()

# ╔═╡ d09a5620-570a-4f16-9e0d-f24fb6ba1fbc
Pkg.status()

# ╔═╡ 933774f1-9a00-4ee0-98e0-980f9e8fb930
Pkg.precompile()

# ╔═╡ 99ea72eb-b944-4d2a-9347-f4d9255ac790
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ c22c5abe-980b-11ee-268c-b16eb65f35f4
# ╠═╡ show_logs = false
asf = ingredients("../src/autostepfinder.jl")

# ╔═╡ d3da2ac5-df51-4d7e-99f5-039595f6f75e
PlutoUI.TableOfContents(title="AutoStepfinder Notebook", indent=true)

# ╔═╡ b2b06a68-dcd5-43e1-9649-8574423e30ec
md"""# Plot image"""

# ╔═╡ 74183602-063d-4aa9-884e-aff77e9606d9
path_free = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/ZFFL62_NAPH3_5e-8M/Microscope/Wide field/no_barium/35/Field1"

# ╔═╡ f94daecf-afda-4b61-a5d1-3ae48c035480
begin
	img = asf.autostepfinder.get_image(path_free, 23)
	heatmap(img)
end

# ╔═╡ 12b432f4-5942-4958-9de8-473a08c364a5
begin
	roi1 = [190,290,300,400]
	evol = asf.autostepfinder.temporal_evolution(path_free, 400, roi1);
	nothing
end

# ╔═╡ 10b8263d-8dc5-4ddc-bd21-b310ac7fe11b
begin 
	imroi = asf.autostepfinder.get_image(path_free, 23)[roi1[1]:roi1[2],roi1[3]:roi1[4]]
	gr()
	#c = cgrad([:red,:yellow,:green], [0.01, 0.9995], categorical = true)

	heatmap(imroi, clim=(1800,2500))
	#display(fig)
end

# ╔═╡ 8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
md"""# Analyse large set of data (free NAPH3)"""

# ╔═╡ bd225196-cea7-490d-b8bf-d515f51a31a3
md""" Select Threshold $(@bind tresH NumberField(0.1:10., default=0.1))"""

# ╔═╡ ccab705b-5193-4618-89f1-ced300100ba3
md""" Check to analyze steps using full ROI: $(@bind zroi CheckBox())"""

# ╔═╡ 680641ef-16b5-4fca-9a57-149bf8c72eb8
md"""## Plot one evolution"""

# ╔═╡ 39a8c1dc-7a7e-4958-a997-039d2eeb20bd
md"""# Analyse Naph3+Ba """

# ╔═╡ 2ba85012-706a-4719-ab87-e0fa91d2f36f
path_ba = "/Users/pabloherrero/DIPC Dropbox/BOLD/Data1/tempo_for_Mikel/commercial_molecules/Photobleaching_1mW_100ms/Field3_long_sequence"

# ╔═╡ 5954a7be-70c1-4d23-b4c2-719f23220212
begin
	imgba = asf.autostepfinder.get_image(path_ba, 23)
	heatmap(imgba,  clim=(1800,2500))
end

# ╔═╡ 4173f292-644b-4a87-b052-c9cfc8bc2765
roi1

# ╔═╡ efc103b6-8695-4daf-9647-2fd1baafc8c8
begin
	roiba = [200,300,200,300]
	evol_ba = asf.autostepfinder.temporal_evolution(path_ba, 400, roi1);
	nothing
end

# ╔═╡ 9d28bc4a-1bcf-44ce-9429-2bc74e17f9ac
begin 
	imroiba = asf.autostepfinder.get_image(path_ba, 2)[roiba[1]:roiba[2],roiba[3]:roiba[4]]
	gr()
	#c = cgrad([:red,:yellow,:green], [0.01, 0.9995], categorical = true)

	heatmap(imroiba, clim=(1800,2500))
	#display(fig)
end

# ╔═╡ 7bd7add8-232a-4734-a272-84456af2496e
md""" Check to analyze steps using full ROI: $(@bind zroiba CheckBox())"""

# ╔═╡ b6a31f37-9600-44f7-80da-a48ad6e9986a


# ╔═╡ 54ec44a2-189e-4b23-8afd-f562378f8b7d
md"""# Compare free/Ba histograms"""

# ╔═╡ 7bb25cd5-0f74-48ce-aeee-3020221f6b9b


# ╔═╡ e477b05f-2fbb-4bb8-a37a-a56c286445dc
md"""# Functions """

# ╔═╡ 7c6dcca9-fd13-4977-8bf7-7f0f5d2b5636
function get_hist_roi(evol, tresH = 0.1, N_iter = 50)
	hist_height = Vector{Any}()
	hist_widths = Vector{Any}()
	n_steps = Vector{Int}()
	
	I, J, _ = size(evol)
	ijyes = reshape([], 0, 2)
	
	for i in 1:I
	    for j in 1:J
			dataX = evol[ i, j, :]
			
			
			S_curve, best_shots, Fits, steptable = asf.autostepfinder.AutoStepMain(dataX, tresH, N_iter)
			
	            if best_shots > 0 && S_curve[best_shots] > tresH 

					LST = size(steptable)[1]
					if LST > 2
						push!(n_steps, size(steptable, 1))
						push!(hist_height, steptable[1:LST-2, 4])
						push!(hist_widths, steptable[1:LST-2, 5])
						
						ijyes = [ijyes;  [i j]]
						
					end
	                    
	            end

			end
	    end
	
	hw = vcat(hist_widths...)
	hh = vcat(hist_height...)
	return hh, hw, n_steps, ijyes
end

# ╔═╡ 78f8ff68-ec80-4988-a240-064b1e2bc93c
if zroi
	hh, hw, n_steps, ijy = get_hist_roi(evol, tresH, 40)
end

# ╔═╡ 4cbbadc4-9948-4a21-ab03-18ed7c73573d
begin 
	p1 = histogram(n_steps, label="Number of steps")
	p2 = histogram(hh, label="Step height")
	p3 = histogram(hw, label="Step width")
	plot(p1, p2, p3, layout=(1,3))
end

# ╔═╡ b1076909-7591-41aa-bada-4747616cec67
md""" Select pixel: $(@bind ijview NumberField(1:size(ijy)[1], default=1))"""

# ╔═╡ 9a6e62bb-289e-4ca7-9f95-cd0caaac7700
begin
	pxi, pxj = ijy[ijview, :]
	asf.autostepfinder.plot_pixel(evol, pxi, pxj, tresH)
end

# ╔═╡ bd1320fb-91d1-448c-be04-f60807e96ba0
if zroiba
	hh_ba, hw_ba, n_steps_ba, ijy_ba = get_hist_roi(evol_ba, tresH, 40)
end

# ╔═╡ afb4dca5-ad8d-481d-8b51-5ba9b27b78dd
begin 
	p1b = histogram(n_steps_ba, label="Number of steps")
	p2b = histogram(hh_ba, label="Step height")
	p3b = histogram(hw_ba, label="Step width")
	plot(p1b, p2b, p3b, layout=(1,3))
end

# ╔═╡ 47b632b2-a17a-4dcf-a68c-edb4a28f1390
md""" Select pixel: $(@bind ijview_ba NumberField(1:size(ijy_ba)[1], default=1))"""

# ╔═╡ 1fb4b44d-2006-4bb2-a97c-4b42dadde372
begin
	pxib, pxjb = ijy_ba[ijview_ba, :]
	asf.LaserLab.plot_pixel(evol_ba, pxib, pxjb, tresH)
end

# ╔═╡ 039f7842-aa3c-4eee-917a-c4d7dd232781
begin
	pc1 = histogram(hh, label="Heights free")
	histogram!(pc1, hh_ba, label = "Heights Barium")
	pc2 = histogram(hw, label="Widths free")
	histogram!(pc2, hw_ba, label="Widths Barium")
	pc3 = histogram(n_steps, label="N free")
	histogram!(pc3, n_steps_ba, label="N Ba")
	
	plot(pc1, pc2, pc3, layout=(1,3))
	
end

# ╔═╡ Cell order:
# ╠═ac949aea-be3f-404d-8a8b-22c25b26ee38
# ╠═980dcf83-1e5c-4559-a42b-b02b0db094b3
# ╠═a03bd74c-dc79-4f88-a142-815a6e94954c
# ╠═ee654cd8-b63a-4e5e-96f8-3104bea4fa60
# ╠═16cbcb03-5dc5-4219-8337-c3024d805338
# ╠═d09a5620-570a-4f16-9e0d-f24fb6ba1fbc
# ╠═933774f1-9a00-4ee0-98e0-980f9e8fb930
# ╠═c22c5abe-980b-11ee-268c-b16eb65f35f4
# ╠═99ea72eb-b944-4d2a-9347-f4d9255ac790
# ╠═d3da2ac5-df51-4d7e-99f5-039595f6f75e
# ╠═b2b06a68-dcd5-43e1-9649-8574423e30ec
# ╠═74183602-063d-4aa9-884e-aff77e9606d9
# ╠═f94daecf-afda-4b61-a5d1-3ae48c035480
# ╠═12b432f4-5942-4958-9de8-473a08c364a5
# ╠═10b8263d-8dc5-4ddc-bd21-b310ac7fe11b
# ╠═8a6bcad1-83f1-44b0-afa0-86f7e09e63c6
# ╠═bd225196-cea7-490d-b8bf-d515f51a31a3
# ╠═4cbbadc4-9948-4a21-ab03-18ed7c73573d
# ╠═ccab705b-5193-4618-89f1-ced300100ba3
# ╠═78f8ff68-ec80-4988-a240-064b1e2bc93c
# ╠═680641ef-16b5-4fca-9a57-149bf8c72eb8
# ╠═b1076909-7591-41aa-bada-4747616cec67
# ╠═9a6e62bb-289e-4ca7-9f95-cd0caaac7700
# ╠═39a8c1dc-7a7e-4958-a997-039d2eeb20bd
# ╠═2ba85012-706a-4719-ab87-e0fa91d2f36f
# ╠═5954a7be-70c1-4d23-b4c2-719f23220212
# ╠═4173f292-644b-4a87-b052-c9cfc8bc2765
# ╠═efc103b6-8695-4daf-9647-2fd1baafc8c8
# ╠═9d28bc4a-1bcf-44ce-9429-2bc74e17f9ac
# ╠═7bd7add8-232a-4734-a272-84456af2496e
# ╠═b6a31f37-9600-44f7-80da-a48ad6e9986a
# ╠═bd1320fb-91d1-448c-be04-f60807e96ba0
# ╠═afb4dca5-ad8d-481d-8b51-5ba9b27b78dd
# ╠═47b632b2-a17a-4dcf-a68c-edb4a28f1390
# ╠═1fb4b44d-2006-4bb2-a97c-4b42dadde372
# ╠═54ec44a2-189e-4b23-8afd-f562378f8b7d
# ╠═039f7842-aa3c-4eee-917a-c4d7dd232781
# ╠═7bb25cd5-0f74-48ce-aeee-3020221f6b9b
# ╠═e477b05f-2fbb-4bb8-a37a-a56c286445dc
# ╠═7c6dcca9-fd13-4977-8bf7-7f0f5d2b5636
