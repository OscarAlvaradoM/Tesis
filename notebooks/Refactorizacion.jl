### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ a5c6a6a8-317a-11ec-25bc-7d7e08c53eb4
begin
	using PlutoUI
	PlutoUI.TableOfContents(aside=true, title="🎣 Método del Volumen Finito 🔥")
end

# ╔═╡ c2edc6c3-dbcc-4050-8992-e160bd9ebabb
begin
	using Plots
	plotly()
end

# ╔═╡ 014b9495-588e-4485-a26f-81807624d7b5
md"""
# Método del Volumen Finito
- Óscar A. Alvarado-Morán
- Oscar A. Esquivel-flores
"""


# ╔═╡ 226a6ac4-3f23-454a-9d7c-737cc316e11c
md"## Mallado"

# ╔═╡ 6eb2ab18-da3f-4433-85f1-ead20aa2fb80
begin
	function set_deltas(dominio::Array)
		(dominio[2:end] - dominio[1:end-1])[2:end-1]
	end
	
	function grid(I::Array, J::Array, K::Array)
		grid = [[i, j, k] for i ∈ I, j ∈ J, k ∈ K]
	end
	
	function get_deltas_faces(faces::Array)
		deltas_faces = [faces[direction][2:end] - faces[direction][1:(end-1)] for direction ∈ 1:3]
		return deltas_faces
	end

	function get_grids(deltas::Array, centers::Array, faces::Array)
		#--------------- Aquí falta ver lo del meshgrid------------------
		grid_deltas = grid(deltas[1], deltas[2], deltas[3])
		grid_coords = grid(centers[1], centers[2], centers[3])
		grid_faces = grid(faces[1], faces[2], faces[3])
		return grid_deltas, grid_coords, grid_faces
	end
end	

# ╔═╡ ca72c120-e0b0-4bdc-8641-1d65aa1fdd9f
# set_volumes_and_lengths(5,5)

# ╔═╡ 7b0d3a85-1eab-4726-a6f1-7b2a782a804d
#set_volumes_and_lengths([5,4],[5,4])

# ╔═╡ b59b0526-5e2a-478f-a16d-4de6572bab5e
function uniform_grid(volume::Real, len::Real)
	volumes = [volume, volume, volume]
	lengths = [len, len, len]
	# Separación entre todos los nodos de cada dimensión
	volume_length = len/volume
	# La frontera "inicial" del arreglo
	start = volume_length/2
	# La frontera "final" del arreglo
	stop = len - start
	
	centers = [collect(start:volume_length:stop) for _ ∈ 1:3]
	
	centers_and_boundaries = [copy(arr) for arr ∈ centers]
 	centers_and_boundaries = [insert!(arr,1,0) for arr ∈ centers_and_boundaries]
 	centers_and_boundaries = [insert!(arr,length(arr)+1,len) for arr ∈ centers_and_boundaries]
	
 	deltas = [length(set_deltas(dom)) != 0 ? set_deltas(dom) : [dom[end]] for dom ∈ centers_and_boundaries]
	
 	faces = [vcat(0,  (coords[1:end-1] + coords[2:end])/2 ,  len) for coords ∈ centers]

 	deltas_faces = get_deltas_faces(faces)
 	
	grid_deltas, grid_coords, grid_faces = get_grids(deltas, centers, faces)
	
 	return volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces
end

# ╔═╡ 8f8169ae-214f-4d2a-8fef-0ba170dcefdf
#uniform_grid(5, 5)

# ╔═╡ 314ec1b1-d930-4f95-b31d-474bd7c8a357
function init_tags(dim::Int, volumes::Array, centers_and_boundaries::Array)
	tags = Dict()
	X, Y, Z = volumes
	for z ∈ 2:Z+1, y ∈ 2:Y+1, x ∈ 2:X+1                   
		t = b = n = s = false
		e = w = :F
		if x == 2 w = Dict()
		elseif x == X+1 e = Dict()
		end
		if dim > 1
			n = s = :F
			if y == 2 n = Dict()
			elseif y == Y+1 s = Dict
			end
			if dim == 3
				t = b = :F
				if z == 2 t = Dict()
				elseif z == Z+1 b = Dict()
				end
			end
		end
		tags["$(x)-$(y)-$(z)"] = Dict(:E => e, :W => w, :N => n, :S => s,
			:T => t, :B => b, "coord" => [centers_and_boundaries[1][x],
				centers_and_boundaries[2][y],
				centers_and_boundaries[3][z]])
	end
	return tags
end

# ╔═╡ f701a90a-1814-4398-a910-7940efdec556
#init_tags(3, unif_grid.volumes, unif_grid.centers_and_boundaries)

# ╔═╡ eb9f3443-0032-41ca-8a55-9a38d8a2e3a9
function init_tags_boundaries(dim::Int, centers_and_boundaries::Array)
	tags_boundaries = Dict()
	X, Y, Z = [length(dom) for dom ∈ centers_and_boundaries]
	t = b = n = s = false
	e = w = true
	if dim > 1
		n = s = true
		if dim == 3 
			t = b = true 
		end
	end
	for z ∈ 1:Z, y ∈ 1:Y, x ∈ 1:X
		# El siguiente cacho de código es para saber si nos encontramos con una frontera
		if x==1 || y==1 || z==1 || x==X || y==Y || z==Z
			var = ""
			if y != 1 && y != Y
				if z != 1 && z != Z
					if x == 1 
						var = :W
						value = w
					elseif x == X 
						var = :E
						value = e
					else continue end
				elseif x != 1 && x != X
					if z == 1 
						var = :T
						value = t
					elseif z == Z 
						var = :B
						value = b
					else continue end
				else continue end
			elseif z != 1 && z != Z
				if x != 1 && x != X
					if y == 1 
						var = :N
						value = n
					elseif y == Y
						var = :S
						value = s
					end
				else continue end
			else continue end
			tags_boundaries["$(x)-$(y)-$(z)"] = Dict("frontera" => Dict(var => value), "coord" => [centers_and_boundaries[1][x],centers_and_boundaries[2][y], centers_and_boundaries[3][z]], "cond" => Dict())
		end
	end
	return tags_boundaries
end

# ╔═╡ fb79c628-780e-40ed-bc3f-9670264bfa3b
#init_tags_boundaries(3, unif_grid.centers_and_boundaries)

# ╔═╡ 4d43abf5-2b9a-4f4e-a060-d2cc7bbf72a9
begin
	function tag_wall(tags::Dict, tags_boundaries::Dict, directions::Array, values::Array, coords::Array, cond_type::Symbol=:D)
		for (idx, key) ∈ enumerate(coords)
			# Primero checamos si NO está en la frontera
			if key ∈ keys(tags)
				tags[key][directions[idx]][cond_type] = values[idx]
			# Si no está fuera de la frontera,suponemos que está en ella
			elseif key ∈ keys(tags_boundaries)
				tags_boundaries[key]["cond"][cond_type] = values[idx]
			end
		end
		return tags, tags_boundaries
	end
	
	function tag_wall(tags::Dict, tags_boundaries::Dict, directions::Array, value::Real, coords::Array, cond_type::Symbol=:D)
	
		for (idx, key) ∈ enumerate(coords)
			# Primero checamos si NO está en la frontera
			if key ∈ keys(tags)
				tags[key][directions[idx]][cond_type] = value
			# Si no está fuera de la frontera,suponemos que está en ella
			elseif key ∈ keys(tags_boundaries)
				tags_boundaries[key]["cond"][cond_type] = value
			end
		end
		return tags, tags_boundaries
	end
	
	function tag_wall(tags::Dict, tags_boundaries::Dict, directions::Array, values::Array, coords::Array, cond_type::Symbol=:D)
	
		for (idx, key) ∈ enumerate(coords)
			# Primero checamos si NO está en la frontera
			if key ∈ keys(tags)
				tags[key][directions[idx]][cond_type] = values[idx]
			# Si no está fuera de la frontera,suponemos que está en ella
			elseif key ∈ keys(tags_boundaries)
				tags_boundaries[key]["cond"][cond_type] = values[idx]
			end
		end
		return tags, tags_boundaries
	end
	
	function tag_wall(tags::Dict, tags_boundaries::Dict, directions::Array, values::Array, cond_type::Symbol=:D)
		for (idx, direction) ∈ enumerate(directions)
			tag_wall(cond_type, tags, tags_boundaries, direction, values[idx])
		end
		return tags, tags_boundaries
	end
	
	function tag_wall(tags::Dict, tags_boundaries::Dict, directions::Array, value::Real, cond_type::Symbol=:D)
		for (idx, direction) ∈ enumerate(directions)
			tag_wall(tags, tags_boundaries, direction, value, cond_type)
		end
		return tags, tags_boundaries
	end
	
	function tag_wall(tags::Dict, tags_boundaries::Dict, direction::Symbol, value::Real, cond_type::Symbol=:D)
		for key ∈ keys(tags)
			if typeof(tags[key][direction]) == Dict
				tags[key][direction][cond_type] = value
			end
		end
		for key ∈ keys(tags_boundaries)
			if get(tags_boundaries[key]["frontera"], direction, false)
				tags_boundaries[key]["cond"][cond_type] = value
			end
		end
		return tags, tags_boundaries
	end
end

# ╔═╡ 43407e34-b7cb-4a4e-9892-59b39ee38af7
#tag_wall(tags, tags_b, :T, 1670, :D)

# ╔═╡ c58ea7b0-daeb-4340-ba1f-cb297ab0ac77
begin
	mutable struct Mesh
		volumes
		lengths
		centers::Array
		centers_and_boundaries::Array
		deltas::Array
		faces::Array
		deltas_faces::Array
		tags
		tags_boundaries
	end 
end

# ╔═╡ 2b833859-dd6a-4dc6-9640-b206cd254413
function set_volumes_and_lengths(volume::Real, len::Real) :: Mesh
	volumes = [volume, 1, 1]
	lengths = [len, len/10, len/10]
	centers = [[(i+0.5)*len/volume for i ∈ 0:volume-1], [len/20], [len/20]]
 	volume_lengths = lengths/volumes

	centers_and_boundaries = [copy(arr) for arr ∈ centers]
 	centers_and_boundaries = [insert!(arr,1,0) for arr ∈ centers_and_boundaries]
 	centers_and_boundaries = [insert!(arr,length(arr)+1,lengths[idx]) for (idx,arr) ∈ enumerate(centers_and_boundaries)]

 	deltas = [length(set_deltas(dom)) != 0 ? set_deltas(dom) : [dom[end]] for dom ∈ centers_and_boundaries]
	
	faces = [vcat(0,  (coords[1:(end-1)] + coords[2:end])/2 ,  lengths[idx]) for (idx,coords) ∈ enumerate(centers)]

 	deltas_faces = get_deltas_faces(faces)
 	grid_deltas, grid_coords, grid_faces = get_grids(deltas, centers, faces)
	
	return volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces
end

# ╔═╡ ed1a1dbd-6aab-48b3-9a06-354cd6adc77e
function set_volumes_and_lengths(volumes::Array, lengths::Array) :: Mesh
	missing_volumes = 3 - length(volumes)
	volumes = vcat(volumes, [1 for _ ∈ 1:missing_volumes])

	missing_lengths = 3 - length(lengths)
	lengths = vcat(lengths, [lengths[1]/10 for _ in 1:missing_lengths])
	centers = [[(j+0.5)*lengths[i]/volumes[i] for j in 0:volumes[i]-1] for i in 1:3]
	
	centers_and_boundaries = [copy(arr) for arr ∈ centers]
	centers_and_boundaries = [insert!(arr,1,0) for arr ∈ centers_and_boundaries]
	centers_and_boundaries = [insert!(arr,length(arr)+1,lengths[idx]) for (idx,arr) ∈ enumerate(centers_and_boundaries)]

	deltas = [length(set_deltas(dom)) != 0 ? set_deltas(dom) : [dom[end]] for dom ∈ centers_and_boundaries]

	faces = [vcat(0,  (coords[1:(end-1)] + coords[2:end])/2 ,  lengths[idx]) for (idx,coords) ∈ enumerate(centers)]

	deltas_faces = get_deltas_faces(faces)
	grid_deltas, grid_coords, grid_faces = get_grids(deltas, centers, faces)
	
	return volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces
end

# ╔═╡ 2bff23a2-ffe8-4288-b473-74ba8d023763
function draw(mesh)
	"""
	Método para graficar la malla. Este método se invoca hasta que se hayan inizializado todas las condiciones
	de frontera.
	"""
	# Graficamos las fronteras, sean o no activas
	dic_colors = Dict(:D => "darkturquoise", :N => "red", :S => "magenta",
		:Off => "white", :I => "gray")
	condiciones = [collect(values(mesh.tags_boundaries[key]["frontera"]))[1] == true ? collect(keys(mesh.tags_boundaries[key]["cond"]))[1] : false for key in keys(mesh.tags_boundaries)]
	
 	colores = [dic_colors[cond] for cond in condiciones]
	
 	# Obtenemos las coordenadas de las fronteras y de los nodos internos.
	# Aquí se pondrán las coordenadas de las fronteras
 	coordenadas = [] 
	# Aquí las coordendas de los nodos internos
 	coord = [] 
 	for i in 1:3
 		push!(coordenadas, [mesh.tags_boundaries[key]["coord"][i] for key in keys(mesh.tags_boundaries)])
 		push!(coord, [mesh.tags[key]["coord"][i] for key in keys(mesh.tags)])
	end
	
	Plots.scatter(coordenadas[1], coordenadas[2], coordenadas[3], markershape=:rect, markersize = 1, alpha = 0.5, color = colores)
 	Plots.scatter!(coord[1], coord[2], coord[3], color = :blue, markersize = 3)
end

# ╔═╡ f5e04936-f01e-487e-831b-b518ddedf103
md"---"

# ╔═╡ 3b71ed0f-53fb-4382-bcfb-786dd85f0ad1
begin
	function get_grid_deltas_centers_and_boundaries(mesh::Mesh, axis::Int=1, orientation::Symbol=:E, reverse::Bool=false)
		centers_and_boundaries = mesh.centers_and_boundaries
		centers = mesh.centers
		deltas_centers_and_boundaries = []
		grid_deltas_centers_and_boundaries = []
		for direction ∈ 1:3
			if direction != axis
				if length(centers_and_boundaries[direction]) == 3
					push!(deltas_centers_and_boundaries, [centers_and_boundaries[direction][2]])
				else
					push!(deltas_centers_and_boundaries, centers[direction])
				end
			else
				push!(deltas_centers_and_boundaries, 
					centers_and_boundaries[direction][2:end] - centers_and_boundaries[direction][1:end-1])
			end
		end
		if orientation ∈ [:E, :S,:B]
			if reverse
				deltas_centers_and_boundaries[axis] = /
				deltas_centers_and_boundaries[axis][1:end-1]
			else
				deltas_centers_and_boundaries[axis] = deltas_centers_and_boundaries[axis][2:end]
			end
		else
			if reverse
				deltas_centers_and_boundaries[axis] = deltas_centers_and_boundaries[axis][2:end]
			else
				deltas_centers_and_boundaries[axis] = deltas_centers_and_boundaries[axis][1:end-1]
			end
		end	
		
		grid_deltas_centers_and_boundaries = grid(deltas_centers_and_boundaries[1],
			deltas_centers_and_boundaries[2], deltas_centers_and_boundaries[3])
	end
	
	function get_deltas_faces_by_axis(mesh::Mesh, axis::Int=1)
		faces = mesh.faces
		
		deltas_faces =  faces[axis][2:end] - faces[axis][1:end-1]
	end
	
	
	function get_array_before_area(mesh::Mesh, direction::Real, extended::Bool=false)
		perpendicular = [i for i ∈ 1:3 if i != direction]
		num_boundaries = mesh.volumes[direction]

		array = [idx ∈ perpendicular ? mesh.deltas_faces[idx] : ones(num_boundaries) for idx ∈ 1:3]
		return array
	end
end

# ╔═╡ 23d190c8-4c51-4f05-b8a4-0573f10d8c11
# begin
# 	grid_centers_and_boundaries = get_grid_deltas_centers_and_boundaries(mesh)
# 	grid_centers_and_boundaries[1,1,1]
# end

# ╔═╡ 912c556a-4ffd-4e49-a35e-ac83c77c112a
# begin
# 	grid_deltas_faces = get_deltas_faces_by_axis(mesh)
# 	grid_deltas_faces
# end

# ╔═╡ b5ed1976-b09a-4a60-83ce-2455d5d63755
md"## Coeficientes"

# ╔═╡ 299fb45f-0579-42a4-8e28-46f8b818ecc9
begin
	mutable struct Coefficients_1d
		mesh::Mesh
		aP::Array{Float64, 3}
		aW::Array{Float64, 3}
		aE::Array{Float64, 3}
		Su::Array{Float64, 3}
		Sp_a::Array{Float64, 3}
		Sp_d::Array{Float64, 3}
		bW::Array{Float64, 3}
		bE::Array{Float64, 3}
		aWW::Array{Float64, 3}
		aEE::Array{Float64, 3}
	end
	
	mutable struct Coefficients_2d
		mesh::Mesh
		aP::Array{Float64, 3}
		aW::Array{Float64, 3}
		aE::Array{Float64, 3}
		Su::Array{Float64, 3}
		Sp_a::Array{Float64, 3}
		Sp_d::Array{Float64, 3}
		bW::Array{Float64, 3}
		bE::Array{Float64, 3}
		aWW::Array{Float64, 3}
		aEE::Array{Float64, 3}
		
		aN::Array{Float64, 3}
		aS::Array{Float64, 3}
		bN::Array{Float64, 3}
		bS::Array{Float64, 3}
		aNN::Array{Float64, 3}
		aSS::Array{Float64, 3}
	end 
	
	mutable struct Coefficients_3d
		mesh::Mesh
		aP::Array{Float64, 3}
		aW::Array{Float64, 3}
		aE::Array{Float64, 3}
		Su::Array{Float64, 3}
		Sp_a::Array{Float64, 3}
		Sp_d::Array{Float64, 3}
		bW::Array{Float64, 3}
		bE::Array{Float64, 3}
		aWW::Array{Float64, 3}
		aEE::Array{Float64 ,3}
		
		aN::Array{Float64, 3}
		aS::Array{Float64, 3}
		bN::Array{Float64, 3}
		bS::Array{Float64, 3}
		aNN::Array{Float64, 3}
		aSS::Array{Float64, 3}
		
		aT::Array{Float64, 3}
		aB::Array{Float64, 3}
		bT::Array{Float64, 3}
		bB::Array{Float64, 3}
		aTT::Array{Float64, 3}
		aBB::Array{Float64, 3}
	end 
end

# ╔═╡ cadfb39c-3487-460b-972e-2af7d348bed9
function init_coefficients(mesh::Mesh)
	dim = sum(mesh.volumes .!= 1)
	vols = mesh.volumes

	aP = zeros(vols[1], vols[2], vols[3])
	aW = zeros(vols[1], vols[2], vols[3])
	aE = zeros(vols[1], vols[2], vols[3])
	Su = zeros(vols[1], vols[2], vols[3])
	Sp_a = zeros(vols[1], vols[2], vols[3])
	Sp_d = zeros(vols[1], vols[2], vols[3])
	bW = zeros(vols[1], vols[2], vols[3]) # Aquí la contribuciones de la serie de Taylor para la condición de frontera
	bE = zeros(vols[1], vols[2], vols[3])
	aWW = zeros(vols[1], vols[2], vols[3]) # Aquí depositaremos lo de la advección a segundo orden
	aEE = zeros(vols[1], vols[2], vols[3])

	if dim > 1
		aN = zeros(vols[1], vols[2], vols[3])
		aS = zeros(vols[1], vols[2], vols[3])
		bN = zeros(vols[1], vols[2], vols[3])
		bS = zeros(vols[1], vols[2], vols[3])
		aNN = zeros(vols[1], vols[2], vols[3])
		aSS = zeros(vols[1], vols[2], vols[3])
	end
	if dim == 3
		aT = zeros(vols[1], vols[2], vols[3])
		aB = zeros(vols[1], vols[2], vols[3])
		bT = zeros(vols[1], vols[2], vols[3])
		bB = zeros(vols[1], vols[2], vols[3])
		aTT = zeros(vols[1], vols[2], vols[3])
		aBB = zeros(vols[1], vols[2], vols[3])
	end
	
	if dim == 1
		coeff = Coefficients_1d(mesh, aP, aW, aE, Su, Sp_a, Sp_d, bW, bE, aWW, aEE)
	elseif dim == 2
		coeff = Coefficients_2d(mesh, aP, aW, aE, Su, Sp_a, Sp_d, bW, bE, aWW, aEE,
		aN, aS, bN, bS, aNN, aSS)
	else
		coeff = Coefficients_3d(mesh, aP, aW, aE, Su, Sp_a, Sp_d, bW, bE, aWW, aEE,
		aN, aS, bN, bS, aNN, aSS, aT, aB, bT, bB, aTT, aBB)
	end
	return coeff
end

# ╔═╡ 3d74fabd-a66c-46f8-8ae3-ecac28587efb
function add_source(coeff, source::Real)
	"""

	"""
	#--------------- Aquí falta ver lo del meshgrid------------------
	x, y, z = grid(coeff.mesh.deltas_faces[1], coeff.mesh.deltas_faces[2],
						  coeff.mesh.deltas_faces[3])
	# Si es una constante
	coeff.Su += source*(x*y*z)
end

# ╔═╡ 1415d40a-91bb-4c26-9d31-23ad698285dd
md"## Difusión
El coefficiente Γ DEBE ser una función, para que con esta construyamos un meshgrid FORZOSO sobre el cual vamos a estar metiendo todo, sobre este podemos meterle el área con multiplicaciones tipo dx\*dy\*dz, el delta y todo lo demás, sólo tenemos que hacer UNITO grid"

# ╔═╡ 4b30e966-8b53-4fa2-bf83-3f904118acc1
function Γ_constant(x::Array, y::Array, z::Array)
	Γ = 1000
	tensor_Γ = ones(length(x), length(y), length(z))
	tensor_Γ = Γ .* tensor_Γ
end

# ╔═╡ a8b43883-37de-4f20-aa2d-97816d13344b
begin
	mutable struct Diffusion
		coeff
		Γ
	end 
end

# ╔═╡ 4b58c8c2-82d7-4fdf-ab9a-6ed248022ed3
begin
	function grid_mult_area(gridmesh, list_3d)
		for (idx_z, dz) ∈ enumerate(list_3d[3])
			for (idx_y, dy) ∈ enumerate(list_3d[2])
				for (idx_x, dx) ∈ enumerate(list_3d[1])
					gridmesh[idx_x, idx_y, idx_z] *= dx*dy*dz
				end
			end
		end
		return gridmesh
	end
	
	function grid_div_delta(gridmesh::Array, list::Array, axis::Int)
		gridmesh_copy = ones(size(gridmesh))
		for idx_z ∈ 1:size(gridmesh)[3]
			for idx_y ∈ 1:size(gridmesh)[2]
				for idx_x ∈ 1:size(gridmesh)[1]
					index_list = [idx_x, idx_y, idx_z] 
					gridmesh[idx_x, idx_y, idx_z] = gridmesh[idx_x, idx_y, idx_z] / list[index_list[axis]]
				end
			end
		end
		return gridmesh
	end
	
	function grid_mult_delta(gridmesh, list::Array, axis::Int)
		for idx_z ∈ 1:size(gridmesh)[3]
			for idx_y ∈ 1:size(gridmesh)[2]
				for idx_x ∈ 1:size(gridmesh)[1]
					index_list = [idx_x, idx_y, idx_z]
					gridmesh[idx_x, idx_y, idx_z] *= list[index_list[axis]]
				end
			end
		end
		return gridmesh
	end
end

# ╔═╡ 0f0e2768-83bc-45f3-9507-cd8d25209bab
begin	
	function get_basic_diffusion_coef(diff::Diffusion, direction::Symbol,
			what_dimension::Int, limit::Int)
		coeff = diff.coeff
		mesh = coeff.mesh
		faces = mesh.faces[what_dimension]

		δ_by_axis = get_deltas_faces_by_axis(mesh, what_dimension)
		array_before_areas = get_array_before_area(mesh, what_dimension)
		Γ = get_Γ_grid(diff, faces, direction, what_dimension)

		# Para obtener los coeficientes de las caras centrales
		diffusion_term = get_interior_terms(Γ, direction, what_dimension,
			array_before_areas, δ_by_axis)
		
		condition_mask, condition_mask_type = get_mask_boundaries(coeff, direction)
		
		
		# Aquí obtenemos Sp (que son las fuentes internas)
		# Parece que estas siempre van. Es la contribución de la frontera que 
		# no va directamente relacionado con la condición de la frontera. 
		# En la matriz del libro, son los extremos de la diagonal principal.
		array_before_boundary_area = copy(array_before_areas)
		
		bound_term_c, bound_term_f, Γ_bound, Γ_bound_next = get_boundary_terms(diff,
			Γ, direction, condition_mask_type, what_dimension,
			array_before_boundary_area, δ_by_axis)
		
		# Esto porque sólo contribuye a volumen central y no a alguno que lo rodee.
		sp = bound_term_c

		# Aquí obtenemos Su (que son las condition_masks de frontera)
		# En la matriz del libro son los extremos de las diagonales que no 
		# son la principal.
		su = get_boundary_condition_mask(diff, Γ, direction, Γ_bound, Γ_bound_next,
			condition_mask, what_dimension, array_before_areas, δ_by_axis)
		
		return diffusion_term, sp, su, bound_term_f
	end
	
	function get_Γ_grid(diff::Diffusion, faces::Array, direction::Symbol,
			what_dimension::Int)
		mesh = diff.coeff.mesh
		centers = mesh.centers
		# Para obtener el arreglo con las Γ's
		λ₁(x, d) = d ∈ [:E,:S,:B] ? x[2:end] : x[1:end-1]
		coords_align = [i!=what_dimension ? var : λ₁(faces, direction) for (i, var) ∈ enumerate(centers)]
		# Esto sí es un meshgrid 
		Γ = diff.Γ(coords_align[1], coords_align[2], coords_align[3])
		
		return Γ
	end

	function get_interior_terms(Γ::Array, direction::Symbol, what_dimension::Int,
				array_before_areas::Array, δ::Array)
		
		# Aquí se obtiene el promedio de las Γ's en la dirección en la que 
		# se está trabajando
		coord_l = [idx!=what_dimension ? (:) : 1:(size(Γ)[idx]-1) for idx ∈ 1:3]
		coord_r = [idx!=what_dimension ? (:) : 2:size(Γ)[idx] for idx ∈ 1:3]
		
		Γ_left = Γ[coord_l[1], coord_l[2], coord_l[3]]
		Γ_right = Γ[coord_r[1], coord_r[2], coord_r[3]]
		
		Γ_mean = 0.5 .* (Γ_left .+ Γ_right)
		
		
		# Le pones un cero en la dirección en donde le estamos mochando un pedazo
		Γ_bar = zeros(size(Γ))
		λ₂(d, idx) = d ∈ [:E,:S,:B] ? (1:(size(Γ_bar)[idx]-1)) : (2:size(Γ_bar)[idx])
		coord₃ = [idx!=what_dimension ? (:) : λ₂(direction, idx)  for idx ∈ 1:3]
		
		Γ_bar[coord₃[1], coord₃[2], coord₃[3]] = Γ_mean
		Γ_area = grid_mult_area(Γ_bar, array_before_areas)
		diff_coef = grid_div_delta(Γ_area, δ, what_dimension)
		
		return diff_coef
	end
	
	# Aquí creamos una máscara que nos dirá qué condición de frontera hay por cada coordenada.
	#### Hay que cambiar dependiendo de qué pase con cada condición de frontera.
	function get_mask_boundaries(coeff, direction::Symbol)
		mesh = coeff.mesh
		tags_boundaries = mesh.tags_boundaries
		condition_mask = []
		condition_mask_type = []
		dict_cond = Dict(:I=>0, :N=>0, :D=>1)
		for tag_b ∈ keys(tags_boundaries)
			if collect(keys(tags_boundaries[tag_b]["frontera"]))[1] == direction
				cond = collect(keys(tags_boundaries[tag_b]["cond"]))[1]
				push!(condition_mask, dict_cond[cond])
				push!(condition_mask_type, cond)
			end
		end
		return condition_mask, condition_mask_type
	end

	function get_boundary_terms(diff, Γ::Array, direction::Symbol,
			condition_mask_type, what_dimension::Int,
			array_before_boundary_area::Array, δ::Array)
		#mesh = diff.coeff.mesh
		ba = [copy(array) for array ∈ array_before_boundary_area]

		λ_area(d, idx) = d ∈ [:E,:S,:B] ? (length(ba[idx]):length(ba[idx])) : (1:1)
		coord_area = [idx!=what_dimension ? (:) : λ_area(direction, idx) for idx ∈ 1:3]
		#[ba[idx][coord] .= 0 for (idx, coord) ∈ enumerate(coord_area)]
		ba = [ba[idx][coord_area[idx]] for idx in 1:3]

		λ(d, idx) = d in [:E,:S,:B] ? (length(Γ[idx]):length(Γ[idx])) : (1:1)
		coord = [idx!=what_dimension ? (:) : λ(direction, idx) for idx ∈ 1:3]
		Γ_bound = Γ[coord[1], coord[2], coord[3]]
		
		λ_1(d, idx) = d ∈ [:E,:S,:B] ? ((size(Γ)[idx]-1):(size(Γ)[idx]-1)) : 2:2
		coord_1 = [idx!=what_dimension ? (:) : λ_1(direction, idx) for idx in 1:3]
		Γ_bound_next = Γ[coord_1[1], coord_1[2], coord_1[3]]

		coef_f, coef_c = get_boundary_coefs(Γ_bound, Γ_bound_next,
			condition_mask_type)

		# Para esta parte habría que agregar las areas correspondientes a las diferentes caras
		area_coef_f = grid_mult_area(coef_f, ba)
		area_coef_c = grid_mult_area(coef_c, ba)
		bound_term_f = grid_div_delta(area_coef_f, δ, what_dimension)
		bound_term_c = grid_div_delta(area_coef_c, δ, what_dimension)

		bound_term_c = to_tensor_number_volumes(bound_term_c, diff, direction,
			what_dimension)
		bound_term_f = to_tensor_number_volumes(bound_term_f, diff, direction,
			what_dimension)
		return bound_term_c, bound_term_f, Γ_bound, Γ_bound_next
	end

	function get_boundary_coefs(Γ_bound, Γ_bound_next, condition_mask_type)
		initial_shape = size(Γ_bound)
		Γ_bound_flatten = vcat(Γ_bound...)
		Γ_bound_next_flatten = vcat(Γ_bound_next...)

		coefs_f = []
		coefs_c = []
		for (idx, condition) ∈ enumerate(condition_mask_type)
			if condition == :D
				coef_c = (9Γ_bound_flatten[idx] - 3Γ_bound_next_flatten[idx])/2
				push!(coefs_c, coef_c)
				coef_f = (Γ_bound_flatten[idx] - Γ_bound_next_flatten[idx]/3)/2
				push!(coefs_f, coef_f)
			elseif condition == :N
				coef_f = 0.
				push!(coefs_f, 0.)
				push!(coefs_c, 0.)
			else
				# Aquí faltaría definir la condición de Robín con un ejemplo, 
				# pero en teoría todo está en el libro
				push!(coefs_f, 0.)
				push!(coefs_c, 0.)
			end
		end
		coefs_f = reshape(coefs_f, initial_shape)
		coefs_c = reshape(coefs_c, initial_shape)
		return coefs_f, coefs_c
	end
	
	function reshape_boundary(array::Array, diff::Diffusion, what_dimension::Int)
		mesh = diff.coeff.mesh
		volumes = mesh.volumes
		
		dims_to_reshape = [idx!=what_dimension ? volumes[idx] : 1 for idx in 1:3]
		array = reshape(array, tuple(dims_to_reshape...))
	end
	
 	function get_boundary_condition_mask(diff::Diffusion, Γ::Array, direction::Symbol,
 			Γ_bound::Array, Γ_bound_next::Array, condition_mask::Array,
 			what_dimension::Int, array_before_areas::Array, δ::Array)
		
 		Γ_cond = (3Γ_bound - Γ_bound_next)/2
 		conds_su = get_mask_boundaries_Su(diff, Γ, direction)
 		conds_su = reshape_boundary(conds_su, diff, what_dimension)
 		ba = copy(array_before_areas)
 		condition_mask = reshape_boundary(condition_mask, diff, what_dimension)

 		λ_area(d, idx) = d ∈ [:E,:S,:B] ? (length(ba[idx]):length(ba[idx])) : (1:1)
 		coord_area = [idx!=what_dimension ? (:) : λ_area(direction,idx) for idx ∈ 1:3]
 		ba = [ba[idx][coord_area[idx]] for idx ∈ 1:3]
 		su = ba

  		λ(d, idx) = d ∈ [:E,:S,:B] ? (length(su[idx])) : 1
  		coord = [idx!=what_dimension ? (:) : λ(direction, idx) for idx ∈ 1:3]
  		tmp_su = [su[idx][coord[idx]] for idx ∈ 1:3]
		
  		tmp_su2 = grid_mult_area(conds_su, tmp_su)
  		div = grid_mult_delta(condition_mask, δ[coord[1],coord[2],coord[3]],
			what_dimension)
  		div[findall(==(0), div)] .= 1
 		su = tmp_su2 ./ div
		
		su = to_tensor_number_volumes(su, diff, direction, what_dimension)
  		return su
 	end


	function get_mask_boundaries_Su(diff::Diffusion, Γ::Array, direction::Symbol)
		mesh = diff.coeff.mesh
		tags_boundaries = mesh.tags_boundaries
		condition_mask = []
		counter = 1
		Γ_vcat = vcat(Γ...)
		for (idx, tag_b) ∈ enumerate(keys(tags_boundaries))
			if collect(keys(tags_boundaries[tag_b]["frontera"]))[1] == direction
				cond = collect(keys(tags_boundaries[tag_b]["cond"]))[1]
				if cond == :I
					push!(condition_mask, 0)
				elseif cond == :N
					push!(condition_mask, tags_boundaries[tag_b]["cond"][cond])
				elseif cond == :D
					gamma_idx = Γ_vcat[counter]
					tag_b_idx = 8tags_boundaries[tag_b]["cond"][cond]/3
					push!(condition_mask, tag_b_idx*gamma_idx)
				counter += 1
				end
			end
		end
		return condition_mask
	end
	
	
	function to_tensor_number_volumes(array::Array, diff::Diffusion,
			direction::Symbol, what_dimension::Int)
		mesh = diff.coeff.mesh
		volumes = mesh.volumes
		
		zero = zeros(volumes[1], volumes[2], volumes[3])
		f(i, d) = d ∈ [:E,:S,:B] ? volumes[i] : 1
		coord = [idx!=what_dimension ? (:) : f(idx, direction) for idx ∈ 1:3]
		
		zero[coord[1], coord[2], coord[3]] = array
		
		return zero
	end
end
	

# ╔═╡ 2be171c0-009c-4006-a5db-e3f3cd0ef490
function get_diffusion_coef(diff::Diffusion, direction::Symbol=:E)
		"""

		"""
		what_dimension, limit = 1,1
		if direction ∈ [:E, :S, :B] limit = -1 end
		if direction ∈ [:N, :S] what_dimension = 2
		elseif direction ∈ [:T, :B] what_dimension = 3
		end
		get_basic_diffusion_coef(diff, direction, what_dimension, limit)
	end


# ╔═╡ e52f6cb9-a17e-49e7-9b93-9855bccf1755
function set_diffusion(coeff, Γ)
	"""

	"""
	mesh = coeff.mesh
	dim = sum(mesh.volumes .!= 1)
	diffusion = Diffusion(coeff, Γ)

	west_diff, sp_w, su_w, bound_term_w = get_diffusion_coef(diffusion, :W)
	east_diff, sp_e, su_e, bound_term_e = get_diffusion_coef(diffusion, :E)
	
	coeff.aW -= west_diff
	coeff.aE -= east_diff
	coeff.bW -= bound_term_w
	coeff.bE -= bound_term_e
	coeff.Sp_d -= sp_e + sp_w
	coeff.Su += su_e + su_w
	coeff.aP -= -east_diff - west_diff
	# Las siguientes dos están al revés porque la frontera sólo afecta al coeficiente anterior 
	# por la forma en que hacemos la proximación con la serie de Taylor. 
	coeff.aE += coeff.bW 
	coeff.aW += coeff.bE
	
	@show coeff.aW
	@show coeff.aE
	@show coeff.Sp_d
	@show coeff.Su
	@show coeff.aP
	
	if dim > 1
		north_diff, sp_n, su_n, bound_term_n = get_diffusion_coef(diffusion, :N)
		south_diff, sp_s, su_s, bound_term_s = get_diffusion_coef(diffusion,  :S)

		coeff.aN -= north_diff
		coeff.aS -= south_diff
		coeff.bN -= bound_term_n
		coeff.bS -= bound_term_s
		coeff.Sp_d -= sp_n + sp_s
		coeff.Su += su_n + su_s
		coeff.aP -= -north_diff - south_diff
		coeff.aS += coeff.bN
		coeff.aN += coeff.bS
	end

	if dim == 3
		top_diff, sp_t, su_t, bound_term_t = get_diffusion_coef(diffusion, :T)
		bottom_diff, sp_b, su_b, bound_term_b = get_diffusion_coef(diffusion, :B)

		coeff.aT -= top_diff
		coeff.aB -= bottom_diff
		coeff.bT -= bound_term_t
		coeff.bB -= bound_term_b
		coeff.Sp_d -= sp_t + sp_b
		coeff.Su += su_t + su_b
		coeff.aP -= -top_diff - bottom_diff
		coeff.aB += coeff.bT
		coeff.aT += coeff.bB
	coeff.aP -= coeff.Sp_d
	end
end

# ╔═╡ c0728e21-240f-4f38-98b7-caa943842afb
begin
	volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces = uniform_grid(5, 0.5)
		
	tags = init_tags(3, volumes, centers_and_boundaries)
	
	tags_b = init_tags_boundaries(3, centers_and_boundaries)
	
	tag_wall(tags, tags_b, [:W, :E, :T, :N, :B], 0, :D)
	tag_wall(tags, tags_b, :S, 100, :D)
	
	mesh = Mesh(volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces, tags, tags_b)
	
	coeff = init_coefficients(mesh)
	
	set_diffusion(coeff, Γ_constant)
	
	coeff
end

# ╔═╡ 3f1cd4d3-9052-4271-9e4e-8965c4e99170
md"## Sistema de ecuaciones"

# ╔═╡ 35fb94be-3647-4156-b82e-b54baa167fba
def __init__(self, coef):
        self.vols = coef.vols
        self.coef = coef
        self.aP = coef.aP
        self.aE, self.aW = coef.aE, coef.aW
        if coef.dim > 1: self.aN, self.aS = coef.aN, coef.aS
        if coef.dim == 3: self.aT, self.aB = coef.aT, coef.aB
        self.Su = coef.Su
        self.N = self.aP.shape[0]*self.aP.shape[1]*self.aP.shape[2]
        self.A, self.b = None, None
        
    def get_diag(self, array, k=0):
        """
        Método para construir una matriz diagonal dado un arreglo 2D y la diagonal en que queremos 
        poner ese arreglo
        """
        if self.N > 50:
            return diags(array, k)
        else:
            return np.diag(array,k)
        
        
    def get_A_matrix(self):
        """
        Método para construir la matriz A, de la ecuación Ax=b, a partir de los coeficientes obtenidos.
        """
        aP = np.ravel(self.aP)
        self.A = self.get_diag(aP)
        if self.coef.dim == 1:
            aW = np.ravel(self.aW)[1:]
            aE = np.ravel(self.aE)[:-1]
            self.A += self.get_diag(aE, k=1) + self.get_diag(aW, k=-1)
            
        elif self.coef.dim == 2:
            aN = np.ravel(self.aN)[1:]
            aS = np.ravel(self.aS)[:-1]
            aW = np.ravel(self.aW)[self.vols[1]:]
            aE = np.ravel(self.aE)[:-self.vols[1]]
            self.A += self.get_diag(aN, k=-1) + self.get_diag(aS, k=1) + \
                        self.get_diag(aE, k=self.vols[1]) + self.get_diag(aW, k=-self.vols[1])
            
        # Esta parte está por verse, hay que encontrar un ejemplo chido en 3D
        elif self.coef.dim == 3:
            aT = np.ravel(self.aT)[1:]
            aB = np.ravel(self.aB)[:-1]
            aN = np.ravel(self.aN)[self.vols[1]:]
            aS = np.ravel(self.aS)[:-self.vols[1]]
            aW = np.ravel(self.aW)[self.vols[1]*self.vols[2]:]
            aE = np.ravel(self.aE)[:-self.vols[1]*self.vols[2]]
            self.A += self.get_diag(aT, k=-1) + self.get_diag(aB, k=1) + \
                        self.get_diag(aN, k=-self.vols[1]) + \
                        self.get_diag(aS, k=self.vols[1]) + \
                        self.get_diag(aW, k=-self.vols[1]*self.vols[2]) +\
                        self.get_diag(aE, k=self.vols[1]*self.vols[2])
        return self.A
    
    
    def get_b_vector(self):
        """
        Método para obtener el vector b, de la ecuación Ax = b, a partir del arreglo Su construido anteriormente
        """
        self.b = np.ravel(self.Su)
        return self.b
    
    
    def solve(self,A,b):
        if len(b)>50:
            return spsolve(A,b)
        else:
            return np.linalg.solve(A,b)
    
    def get_solution(self):
        """
        Método para obtener los valores de la solución al sistema
        """
        A = self.get_A_matrix()
        b = self.get_b_vector()
        return self.solve(A, b)

# ╔═╡ a934ecd4-5d4c-4764-a518-0fb01aeee45d
md"## Ecuación de Poisson 3D"

# ╔═╡ 3e0e018b-a006-45f4-8d7b-e58f944cc972


# ╔═╡ Cell order:
# ╟─014b9495-588e-4485-a26f-81807624d7b5
# ╟─a5c6a6a8-317a-11ec-25bc-7d7e08c53eb4
# ╟─c2edc6c3-dbcc-4050-8992-e160bd9ebabb
# ╟─226a6ac4-3f23-454a-9d7c-737cc316e11c
# ╠═6eb2ab18-da3f-4433-85f1-ead20aa2fb80
# ╠═2b833859-dd6a-4dc6-9640-b206cd254413
# ╠═ca72c120-e0b0-4bdc-8641-1d65aa1fdd9f
# ╠═ed1a1dbd-6aab-48b3-9a06-354cd6adc77e
# ╠═7b0d3a85-1eab-4726-a6f1-7b2a782a804d
# ╠═b59b0526-5e2a-478f-a16d-4de6572bab5e
# ╠═8f8169ae-214f-4d2a-8fef-0ba170dcefdf
# ╠═314ec1b1-d930-4f95-b31d-474bd7c8a357
# ╠═f701a90a-1814-4398-a910-7940efdec556
# ╠═eb9f3443-0032-41ca-8a55-9a38d8a2e3a9
# ╠═fb79c628-780e-40ed-bc3f-9670264bfa3b
# ╠═4d43abf5-2b9a-4f4e-a060-d2cc7bbf72a9
# ╠═43407e34-b7cb-4a4e-9892-59b39ee38af7
# ╠═c58ea7b0-daeb-4340-ba1f-cb297ab0ac77
# ╠═2bff23a2-ffe8-4288-b473-74ba8d023763
# ╟─f5e04936-f01e-487e-831b-b518ddedf103
# ╠═3b71ed0f-53fb-4382-bcfb-786dd85f0ad1
# ╠═23d190c8-4c51-4f05-b8a4-0573f10d8c11
# ╠═912c556a-4ffd-4e49-a35e-ac83c77c112a
# ╟─b5ed1976-b09a-4a60-83ce-2455d5d63755
# ╠═299fb45f-0579-42a4-8e28-46f8b818ecc9
# ╠═cadfb39c-3487-460b-972e-2af7d348bed9
# ╠═e52f6cb9-a17e-49e7-9b93-9855bccf1755
# ╠═3d74fabd-a66c-46f8-8ae3-ecac28587efb
# ╟─1415d40a-91bb-4c26-9d31-23ad698285dd
# ╠═4b30e966-8b53-4fa2-bf83-3f904118acc1
# ╠═a8b43883-37de-4f20-aa2d-97816d13344b
# ╠═2be171c0-009c-4006-a5db-e3f3cd0ef490
# ╠═4b58c8c2-82d7-4fdf-ab9a-6ed248022ed3
# ╠═0f0e2768-83bc-45f3-9507-cd8d25209bab
# ╠═c0728e21-240f-4f38-98b7-caa943842afb
# ╟─3f1cd4d3-9052-4271-9e4e-8965c4e99170
# ╠═35fb94be-3647-4156-b82e-b54baa167fba
# ╟─a934ecd4-5d4c-4764-a518-0fb01aeee45d
# ╠═3e0e018b-a006-45f4-8d7b-e58f944cc972
