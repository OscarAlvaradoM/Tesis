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
		@show x,y,z
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

# ╔═╡ f5e04936-f01e-487e-831b-b518ddedf103
md"---"

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

# ╔═╡ 34bf980d-eebf-4767-907f-df1be41d7da0
function init_mesh(volume::Real, len::Real)
	volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces = uniform_grid(volume, len)
	
	tags = init_tags(3, volumes, centers_and_boundaries)
	
	tags_b = init_tags_boundaries(3, centers_and_boundaries)
	
	tag_wall(tags, tags_b, [:W, :E, :T, :N, :B], 0, :D)
	tag_wall(tags, tags_b, :S, 100, :N)
	
	mesh = Mesh(volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces, tags, tags_b)
end

# ╔═╡ feeda8d2-7aef-46fc-94b3-eb5882c0cfb3
mesh = init_mesh(5, 0.5)

# ╔═╡ 2bff23a2-ffe8-4288-b473-74ba8d023763
function draw(mesh)
	"""
	Método para graficar la malla. Este método se invoca hasta que se hayan inizializado todas las condiciones
	de frontera.
	"""
	# Graficamos las fronteras, sean o no activas
	dic_colors = Dict(:D => "darkturquoise", :N => "red", :S => "magenta", :Off => "white", :I => "gray")
	condiciones = [collect(values(mesh.tags_boundaries[key]["frontera"]))[1] == true ? collect(keys(mesh.tags_boundaries[key]["cond"]))[1] : false for key in keys(mesh.tags_boundaries)]
	
 	colores = [dic_colors[cond] for cond in condiciones]
	
 	# Obtenemos las coordenadas de las fronteras y de los nodos internos.
 	coordenadas = [] # Aquí se pondrán las coordenadas de las fronteras
 	coord = [] # Aquí las coordendas de los nodos internos
 	for i in 1:3
 		push!(coordenadas, [mesh.tags_boundaries[key]["coord"][i] for key in keys(mesh.tags_boundaries)])
 		push!(coord, [mesh.tags[key]["coord"][i] for key in keys(mesh.tags)])
	end
	
	Plots.scatter(coordenadas[1], coordenadas[2], coordenadas[3], markershape=:rect, markersize = 1, alpha = 0.5, color = colores)
 	Plots.scatter!(coord[1], coord[2], coord[3], color = :blue, markersize = 3)
end

# ╔═╡ 93d57deb-f0c7-4df7-bd15-69e95add9d3c
draw(mesh)

# ╔═╡ 8274809e-0d6f-4a94-baac-802f804193ae
length(keys(mesh.tags_boundaries))

# ╔═╡ 3b71ed0f-53fb-4382-bcfb-786dd85f0ad1
begin
	function get_grid_deltas_centers_and_boundaries(centers_and_boundaries::Array, centers::Array, axis::Int=1, orientation::Symbol=:E, reverse::Bool=false)
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
				@show centers_and_boundaries[direction][2:end]
				@show centers_and_boundaries[direction][1:end-1]
				push!(deltas_centers_and_boundaries, 
					centers_and_boundaries[direction][2:end] - centers_and_boundaries[direction][1:end-1])
			end
		end
		@show deltas_centers_and_boundaries
			
		if orientation == :E || orientation == :S || orientation == :B
			if reverse
				deltas_centers_and_boundaries[axis] = /
				deltas_centers_and_boundaries[axis][:end-1]
			else
				deltas_centers_and_boundaries[axis] = deltas_centers_and_boundaries[axis][2:end]
			end
				
		else
			if reverse
				deltas_centers_and_boundaries[axis] = deltas_centers_and_boundaries[axis][2:end]
			else
				deltas_centers_and_boundaries[axis] = deltas_centers_and_boundaries[axis][:end-1]
			end
		end	
		#--------------- Aquí falta ver lo del meshgrid------------------
		grid_deltas_centers_and_boundaries = grid(deltas_centers_and_boundaries[1],
			deltas_centers_and_boundaries[2], deltas_centers_and_boundaries[3])
	end
	
	function get_grid_deltas_faces(centers_and_boundaries::Array, centers::Array, faces::Array, axis::Int=1)
		deltas_faces = []
		for direction ∈ 1:3
			if direction != axis
				if length(centers_and_boundaries[direction]) == 3
					push!(deltas_faces, [centers_and_boundaries[direction][2]])
				else
					push!(deltas_faces, centers[direction])
				end

			else
				push!(deltas_faces, faces[direction][2:end] - faces[direction][1:end-1])
			end
		end
		@show deltas_faces
		#--------------- Aquí falta ver lo del meshgrid------------------
		grid_deltas_faces = grid(deltas_faces[1], deltas_faces[2], deltas_faces[3])
	end
	
	function get_area(mesh::Mesh, direction::Real, extended=false::Boolean)
		perpendicular = [i for i ∈ 1:3 if i != direction]
		num_boundaries = mesh.volumes[direction]

		array = [idx ∈ perpendicular ? mesh.deltas_faces[idx] : ones(num_boundaries) for idx ∈ 1:3]
		if extended
			array = [idx ∈ perpendicular ? mesh.deltas_faces[idx] : ones(num_boundaries+1) for idx in 1:3]
		end
		#--------------- Aquí falta ver lo del meshgrid------------------
		areas_grid = [array[1], array[2], array[3]]
		return areas_grid[perpendicular[1]]*areas_grid[perpendicular[2]]
	end
end

# ╔═╡ 23d190c8-4c51-4f05-b8a4-0573f10d8c11
#get_grid_deltas_centers_and_boundaries(mesh.deltas_centers_and_boundaries, mesh.centers_and_boundaries, mesh.centers)

# ╔═╡ 912c556a-4ffd-4e49-a35e-ac83c77c112a
#get_grid_deltas_faces(mesh.centers_and_boundaries, mesh.centers, mesh.faces)

# ╔═╡ b5ed1976-b09a-4a60-83ce-2455d5d63755
md"## Coeficientes"

# ╔═╡ 299fb45f-0579-42a4-8e28-46f8b818ecc9
begin
	mutable struct Coefficients
		mesh::Mesh
	end 
end

# ╔═╡ cadfb39c-3487-460b-972e-2af7d348bed9
function init_coefficients(coeff::Coefficients)
	dim = coeff.mesh.dim
	mesh = coeff.mesh
	vols = coeff.mesh.volumes

	aP = zeros(vols)
	aW = zeros(vols)
	aE = zeros(vols)
	Su = zeros(vols)
	Sp_a = zeros(vols)
	Sp_d = zeros(vols)
	bW = zeros(vols) # Aquí la contribuciones de la serie de Taylor para la condición de frontera
	bE = zeros(vols)
	aWW = zeros(vols) # Aquí depositaremos lo de la advección a segundo orden
	aEE = zeros(vols)

	if dim > 1
		aN = zeros(vols)
		aS = zeros(vols)
		bN = zeros(vols) 
		bS = zeros(vols)
		aNN = zeros(vols)
		aSS = zeros(vols)
	end
	if dim == 3
		aT = zeros(vols)
		aB = zeros(vols)
		bT = zeros(vols) 
		bB = zeros(vols)
		aTT = zeros(vols)
		aBB = zeros(vols)
	end
end

# ╔═╡ 3d74fabd-a66c-46f8-8ae3-ecac28587efb
begin
	function get_Su(coeff::Coefficients)
		coeff.Su
	end
	    
	function get_Sp_d(coeff::Coefficients)
	    coeff.Sp_d
	end
	
	function get_Sp_a(coeff::Coefficients)
	    coeff.Sp_a
	end
	
	function get_aP(coeff::Coefficients)
	    coeff.aP
	end
	
	function get_aE(coeff::Coefficients)
	    coeff.aE
	end
	
	function get_aW(coeff::Coefficients)
	    coeff.aW
	end
	
	function get_aN(coeff::Coefficients)
	    coeff.aN
	end
	
	function get_aS(coeff::Coefficients)
	    coeff.aS
	end
	
	function get_aT(coeff::Coefficients)
	    coeff.aT
	end
	
	function get_aB(coeff::Coefficients)
	    coeff.aB
	end
	
	function get_bE(coeff::Coefficients)
	    coeff.bE
	end
	
	function get_bW(coeff::Coefficients)
	    coeff.bW
	end
	
	function get_bN(coeff::Coefficients)
	    coeff.bN
	end
	
	function get_bS(coeff::Coefficients)
	    coeff.bS
	end
	
	function get_bT(coeff::Coefficients)
	    coeff.bT
	end
	
	function get_bB(coeff::Coefficients)
	    coeff.bB
	end
	
	function add_source(coeff::Coefficients, source::Real)
	    """
	    
	    """
		#--------------- Aquí falta ver lo del meshgrid------------------
	    x, y, z = [coeff.mesh.deltas_faces[1], coeff.mesh.deltas_faces[2],
	                          coeff.mesh.deltas_faces[3]]
	    coeff.vols = x*y*z
	    # Si es una constante
		coeff.Su += source*coeff.vols
	end
																		
	function add_source(coeff::Coefficients, source)
	    """
	    
	    """ 
	    # Si es una función
		# Aquí se hacen como coordenadas radiales
		crds = sqrt(sum(array([c for c ∈ coeff.mesh.coords if len(c) > 1])^2, 
							  axis = 0))
		coeff.Su += source(crds).reshape(coeff.mesh.volumes)*coeff.vols
	end
end

# ╔═╡ 1415d40a-91bb-4c26-9d31-23ad698285dd
md"## Difusión"

# ╔═╡ a8b43883-37de-4f20-aa2d-97816d13344b
begin
	mutable struct Diffusion
		mesh::Mesh
	end 
end

# ╔═╡ 9cfffb8f-a061-4114-84d8-b376be3a5477
begin
	function set_diffusion(coeff::Coefficients, gamma, velocity_direction::String=None)
	    """
	    
	    """
	    
	    dim = coeff.mesh.dim
	    mesh = coeff.mesh
	    diffusion = Diffusion(mesh, gamma)
		
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
	
	function set_advection(coeff,rho,vel,scheme="upwind1",stgy="b")
	    """
	    
	    """
	    
	    dim = coeff.mesh.dim
	    mesh = coeff.mesh
	    advection = Advection(mesh, rho, vel, scheme)
	    west_c_adv, west_f1_adv, west_f2_adv, west_b1_adv, west_b2_adv, sp_w, su_w, bound_term_w = \
	        get_advection_coef(advection, :W, stgy)
	    east_c_adv, east_f1_adv, east_f2_adv, east_b1_adv, east_b2_adv, sp_e, su_e, bound_term_e = \
	        get_advection_coef(advection, :E, stgy)
	    coeff.aE += east_f1_adv
	    coeff.aW -= west_b1_adv
	    coeff.bW -= bound_term_w
	    coeff.bE -= bound_term_e
	    coeff.aP += -west_c_adv + east_c_adv
	    coeff.Sp_a += sp_w + sp_e
	    coeff.Su -= su_w + su_e
	    coeff.aW += coeff.bW
	    coeff.aE += coeff.bE
				
	    if dim > 1
	        north_c_adv, north_f1_adv, north_f2_adv, north_b1_adv, north_b2_adv, sp_n, su_n, bound_term_n = \
	            advection.get_advection_coef(:N, staggered, velocity_direction, vel,pressure_mesh, stgy)
	        south_c_adv, south_f1_adv, south_f2_adv, south_b1_adv, south_b2_adv, sp_s, su_s, bound_term_s = \
	            advection.get_advection_coef(:S, staggered, velocity_direction, vel,pressure_mesh, stgy)
	        coeff.aN += north_f1_adv
	        coeff.aS -= south_b1_adv
	        coeff.bN -= bound_term_n
	        coeff.bS -= bound_term_s
	        coeff.aP += -north_c_adv + south_c_adv
	        coeff.Sp_a += sp_n + sp_s
	        coeff.Su -= su_n + su_s
	        coeff.aN += coeff.bN
	        coeff.aS += coeff.bS
		end	
	    if dim == 3
	        top_c_adv, top_f1_adv, top_f2_adv, top_b1_adv, top_b2_adv, sp_t, su_t, bound_term_t = \
	            get_advection_coef(advection, :T, stgy)
	        bottom_c_adv, bottom_f1_adv, bottom_f2_adv, bottom_b1_adv, bottom_b2_adv, sp_b, su_b, bound_term_b = \
	            get_advection_coef(advection, :B, stgy)
			
	        coeff.aT += top_f1_adv
	        coeff.aB -= bottom_b1_adv
	        coeff.bT -= bound_term_t
	        coeff.bB -= bound_term_b
	        coeff.aP += -top_c_adv + bottom_c_adv
	        coeff.Sp_a += sp_t + sp_b
	        coeff.Su -= su_t + su_b
	        coeff.aT += coeff.bT
	        coeff.aB += coeff.bB
	    coeff.aP += coeff.Sp_a
		end
	end 
end

# ╔═╡ d44135b0-7825-4742-8829-58e97736fc56
function init_diffusion(diff::Diffusion, Γ::Real)
	""" 

	"""

	diff.Γ = Γ
end

# ╔═╡ 0f0e2768-83bc-45f3-9507-cd8d25209bab
function get_diffusion_coef(diff::Diffusion, direction::Char, stgy::Char='b')
	"""

	"""
	what_dimension, limit = 1,1
	if direction ∈ [:E, :S, :B] limit=-1 end
	if direction ∈ [:N, :S] what_dimension=2
	elseif direction ∈ [:T, :B] what_dimension = 3
	end
        
	get_basic_diffusion_coef(diff, what_dimension, limit)
end

function get_basic_diffusion_coef(diff::Diffusion, what_dimension::Int, limit::Int):
	direction = diff.direction
	mesh = diff.mesh
	faces = diff.mesh.faces[what_dimension]

	#----------------- Aquí vamos de nuevo con lo de Meshgrid----------------
	x, y, z = [mesh.coords[1], mesh.coords[2], mesh.coords[3]]
	δ_d_grid = mesh.get_grid_deltas_faces(axis=what_dimension, orientation=direction)
	δ = δ_d_grid[what_dimension]
	areas = get_area(diff.mesh, direction=what_dimension)
	boundary_area = copy(areas)

	# Para obtener los coeficientes de las caras centrales
	diff = get_interior_terms(diff, x, y, z, direction, faces, what_dimension, limit, areas, δ)

	diff.condition_mask, diff.condition_mask_type = get_mask_boundaries(diff)
	# Aquí obtenemos Sp (que son las fuentes internas)
	# Parece que estas siempre van. Es la contribución de la frontera que no va directamente relacionado
	# con la condición de la frontera. En la matriz del libro, son los extremos de la diagonal principal.
	bound_term_c, bound_term_f = get_boundary_terms(diff, what_dimension, boundary_area, δ)
	sp = bound_term_c # Esto porque sólo contribuye a volumen central y no a alguno que lo rodee.

	# Aquí obtenemos Su (que son las condition_maskes de frontera)
	# En la matriz del libro son los extremos de las diagonales que no son la principal.
	su = get_boundary_condition_mask(diff, what_dimension, limite, boundary_area, δ)
	return diff, sp, su, bound_term_f
end
    
function get_interior_terms(diff::Diffusion, x, y, z, direction, faces, what_dimension, limite, areas, δ):
	λ(x, d) = x[2:end] if d ∈ [:E,:S,:B] else x[:end-1]
	coord = [i!=what_dimension ? var : λ(faces, direction) for (i, var) ∈ enumerate([x,y,z])]
	#----------------- Aquí vamos de nuevo con lo de Meshgrid----------------
	#=
	Γ = diff.Γ(coord[1], coord[2], coord[3])
	coord_l = [(:) if var!=what_dimension else 1:length(Γ[var])-1 for var ∈ 1:3]
	coord_r = [(:) if var!=what_dimension else 2:length(Γ[var]) for var ∈ 1:3]
		
	# Aquí se obtiene el promedio de las Γ's
	Γ_mean = (Γ[coord_l[0], coord_l[1], coord_l[2]] + Γ[coord_r[0], coord_r[1], coord_r[2]])*0.5
	Γ_bar = zeros(Γ.shape)
	=#

	λ_2(d) = d in [:E,:S,:B] ? slice(None,-1) : slice(1,None)
	coord_3 = [var!=what_dimension ? slice(None,None) : λ_2(direction) for var in 1:3]
	Γ_bar[coord_3[1],coord_3[2],coord_3[3]] = Γ_mean
	diff_coef = Γ_bar * areas / δ

	diff.Γ = Γ
	return diff_coef
    
function get_boundary_terms(diff::Diffusion, what_dimension, Γ, boundary_area, δ, stag = False):
	mesh = diff.mesh
	direction = diff.direction
	condition_mask = diff.condition_mask
	condition_mask_type = diff.condition_mask_type
	ba = copy(boundary_area)

	λ_area(d) = d in [:E,:S,:B] ? slice(None,-1) : slice(1,None)
	coord_area = [var!=what_dimension ? slice(None,None) : λ_area(direction) for var in range(3)]
	ba[coord_area[0],coord_area[1],coord_area[2]] = 0

	λ = lambda d: slice(-1,None) if d in ["E","S","B"] else slice(None,1)
	coord = [slice(None,None) if var!=what_dimension else λ(direction) for var in range(3)]
	Γ_bound = Γ[coord[0], coord[1], coord[2]]
	λ_1 = lambda d: slice(-2,-1) if d in ["E","S","B"] else slice(1,2)
	coord_1 = [slice(None,None) if var!=what_dimension else λ_1(direction) for var in range(3)]
	Γ_bound_next = Γ[coord_1[0], coord_1[1], coord_1[2]]
	if stag:
		coef_f, coef_c = self.get_boundary_coefs_staggered(Γ_bound, Γ_bound_next, condition_mask_type)
	else:
		coef_f, coef_c = self.get_boundary_coefs(Γ_bound, Γ_bound_next, condition_mask_type)

	# Para esta parte habría que agregar las areas correspondientes a las diferentes caras
	bound_term_f = (coef_f)* ba / δ
	bound_term_c = (coef_c)* ba / δ

	self.Γ_bound = Γ_bound
	self.Γ_bound_next = Γ_bound_next

	return bound_term_c, bound_term_f
    
function get_boundary_coefs(self, Γ_bound, Γ_bound_next, condition_mask_type):
	initial_shape = Γ_bound.shape
	Γ_bound_flatten = np.ravel(Γ_bound)
	Γ_bound_next_flatten = np.ravel(Γ_bound_next)

	coefs_f = []
	coefs_c = []
	for idx, condition in enumerate(condition_mask_type):
		if condition == "D":
			coef_c = (9*Γ_bound_flatten[idx] - 3*Γ_bound_next_flatten[idx])/2
			coefs_c.append(coef_c)
			coef_f = (Γ_bound_flatten[idx] - Γ_bound_next_flatten[idx]/3)/2
			coefs_f.append(coef_f)
		if condition == "N":
			coef_f = 0
			coefs_f.append(0)
			coefs_c.append(0)
		if condition == "R" or condition == "I":
			# Aquí faltaría definir la condición de Robín con un ejemplo, pero en teoría 
			# todo está en el libro
			coefs_f.append(0)
			coefs_c.append(0)
	coefs_f = np.array(coefs_f).reshape(initial_shape)
	coefs_c = np.array(coefs_c).reshape(initial_shape)

	return coefs_f, coefs_c
        
function get_boundary_coefs_staggered(self, Γ_bound, Γ_bound_next, condition_mask_type):
	initial_shape = Γ_bound.shape
	Γ_bound_flatten = np.ravel(Γ_bound)
	Γ_bound_next_flatten = np.ravel(Γ_bound_next)
	direction = self.direction
	mesh = self._mesh


	coefs_f = []
	coefs_c = []
	for idx, condition in enumerate(condition_mask_type):
		if condition == "D":
			coef_c = (9*Γ_bound_flatten[idx] - 3*Γ_bound_next_flatten[idx])/2
			coefs_c.append(coef_c)
			coef_f = (Γ_bound_flatten[idx] - Γ_bound_next_flatten[idx]/3)/2
			coefs_f.append(coef_f)
		if condition == "N":
			coef_f = 0
			coefs_f.append(0)
			coefs_c.append(0)
		if condition == "R" or condition == "I":
			# Aquí faltaría definir la condición de Robín con un ejemplo, pero en teoría todo 
			# está en el libro
			coefs_f.append(0)
			coefs_c.append(0)  
	# A continuación multiplicamos los valores de la velocidad que ya sabemos de las 
	# condiciones de la frontera
	tags_fronteras = mesh._Mesh__tags_fronteras
	condition_mask = []
	counter = 0
	for tag in tags_fronteras:
		if list(tags_fronteras[tag]["frontera"].keys())[0] == direction:
			cond = list(tags_fronteras[tag]["cond"].keys())[0]
			if cond == "I":
				condition_mask.append(0)
			elif cond == "N":
				# Aquí no es necesariamente 0, pero habría que hacer toda la transformación 
				# como viene en el libro
				condition_mask.append(0)
			elif cond == "D":
				condition_mask.append((tags_fronteras[tag]["cond"][cond]))
			counter += 1

	coefs_f = np.array(condition_mask).reshape(initial_shape)
	coefs_c = np.array(coefs_c).reshape(initial_shape)

	return coefs_f, coefs_c
    
function get_boundary_condition_mask(self, what_dimension, limite, boundary_area, δ):
	direction = self.direction
	Γ_bound = self.Γ_bound
	Γ_bound_next = self.Γ_bound_next
	condition_mask = self.condition_mask
	condition_mask_type = self.condition_mask_type
	Γ_cond = (3*Γ_bound - Γ_bound_next)/2
	conds_su = self.get_mask_boundaries_Su(Γ_cond)
	ba = boundary_area.copy()

	λ_area = lambda d: slice(None,-1) if d in ["E","S","B"] else slice(1,None)
	coord_area = [slice(None,None) if var!=what_dimension else λ_area(direction) for var in range(3)]
	ba[coord_area[0],coord_area[1],coord_area[2]] = 0
	su = ba

	λ = lambda d: slice(-1,None) if d in ["E","S","B"] else slice(None,1)
	coord = [slice(None,None) if var!=what_dimension else λ(direction) for var in range(3)]
	tmp_su = su[coord[0],coord[1],coord[2]]
	su[coord[0],coord[1],coord[2]] = tmp_su*np.array(conds_su).reshape(tmp_su.shape)
	tmp_su = su[coord[0],coord[1],coord[2]]
	div = δ[coord[0],coord[1],coord[2]]*np.array(condition_mask).reshape(tmp_su.shape)
	div = np.where(div == 0., 1, div)
	su[coord[0],coord[1],coord[2]] = tmp_su/div
	return su
    
    # Aquí creamos una máscara que nos dirá qué condición de frontera hay por cada coordenada.
    #### Hay que cambiar dependiendo de qué pase con cada condición de frontera.
function get_mask_boundaries(self, malla):
	direction = self.direction
	tags_fronteras = malla._Mesh__tags_fronteras
	condition_mask = []
	condition_mask_type = []
	dict_cond = {"I":0, "N":0, "D":1}
	for tag in tags_fronteras:
		if list(tags_fronteras[tag]["frontera"].keys())[0] == direction:
			cond = list(tags_fronteras[tag]["cond"].keys())[0]
			condition_mask.append(dict_cond[cond])
			condition_mask_type.append(cond)
	return condition_mask, condition_mask_type
    
function get_mask_boundaries_Su(self, Γ):
	mesh = self._mesh
	direction = self.direction
	tags_fronteras = mesh._Mesh__tags_fronteras
	condition_mask = []
	counter = 0
	for idx, tag in enumerate(tags_fronteras):
		if list(tags_fronteras[tag]["frontera"].keys())[0] == direction:
			cond = list(tags_fronteras[tag]["cond"].keys())[0]
			if cond == "I":
				condition_mask.append(0)
			elif cond == "N":
				condition_mask.append(tags_fronteras[tag]["cond"][cond])
			elif cond == "D":
				condition_mask.append(np.ravel(Γ)[counter]*(8*tags_fronteras[tag]["cond"][cond]/3))
			counter += 1
	return condition_mask

# ╔═╡ aaa69ed2-d5ad-4bc8-825f-1e65b4d36da2
md"## Advección"

# ╔═╡ 07bb4c2d-3e19-4b3d-996b-5d7eba5fbd70


# ╔═╡ 3f1cd4d3-9052-4271-9e4e-8965c4e99170
md"## Sistema de ecuaciones"

# ╔═╡ 35fb94be-3647-4156-b82e-b54baa167fba


# ╔═╡ a934ecd4-5d4c-4764-a518-0fb01aeee45d
md"## Ecuación de Poisson 3D"

# ╔═╡ Cell order:
# ╟─014b9495-588e-4485-a26f-81807624d7b5
# ╟─a5c6a6a8-317a-11ec-25bc-7d7e08c53eb4
# ╟─c2edc6c3-dbcc-4050-8992-e160bd9ebabb
# ╟─226a6ac4-3f23-454a-9d7c-737cc316e11c
# ╟─6eb2ab18-da3f-4433-85f1-ead20aa2fb80
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
# ╟─f5e04936-f01e-487e-831b-b518ddedf103
# ╠═c58ea7b0-daeb-4340-ba1f-cb297ab0ac77
# ╠═34bf980d-eebf-4767-907f-df1be41d7da0
# ╠═feeda8d2-7aef-46fc-94b3-eb5882c0cfb3
# ╠═2bff23a2-ffe8-4288-b473-74ba8d023763
# ╠═93d57deb-f0c7-4df7-bd15-69e95add9d3c
# ╠═8274809e-0d6f-4a94-baac-802f804193ae
# ╠═3b71ed0f-53fb-4382-bcfb-786dd85f0ad1
# ╠═23d190c8-4c51-4f05-b8a4-0573f10d8c11
# ╠═912c556a-4ffd-4e49-a35e-ac83c77c112a
# ╟─b5ed1976-b09a-4a60-83ce-2455d5d63755
# ╠═299fb45f-0579-42a4-8e28-46f8b818ecc9
# ╠═cadfb39c-3487-460b-972e-2af7d348bed9
# ╠═9cfffb8f-a061-4114-84d8-b376be3a5477
# ╠═3d74fabd-a66c-46f8-8ae3-ecac28587efb
# ╟─1415d40a-91bb-4c26-9d31-23ad698285dd
# ╠═a8b43883-37de-4f20-aa2d-97816d13344b
# ╠═d44135b0-7825-4742-8829-58e97736fc56
# ╠═0f0e2768-83bc-45f3-9507-cd8d25209bab
# ╟─aaa69ed2-d5ad-4bc8-825f-1e65b4d36da2
# ╠═07bb4c2d-3e19-4b3d-996b-5d7eba5fbd70
# ╟─3f1cd4d3-9052-4271-9e4e-8965c4e99170
# ╠═35fb94be-3647-4156-b82e-b54baa167fba
# ╟─a934ecd4-5d4c-4764-a518-0fb01aeee45d
