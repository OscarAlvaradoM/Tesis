module FVM

	using Markdown
	using InteractiveUtils
	using LinearAlgebra, SparseArrays
	using Plots
	plotly()



	md"## Mallado"

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
		grid_deltas = grid(deltas[1], deltas[2], deltas[3])
		grid_coords = grid(centers[1], centers[2], centers[3])
		grid_faces = grid(faces[1], faces[2], faces[3])
		return grid_deltas, grid_coords, grid_faces
	end	

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

	mutable struct Mesh
		volumes
		lengths
		centers::Array
		centers_and_boundaries::Array
		deltas::Array
		faces::Array
		deltas_faces::Array
		tags
		tags_boundaries;
	end 

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

	md"## Coeficientes"

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

	function add_source(coeff, source::Real)
		"""

		"""
		x, y, z = grid(coeff.mesh.deltas_faces[1], coeff.mesh.deltas_faces[2],
							  coeff.mesh.deltas_faces[3])
		# Si es una constante
		coeff.Su += source*(x*y*z)
	end

	md"## Difusión"


	mutable struct Diffusion
		coeff
		Γ
	end 

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
		ba = [copy(array) for array ∈ array_before_boundary_area]

		λ_area(d, idx) = d ∈ [:E,:S,:B] ? (length(ba[idx]):length(ba[idx])) : (1:1)
		coord_area = [idx!=what_dimension ? (:) : λ_area(direction, idx) for idx ∈ 1:3]
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
		coeff.aW += coeff.bE;
		
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
			coeff.aN += coeff.bS;
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
		coeff.aP -= coeff.Sp_d;
		end
	end

	md"## Sistema de ecuaciones"

	mutable struct EqSystem
		A::Array
		b::Vector
	end

	function init_eq_system(coeff)
		mesh = coeff.mesh
		volumes = mesh.volumes
		dim = sum(volumes .!= 1)
	    aP = coeff.aP
	    aE, aW = coeff.aE, coeff.aW
		if dim == 1
			A = get_matrix_A(aP, aW, aE)
		elseif dim > 1
			aN, aS = coeff.aN, coeff.aS
			if dim == 2
				A = get_matrix_A(aP, aW, aE, aN, aS)
			elseif dim == 3	
				aT, aB = coeff.aT, coeff.aB
				A = get_matrix_A(aP, aW, aE, aN, aS, aT, aB, volumes)
			end 
		end
	    
	    Su = coeff.Su
		b = get_vector_b(Su)
		
		eq_system = EqSystem(A, b)
	end
	    
	    
	function get_matrix_A(aP::Array, aW::Array, aE::Array)
	    """
	    Método para construir la matriz A, de la ecuación Ax=b, a partir de los coeficientes obtenidos.
	    """
	    aP = [(aP...)...]
	    A = get_diag(aP)
		aW = [(aW...)...][2:end]
		aE = [(aE...)...][1:(end-1)]
		A += get_diag(aE, 1) + get_diag(aW, -1)
		return A
	end
		
	function get_matrix_A(aP::Array, aW::Array, aE::Array, aN::Array, aS::Array)
		aP = [(aP...)...]
	    A = get_diag(aP)
		aN = [(aN...)...][2:end]
		aS = [(aS...)...][1:(end-1)]
		aW = [(aW...)...][volumes[2]:end]
		aE = [(aE...)...][1:(end-volumes[2])]
		A += get_diag(aN, -1) + get_diag(aS, 1) + get_diag(aE, volumes[2]) + get_diag(aW, -volumes[2])
		return A
	end

	function get_matrix_A(aP::Array, aW::Array, aE::Array, aN::Array, aS::Array,
			aT::Array, aB::Array, volumes::Array)
		aP = [(aP...)...]
	    A = get_diag(aP)
		aT = [(aT...)...][2:end]
		aB = [(aB...)...][1:(end-1)]
		
		aN = [(aN...)...][(volumes[2]+1):end]
		aS = [(aS...)...][1:(end-volumes[2])]
		
		aW = [(aW...)...][((volumes[2]*volumes[3])+1):end]
		aE = [(aE...)...][1:(end-volumes[2]*volumes[3])]
		A += get_diag(aT, -1) + get_diag(aB, 1) + get_diag(aN, -volumes[2]) + get_diag(aS, volumes[2]) + get_diag(aW, -volumes[2]*volumes[3]) + get_diag(aE, volumes[2]*volumes[3])
	    return A
	end
			
	function get_diag(array, k=0)
	    """
	    Método para construir una matriz diagonal dado un arreglo 2D y 
		la diagonal en que queremos poner ese arreglo
	    """
		N = size(array)[1]
				
	    if N > 50
		return spdiagm(k => array)
	    else
		return diagm(k => array)
		end
	end


	function get_vector_b(Su)
	    """
	    Método para obtener el vector b, de la ecuación Ax = b, a partir del arreglo Su construido anteriormente
	    """
	    b = [(Su...)...]
	    return b
	end

	function solve(A,b)
	    return A\b
	end

	function get_solution(eq_system::EqSystem)
	    """
	    Método para obtener los valores de la solución al sistema
	    """
	    A = eq_system.A
	    b = eq_system.b
	    return solve(A, b)
	end
	
	function plot_solution(solution::Array, mesh::Mesh)
		(X,Y,Z) = mesh.centers
		coordinates = [[x,y,z] for x ∈ X, y ∈ Y,  z ∈ Z]
		x = [[coordinate[1] for coordinate ∈ coordinates]...]
		y = [[coordinate[2] for coordinate ∈ coordinates]...]
		z = [[coordinate[3] for coordinate ∈ coordinates]...]
		scatter(x,y,z, marker_z=solution, color=:plasma)
    end
end
