### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ a5c6a6a8-317a-11ec-25bc-7d7e08c53eb4
begin
	using PlutoUI
	PlutoUI.TableOfContents(aside=true, title="🎣 Método del Volumen Finito 🔥")
end

# ╔═╡ 014b9495-588e-4485-a26f-81807624d7b5
md"""
# Método del Volumen Finito
- Óscar A. Alvarado-Morán
- Oscar A. Esquivel-flores
"""


# ╔═╡ 226a6ac4-3f23-454a-9d7c-737cc316e11c
md"## Mallado"

# ╔═╡ de5fb957-0121-4830-974a-4587259da8fb
begin
	mutable struct Mesh
		dim::Int
	end 
end

# ╔═╡ ed1a1dbd-6aab-48b3-9a06-354cd6adc77e
begin
	function set_volumes_and_lengths(mesh::Mesh, volumes::Real, lengts::Real)
		mesh.volumes = (volumes, 1, 1)
		mesh.lengths = [lengths, lengths/10, lengths/10]
	end
	
	function set_volumes_and_lengths(mesh::Mesh, volumes::Array, lengts::Array)
		missing_volumes = 3 - length(volumes)
		mesh.volumes = vcat(volumes, [1 for _ ∈ 1:missing_volumes])
		
		missing_lengths = 3 - length(lengths)
		mesh.lengths = vcat(lengths, [lengths[1]/10 for _ in 1:missing_lengths])
	end
end

# ╔═╡ c58ddc15-3824-4ef4-b371-1b3ce6fb4fe3
function set_deltas(dominio::Array)
	(dominio[2:end] - dominio[:end-1])[2:end-1]
end

# ╔═╡ b59b0526-5e2a-478f-a16d-4de6572bab5e
function unifor_grid(mesh::Mesh)
	# Separación entre todos los nodos de cada dimensión
	volume_length = mesh.lengths/mesh.volumes
	# La frontera "inicial" del arreglo
	start = volume_length/2
	# La frontera "final" del arreglo
	stop = mesh.lengths - volume_length/2
	
	mesh.centers = /
	[LinRange(strt, stp, vol) for (strt, stp, vol) ∈ zip(start, stop, mesh.volumes)]
	
	mesh.centers_and_boundaries = /
	[insert!(arr,1,length(arr)) for (idx, arr) ∈ enumerate(mesh.centers)]
	
	mesh.centers_and_boundaries = /
	[insert!(arr,1,volume_length[idx]) for (idx, arr) ∈ enumerate(mesh.centers)]
	
	mesh.deltas = [length(set_deltas(dom)) != 0 ? set_deltas(dom) : [dom[end-1]] for dom ∈ mesh.centers_and_boundaries]
	
	mesh.centers_and_boundaries = [[dom] for dom ∈ centers_and_boundaries]
	
	mesh.faces = [vcat([mesh.centers_and_boundaries[idx][1]], 
			[(coords[1:end] + coords[2:end])/2] , 
			[mesh.centers_and_boundaries[idx][end]]
			) for (idx, coords) ∈ enumerate(mesh.centers)]
	
	get_deltas_faces(mesh)
	get_grids(mesh)
end

# ╔═╡ eb9f3443-0032-41ca-8a55-9a38d8a2e3a9
begin
	function init_tags(mesh::Mesh)
		
		mesh.tags = Dict()
		X, Y, Z = mesh.volumes
		for z ∈ 2:Z, y ∈ 2:Y, x ∈ 2:X                   
			t = b = n = s = "Off"
			e = w = "F"
			if x == 2 w = Dict()
			elseif x == X e = Dict()
			end
			if mesh.dim > 1
				n = s = "F"
				if y == 2 n = Dict()
				elseif y == Y s = Dict
				end
				if mesh.dim == 3
					t = b = "F"
					if z == 2 t = Dict()
					elseif z == Z b = Dict()
					end
				end
			end
			mesh.tags["$(x)-$(y)-$(z)"] = Dict("E" => e, "W" => w, "N" => n, "S" => s,
				"T" => t, "B" => b, "coord" => [mesh.centers_and_boundaries[1][x],
					mesh.centers_and_boundaries[2][y],
					mesh.centers_and_boundaries[3][z]])
		end
	end

	
	function init_tags_boundaries(mesh::Mesh)

		mesh.tags_boundaries = Dict()
		X, Y, Z = [length(dom) for dom ∈ mesh.centers_and_boundaries]
		for z ∈ 1:Z, y ∈ 1:Y, x ∈ 1:X
			t = b = n = s = "Off"
			e = w = "ON"
			if mesh.dim > 1
				n = s = "ON"
				if mesh.dim == 3 t = b = "ON" end
			end
			# El siguiente cacho de código es para saber si nos encontramos con una frontera
			if x==1 | y==1 | z==1 | x==X | y==Y | z==Z
				var = ""
				
				if y != 1 & y != Y
					if z != 1 & z != Z
						if x == 1 var = "W"; value = w
						elseif x == X var = "E"; value = e
						else continue
						end
					elseif x != 1 & x != X
						if z == 1 var = "T"; value = t
						elseif z == Z var = "B"; value = b
						else continue
						end
					else continue
					end
					mesh.tags_boundaries["$(x)-$(y)-$(z)"] = /
					Dict("frontera" => Dict(var => value), 
						"coord" => [mesh.centers_and_boundaries[1][x],
							mesh.centers_and_boundaries[2][y],
							mesh.centers_and_boundaries[3][z]], 
						"cond" => Dict())
					
				elseif z != 1 & z != Z
					if x != 1 & x != X
						if y == 1 var = "N"; value = n
						elseif y == Y var = "S"; value = s
						mesh.tags_boundaries["$(x)-$(y)-$(z)"] = /
							Dict("frontera" => Dict(var => value), 
								"coord" => [mesh.centers_and_boundaries[1][x],
									mesh.centers_and_boundaries[2][y],
									mesh.centers_and_boundaries[3][z]], 
								"cond" => Dict())
						end
					end
				else continue
				end
			end
		end
	end
end

# ╔═╡ b8e3439d-5797-42dd-a6ee-44a27d913408
function init_mesh(mesh::Mesh, volumes::Real, lengths::Real)
	mesh.volumes = [1,1,1]
	mesh.lengths = [0.01, 0.01, 0.01]

	# --Valores por functionecto para posición y separación del mallado--

	# Coordenadas de los centros de los volúmenes
	mesh.centers = [[0.005] for _ ∈ 1:3]
	# Coordenadas de los centros y además las orillas.
	mesh.centers_and_boundaries = [[[0]+[mesh.centers[i][1]]+[mesh.lengths[i]]] for i ∈ 1:3]
	# Diferencia entre los centros de los volúmenes
	deltas = [[0] for _ ∈ 1:3]
	# --Hasta aquí tenemos un cubito--

	set_volumes_and_lengths(volumes, length)
	
	#---------------------------------------------------------------
	if volumes & lengths
		uniform_grid(mesh)
		init_tags(mesh)
		init_tags_boundaries(mesh)
	end
end

# ╔═╡ 1984a707-6cd3-43f2-8b06-4cc99781a38c
begin
	function tag_wall(mesh::Mesh, direction::Real, tag::Symbol, value::Real)

        for key ∈ keys(mesh.tags)
            if typeof(mesh.tags[key][direction]) == Dict
                mesh.tags[key][direction][tag] = value
			end
		end
        for key ∈ keys(mesh.tags_boundaries)
            if mesh.tags_boundaries[key]["frontera"][direction] == "ON"
                mesh.tags_boundaries[key]["cond"][tag] = value
			end
		end
	end

    function tag_wall_dirichlet(mesh::Mesh, directions::Array, value::Real, coords::Array)
        
		for (idx, key) ∈ enumerate(coords)
			if key ∈ keys(mesh.tags)
				mesh.tags[key][directions[idx]]["D"] = value[idx]
			elseif key ∈ keys(mesh.tags_boundaries)
				mesh.tags_boundaries[key]["cond"]["D"] = value[idx]
			end
		end
	end
	
	function tag_wall_dirichlet(mesh::Mesh, directions::Array, value::Real)
  
		for (idx, direction) ∈ enumerate(directions)
			tag_wall(direction, "D", value[idx])
		end
	end
	
	function tag_wall_dirichlet(mesh::Mesh, direction::Real, value::Real)
		tag_wall(direction, "D", value)
	end
                    
    function tag_wall_neumann(mesh::Mesh, directions::Array, value::Real, coords::Array)
      
		for (idx, key) ∈ enumerate(coords)
			if key ∈ keys(mesh.tags)
				mesh.tags[key][direction[idx]]["N"] = value[idx]
			elseif key ∈ keys(mesh.tags_boundaries)
				mesh.tags_boundaries[key]["cond"]["N"] = value[idx]
			end
		end
	end
				
	function tag_wall_neumann(mesh::Mesh, directions::Array, value::Real)
		
		for (idx, direction) ∈ enumerate(directions)
			tag_wall(direction, "N", value[idx])
		end
	end
	
	function tag_wall_neumann(mesh::Mesh, direction::Real, value::Real)
			tag_wall(direction, "N", value)
	end
	
	function tag_wall_insulated(mesh::Mesh, directions::Array, coords::Array)
        
		for (idx, key) ∈ enumerate(coords)
			if key ∈ kes(mesh.tags)
				mesh.tags[key][directions[idx]]["I"] = None
			elseif key ∈ keys(mesh.tags_boundaries)
				mesh.tags_boundaries[key]["cond"]["I"] = None
			end
		end
	end
	
	function tag_wall_insulated(mesh::Mesh, directions::Array)
            
		for (idx, direction) in enumerate(directions)
			mesh.tag_wall(direction, "I", None)
		end
	end
	
    function tag_wall_insulated(mesh::Mesh, direction::Real)
		tag_wall(direction, "I", None)
	end
end

# ╔═╡ 05bee847-1a28-49a4-89a6-9acdfcd15d35
begin
	function set_centers_and_boundaries(mesh::Mesh, centers_and_boundaries::Array)
		
		set_principal_things(mesh, centers_and_boundaries)
		set_faces(mesh)
		set_secondary_things()
	end
	
    function set_centers_and_boundaries(mesh::Mesh, centers_and_boundaries::Array, faces::Array)
		
		set_principal_things(mesh, centers_and_boundaries)
		set_faces(mesh, faces)
		set_secondary_things()
	end
       
	function set_principal_things(mesh::Mesh, centers_and_boundaries::Array)
        
        # Tendría que ser una tupla de tuplas/listas/arreglos para que sea válido.
        if typeof(centers_and_boundaries) == Array{Real}
            tmp_array = vcat(centers_and_boundaries, 
				mesh.centers_and_boundaries[2], 
				mesh.centers_and_boundaries[3])
            centers_and_boundaries = tmp_array
		end
        # Asigna los atributos de la mesh correspondientes    
        mesh.centers_and_boundaries = [centers_and_boundaries[i] for i in 1:3]
		mesh.centers = [centers_and_boundaries[i][2:end-1] for i ∈ 1:3]
		mesh.lengths = [centers_and_boundaries[i][end-1] for i ∈ 1:3]
		mesh.volumes = [length(centers_and_boundaries[i][2:end-1]) for i ∈ 1:3]       
		mesh.deltas = [length(set_deltas(dom)) != 0 ? set_deltas(dom) : dom[end] for dom ∈ mesh.centers_and_boundaries]
	end
	
	function set_secondary_things(mesh::Mesh)	
		get_deltas_faces(mesh)
        init_tags(mesh)
        init_tags_boundaries(mesh)
        get_grids(mesh)
	end
        
	function set_faces(mesh::Mesh, faces::Array{Real})
		
        # Si me está pasando una lista (o sea, es de una dimensión)
		mesh.faces = hcat(centers_and_boundaries[1],
					faces,
					centers_and_boundaries[end])
	end
			
	function set_faces(mesh::Mesh, faces::Array{Array})
		
		# Suponemos aquí que nos está pasando una lista de listas (o tupla de tuplas)
		for (idx, face_1dim) ∈ enumerate(faces)
			faces[idx] = hcat(centers_and_boundaries[idx][1], 
				face_1dim, 
				dominio[idx][end])
		end
		mesh.faces = faces
	end
	
	function set_faces(mesh::Mesh)
	
		mesh.faces = [hcat(mesh.centers_and_boundaries[idx][0],
				(center[:end-1] + center[2:end])/2,
				mesh.centers_and_boundaries[idx][end])
			for (idx, center) ∈ enumerate(mesh.centers)]
	end
end


# ╔═╡ a4225cfd-7b75-4d9f-80d0-70a9aa4a007b
[i+j+k for i in 1:3, j in 3:6, k in 6:9]

# ╔═╡ 3b71ed0f-53fb-4382-bcfb-786dd85f0ad1
begin
	function get_grid_deltas_centers_and_boundaries(mesh::Mesh, axis=0::Int, orientation='E'::Char, reverse=false::Bool)
		deltas_centers_and_boundaries = []
		mesh.grid_deltas_centers_and_boundaries = []
		for direction ∈ 1:3
			if direction != axis
				if length(centers_and_boundaries[direction]) == 3
					push!(deltas_centers_and_boundaries, mesh.centers_and_boundaries[direction][2])
				else
					push!(deltas_centers_and_boundaries, mesh.centers[direction])
				end
			else
				push!(deltas_centers_and_boundaries, 
					mesh.centers_and_boundaries[direction][2:end] /
					- mesh.centers_and_boundaries[direction][:end-1])
			end
		end
			
		if orientation == 'E' | orientation == 'S' | orientation == 'B'
			if reverse
				deltas_centers_and_boundaries[axis] = /
				deltas_centers_and_boundaries[axis][:end-1]
			else
				deltas_centers_and_boundaries[axis] = /
				deltas_centers_and_boundaries[axis][2:]
			end
				
		else
			if reverse:
				deltas_centers_and_boundaries[axis] = /
				deltas_centers_and_boundaries[axis][2:end]
			else
				deltas_centers_and_boundaries[axis] = /
				deltas_centers_and_boundaries[axis][:end-1]
			end
		end	
		#--------------- Aquí falta ver lo del meshgrid------------------
		mesh.grid_deltas_centers_and_boundaries = [deltas_centers_and_boundaries[1],
			deltas_centers_and_boundaries[2], deltas_centers_and_boundaries[3]]
		return mesh.grid_deltas_centers_and_boundaries
	end
								

	function get_grid_deltas_faces(mesh::Mesh, axis=0::Int, orientation='E'::Char, reverse=false::Bool):
		deltas_faces = []
		mesh.grid_deltas_faces = []
		for direction ∈ 1:3
			if direction != axis
				if length(centers_and_boundaries[direction]) == 3
					deltas_faces.append(centers_and_boundaries[direction][2])
				else
					push!(deltas_faces, mesh.centers[direction])
				end
				
			else
				push!(deltas_faces, mesh.faces[direction][2:end] - meshfaces[direction][:end-1])
			end
		end
		#--------------- Aquí falta ver lo del meshgrid------------------
		mesh.grid_deltas_faces = [deltas_faces[1], deltas_faces[2], deltas_faces[3]]
		return mesh.grid_deltas_faces


	function get_deltas_faces(mesh):
		mesh.deltas_faces = []
		for direction in 1:3
				mesh.deltas_faces.append(mesh.faces[direction][2:end] - meshfaces[direction][:end-1])

	function get_grids(mesh):
		#--------------- Aquí falta ver lo del meshgrid------------------
		mesh.grid_deltas = [mesh.deltas[1], mesh.deltas[2], mesh.deltas[3]]
		mesh.grid_coords = [mesh.centers[1], mesh.centers[2], mesh.centers[3]]
		mesh.grid_faces = [mesh.faces[1], mesh.faces[2], mesh.faces[3]]
end

# ╔═╡ 2bff23a2-ffe8-4288-b473-74ba8d023763
function draw(mesh):
        """
        Método para graficar la malla. Este método se invoca hasta que se hayan inizializado todas las condiciones
        de frontera.
        """
        # Graficamos las fronteras, sean o no activas
        dic_colors = {"D": "darkturquoise", "N": "red", "S": "magenta", "Off": "white", "I": "gray"}
        condiciones = [list(mesh.__tags_fronteras[key]["cond"].keys())[0] \
                       if list(mesh.__tags_fronteras[key]["frontera"].values())[0] == "ON" \
                       else "Off" for key in list(mesh.__tags_fronteras.keys())]
        colores = [dic_colors[cond] for cond in condiciones]
        # Obtenemos las coordenadas de las fronteras y de los nodos internos.
        coordenadas = [] # Aquí se pondrán las coordenadas de las fronteras
        coord = [] # Aquí las coordendas de los nodos internos
        for i in range(3):
            coordenadas.append([mesh.__tags_fronteras[key]["coord"][i] \
                                for key in list(mesh.__tags_fronteras.keys())])
            coord.append([mesh.__tags[key]["coord"][i] for key in list(mesh.__tags.keys())])
        fig = go.Figure(data = go.Scatter3d(x = coordenadas[0], y = coordenadas[1], z = coordenadas[2],
                                          mode = 'markers', marker = dict(color = colores, 
                                                                          symbol = "square", size = 2)))
        fig.add_trace(go.Scatter3d(x = coord[0], y = coord[1], z = coord[2],
                                              mode = 'markers', marker = dict(color = "blue", size = 5)))
        fig.show()

# ╔═╡ 5b4d3eeb-125e-4fe5-bd7a-2af50e3de107
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

# ╔═╡ 9cfffb8f-a061-4114-84d8-b376be3a5477
begin
	function set_diffusion(coeff::Coefficients, gamma, velocity_direction::String=None)
	    """
	    
	    """
	    
	    dim = coeff.mesh.dim
	    mesh = coeff.mesh
	    diffusion = Diffusion(mesh, gamma)
		
	    west_diff, sp_w, su_w, bound_term_w = get_diffusion_coef(diffusion,"W")
	    east_diff, sp_e, su_e, bound_term_e = get_diffusion_coef(diffusion,"E")
		
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
	        north_diff, sp_n, su_n, bound_term_n = get_diffusion_coef(diffusion,"N")
	        south_diff, sp_s, su_s, bound_term_s = get_diffusion_coef(diffusion, "S")
			
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
	        top_diff, sp_t, su_t, bound_term_t = get_diffusion_coef(diffusion, "T")
	        bottom_diff, sp_b, su_b, bound_term_b = get_diffusion_coef(diffusion, "B")
			
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
	        get_advection_coef(advection, "W", stgy)
	    east_c_adv, east_f1_adv, east_f2_adv, east_b1_adv, east_b2_adv, sp_e, su_e, bound_term_e = \
	        get_advection_coef(advection, "E", stgy)
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
	            advection.get_advection_coef("N", staggered, velocity_direction, vel,pressure_mesh, stgy)
	        south_c_adv, south_f1_adv, south_f2_adv, south_b1_adv, south_b2_adv, sp_s, su_s, bound_term_s = \
	            advection.get_advection_coef("S", staggered, velocity_direction, vel,pressure_mesh, stgy)
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
	            get_advection_coef(advection, "T", stgy)
	        bottom_c_adv, bottom_f1_adv, bottom_f2_adv, bottom_b1_adv, bottom_b2_adv, sp_b, su_b, bound_term_b = \
	            get_advection_coef(advection, "B", stgy)
			
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
	    x, y, z = np.meshgrid(coeff.mesh.deltas_faces[1],coeff.mesh.deltas_faces[2],
	                          coeff.mesh.deltas_faces[3])
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
# ╟─226a6ac4-3f23-454a-9d7c-737cc316e11c
# ╠═de5fb957-0121-4830-974a-4587259da8fb
# ╠═b8e3439d-5797-42dd-a6ee-44a27d913408
# ╠═ed1a1dbd-6aab-48b3-9a06-354cd6adc77e
# ╠═b59b0526-5e2a-478f-a16d-4de6572bab5e
# ╠═c58ddc15-3824-4ef4-b371-1b3ce6fb4fe3
# ╠═eb9f3443-0032-41ca-8a55-9a38d8a2e3a9
# ╠═1984a707-6cd3-43f2-8b06-4cc99781a38c
# ╠═05bee847-1a28-49a4-89a6-9acdfcd15d35
# ╠═a4225cfd-7b75-4d9f-80d0-70a9aa4a007b
# ╠═3b71ed0f-53fb-4382-bcfb-786dd85f0ad1
# ╠═2bff23a2-ffe8-4288-b473-74ba8d023763
# ╠═5b4d3eeb-125e-4fe5-bd7a-2af50e3de107
# ╟─b5ed1976-b09a-4a60-83ce-2455d5d63755
# ╠═299fb45f-0579-42a4-8e28-46f8b818ecc9
# ╠═cadfb39c-3487-460b-972e-2af7d348bed9
# ╠═9cfffb8f-a061-4114-84d8-b376be3a5477
# ╠═3d74fabd-a66c-46f8-8ae3-ecac28587efb
# ╟─1415d40a-91bb-4c26-9d31-23ad698285dd
# ╠═a8b43883-37de-4f20-aa2d-97816d13344b
# ╟─aaa69ed2-d5ad-4bc8-825f-1e65b4d36da2
# ╠═07bb4c2d-3e19-4b3d-996b-5d7eba5fbd70
# ╟─3f1cd4d3-9052-4271-9e4e-8965c4e99170
# ╠═35fb94be-3647-4156-b82e-b54baa167fba
# ╟─a934ecd4-5d4c-4764-a518-0fb01aeee45d
