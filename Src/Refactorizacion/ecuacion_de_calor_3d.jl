# Para ver si estamos haciendo bien las cosas en 3D
include("FVM.jl")

vols = 10
dim = 3
# Primero todo lo involucrado con el mallado
volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces = FVM.uniform_grid(vols, 10)
tags = FVM.init_tags(dim, volumes, centers_and_boundaries)
tags_b = FVM.init_tags_boundaries(dim, centers_and_boundaries)

FVM.tag_wall(tags, tags_b, [:W, :E, :T, :N, :B], 0, :D)
FVM.tag_wall(tags, tags_b, :S, 100, :D)

mesh = FVM.Mesh(volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces, tags, tags_b);

# Luego todo lo involucrado con los cálculos del FVM
#Γ = 1000
#tensor_Γ = ones(length(x), length(y), length(z))
#tensor_Γ = Γ .* tensor_Γ
function Γ_constant(x::Array, y::Array, z::Array)
    Γ = 1000
    tensor_Γ = ones(length(x), length(y), length(z))
    tensor_Γ = Γ .* tensor_Γ
end
# Le agregamos la difusión
coeff = FVM.init_coefficients(mesh)
FVM.set_diffusion(coeff, Γ_constant);

# Obtenemos la solución
equation_system = FVM.init_eq_system(coeff)
solution = FVM.get_solution(equation_system);

FVM.plot_solution(solution, mesh)


