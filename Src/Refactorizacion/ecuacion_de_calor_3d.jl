# Para ver si estamos haciendo bien las cosas en 3D
include("FVM_2.jl")

vols = 10
dims = 3
# Primero todo lo involucrado con el mallado
mesh = FVM.initialize_mesh(dims, vols, 0.1)

# Luego las condiciones de frontera
FVM.tag_wall(mesh, [:W, :E, :S, :N, :B], 0, :D)
FVM.tag_wall(mesh, :T, 100, :D)

# Luego todo lo involucrado con los cálculos del FVM
function Γ_constant(x::Array, y::Array, z::Array)
    Γ = 1000
    tensor_Γ = ones(length(x), length(y), length(z))
    Γ .* tensor_Γ
end
# Le agregamos la difusión
coeff = FVM.init_coefficients(mesh)
FVM.set_diffusion(coeff, Γ_constant);

# Obtenemos la solución
equation_system = FVM.init_eq_system(coeff)
solution = FVM.get_solution(equation_system);

FVM.plot_solution(solution, mesh)