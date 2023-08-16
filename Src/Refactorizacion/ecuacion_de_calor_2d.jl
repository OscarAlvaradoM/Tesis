# Para ver si estamos haciendo bien las cosas en 3D
include("FVM.jl")

vols = 3
dims = 2
# Primero todo lo involucrado con el mallado
mesh = FVM.initialize_mesh(dims, vols, 0.1)

# Luego las condiciones de frontera
FVM.tag_wall(mesh, [:E, :W, :N], 0, :D)
FVM.tag_wall(mesh, :S, 100, :D)

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
# print(equation_system.A)
# print(equation_system.b)

#print(coeff)
FVM.plot_solution(solution, mesh)