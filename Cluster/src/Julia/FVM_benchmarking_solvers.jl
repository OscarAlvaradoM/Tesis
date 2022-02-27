using Statistics, CSV, Tables

#include("../notebooks/FVM.jl")
include("FVM.jl")

include("solvers.jl")

function Γ_constant(x::Array, y::Array, z::Array)
    Γ = 1000
    tensor_Γ = ones(length(x), length(y), length(z))
    tensor_Γ = Γ .* tensor_Γ
end

function getsystem(vol)
    # Construcción de la malla.
    volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces = FVM.uniform_grid(vol, vol/10)
    tags = FVM.init_tags(3, volumes, centers_and_boundaries)
    tags_b = FVM.init_tags_boundaries(3, centers_and_boundaries)
    FVM.tag_wall(tags, tags_b, [:W, :E, :T, :N, :B], 0, :D)
    FVM.tag_wall(tags, tags_b, :S, 100, :D)
    mesh = FVM.Mesh(volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces, tags, tags_b)
    # Generando los coeficientes
    coeff = FVM.init_coefficients(mesh)
    FVM.set_diffusion(coeff, Γ_constant);
    # Construyendo el sistema de ecuaciones
    equation_system = FVM.init_eq_system(coeff)
    A = equation_system.A
    b = equation_system.b
    return A, b
end

function iteratevolumes(volumes::Array{Int, 1}, times::Int)
    for volume ∈ volumes
        A, b = getsystem(volume)
        solvers = [Solvers.vectorizedjacobi, Solvers.paralleljacobi, 
        Solvers.vectorizedgaussseidel, Solvers.parallelgaussseidel, 
        Solver.vectorizedsor, Solver.parallelsor, 
        Solver.gmres, Solver.parallelgmres, 
        Solver.gmresreiniciado, Solver.parallelgmresreiniciado, 
        Solver.gmresprecondicionado, Solver.parallelgmresprecondicionado,
        Solver.gmresprecondicionadoreiniciado , Solver.parallelgmresprecondicionadoreiniciado]
        names = ["Jacobi", "Gauss-Seidel", "SOR", "Jacobi", "Gauss-Seidel", "SOR"]
        idx_1 = 1
        idx_2 = 1
        for solver ∈ solvers
            println("Comencé el de $volume volúmenes con el solver $solver")
            if "precondicionado" in string(solver)
            	if "parallel" not in string(solver)
            		list_of_time_statistics = list_of_statistics(solver, names[idx_1], A, b, volume, times)
            		idx_1+=1
	            	writetofile(solver, names[idx_1], list_of_time_statistics')            		
            	else
            		list_of_time_statistics = list_of_statistics(solver, names[idx_2], A, b, volume, times)
            		idx_2+=1
	            	writetofile(solver, names[idx_2], list_of_time_statistics')
            	end
            else
	            list_of_time_statistics = list_of_statistics(solver, A, b, volume, times)
	            writetofile(solver, list_of_time_statistics')
	        end
	    @show list_of_time_statistics
        end
        println("Terminé el de $volume volúmenes")
    end
end

function list_of_statistics(solver, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, volume::Int, times::Int)
    time_list, μ, σ, iterations_list, μ_iterations, σ_iterations = getstatistics(solver, A, b, times)
    time_statistics = [volume, volume^3, μ, σ, μ_iterations, σ_iterations]
    push!(time_statistics, time_list...)
    push!(time_statistics, iterations_list...)
    return time_statistics
end

function list_of_statistics(solver, name::String, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, volume::Int, times::Int)
    time_list, μ, σ, iterations_list, μ_iterations, σ_iterations = getstatistics(solver, name, A, b, times)
    time_statistics = [volume, volume^3, μ, σ, μ_iterations, σ_iterations]
    push!(time_statistics, time_list...)
    push!(time_statistics, iterations_list...)
    return time_statistics
end

function getstatistics(solver, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, times::Int)
    time_list = []
    iterations_list = []
    for i ∈ 1:(times+1)
	println("Estoy en la medición $i del tiempo")
        t, iterations = measuretime(solver, A, b)
        push!(time_list, t)
        push!(iterations_list, iterations)
	@show t
	@show iterations
    end
    times_without_compiling = time_list[2:end]
    iterations_without_compiling = iterations_list[2:end]
    μ = mean(times_without_compiling)
    σ = std(times_without_compiling)
    μ_iterations = mean(iterations_without_compiling)
    σ_iterations = std(iterations_without_compiling)
    
    return time_list, μ, σ, iterations_without_compiling, μ_iterations, σ_iterations
end

function getstatistics(solver, name::String, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, times::Int)
    time_list = []
    iterations_list = []
    for i ∈ 1:(times+1)
	println("Estoy en la medición $i del tiempo")
        t, iterations = measuretime(solver, name, A, b)
        push!(time_list, t)
        push!(iterations_list, iterations)
	@show t
	@show iterations
    end
    times_without_compiling = time_list[2:end]
    iterations_without_compiling = iterations_list[2:end]
    μ = mean(times_without_compiling)
    σ = std(times_without_compiling)
    μ_iterations = mean(iterations_without_compiling)
    σ_iterations = std(iterations_without_compiling)
    
    return time_list, μ, σ, iterations_without_compiling, μ_iterations, σ_iterations
end

function measuretime(solver, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64})
    start_time = time()
    iterations = solver(A,b)
    finish_time = time()
    t = finish_time - start_time
    return t, iterations
end

function measuretime(solver, name, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64})
    start_time = time()
    iterations = solver(A, b, name)
    finish_time = time()
    t = finish_time - start_time
    return t, iterations
end

function writetofile(solver, list)
    file_name = string(solver)
    CSV.write("../../benchmarking/solvers/$(file_name).csv", Tables.table(list), delim = ',',append=true)
    #CSV.write("benchmarking/$(file_name).csv", Tables.table(list), delim = ',',append=true)
end

function writetofile(solver, name, list)
    file_name = string(solver)
    CSV.write("../../benchmarking/solvers/$(file_name)_$(name).csv", Tables.table(list), delim = ',',append=true)
end

volumes = [10,20,30,40,50,60,80,100]
times = 10
iteratevolumes(volumes, times)
