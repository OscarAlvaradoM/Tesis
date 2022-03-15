using BenchmarkTools, Statistics, CSV, Tables

include("FVM.jl")

function benchmark_mesh(vols::Int)
    volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces = FVM.uniform_grid(vols, 10)
    tags = FVM.init_tags(3, volumes, centers_and_boundaries)
    tags_b = FVM.init_tags_boundaries(3, centers_and_boundaries)

    FVM.tag_wall(tags, tags_b, [:W, :E, :T, :N, :B], 0, :D)
    FVM.tag_wall(tags, tags_b, :S, 100, :D)

    mesh = FVM.Mesh(volumes, lengths, centers, centers_and_boundaries, deltas, faces, deltas_faces, tags, tags_b);
end

function Γ_constant(x::Array, y::Array, z::Array)
    Γ = 1000
    tensor_Γ = ones(length(x), length(y), length(z))
    tensor_Γ = Γ .* tensor_Γ
end

function benchmark_set_boundary_conditions(mesh::FVM.Mesh)
    coeff = FVM.init_coefficients(mesh)
    FVM.set_diffusion(coeff, Γ_constant);
    return coeff
end

function benchmark_solutions(coeff::FVM.Coefficients_3d)
    equation_system = FVM.init_eq_system(coeff)
    solution = FVM.get_solution(equation_system);
end

function iteratevolumes(volumes::Array{Int, 1}, times::Int)
    for volume ∈ volumes
        mesh = benchmark_mesh(volume)
        coefficients = benchmark_set_boundary_conditions(mesh)
        functions = [benchmark_mesh, benchmark_set_boundary_conditions, benchmark_solutions]
        args = [volume, mesh, coefficients]
        
        for (f, arg) ∈ zip(functions, args)
            println("Comencé el de $volume volúmenes con la función $f")
            list_of_time_statistics = list_of_statistics(f, arg, volume, times)
            writetofile(f, list_of_time_statistics')
        end
        println("Terminé el de $volume volúmenes")
    end
end

function list_of_statistics(f, arg, volume::Int, times::Int)
    time_list, μ, σ = getstatistics(f, arg, times)
    time_statistics = [volume, volume^3, μ, σ]
    push!(time_statistics, time_list...)
    return time_statistics
end

function getstatistics(f, arg, times)
    time_list = []
    for epoch ∈ 1:(times+1)
        t = measuretime(f, arg)
        push!(time_list, t)
    end
    times_without_compiling = time_list[2:end]
    μ = mean(times_without_compiling)
    σ = std(times_without_compiling)
    
    return time_list, μ, σ
end

function measuretime(f, arg)
    start_time = time()
    f(arg)
    finish_time = time()
    t = finish_time - start_time
    return t
end

function writetofile(f, list)
    file_name = string(f)
    CSV.write("../benchmarking/$(file_name).csv", Tables.table(list), delim = ',',append=true)
end

volumes = [10,20,30,40,50,60,80,100]
times = 10
iteratevolumes(volumes, times)


