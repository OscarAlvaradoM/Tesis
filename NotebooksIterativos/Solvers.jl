using LinearAlgebra, SparseArrays, BenchmarkTools

function vectorizedjacobi(A, b, ϵ = 1e-5)
    # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= D⁻¹(L+U)x⁽ᵏ⁻¹⁾+D⁻¹b
    # Con D := La matriz diagonal A como matriz cuadrada,
    # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
    # Definimos a N como la suma de estas dos últimas matrices.
    dA = diag(A)
    N = - triu(A,1) - tril(A,-1)
    x = 0*b;
    normres = [];
    while norm(A*x - b) > ϵ
        x = (N*x + b)./dA
    end
    return x
end

function paralleljacobi(A, b, ϵ = 1e-5)
    dA = diag(A);
    dA = CuArray{Float64}(dA)
    N = - triu(A,1) - tril(A,-1)
    N = CuSparseMatrixCSR{Float64}(N)
    
    A = CuSparseMatrixCSR{Float64}(A)
    b = CuArray{Float64}(b)
    x = CUDA.zeros(Float64, size(b))
    while norm(A*x - b) > ϵ
        x = (N*x + b)./dA
    end
    return x
end

function vectorizedgaussseidel(A, b, ϵ = 1e-5)
    # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
    # Con D := La matriz diagonal A como matriz cuadrada,
    # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
    U = - triu(A,1) 
    L = - tril(A,-1)
    L₀ = tril(A,0)
    x = 0*b
    normres = [];
    while norm(A*x - b) > ϵ
        x = L₀\(U*x + b)
    end
    return x
end

function parallelgaussseidel(A, b, ϵ = 1e-5)
    # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
    # Con D := La matriz diagonal A como matriz cuadrada,
    # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
    U = - triu(A,1) 
    L₀ = tril(A,0)
    U = CuSparseMatrixCSR{Float64}(U)
    L₀ = CuSparseMatrixCSR{Float64}(L₀)

    A = CuSparseMatrixCSR{Float64}(A)
    b = CuArray{Float64}(b)
    x = CUDA.zeros(Float64, size(b))
    while norm(A*x - b) > ϵ
        x = L₀\(U*x + b)
    end
    return x
end

function vectorizedsor(A, b, ω = 0.49, ϵ = 1e-5)
    # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-ωL)⁻¹(ωU+(1-ω)D)x⁽ᵏ⁻¹⁾+ω(D-ωL)⁻¹b
    # Con D := La matriz diagonal A como matriz cuadrada,
    # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
    # ω es justo un argumento de sobrerelajación, que sirve para hacer más rápido al método. Hay que escogerlo 
    # ente (0,1)
    D = Diagonal(A)
    U = - triu(A,1) 
    L = - tril(A,-1)
    x = 0*b;
    counter = 0
    while norm(A*x - b) > ϵ
        x = (D-ωL)\((ω*U+(1-ω)*D)*x + ω*b)
    counter+=1
    end
    println("Terminé en $counter iteraciones")
    return x
end

function parallelsor(A, b, ω, ϵ = 1e-5)
    # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
    # Con D := La matriz diagonal A como matriz cuadrada,
    # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
    D = Diagonal(A)
    D = CuArray{Float64}(D)
    U = - triu(A,1) 
    L = - tril(A,-1)
    U = CuSparseMatrixCSR{Float64}(U)
    L = CuSparseMatrixCSR{Float64}(L)

    A = CuSparseMatrixCSR{Float64}(A)
    b = CuArray{Float64}(b)
    x = CUDA.zeros(Float64, size(b))
    while norm(A*x - b) > ϵ
        x = (D-ωL)\((ω*U+(1-ω)*D)*x + ω*b)
    end
    return x
end

# # Algoritmo sacado directo del pseudocódigo
# function conjugategradient(A, b, ϵ = 1e-5)
#     x = rand(length(b))
#     r = b - A*x
#     ρ = [0, r'*r]
#     res = b - A*x
#     p = Array{Float64}(undef, size(r,1), 1)
    
#     counter = 0
#     while  norm(res) > ϵ
#         # Aquí iría la condición más chida de precondicionamiento
#         ρ[1] = ρ[2] # Utilizamos sólo dos actualizaciones hacia atrás, por eso esto sí es una lista.
#         ρ[2] = r'*r # Aquí cambia cuando es precondicionada
#         if counter == 0
#             p = r
#         else
#             β = ρ[2]/ρ[1]
#             p = r + β*p
#         end
#         q = A*p
#         α = ρ[2]/(p'*q)
#         x = x + α*p
#         r = r - α*q
#         res = b - A*x
#         counter+=1
#     end
#     println("Se necesitaron $counter iteraciones para obtener la solución con ϵ = $ϵ")
#     return x
# end

# # Algoritmo sacado directo del pseudocódigo
# function parallelconjugategradient(A, b, ϵ = 1e-5)
#     A = CuSparseMatrixCSR(A)
#     b = CuArray{Float32}(b)
#     x = rand(length(b))
#     x = CuArray{Float32}(x)
    
#     r = b - A*x
#     d = r'*r
#     ρ = [0.0f0, d]
#     res = b - A*x
#     p = CuArray{Float32}(undef, size(r,1), 1)
    
#     counter = 0
#     while  norm(res) > ϵ
#         # Aquí iría la condición más chida de precondicionamiento
#         ρ[1] = ρ[2] # Utilizamos sólo dos actualizaciones hacia atrás, por eso esto sí es una lista.
#         ρ[2] = r'*r # Aquí cambia cuando es precondicionada
#         if counter == 0
#             p = r
#         else
#             β = ρ[2]/ρ[1]
#             p = r + β*p
#         end
#         q = A*p
#         α = ρ[2]/(p'*q)
#         x = x + α*p
#         r = r - α*q
#         res = b - A*x
#         counter+=1
#     end
#     println("Se necesitaron $counter iteraciones para obtener la solución con ϵ = $ϵ")
#     return x
# end
