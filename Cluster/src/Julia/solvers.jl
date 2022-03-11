module Solvers
    using LinearAlgebra, SparseArrays, BenchmarkTools
    using CUDA, CUDA.CUSPARSE

    function jacobi(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
        # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= D⁻¹(L+U)x⁽ᵏ⁻¹⁾+D⁻¹b
        # Con D := La matriz diagonal A como matriz cuadrada,
        # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
        # Definimos a N como la suma de estas dos últimas matrices.
        d = diag(A)
        N = - triu(A,1) - tril(A,-1)
        x = rand(length(b))
        counter = 0
        while norm(A*x - b) > ϵ
            x = (N*x + b) ./ d
            counter+=1
        end
        return counter
    end

    function paralleljacobi(A, b, ϵ = 1e-5)
        d = diag(A)
        d_d = CuArray{Float64}(d)
        N = - triu(A,1) - tril(A,-1)
        d_N = CuSparseMatrixCSR{Float64}(N)

        d_A = CuSparseMatrixCSR{Float64}(A)
        d_b = CuArray{Float64}(b)
        d_x = CUDA.rand(Float64, size(b))
        counter = 0
        while norm(d_A*d_x - d_b) > ϵ
            d_x = (d_N*d_x + d_b) ./ d_d
            counter+=1
        end
        return counter
    end


    function jacobiparallelsolver(A, b, x, d, N, ϵ = 1e-5)
        while norm(A*x - b) > ϵ
            x = (N*x + b) ./ d
        end
        return x
    end


    function gaussseidel(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
        # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
        # Con D := La matriz diagonal A como matriz cuadrada,
        # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
        U = - triu(A,1) 
        L₀ = tril(A,0)
        x = rand(length(b))
        counter = 0
        while norm(A*x - b) > ϵ
            x = L₀\(U*x + b)
            counter+=1
        end
        return counter
    end

    function parallelgaussseidel(A, b, ϵ = 1e-5)
        # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
        # Con D := La matriz diagonal A como matriz cuadrada,
        # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
        U = - triu(A,1)
        L₀ = tril(A,0)
        d_U = CuSparseMatrixCSC{Float64}(U)

        d_A = CuSparseMatrixCSC{Float64}(A)
        d_b = CuArray{Float64}(b)
        d_x = CUDA.rand(Float64, size(b))

        counter = 0
        while norm(d_A*d_x - d_b) > ϵ
            new_b = Array(d_U*d_x + d_b)
            d_x = CuArray{Float64}(L₀\new_b)
            counter+=1
        end
        return counter
    end


    function sor(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ω = 1.79::Float64, ϵ = 1e-5)
        # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-ωL)⁻¹(ωU+(1-ω)D)x⁽ᵏ⁻¹⁾+ω(D-ωL)⁻¹b
        # Con D := La matriz diagonal A como matriz cuadrada,
        # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
        # ω es justo un argumento de sobrerelajación, que sirve para hacer más rápido al método. Hay que escogerlo 
        # entre (0,1)
        D = Diagonal(A)
        U = - triu(A,1) 
        L = - tril(A,-1)
        x = rand(length(b))

        L′ = D-ω*L
        U′ = ω*U+(1.0-ω)*D
        b′ = ω*b
        counter = 0
        while norm(A*x - b) > ϵ
            x = L′ \ (U′*x + b′)
            counter+=1
        end
        return counter
    end

    function parallelsor(A, b, ω = 1.79, ϵ = 1e-5)
        # Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
        # Con D := La matriz diagonal A como matriz cuadrada,
        # -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
        D = Diagonal(A)
        U = - triu(A,1)
        L = - tril(A,-1)

        L′ = D-ω*L
        U′ = ω*U+(1.0-ω)*D
        b′ = ω*b

        d_U′ = CuSparseMatrixCSC{Float64}(U′)
        d_b′ = CuArray{Float64}(b′)

        d_A = CuSparseMatrixCSC{Float64}(A)
        d_b = CuArray{Float64}(b)
        d_x = CUDA.rand(Float64, size(b))

        counter = 0
        while norm(d_A*d_x - d_b) > ϵ
            new_b = Array(d_U′*d_x + d_b′)
            d_x = CuArray{Float64}(L′\new_b)
            counter+=1
        end
        return counter
    end


    function leastsquares(H, r)
        r′ = zeros(size(H)[1])
        r′[1] = norm(r)
        x = H \ r′
    end

    # Algoritmo sacado directo del pseudocódigo
    function gmres(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
        x₀ = rand(length(b))
        residual₀ = b - A*x₀
        q = [residual₀ / norm(residual₀)]

        k = 1
        x = x₀
        H = zeros(2,1)
        residual = residual₀
        counter = 0
        while norm(residual) > ϵ
            y = A*q[k]
            for j ∈ 1:k
                H[j,k] = q[j]' * y
                y -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(y)
            push!(q, y/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            if k % 10 == 0
                c = leastsquares(H, residual₀)
                Q = hcat(q...)
                x = Q*c + x₀
                residual = A*x - b
            end
            k += 1
            counter+=1
        end
        return counter
    end

    # Algoritmo sacado directo del pseudocódigo
    function parallelgmres(A, b, ϵ = 1e-5)
        d_b = CuArray{Float64}(b)
        d_A = CuSparseMatrixCSR{Float64}(A)
        d_x₀ = CUDA.rand(Float64, size(b))

        d_residual₀ = d_b - d_A*d_x₀
        q = [residual₀ / norm(residual₀)]

        k = 1
        d_x = d_x₀
        counter = 0
        H = zeros(2,1)
        d_residual = d_residual₀
        while norm(d_residual) > ϵ
            d_y = d_A*q[k]
            for j ∈ 1:k
                H[j,k] = q[j]' * d_y
                d_y -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(d_y)
            push!(q, d_y/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            if k % 10 == 0
                d_H = CuArray{Float64}(H)
                d_c = leastsquares(d_H, d_residual₀)
                Q = hcat(q...)
                #c = CuArray{Float64}(c)
                d_Q = CuSparseMatrixCSR{Float64}(Q)
                x = d_Q*d_c + d_x₀
                d_residual = d_A*d_x - d_b
            end
            k += 1
            counter+=1
        end
        return counter
    end

    function reiniciarvariables(x::Vector, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64})
        k = 1
        x₀ = x
        r = b - A*x₀
        H = zeros(2,1)
        q = [r / norm(r)]
        return x₀, r, q, k, H
    end


    function obtenermodulo(dims)
        if dims == 10
            modulo = 15
        elseif dims == 20
            modulo = 20
        elseif dims == 30
            modulo = 35
        elseif dims == 40
            modulo = 25
        elseif dims == 50
            modulo = 15
        elseif dims == 60
            modulo = 15
        elseif dims == 80
            modulo = 30
        elseif dims == 100
            modulo = 35
        end
        return modulo
    end


    function gmresreiniciado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
        x₀ = rand(length(b))
        residual₀ = b - A*x₀
        q = [residual₀ / norm(residual₀)]

        k = 1
        x = x₀
        H = zeros(2,1)
        residual = residual₀
        counter = 0
        modulo = obtenermodulo(cbrt(length(b)))
        while norm(residual) > ϵ
            y = A*q[k]
            for j ∈ 1:k
                H[j,k] = q[j]' * y
                y -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(y)
            push!(q, y/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            k += 1
            if k % modulo == 0 
                c = leastsquares(H, residual₀)
                Q = hcat(q...)
                x = Q*c + x₀
                residual = A*x - b
                x₀, residual₀, q, k, H = reiniciarvariables(x, A, b)
            end
            counter+=1
        end
        return counter
    end


    function reiniciarvariables(x, A, b)
        k = 1
        x₀ = x
        r = b - A*x₀
        H = zeros(2,1)
        q = [r / norm(r)]
        return x₀, r, q, k, H
    end


    function parallelgmresreiniciado(A, b, ϵ = 1e-5)
        d_b = CuArray{Float64}(b)
        d_A = CuSparseMatrixCSR{Float64}(A)
        d_x₀ = CUDA.rand(Float64, size(b))

        d_residual₀ = d_b - d_A*d_x₀
        q = [d_residual₀ / norm(d_residual₀)]

        k = 1
        d_x = d_x₀
        counter = 0
        H = zeros(2,1)
        d_residual = d_residual₀
        modulo = obtenermodulo(cbrt(length(b)))
        while norm(d_residual) > ϵ
            d_y = d_A*q[k]
            for j ∈ 1:k
                H[j,k] = q[j]' * d_y
                d_y -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(d_y)
            push!(q, d_y/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            k += 1
            if k % modulo == 0
                d_H = CuArray{Float64}(H)
                d_c = leastsquares(d_H, d_residual₀)
                Q = hcat(q...)
                #c = CuArray{Float64}(c)
                d_Q = CuSparseMatrixCSR{Float64}(Q)
                d_x = d_Q*d_c + d_x₀
                residual = d_A*d_x - d_b
                d_x₀, d_residual₀, q, k, d_H = reiniciarvariables(d_x, d_A, d_b)
            end
            counter+=1
        end
        return counter
    end

    function precondition(name, A, b)
        if name == "Jacobi"
            M_jacobi = Diagonal(A)
            return [M_jacobi]

        elseif name == "Gauss-Seidel"
            D = Diagonal(A)
            U = triu(A,1) 
            L = tril(A,-1)
            M_gauss_seidel_1 = I+(L*inv(D))
            M_gauss_seidel_2 = D+U
            return [M_gauss_seidel_1, M_gauss_seidel_2]

        elseif name == "SOR"
            α = 1.8
            D = Diagonal(A)
            U = triu(A,1) 
            L = tril(A,-1)
            M_sor_1 = I+(α*L*inv(D))
            M_sor_2 = D+α*U
            return [M_sor_1, M_sor_2]
        end
    end


    function gmresprecondicionado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, precondition_name::String, ϵ = 1e-5)
        Ms = precondition(precondition_name, A, b)

        x₀ = rand(length(b))
        residual₀ = (b - A*x₀)
        [residual₀ = M\residual₀ for M ∈ Ms]
        q = [residual₀ / norm(residual₀)]

        k = 1
        x = x₀
        H = zeros(2,1)
        residual = residual₀
        counter = 0
        while norm(residual) > ϵ
            ω = A*q[k]
            [ω = M\ω for M ∈ Ms]
            for j ∈ 1:k
                H[j,k] = q[j]' * ω
                ω -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(ω)
            push!(q, ω/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            if k % 10 == 0
                c = leastsquares(H, residual₀)
                Q = hcat(q...)
                x = Q*c + x₀
                residual = A*x - b
                [residual = M\residual for M ∈ Ms]
            end
            k += 1
            counter+=1
        end
        return counter
    end


    function parallelgmresprecondicionado(A, b, precondition_name::String, ϵ = 1e-5)
        Ms = precondition(precondition_name)
        d_Ms = [CuArray{Float64}(M) for M ∈ Ms]

        d_b = CuArray{Float64}(b)
        d_A = CuSparseMatrixCSR{Float64}(A)
        d_x₀ = CUDA.rand(Float64, size(b))

        d_residual₀ = (d_b - d_A*d_x₀)
        [d_residual₀ = d_M\d_residual₀ for d_M ∈ d_Ms]
        q = [d_residual₀ / norm(d_residual₀)]

        k = 1
        d_x = d_x₀
        counter = 0
        H = zeros(2,1)
        d_residual = d_residual₀
        while norm(d_residual) > ϵ
            d_ω = d_A*q[k]
            [d_ω = d_M\d_ω for d_M ∈ d_Ms]
            for j ∈ 1:k
                H[j,k] = q[j]' * d_ω
                d_ω -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(d_ω)
            push!(q, d_ω/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            if k % 10 == 0
                d_H = CuArray{Float64}(H)
                d_c = leastsquares(d_H, d_residual₀)
                Q = hcat(q...)
                #c = CuArray{Float64}(c)
                d_Q = CuSparseMatrixCSR{Float64}(Q)
                d_x = d_Q*d_c + d_x₀
                d_residual = d_A*d_x - d_b
                [d_residual = d_M\d_residual for d_M ∈ d_Ms]
            end
            k += 1
            counter+=1
        end
        return counter
    end

    function obtenermoduloprecondicionado(precondition_name, dims)
        if precondition_name == "Jacobi"
            if dims == 10
                modulo = 15
            elseif dims == 20
                modulo = 20
            elseif dims == 30
                modulo = 35
            elseif dims == 40
                modulo = 25
            elseif dims == 50
                modulo = 15
            elseif dims == 60
                modulo = 15
            elseif dims == 80
                modulo = 30
            elseif dims == 100
                modulo = 35
            end

        elseif precondition_name == "Gauss-Seidel"
            if dims == 10
                modulo = 20
            elseif dims == 20
                modulo = 30
            elseif dims == 30
                modulo = 25
            elseif dims == 40
                modulo = 20
            elseif dims == 50
                modulo = 20
            elseif dims == 60
                modulo = 25
            elseif dims == 80
                modulo = 30
            elseif dims == 100
                modulo = 25
            end		
        elseif precondition_name == "SOR"
            if dims == 10
                modulo = 20
            elseif dims == 20
                modulo = 15
            elseif dims == 30
                modulo = 25
            elseif dims == 40
                modulo = 30
            elseif dims == 50
                modulo = 15
            elseif dims == 60
                modulo = 15
            elseif dims == 80
                modulo = 30
            elseif dims == 100
                modulo = 20
            end
        end
        return modulo
    end

    function reiniciarvariablesprecondicionado(x::Vector, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, Ms::Vector)
        k = 1
        x₀ = x
        H = zeros(2,1) 
        residual₀ = b - A*x₀
        [residual₀ = M\residual₀ for M ∈ Ms]
        q = [residual₀ / norm(residual₀)]
        return x₀, residual₀, q, k, H
    end

    function gmresprecondicionadoreiniciado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, precondition_name::String, ϵ = 1e-5)
        x₀ = rand(length(b))
        Ms = precondition(precondition_name, A, b)
        residual₀ = b - A*x₀
        [residual₀ = M\residual₀ for M ∈ Ms]
        q = [residual₀ / norm(residual₀)]

        k = 1
        x = x₀
        H = zeros(2,1)
        residual = residual₀
        counter = 0
        modulo = obtenermoduloprecondicionado(precondition_name, cbrt(length(b)))
        while norm(residual) > ϵ
            ω = A*q[k]
            [ω = M\ω for M ∈ Ms]
            for j ∈ 1:k
                H[j,k] = q[j]' * ω
                ω -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(ω)
            push!(q, ω/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            k += 1
            if k % modulo == 0
                c = leastsquares(H, residual₀)
                Q = hcat(q...)
                x = Q*c + x₀
                residual = A*x - b
                [residual = M\residual for M ∈ Ms]
                x₀, residual₀, q, k, H = reiniciarvariablesprecondicionado(x, A, b, Ms)
            end
            counter+=1
        end
        return counter
    end

    function reiniciarvariablesprecondicionado(x, A, b, Ms::Vector)
        k = 1
        x₀ = x
        H = zeros(2,1)

        residual₀ = b - A*x₀
        [residual₀ = M\residual₀ for M ∈ Ms]
        q = [residual₀ / norm(residual₀)]
        return x₀, residual₀, q, k, H
    end

    function parallelgmresprecondicionadoreiniciado(A, b, precondition_name::String, ϵ = 1e-5)
        d_b = CuArray{Float64}(b)
        d_A = CuSparseMatrixCSR{Float64}(A)
        d_x₀ = CUDA.rand(Float64, size(b))

        Ms = precondition(precondition_name)
        d_Ms = [CuArray{Float64}(M) for M ∈ Ms]
        d_residual₀ = d_b - d_A*d_x₀
        [d_residual₀ = d_M\d_residual₀ for d_M ∈ d_Ms]
        q = [d_residual₀ / norm(d_residual₀)]

        k = 1
        d_x = d_x₀
        counter = 0
        H = zeros(2,1)
        d_residual = d_residual₀
        modulo = obtenermoduloprecondicionado(precondition_name, cbrt(length(b)))
        while norm(d_residual) > ϵ
            d_ω = d_A*q[k]
            [d_ω = d_M\d_ω for d_M ∈ d_Ms]
            for j ∈ 1:k
                H[j,k] = q[j]' * d_ω
                d_ω -= H[j,k]*q[j]
            end
            H[k+1,k] = norm(d_ω)
            push!(q, d_ω/H[k+1,k])
            H = vcat(H, zeros(1, size(H)[2]))
            H = hcat(H, zeros(size(H)[1], 1))
            k += 1
            if k % modulo == 0
                d_H = CuArray{Float64}(H)
                d_c = leastsquares(d_H, d_residual₀)
                Q = hcat(q...)
                #c = CuArray{Float64}(c)
                d_Q = CuSparseMatrixCSR{Float64}(Q)
                d_x = d_Q*d_c + d_x₀
                d_residual = d_A*d_x - d_b
                [d_residual = d_M\d_residual for d_M ∈ d_Ms]
                d_x₀, d_residual₀, q, k, d_H = reiniciarvariablesprecondicionado(d_x, d_A, d_b, d_Ms)
            end
            counter+=1
        end
        return counter
    end
end
