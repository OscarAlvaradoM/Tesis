module Solvers
	using LinearAlgebra, SparseArrays, BenchmarkTools
	using CUDA, CUDA.CUSPARSE

	function vectorizedjacobi(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
		# Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= D⁻¹(L+U)x⁽ᵏ⁻¹⁾+D⁻¹b
		# Con D := La matriz diagonal A como matriz cuadrada,
		# -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
		# Definimos a N como la suma de estas dos últimas matrices.
		D = Diagonal(A)
		N = - triu(A,1) - tril(A,-1)
		x = rand(length(b))
		while norm(A*x - b) > ϵ
		    x = D\(N*x + b)
		end
		return Vector(x)
	end

	function paralleljacobi(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
		D = diag(A);
		D = CuArray{Float64}(D)
		N = - triu(A,1) - tril(A,-1)
		N = CuSparseMatrixCSR{Float64}(N)
		
		A = CuSparseMatrixCSR{Float64}(A)
		b = CuArray{Float64}(b)
		x = CUDA.rand{Float64}(length(b))
		while norm(A*x - b) > ϵ
		    x = D\(N*x + b)
		end
		return x
	end

	function vectorizedgaussseidel(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
		# Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
		# Con D := La matriz diagonal A como matriz cuadrada,
		# -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
		U = - triu(A,1) 
		L₀ = tril(A,0)
		x = rand(length(b))
		while norm(A*x - b) > ϵ
		    x = L₀\(U*x + b)
		end
		return x
	end

	function parallelgaussseidel(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
		# Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-L)⁻¹(Ux⁽ᵏ⁻¹⁾+b)
		# Con D := La matriz diagonal A como matriz cuadrada,
		# -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
		U = - triu(A,1) 
		L₀ = tril(A,0)
		U = CuSparseMatrixCSR{Float64}(U)
		L₀ = CuSparseMatrixCSR{Float64}(L₀)

		A = CuSparseMatrixCSR{Float64}(A)
		b = CuArray{Float64}(b)
		x = CUDA.rand{Float64}(length(b))
		while norm(A*x - b) > ϵ
		    x = L₀\(U*x + b)
		end
		return x
	end

	function vectorizedsor(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ω = 1.79::Float64, ϵ = 1e-5)
		# Según el libro, viene siendo así el método matricial: x⁽ᵏ⁾= (D-ωL)⁻¹(ωU+(1-ω)D)x⁽ᵏ⁻¹⁾+ω(D-ωL)⁻¹b
		# Con D := La matriz diagonal A como matriz cuadrada,
		# -L y -U las matrices estríctamente Inferior (Lower) y Superior (Upper) de A, respectivamente
		# ω es justo un argumento de sobrerelajación, que sirve para hacer más rápido al método. Hay que escogerlo 
		# ente (0,1)
		D = Diagonal(A)
		U = - triu(A,1) 
		L = - tril(A,-1)
		x = rand(length(b))
		
		L′ = D-ω*L
		U′ = ω*U+(1.0-ω)*D
		b′ = ω*b
		while norm(A*x - b) > ϵ
		    x = L′ \ (U′*x + b′)
		end
		return x
	end

	function parallelsor(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ω = 1.79::Float64, ϵ = 1e-5)
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
		x = CUDA.rand{Float64}(length(b))
		
		L′ = D-ω*L
		U′ = ω*U+(1.0-ω)*D
		b′ = ω*b
		while norm(A*x - b) > ϵ
		    x = L′ \ (U′*x + b′)
		end
		return x
	end

	function leastsquares(H::Matrix{Float64}, r)
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
		end
		return x
	end

	# Algoritmo sacado directo del pseudocódigo
	function parallelgmres(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
		b = CuArray{Float64}(b)
		A = CuSparseMatrixCSR{Float64}(A)
		x₀ = CUDA.rand{Float64}(length(b))
		
		residual₀ = b - A*x₀
		q = [residual₀ / norm(residual₀)]
		
		k = 1
		x = x₀
		H = zeros(2,1)
		residual = residual₀
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
		        c = CuArray{Float64}(c)
		        Q = CuSparseMatrixCSR{Float64}(Q)
		        x = Q*c + x₀
		        residual = A*x - b
		    end
		    k += 1
		end
		return x
	end

	function reiniciarvariables(x::Vector, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64})
		k = 1
		x₀ = x
		r = b - A*x₀
		H = zeros(2,1)
		q = [r / norm(r)]
		return x₀, r, q, k, H
	end
		
	function gmresreiniciado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
		x₀ = rand(length(b))
		residual₀ = b - A*x₀
		q = [residual₀ / norm(residual₀)]
		
		k = 1
		x = x₀
		H = zeros(2,1)
		residual = residual₀
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
		    #if k % 25 == 0 # cuando tenemos 30³ nuestra matriz ( Q = 5,400,000)
		    if k % 25 == 0 # cuando tenemos 40³ nuestra matriz ( Q = 10,240,000)
		        c = leastsquares(H, residual₀)
		        Q = hcat(q...)
		        x = Q*c + x₀
		        residual = A*x - b
		        x₀, residual₀, q, k, H = reiniciarvariables(x, A, b)
		    end
		end
		return x
	end

	function reiniciarvariables(x, A, b)
		k = 1
		x₀ = x
		r = b - A*x₀
		H = zeros(2,1)
		q = [r / norm(r)]
		return x₀, r, q, k, H
	end
		
	function parallelgmresreiniciado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, ϵ = 1e-5)
		b = CuArray{Float64}(b)
		A = CuSparseMatrixCSR{Float64}(A)
		x₀ = CUDA.rand{Float64}(length(b))
		
		residual₀ = b - A*x₀
		q = [residual₀ / norm(residual₀)]
		
		k = 1
		x = x₀
		H = zeros(2,1)
		residual = residual₀
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
		    #if k % 25 == 0 # cuando tenemos 30³ nuestra matriz ( Q = 5,400,000)
		    if k % 25 == 0 # cuando tenemos 40³ nuestra matriz ( Q = 10,240,000)
		        c = leastsquares(H, residual₀)
		        Q = hcat(q...)
		        c = CuArray{Float64}(c)
		        Q = CuSparseMatrixCSR{Float64}(Q)
		        x = Q*c + x₀
		        residual = A*x - b
		        x₀, residual₀, q, k, H = reiniciarvariables(x, A, b)
		    end
		end
		return x
	end

	function precondition(name)
		if name == "Jacobi"
			M_jacobi = Diagonal(A)
			return [M_jacobi]
		
		elif name == "Gauss-Seidel"
			D = Diagonal(A)
			U = triu(A,1) 
			L = tril(A,-1)
			M_gauss_seidel_1 = I+(L*inv(D))
			M_gauss_seidel_2 = D+U
			return [M_gauss_seidel_1, M_gauss_seidel_2]
		
		elif name == "SOR"
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
		Ms = precondition(precondition_name)
		
		x₀ = rand(length(b))
		residual₀ = (b - A*x₀)
		[residual₀ = M\residual₀ for M ∈ Ms]
		q = [residual₀ / norm(residual₀)]
		
		k = 1
		x = x₀
		H = zeros(2,1)
		residual = residual₀
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
		end
		return x
	end

	function parallelgmresprecondicionado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, precondition_name::String, ϵ = 1e-5)
		Ms = precondition(precondition_name)
		
		b = CuArray{Float64}(b)
		A = CuSparseMatrixCSR{Float64}(A)
		x₀ = CUDA.rand{Float64}(length(b))
		
		residual₀ = (b - A*x₀)
		[residual₀ = M\residual₀ for M ∈ Ms]
		q = [residual₀ / norm(residual₀)]
		
		k = 1
		x = x₀
		H = zeros(2,1)
		residual = residual₀
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
		        c = CuArray{Float64}(c)
		        Q = CuSparseMatrixCSR{Float64}(Q)
		        x = Q*c + x₀
		        residual = A*x - b
		        [residual = M\residual for M ∈ Ms]
		    end
		    k += 1
		end
		return x
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

	# Algoritmo sacado directo del pseudocódigo
	function gmresprecondicionadoreiniciado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, precondition_name::String, ϵ = 1e-5)
		x₀ = rand(length(b))
		Ms = precondition(precondition_name)
		residual₀ = b - A*x₀
		[residual₀ = M\residual₀ for M ∈ Ms]
		q = [residual₀ / norm(residual₀)]
		
		k = 1
		x = x₀
		H = zeros(2,1)
		residual = residual₀
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
		    # if k % 20 == 0 # Para 40³ sirvió chido el módulo 20 para el de Jacobi
		    # if k % 25 == 0 # Para 40³ sirvió chido el módulo 25 para el de Gauss-Seidel
		    if k % 30 == 0 # Para 40³ sirvió chido el módulo 30 para el de SOR
		        c = leastsquares(H, residual₀)
		        Q = hcat(q...)
		        c = CuArray{Float64}(c)
		        Q = CuSparseMatrixCSR{Float64}(Q)
		        x = Q*c + x₀
		        residual = A*x - b
		        [residual = M\residual for M ∈ Ms]
		        x₀, residual₀, q, k, H = reiniciarvariablesprecondicionado(x, A, b, Ms)
		    end
		end
		return x
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

	# Algoritmo sacado directo del pseudocódigo
	function parallelgmresprecondicionadoreiniciado(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, precondition_name::String, ϵ = 1e-5)
		x₀ = rand(length(b))
		Ms = precondition(precondition_name)
		residual₀ = b - A*x₀
		[residual₀ = M\residual₀ for M ∈ Ms]
		q = [residual₀ / norm(residual₀)]
		
		k = 1
		x = x₀
		H = zeros(2,1)
		residual = residual₀
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
		    # if k % 20 == 0 # Para 40³ sirvió chido el módulo 20 para el de Jacobi
		    # if k % 25 == 0 # Para 40³ sirvió chido el módulo 25 para el de Gauss-Seidel
		    if k % 30 == 0 # Para 40³ sirvió chido el módulo 30 para el de SOR
		        c = leastsquares(H, residual₀)
		        Q = hcat(q...)
		        x = Q*c + x₀
		        residual = A*x - b
		        [residual = M\residual for M ∈ Ms]
		        x₀, residual₀, q, k, H = reiniciarvariablesprecondicionado(x, A, b, Ms)
		    end
		end
		return x
	end
end