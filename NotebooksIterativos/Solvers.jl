module Solvers
	using LinearAlgebra, SparseArrays, BenchmarkTools

	function jacobi(A, b, ϵ = 1e-5)
		n = size(A)[1]
		x = b
		x_bar = zeros(n)

		counter = 0
		while norm(A*x - b) > ϵ
		    for i ∈ 1:n
		        x_bar[i] = 0
		        for j ∈ 1:n
		            if j != i
		                x_bar[i] = x_bar[i] + A[i,j]*x[j]
		            end
		        end
		        x_bar[i] = (b[i] - x_bar[i])/A[i,i]
		    end
		    x = x_bar
		    counter+=1
		end
		return counter
	end

	function vectorizedjacobi(A, b, ϵ = 1e-5)
		D = diag(A)
		N = - triu(A,1) - tril(A,-1)
		x = 0*b;
		normres = [];
		counter = 0
		while norm(A*x - b) > ϵ
		    println(norm(A*x - b))
                    x = (N*x + b)./D
                    counter+=1
                    println(counter)
		end
		return counter
	end

	function paralleljacobi(A, b, ϵ = 1e-5)
		dA = diag(A);
		dA = CuArray(dA)
		N = - triu(A,1) - tril(A,-1)
		N = CuSparseMatrixCSR(N)
		
		A = CuSparseMatrixCSR(A)
		b = CuArray(b)
		x = CUDA.zeros(Float32, size(b))
		counter = 0
		while norm(A*x - b) > ϵ
		    x = (N*x + b)./dA
		    counter+=1
		end
		return counter
	end

	# function gauss_seidel(A, b, ϵ = 1e-5)
	#     n = size(A)[1]
	#     x = b
		
	#     counter = 0
	#     while norm(A*x - b) > ϵ
	#         x_next = zeros(n)
	#         for i = 1:n
	#             σ = 0
	#             for j = 1:i-1
	#                 if j <= i
	#                     σ = σ + A[i,j]*x_next[j]
	#                 end
	#             end
	#             for j = i+1:n
	#                 σ = σ + A[i,j]*x[j]
	#             end
	#             x_next[i] = (b[i] - σ)/A[i,i]
	#         end
	#         x = x_next
	#         counter+=1
	#     end
	#     println("Se necesitaron $counter iteraciones para obtener la solución con ϵ = $ϵ")
	#     return x
	# end 

	# A, b = get_A_b(25_000)
	# t = @benchmark gauss_seidel(A,b,1e-11)

	# Algoritmo sacado directo del pseudocódigo
	function conjugategradient(A, b, ϵ = 1e-5)
		#x = rand(length(b))
		x = b
		r = b - A*x
		ρ = [0, r'*r]
		res = b - A*x
		p = Array{Float64}(undef, size(r,1), 1)
		
		counter = 0
		while norm(res) > ϵ
		    println(norm(res))
		    # Aquí iría la condición más chida de precondicionamiento
		    ρ[1] = ρ[2] # Utilizamos sólo dos actualizaciones hacia atrás, por eso esto sí es una lista.
		    ρ[2] = r'*r # Aquí cambia cuando es precondicionada
		    if counter == 0
		        p = r
		    else
		        β = ρ[2]/ρ[1]
		        p = r + β*p
		    end
		    q = A*p
		    α = ρ[2]/(p'*q)
		    x = x + α*p
		    r = r - α*q
		    res = b - A*x
		    counter+=1
		    println(counter)
		end
		return counter
	end

	# Algoritmo sacado directo del pseudocódigo
	function parallelconjugategradient(A, b, ϵ = 1e-5)
		A = CuSparseMatrixCSR(A)
		b = CuArray{Float32}(b)
		#x = rand(length(b))
		#x = CuArray{Float32}(x)
		x = b
		
		r = b - A*x
		d = r'*r
		ρ = [0.0f0, d]
		res = b - A*x
		p = CuArray{Float32}(undef, size(r,1), 1)
		
		counter = 0
		while  norm(res) > ϵ
		    # Aquí iría la condición más chida de precondicionamiento
		    ρ[1] = ρ[2] # Utilizamos sólo dos actualizaciones hacia atrás, por eso esto sí es una lista.
		    ρ[2] = r'*r # Aquí cambia cuando es precondicionada
		    if counter == 0
		        p = r
		    else
		        β = ρ[2]/ρ[1]
		        p = r + β*p
		    end
		    q = A*p
		    α = ρ[2]/(p'*q)
		    x = x + α*p
		    r = r - α*q
		    res = b - A*x
		    counter+=1
		end
		return counter
	end
end
