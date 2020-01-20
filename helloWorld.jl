a = 2 + 2
println(a)

b = a ÷ 3

α = 0.5
Vf(u) = α*u
Vf(2)
sin(2π)

@code_llvm 2*5

2^5

@code_native 2^-5 # should give error, but does not

function test1()
    a = zeros(3)
    @inbounds for i = 1:4
        println(a[i])
        a[i] = i
        println(a[i])
    end
end


test1()

b = Vector{Int}(undef, 3)
b[1] = 1
b[2] = 2
b[3] = 3
typeof(b[2])

c = Vector{Any}(undef, 3)
c[1] = 123
c[2] = "hello"
c[3] = "world"

println(c)

d = Vector{Union{Float64, Int}}(undef, 3)
d[1] = 1.0
d[2] = 2
d[3] = "banana"

typeof(d[1])
typeof(d[2])
typeof(d[3])

@code_warntype 2^5

function expo(x,y)
    if y>0
        return x^y
    else
        x = convert(Float64,x)
        return x^y
    end
end

@code_warntype expo(2,3)

using Traceur
@trace expo(3,4)

function bad_idea()
    aa = 3.0
    for i in 1:20
        aa += i
    end
end


@time bad_idea()

aaa = 3
function test2()
    a = 3.0
    @time for i = 1:20
        a += 1
    end
end
test2()

a = [1,2,3,4,5]
b =  map((x)->x^2, a)
println(b)

A = 1:5
B = [1 2
    3 4
    5 6
    7 8
    9 10]

C = broadcast(+, A, B)
println(C)

println(sin.(A))
D = similar(C)

x = [1,2,3,4,5]
println(map((x)-> x^2, x))

println(sin.(x))

function testLoops(n)
    b = rand(n,n)
    c = 0
    @time for i in 1:n, j in 1:n
        c += b[i,j]
    end

    @time for j in 1:n, i in 1:n
        c += b[i,j]
    end

    bidx = eachindex(b)
    for ij in bidx
        c += b[ij]
    end
end

testLoops(1000)

A = [1 2 3; 4 5 6; 7 8 9]
for i in 1:2, j in 1:3
    println(A[i,j])
end

println(A)
B = A[1:3, 1:3]
C = A
B[1,1] = 1000
println(B)
println(A)
C[1,1] = 1000
println(C)
println(A)

E = view(A, 1:3, 1:3)
E[1,1] = 1001
println(E)
println(A)

M1 = [1; 2; 3]
M2 = [2; 2; 3]
transpose(M1)*M2
typeof(transpose(M1))

M1+M2

A = rand(4,4)
B = rand(4,4)
b = [1; 2; 3; 4]
c = 1:4
# A x = c --> x = A\c
A\c
using LinearAlgebra
I4x4 = Matrix{Float64}(I,4,4)
println(A\I4x4)

A = rand(3,3)
B = rand(3,3)
C = rand(3,3)
i = 1
@time C[i,:] = A[i,:] .* B[i,:]
@time C[i,:] .= view(A,i,:) .* view(B,i,:)

using SparseArrays
A = sparse([1;2;3], [2;2;1], [3;4;5])
println(A)
println(Array(A))

using LinearAlgebra
A = Tridiagonal(2:5, 1:5, 1:4)
print(A)
fieldnames(typeof(A))
typeof(A)

A = rand(5,5)
λ = 2
A - λ*I

b = 1:5
A = rand(5,5)
# A = QR


Q,R = qr(A)
Q'*Q
# A^-1 = R^-1 QT
@time x1 = A\b
@time x2 = inv(R)*Q'*b
q = qr(A)
@time x3 = q\b
println(x1)
println(x2)
println(x3)

using Random
rng = MersenneTwister(1234)
randnumbers = randn(rng, ComplexF64, (2,3))
println(randnumbers)
