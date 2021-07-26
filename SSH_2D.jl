using LinearAlgebra
using PyPlot
##################
#     |   |
#  --A == B--
#    || / ||
#  --C == D--
#    |   |
##################
# open boundary condition
N = 10 # the number of cells
#hopping terms 
w = 1 # inter e.g. A -> outer
v = 1/5 # intra e.g. A -> B
u = 0.6 # B -> C
t = 0 # A -> D

function Hamiltonian(N::Int, w::Number, v::Number, u::Number, t::Number)
    n = 4*N   
    hopping1 = zeros(ComplexF64,n^2-1)
    for i = 1:n^2-1
        if i % n == 0
            hopping1[i] = 0
        elseif i % n % 2 == 0
            hopping1[i] = w
        else
            hopping1[i] = v
        end
    end

    hopping2 = zeros(ComplexF64,n^2-n)
    for i = 1:n-1
        if i % 2 == 0
            hopping2[1+(i-1)*n:i*n] .= w
        else
            hopping2[1+(i-1)*n:i*n] .= v
        end
    end

    hopping3 = zeros(ComplexF64,n^2-n+1)
    for i = 1:n^2-n+1
        if i % 2 == 0
            hopping3[i] = u
        end
    end
    for i = 1:(n-1)
        if i % 2 == 0
            hopping3[i*n-n+1:i*n] .= 0
        end 
    end

    hopping4 = zeros(ComplexF64,n^2-n-1)
    for i in 1:n^2-n-1
        if i % 2 != 0
            hopping4[i] = t
        end
    end
    for i in 1:(n-1)
        if i % 2 == 0
            hopping4[i*n-n+1:i*n] .= 0
        end
    end
    H = diagm(1 => hopping1, -1=> hopping1, n => hopping2, -n => hopping2, n-1 => hopping3, -n+1 => hopping3, n-3 => hopping4, -n+3 => hopping4)
end
H = Hamiltonian(N, w, v, u, t)
E,V = eigen(H)
# bandgap
figure()
plot(real(E),".")
grid()
# field
V1 = V[:,4878]
V1 = reshape(V1,4*N,4*N)
figure();imshow(abs.(V1));colorbar()

# One direction is the periodic boundary condition
function Hamiltonian(N::Int, ky::Number, w::Number, v::Number)
    n = 4*N
    hopping1 = zeros(n-1)
    for i = 1:n-1
        if i % (n/2) == 0
            hopping1[i] = 0
        elseif i % n % 2 == 0
            hopping1[i] = w
        else
            hopping1[i] = v
        end
    end
    hopping2 = ones(ComplexF64,Int(n/2))*(w+v*exp(im*ky))
    H = diagm(1 => hopping1, -1=> hopping1, Int(n/2) => hopping2, -Int(n/2) => conj(hopping2)) # Hamiltonian
end
N = 40
w = 3 #inter
v = 1 #intra
kym = range(-pi,pi,length=101)
E = zeros(ComplexF64,length(kym),4*N)
V = zeros(ComplexF64,length(kym),4*N,4*N)
i = 1
for ky in kym
    H = Hamiltonian(N, ky, w, v)
    E[i,:] .= eigvals(H)    
    V[i,:,:] = eigvecs(H)
    i = i+1
end
ind_k = 2
ind_E = N
V1 = V[ind_k,:,ind_E]
V1 = reshape(V1,Int(2*N),2)
figure();imshow(abs.(transpose(V1)))
Ei = round(real(E[ind_E]),digits=2); ki = round(kym[ind_k]/pi,digits=2)
title("E=$Ei, k=$ki π")
figure()
plot(kym,E,"k",linewidth=0.5)
xticks([-pi,0,pi],[L"-π",L"0",L"π"])
grid()
xlim([-pi,pi])
xlabel(L"k_y")
ylabel(L"E")

# periodic boundary condition
w = 3; v = 1
Γ = [0;0]
X = [pi;0]
M = [pi;pi]
kx0 = [range(0,pi,length=10);range(1,1,length=10)*pi;range(pi,0,length=10)]
ky0 = [range(0,0,length=10);range(0,pi,length=10);range(pi,0,length=10)]
k = hcat(kx0,ky0)
E = zeros(ComplexF64,length(kx0),4)
i = 1
for (kx, ky) in zip(kx0,ky0)
    H = [0 w+v*exp(-im*kx) w+v*exp(-im*ky) 0;
        w+v*exp(im*kx) 0 0 w+v*exp(-im*ky);
        w+v*exp(im*ky) 0 0 w+v*exp(-im*kx);
        0 w+v*exp(im*ky) w+v*exp(im*kx) 0] 
    E[i,:] .= eigvals(H)    
    i = i+1
end
figure();plot(real.(E))
