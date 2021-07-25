using LinearAlgebra
using PyPlot
# open boundary condition
Nx = 20
n = 4*Nx
w = 3 #inter
v = 1 #intra
hopping1 = zeros(n^2-1)
for i = 1:n^2-1
    if i % n == 0
        hopping1[i] = 0
    elseif i % n % 2 == 0
        hopping1[i] = w
    else
        hopping1[i] = v
    end
end

hopping2 = zeros(n^2-n)
for i = 1:n-1
    if i % 2 == 0
        hopping2[1+(i-1)*n:i*n] .= w
    else
        hopping2[1+(i-1)*n:i*n] .= v
    end
end
# Hamiltonian
H = diagm(1 => hopping1, -1=> hopping1, n => hopping2, -n => hopping2)
E,V = eigen(H)
# bandgap
figure()
plot(E,".")
grid()
# field
V1 = V[:,1600]
V1 = reshape(V1,n,n)
figure();imshow(V1)

# One direction is the periodic boundary condition
Nx = 40
n = 4*Nx
w = 3 #inter
v = 1 #intra
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
kxm = range(-pi,pi,length=101)
E = zeros(ComplexF64,length(kxm),n)
V = zeros(ComplexF64,length(kxm),n,n)
i = 1
for kx in kxm
    hopping2 = ones(ComplexF64,Int(n/2))*(w+v*exp(im*kx))
    H = diagm(1 => hopping1, -1=> hopping1, Int(n/2) => hopping2, -Int(n/2) => conj(hopping2)) # Hamiltonian
    E[i,:] .= eigvals(H)    
    V[i,:,:] = eigvecs(H)
    i = i+1
end
ind_k = 2
ind_E = Nx
V1 = V[ind_k,:,ind_E]
V1 = reshape(V1,Int(n/2),2)
figure();imshow(abs.(transpose(V1)))
Ei = round(real(E[ind_E]),digits=2); ki = round(kxm[ind_k]/pi,digits=2)
title("E=$Ei, k=$ki π")
figure()
plot(kxm,E,"k",linewidth=0.5)
xticks([-pi,0,pi],[L"-π",L"0",L"π"])
grid()
xlim([-pi,pi])
xlabel(L"k_y")
ylabel(L"E")

# periodic boundary condition
w = 3; v = 1
kx0 = [range(0,pi,length=10);range(1,1,length=10)*pi;range(pi,0,length=10)]
ky0 = [range(0,0,length=10);range(0,pi,length=10);range(pi,0,length=10)]
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
