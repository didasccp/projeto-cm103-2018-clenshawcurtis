function quad_cc(f, n)
    x = zeros(1,n+1);
    c = zeros(1,n+1);
    w = zeros(1,n+1);
    b = zeros(1,div(n,2)+1);
    for i=0:n
        x[i+1] = cos((i)*pi/n);
        if mod(i,n)==0
            c[i+1] = 1;
        else
            c[i+1] = 2;
        end
    end
    for i=0:div(n,2)
        if i==n/2
            b[i+1]=1;
        else
            b[i+1]=2;
        end
    end
    for i=0:n
        soma=0
        for j=1:div(n,2)
            soma = soma + b[j+1]*cos(2*j*(i)*pi/n)/(4*j^2 -1);
        end
        w[i+1] = c[i+1]*(1-soma)/n
    end
    soma=0;
    for i=0:n
        soma = soma + f(x[i+1])*w[i+1];
    end
    return soma
end

function quad_gl(f, n)
    if n==2
        a = sqrt(3)/3
        return f(-a)+f(a);        
    elseif n==3
        a = sqrt(3/5)
        b = 5/9
        return f(-a)*b + f(0)*8/9 + f(a)*b        
    elseif n == 4
        u = 2*sqrt(6/5);
        a1 = sqrt((3-u)/7)
        a2 = sqrt((3+u)/7)
        b1 = (18+ sqrt(30))/36
        b2 = (18- sqrt(30))/36
        return f(-a1)*b1 + f(-a2)*b2 + f(a2)*b2 + f(a1)*b1;        
    elseif n== 5
        u = 2sqrt(10/7)
        a1 = sqrt(5-u)/3
        a2 = sqrt(5+u)/3
        c = 128/225
        b1 = (322+13sqrt(70))/900
        b2 = (322-13sqrt(70))/900
        return f(-a1)*b1 + f(a1)*b1 + c*f(0) + f(-a2)*b2 + f(a2)*b2
    elseif n==6
        a1 = 0.9324695142
        a2 = 0.6612093864
        a3 = 0.2386191860
        b1 = 0.1713244923
        b2 = 0.3607615730
        b3 = 0.4679139345
        return (f(-a1)+f(a1))*b1 +b2*(f(-a2)+f(a2))+ b3*(f(-a3)+f(a3));
    elseif n==7
        a1 = 0.9491079123
        a2 = 0.7415311855
        a3 = 0.4058451513
        b1 = 0.1294849661
        b2 = 0.2797053914
        b3 = 0.3818300505
        c = 0.4179591836
        return b1*(f(-a1)+f(a1)) + b2*(f(-a2)+f(a2)) + b3*(f(-a3)+f(a3)) + c*f(0)
    elseif n==8
        a1 = 0.9602898564
        a2 = 0.7966664774
        a3 = 0.5255324099
        a4 = 0.1834346424
        b1 = 0.1012285362
        b2 = 0.2223810344
        b3 = 0.3137066458
        b4 = 0.3626837833
        return b1*(f(-a1)+f(a1)) + b2*(f(-a2)+f(a2)) + b3*(f(-a3)+f(a3)) + b4*(f(-a4)+f(a4))
    end    
    return 0;
end

function simpsonrep(f, m)
    a=-1;
    b=1;
    if m%2==0
        m=m+1
    end
    h=(b-a)/(m-1)
    v=zeros(m,1)
    for i=1:m
        v[i]=a+h*(i-1)
    end
    soma1=0;
    soma2=0;
    for i=2:2:m-1
        soma1=soma1+f(v[i])
        soma2=soma2+f(v[i+1])
    end
    soma2=soma2-f(v[m])
    return h*(f(v[1])+4*soma1+2*soma2+f(v[m]))/3 
end

function trap_rep(f, n)
    h = 2/n
    soma = f(-1)+f(1)
    for i=1:n-1
        soma = soma + 2*f(-1 + i*h)
    end
    return h*soma/2
end

function f1(x)
    if x == 0
        return 0;
    else
        return sin(1/x)*3*x^2 - x*cos(1/x) + x^2
    end
end

f2(x) = 2x*exp(x^2) + 1

f3(x) = exp(x+4)*sin(x)

f4(x) = sqrt(x+1)

f5(x) = 1000x^5 + 50x^4 - 95(x^3) /3 - 23(x^2)/2 + 6x

f6(x) = 500x^5 - 3625(x^3)/3 + 1296x + 100

f7(x) = 5000x^4 + 200x^3 -95x^2 -23x + 6

f8(x) = log(x+2)

f9(x) = sin(x)^5 + 1

f10(x) = 1/sqrt(1.21 -x^2)

F = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10]

g1(x) = sin(1/x) * x^3 + (x^3)/3

g2(x) = exp(x^2) + x

g3(x) = exp(x+4)*(sin(x) - cos(x))/2

g4(x) = (2/3)*sqrt(x+1)^3

g5(x) = (1000/6)x^6 + 10x^5 -(95/12)x^4 - (23/6)x^3 + 3x^2

g6(x) = (500/6)x^6 - (3625/12)x^4 + 648x^2 + 100x

g7(x) = f5(x)

g8(x) = (x+2)*log(x+2) - x

g9(x) = -(cos(x)^5)/5 +(2cos(x)^3)/3 - cos(x) + x

g10(x) = asin(x/(1.1))

G = [g1 g2 g3 g4 g5 g6 g7 g8 g9 g10]

using Plots
Ec=zeros(10,7)
Eg=zeros(10,7)
for i=1:10
    f=F[i];
    g=G[i];
    for n=2:8
        A = quad_gl(f, n)
        B = quad_cc(f, n)
        C = g(1) - g(-1)
        Eg[i, n-1] = abs((A-C)/C)
        Ec[i, n-1] = abs((B-C)/C)
    end
end
N = zeros(7)
for i=2:8
    N[i-1] = i
end

EcM = zeros(1, 7)
EgM = zeros(1, 7)

for i=1:7
    EcM[1, i] = mean(Ec[:, i])
    EgM[1, i] = mean(Eg[:, i])
end

scatter(N, EgM[1,:], c=:red, lab=:"Gaussiana")
scatter!(N, EcM[1,:], yaxis = :log, c=:blue, lab=:"Clenshaw-Curtis")
title!("Gaussiana x Clenshaw-Curtis")
name = @sprintf("bom dia")
png(name)

u=100;

Ec2 = zeros(10, u)
Es = zeros(10, u)
Et = zeros(10, u)

for i=1:10
    f = F[i]
    g = G[i]
    for n=2:u+1
        A = trap_rep(f, n)
        B = quad_cc(f, n)
        D = simpsonrep(f, n)
        C = g(1) - g(-1)
        Ec2[i, n-1] = abs((B-C)/C)
        Es[i, n-1] = abs((D-C)/C)
        Et[i, n-1] = abs((A-C)/C)
    end
end

Ec2M = zeros(1, u)
EsM = zeros(1, u)
EtM = zeros(1, u)

for i=1:u
    Ec2M[1, i] = mean(Ec2[:,i])
    EsM[1, i] = mean(Es[:,i])
    EtM[1, i] = mean(Et[:,i])
end

M = zeros(u)
for i=1:u
    M[i] = i+1
end

scatter(M, Ec2M[1,:], c=:red, lab=:"Clenshaw-Curtis")
scatter!(M, EsM[1,:], yaxis = :log, c=:blue, lab=:"Simpson Rep.")
scatter!(M, EtM[1,:], yaxis = :log, c=:green, lab=:"Trapezio Rep.")
title!("Clenshaw-Curtis x Simpson e Trapezio")
name = @sprintf("bom dia 2")
png(name)