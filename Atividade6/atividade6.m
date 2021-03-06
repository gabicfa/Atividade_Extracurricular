 function atividade6(R, V)


r = norm(R);
v = norm(V);
vr= sum(V.*R)/r;

H = cross(R,V);
h = norm(H);
disp('h(km^2/s) = ') 
disp(H);

Hz = H(3);
i = acos(Hz/h);
i = radtodeg(i);
disp('i(graus) = ')
disp(i);

N = cross([0 0 1], H);
n = norm(N);
Nx = N(1,1);
Ny = N(1,2);
if(Ny>= 0)
    omega = acos(Nx/n);
else
    omega = 2*pi - acos(Nx/n);
end
omega = radtodeg(omega);
disp('W(graus) = ')
disp(omega)

u = 3.9860040*(10^5);
E = (1/u) * ((v^2-(u/r))*R- (r*vr*V));
e = norm(E);
disp('e = ')
disp(e);

Ez = E(3);
if Ez>=0
    disp('Entrei aqui no primeiro')
%     w = acos((sum(N.*E))/(n*e));
    w = acos((dot(N, E))/(n*e));
else
    disp('Entrei aqui no segundo')
    w = 2*pi - acos((sum(N.*E))/(n*e));
end
w = radtodeg(w);
disp('w(graus) = ')
disp(w);

if vr >=0
    O = acos((sum(E.*R))/(e*r));
else
    O = 2*pi - acos((sum(E.*R))/(e*r));
end
O = radtodeg(O);
disp('O(graus) = ')
disp(O);

rp = (h^2/u)*(1/(1+e*cos(0)));
disp('rp(km) = ')
disp(rp);

ra = (h^2/u)*(1/(1+e*cos(pi)));
disp('ra(km) = ')
disp(ra);

a = (1/2)*(rp+ra);
disp('a(km) = ')
disp(a);

T= (((2*pi)/(sqrt(u)))* (a^(3/2)))/3600;
disp('T(hr) = ')
disp(T);
end