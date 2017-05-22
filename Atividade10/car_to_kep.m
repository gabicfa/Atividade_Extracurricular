function elem_orb = car_to_kep(R,V)

global mu
global deg
deg = pi/180;
mu = 398600;
eps = 1.e-10;

r= norm(R);
v= norm(V);

vr = dot(R,V)/r;

H = cross(R,V);
h = norm(H);

incl = acos(H(3)/h); % incl = i;

N =cross([0 0 1],H);
n =norm(N);

if n ~= 0
    RA = acos(N(1)/n); % RA = omega
    if N(2) < 0
        RA = 2*pi-RA;
    end
else
    RA = 0;
end
E = 1/mu*((v^2 - mu/r)*R -r*vr*V);
e = norm(E); % e= e

if n ~= 0
    if e > eps
        w = acos(dot(N,E)/n/e); % w=w
        if E(3)< 0
            w = 2*pi -w;
        end
    else
        w=0;
    end
else
    w = 0;
end

if e > eps
    TA = acos(dot(E,R)/e/r); % TA = theta
    if vr < 0
        TA = 2*pi -TA;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0
        TA = acos(dot(N,R)/n/r);
    else
        TA = 2*pi - acos(dot(N,R)/n/r);
    end
end

% a = h^2/mu/(1-e^2);

elem_orb = [h e RA/deg incl/deg w/deg TA/deg];