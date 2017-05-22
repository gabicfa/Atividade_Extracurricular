function [R, V] = kep_to_car(elem_orb)

global mu
global deg
mu = 398600;
deg = pi/180;

h = elem_orb(1); %momento angular
e = elem_orb(2); %ecentricidade
RA = elem_orb(3)*deg; %ascens?o direita do n? ascendente(rad)
incl = elem_orb(4)*deg; %inclinacao da orbita(rad)
w = elem_orb(5)*deg; %argumento do perigeu (rad)
TA = elem_orb(6)*deg; %anomalia verdadeira(rad)

%rp - position vector in the perifocal frame (km)
%vp - velocity vector in the perifocal frame (km/s)
rp = ((h^2/mu) * (1/(1 + e*cos(TA)))) *((cos(TA)*[1;0;0] + sin(TA)*[0;1;0]));
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

R3_W = [[ cos(RA)  sin(RA)  0]
        [ -sin(RA)  cos(RA) 0]
        [    0        0     1]];

R1_i = [[1       0          0]
        [0   cos(incl)  sin(incl)]
        [0  -sin(incl)  cos(incl)]];
    
R3_w = [[ cos(w)  sin(w)  0]
        [-sin(w)  cos(w)  0]
        [   0       0     1]];
    
Q_pX = R3_W'*R1_i'*R3_w';

r = Q_pX*rp;
v = Q_pX*vp;

R = r';
V = v';

        