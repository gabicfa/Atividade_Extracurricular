function [rrel_xyz, vrel_xyz, arel_xyz] = relative_motion(ra,va,rb,vb)

global mu

aA = -mu*(ra/(norm(ra)^3));
aB = -mu*(rb/(norm(rb)^3));

i = (ra/norm(ra));

hA = cross(ra,va);

k = (hA/norm(hA));

j = cross(k,i);

display(norm(ra)^2);
W = hA/((norm(ra)^2));

W_ponto = (-(2*(dot(ra,va)))/(norm(ra)^2))*W;

rrel = rb-ra;

vrel = vb-va-(cross(W,rrel));

arel = aB-aA-(cross(W_ponto,rrel))-(cross(W,cross(W,rrel)))-2*(cross(W,vrel));

Q = [[i(1) i(2) i(3)]
    [ j(1) j(2) j(3)]
    [ k(1) k(2) k(3)]];

rrel_xyz = Q*rrel';
vrel_xyz = Q*vrel';
arel_xyz = Q*arel';


