function atividade7()

    ut = 3.9860040*10^5;
    Rt = 6378.16;
    a = 1.5*Rt;
    e = 0.1; 
    
    i = 30;
    omega = 45;
    w = 60;
    
    i = i*pi/180;
    omega = omega*pi/180;
    w = w*pi/180;
    n = sqrt(ut/a^3);
    
    T= [16 1 5];
    D= [22 6 1962];
    
    t=[2 15 0];
    d=[23 6 1962];
    
    Delta_JD = juliano(T, D, t, d);
    display(Delta_JD)
    display(n)
    M = n*Delta_JD;
    display(M)

    E0 = 0;
    iter = 0;
    itermax = 100;

    Er = E0;
    es = 0.05;
    
    E = newton(Er,M,e,iter, itermax,es);
    
    r = a*(1 - e*cos(E));
    X = a*(cos(E) - e);
    Y = a*sqrt(1 - e^2)*sin(E);
    Z = 0;
    Vx = -n*(a^2)*sin(E)/r;
    Vy = n*(a^2)*sqrt(1 - e^2)*cos(E)/r;
    Vz = 0;
    
    Pos = [X;Y;Z]
    V = [Vx;Vy;Vz]
    disp(omega);
    disp(w);
    disp(i);
    Rot= rotation(omega,i,w)

    Pos_0 = Rot*Pos;
    V_0 = Rot*V;
    
    wt = 7.29*10^-5;
    Pos_0(1) = Pos_0(1)+wt*V_0(2);
    Pos_0(2) = Pos_0(2)-wt*V_0(1);
    
    display(Pos_0);
    display(V_0);
    
    t0 = 0;
    tf = 86200;
    
    options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6); 
    
    [~,valores_saida] = ode45( @odefun,[t0 tf],[Pos_0(1) Pos_0(2) Pos_0(3) V_0(1) V_0(2) V_0(3)],options);
    [lon,lat,~] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));
    
%     figure(1);
%     plot3(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3),'r');
%     title('Movimento do corpo');
%     xlabel('Posicao em x [km]');
%     ylabel('Posicao em y [km]');
%     zlabel('Posicao em z [km]');
%     grid on;
%     axis([-2 2 -2 2 -1 1]*10^4);
%     
%     figure(2);
%     geoshow('landareas.shp', 'FaceColor', [0.8 1 0.8]);
%     hold on;
%     plot(lon*180/pi, lat*180/pi, '.k');
%     axis equal;
    
    
end

function[Delta_JD] = juliano(T, D, t, d)
        UT1=T(1)+(T(2)/60)+(T(3)/3600);
        ut2=t(1)+(t(2)/60)+(t(3)/3600);

        mc1 = floor(7*(D(3)+floor((D(2)+9)/12)/4));
        Jo1 = 367*D(3)-mc1+floor(257*D(2)/9)+D(1)+1721013.5;

        mc2 = floor(7*(d(3)+floor((d(2)+9)/12)/4));
        Jo2 = 367*d(3)-mc2+floor(257*d(2)/9)+d(1)+1721013.5;

        JD1 = Jo1 +(UT1/24);
        JD2 = Jo2 +(ut2/24);

        Delta_JD = (JD2 - JD1)*86400;
    end
function [E] = newton(Er,M,e,iter, itermax,es)
    while(1)
        Erold = Er;
        x1 = Erold - e*sin(Erold)-M;
        x2 = 1 - e*cos(Erold);
        Er = Erold - (x1/x2)
        iter = iter + 1;
        if Er ~= 0
            ea = abs((Er - Erold)/Er)*100;
        end
        if ea < es || iter >= itermax
            E = Er
            break;
        end
    end
end
function Rot = rotation(omega,i,w)
    Rz_omega = [cos(-omega) sin(-omega) 0;-sin(-omega) cos(-omega) 0;0 0 1];
    Rx_i = [1 0 0;0 cos(-i) sin(-i);0 -sin(-i) cos(-i)];
    Rz_w = [cos(-w) sin(-w) 0;-sin(-w) cos(-w) 0;0 0 1];

    Rot = Rz_omega*Rx_i;
    Rot = Rot*Rz_w;
end
    
function [G] = odefun(~,I)
    
    G = zeros(6,1);
    
    w = 7.29*10^-5;
    ut = 3.9860040*(10^5);
    r = sqrt(I(1)^2 + I(2)^2 + I(3)^2);
    
    G(1) = I(4);
    G(2) = I(5);
    G(3) = I(6);
    
    G(4) = -ut*I(1)/(r^3)+(w^2)*I(1)+2*w*I(5);
    G(5) = -ut*I(2)/(r^3)+(w^2)*I(2)-2*w*I(4);
    G(6) = -ut*I(3)/(r^3);

end
