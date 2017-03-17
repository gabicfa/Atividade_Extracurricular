function atividade5()

    ut = 3.9860040*10^5;
    hora = 2;
    min = 15;
    seg = 10;
    altitude = 3000;
    Rt = 6378;

    a = Rt + altitude;
    e = 0.1;
    
    i = 30;
    omega = 45;
    w = 60;
    
    i = i*pi/180;
    omega = omega*pi/180;
    w = w*pi/180;

    T = (hora + (min/60) + (seg/3600))*3600;
    n = sqrt(ut/a^3);
    M = n*T;

    E0 = 0;
    iter = 0;
    itermax = 100;

    Er = E0;
    es = 0.05;

    while(1)
        Erold = Er;
        x1 = Erold - e*sin(Erold)-M;
        x2 = 1 - e*cos(Erold);
        Er = Erold - (x1/x2);
        iter = iter + 1;
        if Er ~= 0
            ea = abs((Er - Erold)/Er)*100;
        end
        if ea < es || iter >= itermax
            E = Er;
            break;
        end
    end
    
    r = a*(1 - e*cos(E));
    X = a*(cos(E) - e);
    Y = a*sqrt(1 - e^2)*sin(E);
    Z = 0;
    Vx = -n*(a^2)*sin(E)/r;
    Vy = n*(a^2)*sqrt(1 - e^2)*cos(E)/r;
    Vz = 0;

    Rz_omega = [cos(-omega) sin(-omega) 0;-sin(-omega) cos(-omega) 0;0 0 1];
    Rx_i = [1 0 0;0 cos(-i) sin(-i);0 -sin(-i) cos(-i)];
    Rz_w = [cos(-w) sin(-w) 0;-sin(-w) cos(-w) 0;0 0 1];

    Rot = Rz_omega*Rx_i;
    Rot = Rot*Rz_w;

    Pos = [X;Y;Z];
    V = [Vx;Vy;Vz];

    Pos_0 = Rot*Pos;
    V_0 = Rot*V;

    t0 = 0;
    tf = 86200;
    
    Xl = 3.714209622133673 * (10^5);
    Yl = 1.069961860819459 * (10^5);
    Zl = 0.066447602321645 * (10^5);
    
    Vxl = -0.298536192639590;
    Vyl = -0.891257374756870;
    Vzl = 0.092657657414370;
    
    options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6);
    
    
    Iniciais = [Pos_0(1) Pos_0(2) Pos_0(3) V_0(1) V_0(2) V_0(3) Xl Yl Zl Vxl Vyl Vzl];
    
    [~,valores_saida] = ode45(@odefun,[t0 tf],Iniciais,options);
    [lon,lat,~] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));
    
    figure(1);
    plot3(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3),'r');
    title('Movimento do corpo');
    xlabel('Posicao em x [km]');
    ylabel('Posicao em y [km]');
    zlabel('Posicao em z [km]');
    grid on;
    axis([-1 1 -1 1 -0.5 0.5]*10^4);
    
    figure(2);
    geoshow('landareas.shp', 'FaceColor', [0.8 1 0.8]);
    hold on;
    plot(lon*180/pi, lat*180/pi, '.k');
    axis equal;
    
    
end


function [G] = odefun(~,I)

%     I[1] = xs;
%     I[2] = ys;
%     I[3] = zs;

%     I[4] = vxs;
%     I[5] = vys;
%     I[6] = vzs;

%     I[7] = xl;
%     I[8] = yl;
%     I[9] = zl;

%     I[10] = vxl;
%     I[11] = vyl;
%     I[12] = vzl;

    
    G = zeros(12,1);
    Rt = 6378;
    w = 7.29*10^-5;
    ut = 3.9860040*(10^5);
    ul = 4.9028*(10^3);
    
    rs = sqrt(I(1)^2 + I(2)^2 + I(3)^2);
    rl = sqrt(I(7)^2 + I(8)^2 + I(9)^2);
    rls = sqrt((I(7) - I(1))^2 + (I(8) - I(2))^2 + (I(9) - I(3))^2);
    
    J2 = 0.001082;
    J3 = -0.0000025323;
    
    factor = (1 + (J2*(Rt/rs)^2)*((3/2)*(1 - (5*(I(3)^2)/(rs^2)))) + (J3*(Rt/rs)^3)*((5/2)*(3 - (7*(I(3)^2)/(rs^2))))*(I(3)/rs));
    factor_z = (ut/(rs^2))*J3*((Rt/rs)^3)*(3/2);
  
    G(1) = I(4) + w*I(2);
    G(2) = I(5) - w*I(1);
    G(3) = I(6);
    
    G(7) = I(10) + w*I(8);
    G(8) = I(11) - w*I(7); 
    G(9) = I(12);
    
    G(4) = (-ut*I(1)/(rs^3))*factor + ul*((I(7) - I(1))/(rls^3) - I(7)/(rl^3));
    G(5) = (-ut*I(2)/(rs^3))*factor + ul*((I(8) - I(2))/(rls^3) - I(8)/(rl^3));
    G(6) = (-ut*I(3)/(rs^3))*factor +factor_z + ul*((I(9) - I(3))/(rls^3) - I(9)/(rl^3));
    
    G(10) = -(ut + ul)*I(7)/(rl^3);
    G(11) = -(ut + ul)*I(8)/(rl^3);
    G(12) = -(ut + ul)*I(9)/(rl^3);
end