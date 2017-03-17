%Atividade 3 - Iniciação cientifica

function saida = atividade2(x0,y0,z0,vx0,vy0,vz0,tf);
    close all;
    u = 3.98600440*(10^5);     
    Rt = 6378;
    options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6); 

    [t, X] = ode45( @simula, [0 tf], [x0 y0 z0 vx0 vy0 vz0], options);
    [lon,lat,r] = cart2sph(X(:,1),X(:,2),X(:,3));
    
    [t, Xc] = ode45( @simula_corrigido, [0 tf], [x0 y0 z0 vx0 vy0 vz0], options);
    [lonc,latc,rc] = cart2sph(Xc(:,1),Xc(:,2),Xc(:,3));
    
    %Figura para a posiçao X(t) X Y(t) X Z(t)
    
    figure(1);
    plot3( X(:,1),X(:,2),X(:,3), 'r');
    grid on;
    
    title('Órbita - integração numérica');
    xlabel('Posição em X [km]');
    ylabel('Posição em Y [km]');
    zlabel('Posição em Z [km]');
    
    figure(2);
    geoshow('landareas.shp', 'FaceColor', [0.8 1 0.8]);
    hold on;
    plot(lon*180/pi, lat*180/pi, '.k');
    axis equal;
    
    figure(3);
    plot3( Xc(:,1),Xc(:,2),Xc(:,3), 'r');
    grid on;
    
    title('Órbita corrigida - integração numérica');
    xlabel('Posição em X [km]');
    ylabel('Posição em Y [km]');
    zlabel('Posição em Z [km]');    
    
    figure(4);
    geoshow('landareas.shp', 'FaceColor', [0.8 1 0.8]);
    hold on;
    plot(lonc*180/pi, latc*180/pi, '.k');
    axis equal;

end

function X_ponto = simula(t,X);

    X_ponto = zeros(6,1);
    
    u = 3.9860044*10^5;     
    Rt = 6378;      
    
    X_ponto(1) = X(4);
    X_ponto(2) = X(5);
    X_ponto(3) = X(6);
    X_ponto(4) = (-u*X(1))/((Rt + sqrt(X(1)^2 + X(2)^2 + X(3)^2))^3);
    X_ponto(5) = (-u*X(2))/((Rt + sqrt(X(1)^2 + X(2)^2 + X(3)^2))^3);
    X_ponto(6) = (-u*X(3))/((Rt + sqrt(X(1)^2 + X(2)^2 + X(3)^2))^3);
    
end

function X_ponto = simula_corrigido(t,X);

    X_ponto = zeros(6,1);
    
    u = 3.9860044*10^5;     
    Rt = 6378;      
    w = 7.29*10^-5;
    
    r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);
    
    X_ponto(1) = X(4) + w*X(2);
    X_ponto(2) = X(5) - w*X(1);
    X_ponto(3) = X(6);
    
    X_ponto(4) = (-u*X(1))/((Rt + r)^3);
    X_ponto(5) = (-u*X(2))/((Rt + r)^3);
    X_ponto(6) = (-u*X(3))/((Rt + r)^3);
    
end
