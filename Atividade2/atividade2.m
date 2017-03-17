%Atividade 2 - Iniciação cientifica

function saida = atividade2(x1_0,x2_0,x3_0,x4_0,tf);
    close all;
    u = 3.98600440*(10^5);     
    Rt = 6378;

    [ t, X] = ode45( @simula, [0 tf], [x1_0 x2_0 x3_0 x4_0]);
    
    %Calculo da altitude para cada instante
    
    for i =1:length(X)
        X(i,5) = sqrt(X(i,1)^2 + X(i,2)^2);
    end
    
    %Figura para a posiçao X(t) X Y(t)
    
    figure(1);
    plot( X(:,1),X(:,2), 'b');
    grid on;
    
    title('Órbita "numérica"');
    xlabel('Posição em X [km]');
    ylabel('Posição em Y [km]');
    
    %Figura para altitude H(t)
    figure(2);
    plot( t, X(:,5), 'r');
    grid on;
    
    title('Altitude X Tempo');
    xlabel('Tempo [s]');
    ylabel('Altitude [km]');
    
    %Valores das altitudes maximas e minimas
    rmax = max(X(:,5));
    rmin = min(X(:,5));
    
    %Localizando o indice onde se encontram os valores máximos e minimos
    %de altitude
    [lin_max, col_max] = find(X == rmax);
    [lin_min, col_min] = find(X == rmin);
    
    %Cálculo dos angulos a partir dos valores respectivos aos indices
    %achados
    teta_max = atan2(X(lin_max,2),X(lin_max,1));
    teta_min = atan2(X(lin_min,2),X(lin_min,1));
    %Somando o raio da Terra aos valores máximos e mínimos
    r1 = Rt + rmin;
    r2 = Rt + rmax;
    
    %Calculo da Excentricidade (e)
    e = (r2 - r1)/(r1*cos(teta_min) - r2*cos(teta_max));
    disp(abs(e));
    %Calculo do momento angular (h)
    h = sqrt(r1*u*(1 + e*cos(teta_min)));
    
    %Calculo do perigeu (pe) e do apogeu (ap)
    rp = (h^2)/(u*(1 + e));
    ra = (h^2)/(u*(1 - e));
    pe = rp - Rt;
    ap = ra - Rt;
    
    %Calculo do semi-eixo maior (a)
    a = (ap + pe)/2;
    
    %Calculo da velocidadeangular média (n)
    n = sqrt(u/(a^3));

    Xa = [];
    Ya = [];
    Xa_ponto = [];
    Ya_ponto = [];
    
    for E = 0:0.01:2*pi + 0.1
        x = a*(cos(E) - e);
        y = a*(sqrt(1 - e^2))*sin(E);
        r = (h^2)/(u*(1 + e*cos(E)));
        xa_ponto = -n*(a^2)*sin(E)/r;
        ya_ponto = n*(a^2)*cos(E)/r;
        
        Xa(length(Xa) + 1) = x;
        Ya(length(Ya) + 1) = y;
    end
    
    figure(3);
    plot(Xa,Ya, 'g');
    grid on;
    hold on;
    plot( X(:,1),X(:,2), 'r');
    axis equal;
    title('Orbita Numerica X Orbita Analitica');
    xlabel('Coordenada em x [km]');
    ylabel('Coordenada em Y [km]');
    lgd = legend('Kepler', 'Integrado','Location','northeast');
    title(lgd,'Métodos para calculo');
end

function X_ponto = simula(t,X);

    X_ponto = zeros(4,1);
    
    u = 3.9860044*10^5;     
    Rt = 6378;      
    
    X_ponto(1) = X(3);
    X_ponto(2) = X(4);
    X_ponto(3) = (-u*X(1))/((Rt + sqrt(X(1)^2 + X(2)^2))^3);
    X_ponto(4) = (-u*X(2))/((Rt + sqrt(X(1)^2 + X(2)^2))^3);
    
end


