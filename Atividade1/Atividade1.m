%Atividade 1

function saida = atividade1( x1_0, x2_0, tf);
    
    function b = yval(tf)
        [ts,Xs] = ode45( @simula,[0,tf],[x1_0 x2_0]);         
         b = Xs(end,2);                         
        
    end

    [ t, X] = ode45( @simula, [0 tf], [x1_0 x2_0]);
     
    valor_final_x = X(end,1);
    disp(valor_final_x);
    
    x0 = fzero(@yval, 1);
    disp(x0);
    %Gráficos com subplot
    
    figure(1);
    s(1) = subplot( 1, 2, 1);
    plot( t, X( :,1), 'r');
    
    
    s(2) = subplot( 1, 2, 2);
    plot ( t, X( :,2), 'b');
    
    %Nomeando gráficos e seus eixos
    
    title(s(1),'Posiçao X Tempo');
    xlabel(s(1), 'Tempo[s]');
    ylabel(s(1), 'Posição [km]');
    
    title(s(2),'Velocidade X Tempo');
    xlabel(s(2), 'Tempo[s]');
    ylabel(s(2), 'Velocidade [km/s]');
    
  
    %Gráficos separados 
    
    figure(2);
    grid on;
    hold on;
    
    [ t, X] = ode45( @simula, [0 210], [10 1]);
    plot( t, X(:,1), 'r');
    [ t, X_2] = ode45( @simula, [0 210],[10 1.1]);
    plot( t, X_2(:,1), 'g');
    [ t, X_3] = ode45( @simula, [0 210],[10 1.4]);
    plot( t, X_3(:,1), 'b');
    
    title('Posições X Tempo');
    xlabel('Tempo [s]');
    ylabel('Posição [km]');
    legend('vxo_1 = 1km/s','vxo_2 = 1.1 km/s','vxo_3 = 1.4 km/s');
    
    figure(3);
    grid on;
    hold on;
    
    plot( t, X(:,2), 'r');
    plot( t, X_2(:,2), 'g');
    plot( t, X_3(:,2), 'b');
    
    title('Velocidades X Tempo');
    xlabel('Tempo [s]');
    ylabel('Velocidade [km/s]');
    legend('vxo_1 = 1km/s','vxo_2 = 1.1 km/s','vxo_3 = 1.4 km/s');
   
end



function X_ponto = simula(t,X);

    X_ponto = zeros(2,1);
    g = 9.81*0.001;  %aceleração da gravidade, sendo igual a x'', em km/s^2
    
    X_ponto(1) = X(2);
    X_ponto(2) = -g;
    
end



