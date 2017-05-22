function varargout = IC(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IC_OpeningFcn, ...
                   'gui_OutputFcn',  @IC_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function varargout = IC_OutputFcn(~, eventdata, handles) 
varargout{1} = handles.output;

function IC_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for IC
handles.output = hObject;
set(handles.plot3D,'Xlim',[-2.5 2.5]*10^4,'Ylim',[-2.5 2.5]*10^4,'Zlim',[-1.5 1.5]*10^4);
set(handles.plotMap,'Xlim',[-180 180]);

%Plot da figura central!
Earth(handles);
geoshow(handles.plotMap,'landareas.shp','FaceColor',[0.8 1 0.8]);
axis(handles.plot3D,[-2.5 2.5 -2.5 2.5 -1.5 1.5]*10^4);
axis(handles.plotMap,[-180 180 -90 90]);

hold(handles.plotMap, 'on');
set(handles.plotMap, 'color', 'b','Xcolor','w','Ycolor','w','Zcolor','w');

hold(handles.plot3D, 'on');
grid(handles.plot3D,'on');

set(handles.plot3D,'GridLineStyle','--');
global cor_eixo
cor_eixo = [0.9 0.4 0.74];

set(handles.plot3D,'color','k','Xcolor','w','Ycolor','w','Zcolor','w');

guidata(hObject, handles);

global T
global t
global D
global d
global Dp
global Tp
global sim_check 
sim_check = 0;
d = [1 1 1900];
D = [1 1 1900];
T = [0 0 0];
t = [0 0 0];
Dp = [1 1 1900];
Tp = [0 0 0];

global inercial
inercial = 0;

global pertubacao
pertubacao =0;

global Xl Yl Zl Vxl Vyl Vzl
Xl = 3.714209622133673 * (10^5);
Yl = 1.069961860819459 * (10^5);
Zl = 0.066447602321645 * (10^5);

Vxl = -0.298536192639590;
Vyl = -0.891257374756870;
Vzl = 0.092657657414370;

% UIWAIT makes IC wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.edit_e,'enable','off');
set(handles.edit_i,'enable','off');
set(handles.edit_w,'enable','off');
set(handles.edit_omega,'enable','off');
set(handles.edit_theta,'enable','off');
set(handles.edit_rp,'enable','off');
set(handles.edit_ra,'enable','off');
set(handles.edit_a,'enable','off');
set(handles.simular_K,'enable','off','String','Simular');

function simular_K_Callback(hObject, eventdata, handles)
global e i w omega X0 Y0 Z0 Vx0 Vy0 Vz0 Dp D d Tp T t Th theta ra rp a sim_check Xl Yl Zl Vxl Vyl Vzl
sim_check = 1;

[Pos_0,V_0] = ValoresIniciais(e,i,omega,w,Tp,Dp,T,D);

Delta_JD = juliano(T, D, t, d);
Delta_JD_per = juliano(Tp, Dp, T, D);

if Delta_JD == 0
    errordlg('O instante inicial n?o pode ser id?ntico ao instante final.','Erro');
    return
end

% if Delta_JD_per == 0
%     errordlg('O instante de passagem pelo perigeu n?o pode ser id?ntico ao instante de refer?ncia.','Erro');
%     return
% end

if Delta_JD < 0 
    errordlg('O instante inicial n?o pode ser maior que o instante final.','Erro');
    return
end

if Delta_JD_per < 0 
    errordlg('O instante de passagem pelo perigeu n?o pode ser maior que o instante de refer?ncia.','Erro');
    return
end


%Aplica??o do ode45 para os valores iniciais calculados
global pertubacao
if pertubacao == 0
   options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6); 
    [tempo,valores_saida] = ode45( @odefun,[0 Delta_JD],[Pos_0(1) Pos_0(2) Pos_0(3) V_0(1) V_0(2) V_0(3)],options);
    [lon,lat,r] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3)); 
end
if pertubacao == 1
    options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6);
 
    Iniciais = [Pos_0(1) Pos_0(2) Pos_0(3) V_0(1) V_0(2) V_0(3) Xl Yl Zl Vxl Vyl Vzl];
    [~,valores_saida] = ode45(@odefun_pert,[0 Delta_JD],Iniciais,options);
    [lon,lat,~] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));
    display(valores_saida(size(valores_saida,1),:));

end



global p
delete(p);
p = plot3(handles.plot3D,valores_saida(:,1),valores_saida(:,2),valores_saida(:,3), 'r');

global p2
delete(p2);
p2 = plot(handles.plotMap,lon*180/pi,lat*180/pi,'.k');
axis(handles.plot3D,'equal');
rotate3d(handles.plot3D,'on');

%Seta os valores em suas respectivas caixas de texto
X0 = Pos_0(1);
Y0 = Pos_0(2);
Z0 = Pos_0(3);

Vx0 = V_0(1);
Vy0 = V_0(2);
Vz0 = V_0(3);

r = norm(Pos_0);
v = norm(V_0);
vr= dot(V_0,Pos_0)/r;

H = cross(Pos_0,V_0);
h = norm(H);

Hz = H(3);

N = cross([0 0 1], H);
n = norm(N);
Nx = N(1);
Ny = N(2);

u = 3.9860040*(10^5);
E = (1/u) * ((v^2 - (u/r))*Pos_0 - (r*vr*V_0));

if vr >= 0
    theta = acos(dot(E,Pos_0)/(e*r));
else
    theta = 2*pi - acos(dot(E,Pos_0)/(e*r));
end

theta = radtodeg(theta);

rp = (h^2/u)*(1/(1+e*cos(0)));


ra = (h^2/u)*(1/(1+e*cos(pi)));

a = (1/2)*(rp+ra);

Th = (((2*pi)/(sqrt(u)))* (a^(3/2)))/3600;

set(handles.edit_X0,'String',X0);
set(handles.edit_Y0,'String',Y0);
set(handles.edit_Z0,'String',Z0);

set(handles.edit_Vx0,'String',Vx0);
set(handles.edit_Vy0,'String',Vy0);
set(handles.edit_Vz0,'String',Vz0);

set(handles.edit_theta,'String',theta);
set(handles.edit_rp,'String',rp);
set(handles.edit_ra,'String',ra);
set(handles.edit_a,'String',a);

function simular_C_Callback(hObject, eventdata, handles)

global e i w omega X0 Y0 Z0 Vx0 Vy0 Vz0 Dp D d Tp T t Th ra rp theta sim_check Xl Yl Zl Vxl Vyl Vzl
sim_check = 1;

X0 = eval(get(handles.edit_X0,'String'));
Y0 = eval(get(handles.edit_Y0,'String'));
Z0 = eval(get(handles.edit_Z0,'String'));
R= [X0 Y0 Z0];
V =[Vx0 Vy0 Vz0];

r = norm(R);
v = norm(V);
vr= dot(V,R)/r;

H = cross(R,V);
h = norm(H);

Hz = H(3);

i = acos(Hz/h);
i = radtodeg(i);

N = cross([0 0 1], H);
n = norm(N);
Nx = N(1);
Ny = N(2);

if(Ny >= 0)
    omega = acos(Nx/n);
else
    omega = 2*pi - acos(Nx/n);
end
omega = radtodeg(omega);

u = 3.9860040*(10^5);
E = (1/u) * ((v^2 - (u/r))*R- (r*vr*V));
e = norm(E);

Ez = E(3);
if Ez >= 0
    w = acos((dot(N,E))/(n*e));
else
    w = 2*pi - acos(dot(N,E)/(n*e));
end

w = radtodeg(w);

if vr >= 0
    theta = acos(dot(E,R)/(e*r));
else
    theta = 2*pi - acos(dot(E,R)/(e*r));
end

theta = radtodeg(theta);

rp = (h^2/u)*(1/(1+e*cos(0)));

ra = (h^2/u)*(1/(1+e*cos(pi)));

a = (1/2)*(rp+ra);

Th = (((2*pi)/(sqrt(u)))* (a^(3/2)))/3600;


Delta_JD = juliano(T, D, t, d);
Delta_JD_per = juliano(Tp, Dp, T, D);

if Delta_JD == 0
    errordlg('O instante inicial n?o pode ser id?ntico ao instante final.','Erro');
    return
end

if Delta_JD_per == 0
    errordlg('O instante de passagem pelo perigeu n?o pode ser id?ntico ao instante de refer?ncia.','Erro');
    return
end

if Delta_JD < 0 
    errordlg('O instante inicial n?o pode ser maior que o instante final.','Erro');
    return
end

if Delta_JD_per < 0 
    errordlg('O instante de passagem pelo perigeu n?o pode ser maior que o instante de refer?ncia.','Erro');
    return
end
global pertubacao
if pertubacao == 0
    options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6); 
    [tempo,valores_saida] = ode45( @odefun,[0 Delta_JD],[X0 Y0 Z0 Vx0 Vy0 Vz0],options);
    [lon,lat,r] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));
end
if pertubacao == 1
    options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6);
    
    Iniciais = [X0 Y0 Z0 Vx0 Vy0 Vz0 Xl Yl Zl Vxl Vyl Vzl];
    [~,valores_saida] = ode45(@odefun_pert,[0 Delta_JD],Iniciais,options);
    [lon,lat,~] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));
    display(valores_saida(size(valores_saida,1),:));

end



global p
delete(p);
p = plot3(handles.plot3D,valores_saida(:,1),valores_saida(:,2),valores_saida(:,3), 'r');

global p2
delete(p2);
p2 = plot(handles.plotMap,lon*180/pi,lat*180/pi,'.k');
%axis equal;
axis(handles.plot3D,'equal');
rotate3d(handles.plot3D,'on');

set(handles.edit_w,'String',w);
set(handles.edit_e,'String',e);
set(handles.edit_ra,'String',ra);
set(handles.edit_a,'String',a);
set(handles.edit_rp,'String',rp);
set(handles.edit_theta,'String',theta);
set(handles.edit_i,'String',i);
set(handles.edit_omega,'String',omega);

function Rot = rotation(omega,i,w)
    Rz_omega = [cosd(-omega) sind(-omega) 0;-sind(-omega) cosd(-omega) 0;0 0 1];
    Rx_i = [1 0 0;0 cosd(-i) sind(-i);0 -sind(-i) cosd(-i)];
    Rz_w = [cosd(-w) sind(-w) 0;-sind(-w) cosd(-w) 0;0 0 1];

    Rot = Rz_omega*Rx_i;
    Rot = Rot*Rz_w;
    
function [G] = odefun(~,I)
    global inercial
    G = zeros(6,1);

    w = 7.29*10^-5;
    ut = 3.9860040*(10^5);
    r = sqrt(I(1)^2 + I(2)^2 + I(3)^2);

    G(1) = I(4);
    G(2) = I(5);
    G(3) = I(6);
    
    if inercial == 1
        G(4) = -ut*I(1)/(r^3)+(w^2)*I(1)+2*w*I(5);
        G(5) = -ut*I(2)/(r^3)+(w^2)*I(2)-2*w*I(4);
        G(6) = -ut*I(3)/(r^3);    
    
    else
        G(4) = -ut*I(1)/(r^3);
        G(5) = -ut*I(2)/(r^3);
        G(6) = -ut*I(3)/(r^3); 
    end
function [G] = odefun_pert(~,I)

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
    
    global inercial;
    if (inercial == 1)
        G(1) = I(4) + w*I(2)+2*w*I(5);
        G(2) = I(5) - w*I(1)-2*w*I(4);
    else
        G(1) = I(4);
        G(2) = I(5);
    end 
    
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

function Earth(handles)
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
global imagem
imagem = 'https://www.evl.uic.edu/pape/data/Earth/2048/BigEarth.jpg';
erad    = 6718; % equatorial radius (km)
prad    = 6718; % polar radius (km)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
globe = surf(handles.plot3D,x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata = imread(imagem);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

function [Pos_0,V_0] = ValoresIniciais(e,i,omega,w,Tp,Dp,T,D)

global ut 
ut = 3.9860040*(10^5);
global Rt 
Rt = 6378.16;

a = 1.5*Rt;
n = sqrt(ut/a^3);

Delta_JD_per = juliano(Tp, Dp, T, D);
M = n*Delta_JD_per;

iter = 0;
Er = 0;
ea = 100;

while(1)
    Erold = Er;
    x1 = Erold - e*sin(Erold)-M;
    x2 = 1 - e*cos(Erold);
    Er = Erold - (x1/x2);
    iter = iter + 1;
    if Er ~= 0
        ea = abs((Er - Erold)/Er)*100;
    end
    if ea < 0.005 || iter >= 100
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

Pos = [X;Y;Z];
V = [Vx;Vy;Vz];


Rot= rotation(omega,i,w);

Pos_0 = Rot*Pos;
V_0 = Rot*V;

%Correcao dos valores para o modelo inercial
global wt 
wt = 7.29*10^-5;
global inercial;


function [Delta_JD] = juliano(T, D, t, d)

    UT1 = T(1)+(T(2)/60)+(T(3)/3600);
    ut2 = t(1)+(t(2)/60)+(t(3)/3600);

    mc1 = floor(7*(D(3)+floor((D(2)+9)/12)/4));
    Jo1 = 367*D(3)-mc1+floor(257*D(2)/9)+D(1)+1721013.5;

    mc2 = floor(7*(d(3)+floor((d(2)+9)/12)/4));
    Jo2 = 367*d(3)-mc2+floor(257*d(2)/9)+d(1)+1721013.5;

    JD1 = Jo1 +(UT1/24);
    JD2 = Jo2 +(ut2/24);

    Delta_JD = (JD2 - JD1)*86400;

function pertub_Callback(hObject, eventdata, handles)
global pertubacao;
if (get(hObject,'Value') == get(hObject,'Max'))
	pertubacao =1;
else
	pertubacao = 0;
end
function inerc_Callback(hObject, eventdata, handles)
global inercial;
if (get(hObject,'Value') == get(hObject,'Max'))
	inercial = 1;
else
	inercial = 0;
end

function Save_button_Callback(hObject, eventdata, handles)
    global i e w omega T t D d Tp Dp inercial
    %Organiza os dados em uma struct
    
    [Pos_0,V_0] = ValoresIniciais(e,i,omega,w,Tp,Dp,T,D);
    Configuracao.inercial = inercial;
    Configuracao.dia0 = D(1);
    Configuracao.mes0 = D(2);
    Configuracao.ano0 = D(3);
    Configuracao.diaF = d(1);
    Configuracao.mesF = d(2);
    Configuracao.anoF = d(3);
    Configuracao.hora0 = T(1);
    Configuracao.min0 = T(2);
    Configuracao.seg0 = T(3);
    Configuracao.horaF = t(1);
    Configuracao.minF = t(2);
    Configuracao.segF = t(3);
    Configuracao.diaP = Dp(1);
    Configuracao.mesP = Dp(2);
    Configuracao.anoP = Dp(3);
    Configuracao.horaP = Tp(1);
    Configuracao.minP = Tp(2);
    Configuracao.segP = Tp(3);
    Configuracao.e = e;
    Configuracao.i = i;
    Configuracao.w = w;
    Configuracao.omega = omega;
    Configuracao.x0 = Pos_0(1);
    Configuracao.y0 = Pos_0(2);
    Configuracao.z0 = Pos_0(3);
    Configuracao.vx0 = V_0(1);
    Configuracao.vy0 = V_0(2);
    Configuracao.vz0 = V_0(3);
    %Pega o horario e a data para gerar o arquivo
    c = clock;
    nome = strcat('SetKepler_(',num2str(c(3)),'_',num2str(c(2)),'_',num2str(c(1)),')_(',num2str(c(4)),'_',num2str(c(5)),'_',num2str(c(6)),').mat');
    save(nome,'Configuracao');
    msgbox('Configuracao Salva');    
function Load_button_Callback(hObject, eventdata, handles)
[filename pathname] = uigetfile({'*.MAT'},'File Selector');
fullpathname = strcat(pathname,filename);
load(fullpathname);
global i e w omega T t D d Tp Dp inercial X0 Y0 Z0 Vx0 Vy0 Vz0

X0 = Configuracao.x0;
Y0 = Configuracao.x0;
Z0 = Configuracao.x0;
Vx0 = Configuracao.vx0;
Vy0 = Configuracao.vy0;
Vz0 = Configuracao.vz0;

D = [Configuracao.dia0,Configuracao.mes0, Configuracao.ano0];
d = [Configuracao.diaF,Configuracao.mesF, Configuracao.anoF];
T = [Configuracao.hora0,Configuracao.min0, Configuracao.seg0];
t = [Configuracao.horaF,Configuracao.minF, Configuracao.segF];

Tp = [Configuracao.horaP, Configuracao.minP, Configuracao.segP];
Dp = [Configuracao.diaP, Configuracao.mesP, Configuracao.anoP];

i = Configuracao.i;
e = Configuracao.e;
w = Configuracao.w;
omega = Configuracao.omega;

inercial = Configuracao.inercial;

set(handles.edit_X0,'String',Configuracao.x0);
set(handles.edit_Y0,'String',Configuracao.y0);
set(handles.edit_Z0,'String',Configuracao.z0);
set(handles.edit_Vx0,'String',Configuracao.vx0);
set(handles.edit_Vy0,'String',Configuracao.vy0);
set(handles.edit_Vz0,'String',Configuracao.vz0);
set(handles.edit_e,'String',Configuracao.e);
set(handles.edit_i,'String',Configuracao.i);
set(handles.edit_w,'String',Configuracao.w);
set(handles.edit_omega,'String',Configuracao.omega);

set(handles.Dia_per,'Value',Configuracao.diaP);
set(handles.Mes_per,'Value',Configuracao.mesP);

set(handles.Dia_0,'Value',Configuracao.dia0);
set(handles.Mes_0,'Value',Configuracao.mes0);
% set(handles.Ano_0,'Value',Configuracao.anoF);


set(handles.Dia_F,'Value',Configuracao.diaF);
set(handles.Mes_F,'Value',Configuracao.mesF);
% set(handles.Ano_F,'Value',Configuracao.anoF);

set(handles.Hora_per,'Value',Configuracao.horaP+1);
set(handles.Min_per,'Value',Configuracao.minP+1);
set(handles.Seg_per,'Value',Configuracao.segP+1);

set(handles.Hora_0,'Value',Configuracao.hora0+1);
set(handles.Min_0,'Value',Configuracao.min0+1);
set(handles.Seg_0,'Value',Configuracao.seg0+1);

set(handles.Hora_F,'Value',Configuracao.horaF+1);
set(handles.Min_F,'Value',Configuracao.minF+1);
set(handles.Seg_F,'Value',Configuracao.segF+1);
% if inercial == 1
%     set(handles.inerc,'Value','Max');
% end
function togglebutton1_Callback(hObject, eventdata, handles)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
	set(handles.text1,'String','Output');
    set(handles.text2,'String','Input');

    setCar_OFF(handles);
    setKep_ON(handles);
    
elseif button_state == get(hObject,'Min')
	set(handles.text2,'String','Output');
    set(handles.text1,'String','Input');
    
    setCar_ON(handles);
    setKep_OFF(handles);
end
function setKep_ON(handles)
    set(handles.edit_e,'enable','on');
    set(handles.edit_i,'enable','on');
    set(handles.edit_w,'enable','on');
    set(handles.edit_omega,'enable','on');
    set(handles.simular_K,'enable','on','String','Simular');
function setKep_OFF(handles)
    set(handles.edit_e,'enable','off');
    set(handles.edit_i,'enable','off');
    set(handles.edit_w,'enable','off');
    set(handles.edit_omega,'enable','off');
    set(handles.simular_K,'enable','off','String','Simular');
function setCar_OFF(handles)
    set(handles.edit_X0,'enable','off');
    set(handles.edit_Y0,'enable','off');
    set(handles.edit_Z0,'enable','off');
    set(handles.edit_Vx0,'enable','off');
    set(handles.edit_Vy0,'enable','off');
    set(handles.edit_Vz0,'enable','off');
    set(handles.simular_C,'enable','off','String','Simular');
function setCar_ON(handles)
    set(handles.edit_X0,'enable','on');
    set(handles.edit_Y0,'enable','on');
    set(handles.edit_Z0,'enable','on');
    set(handles.edit_Vx0,'enable','on');
    set(handles.edit_Vy0,'enable','on');
    set(handles.edit_Vz0,'enable','on');
    set(handles.simular_C,'enable','on','String','Simular');
    
function edit_e_Callback(hObject, eventdata, handles)
global e
e = eval(get(gcbo,'string'));
function edit_e_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_i_Callback(hObject, eventdata, handles)
global i
i = eval(get(gcbo,'string'));
function edit_i_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_w_Callback(hObject, eventdata, handles)
global w
w = eval(get(gcbo,'string'));
function edit_w_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_omega_Callback(hObject, eventdata, handles)
global omega
omega = eval(get(gcbo,'string'));
function edit_omega_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_theta_Callback(hObject, eventdata, handles)
global theta
theta = eval(get(gcbo,'string'));
function edit_theta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_rp_Callback(hObject, eventdata, handles)
global rp
rp = eval(get(gcbo,'string'));
function edit_rp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ra_Callback(hObject, eventdata, handles)
global ra
ra = eval(get(gcbo,'string'));
function edit_ra_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_a_Callback(hObject, eventdata, ~)
global a
a = eval(get(gcbo,'string'));
function edit_a_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_per_Callback(hObject, eventdata, handles)
global per
per = eval(get(gcbo,'string'));
function edit_per_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_X0_Callback(hObject, eventdata, handles)
global X0
X0 = eval(get(gcbo,'string'));
function edit_X0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Y0_Callback(hObject, eventdata, handles)
global Y0
Y0 = eval(get(gcbo,'string'));
function edit_Y0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Z0_Callback(hObject, eventdata, handles)
global Z0
Z0 = eval(get(gcbo,'string'));
function edit_Z0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vx0_Callback(hObject, eventdata, handles)
global Vx0
Vx0 = eval(get(gcbo,'string'));
function edit_Vx0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vy0_Callback(hObject, eventdata, handles)
global Vy0
Vy0 = eval(get(gcbo,'string'));
function edit_Vy0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vz0_Callback(hObject, eventdata, handles)
global Vz0
Vz0 = eval(get(gcbo,'string'));
function edit_Vz0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Dia_0_Callback(hObject, eventdata, handles)
global D
s = get(hObject,'String');
i = get(hObject,'Value');
D(1) = str2num(s{i});
function Dia_0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Mes_0_Callback(hObject, eventdata, handles)
global D
s = get(hObject,'String');
i = get(hObject,'Value');
D(2) = str2num(s{i});
function Mes_0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ano_0_Callback(hObject, eventdata, handles)
global D 
s = get(hObject,'String');
i = get(hObject,'Value');
D(3) = str2num(s{i});
function Ano_0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Hora_0_Callback(hObject, eventdata, handles)
global T
s = get(hObject,'String');
i = get(hObject,'Value');
T(1) = str2num(s{i});
function Hora_0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Min_0_Callback(hObject, eventdata, handles)
global T
s = get(hObject,'String');
i = get(hObject,'Value');
T(2) = str2num(s{i});
function Min_0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Seg_0_Callback(hObject, eventdata, handles)
global T
s = get(hObject,'String');
i = get(hObject,'Value');
T(3) = str2num(s{i});
function Seg_0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Dia_F_Callback(hObject, eventdata, handles)
global d
s = get(hObject,'String');
i = get(hObject,'Value');
d(1) = str2num(s{i});
function Dia_F_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Mes_F_Callback(hObject, eventdata, handles)
global d
s = get(hObject,'String');
i = get(hObject,'Value');
d(2) = str2num(s{i});
function Mes_F_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ano_F_Callback(hObject, eventdata, handles)
global d
s = get(hObject,'String');
i = get(hObject,'Value');
d(3) = str2num(s{i});
function Ano_F_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Hora_F_Callback(hObject, eventdata, handles)
global t
s = get(hObject,'String');
i = get(hObject,'Value');
t(1) = str2num(s{i});
function Hora_F_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Min_F_Callback(hObject, eventdata, handles)

global t
s = get(hObject,'String');
i = get(hObject,'Value');
t(2) = str2num(s{i});
function Min_F_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Seg_F_Callback(hObject, eventdata, handles)
global t
s = get(hObject,'String');
i = get(hObject,'Value');
t(3) = str2num(s{i});
function Seg_F_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Hora_per.
function Hora_per_Callback(hObject, eventdata, handles)

global Tp
s = get(hObject,'String');
i = get(hObject,'Value');
Tp(1) = str2num(s{i});
function Hora_per_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Min_per.
function Min_per_Callback(hObject, eventdata, handles)
global Tp
s = get(hObject,'String');
i = get(hObject,'Value');
Tp(2) = str2num(s{i});
function Min_per_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Seg_per.
function Seg_per_Callback(hObject, eventdata, handles)
global Tp
s = get(hObject,'String');
i = get(hObject,'Value');
Tp(3) = str2num(s{i});
function Seg_per_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Dia_per.
function Dia_per_Callback(hObject, eventdata, handles)
global Dp
s = get(hObject,'String');
i = get(hObject,'Value');
Dp(1) = str2num(s{i});
function Dia_per_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Mes_per.
function Mes_per_Callback(hObject, eventdata, handles)
global Dp
s = get(hObject,'String');
i = get(hObject,'Value');
Dp(2) = str2num(s{i});
function Mes_per_CreateFcn(hObject, eventdata, handles)


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Ano_per.
function Ano_per_Callback(hObject, eventdata, handles)
global Dp
s = get(hObject,'String');
i = get(hObject,'Value');
Dp(3) = str2num(s{i});
function Ano_per_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function animate_Callback(hObject, eventdata, handles)

global e i w omega Dp D d Tp T t imagem sim_check

if sim_check == 0 
    errordlg('Voc? precisa simular algo antes de poder gerar uma anima??o!','Erro');
    return;
end

e = eval(get(handles.edit_e,'String'));
i = eval(get(handles.edit_i,'String'));
w = eval(get(handles.edit_w,'String'));
omega = eval(get(handles.edit_omega,'String'));

[Pos_0,V_0] = ValoresIniciais(e,i,omega,w,Tp,Dp,T,D);

Delta_JD = juliano(T, D, t, d);
Delta_JD_per = juliano(Tp, Dp, T, D);

if Delta_JD == 0
    errordlg('O instante inicial n?o pode ser id?ntico ao instante final.','Erro');
    return
end

if Delta_JD_per == 0
    errordlg('O instante de passagem pelo perigeu n?o pode ser id?ntico ao instante de refer?ncia.','Erro');
    return
end

if Delta_JD < 0 
    errordlg('O instante inicial n?o pode ser maior que o instante final.','Erro');
    return
end

if Delta_JD_per < 0 
    errordlg('O instante de passagem pelo perigeu n?o pode ser maior que o instante de refer?ncia.','Erro');
    return
end


%Aplica??o do ode45 para os valores iniciais calculados
options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6); 

[tempo,valores_saida] = ode45( @odefun,[0 Delta_JD],[Pos_0(1) Pos_0(2) Pos_0(3) V_0(1) V_0(2) V_0(3)],options);

animacao = figure(1);

ax = gca;
figura = gcf;

npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; 

erad    = 6718; % equatorial radius (km)
prad    = 6718; % polar radius (km)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata = imread(imagem);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

title('Movimento do corpo','Color','w','FontSize',24,'FontName','Calibri');
xlabel('Posi??o em x [km]');
ylabel('Posi??o em y [km]');
zlabel('Posi??o em z [km]');
grid on;
set(figura,'color','black');
ax.LineWidth = 1.5;
ax.GridLineStyle = '--';
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';

hold on;

axis equal;
rotate3d on;
comet3(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));
