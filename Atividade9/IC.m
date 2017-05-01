function varargout = IC(varargin)
% IC MATLAB code for IC.fig
%      IC, by itself, creates a new IC or raises the existing
%      singleton*.
%
%      H = IC returns the handle to a new IC or the handle to
%      the existing singleton*.
%
%      IC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IC.M with the given text1 arguments.
%
%      IC('Property','Value',...) creates a new IC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IC

% Last Modified by GUIDE v2.5 26-Apr-2017 08:19:56

% Begin initialization code - DO NOT EDIT
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
% End initialization code - DO NOT EDIT

% --- Executes just before IC is made visible.
function IC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IC (see VARARGIN)

% Choose default command line output for IC
handles.output = hObject;
set(handles.plot3D,'Xlim',[-2.5 2.5]*10^4,'Ylim',[-2.5 2.5]*10^4,'Zlim',[-1.5 1.5]*10^4);
%Plot da figura central!
Earth(handles);
axis(handles.plot3D,[-2.5 2.5 -2.5 2.5 -1.5 1.5]*10^4);
hold on;
grid on;
set(handles.plot3D,'GridLineStyle','--');
global cor_eixo
cor_eixo = [0.9 0.4 0.74];

set(handles.plot3D,'color','k','Xcolor','w','Ycolor','w','Zcolor','w');
guidata(hObject, handles);

global T
global t
global D
global d
d = [1 1 2000];
D = [1 1 2000];
T = [1 1 1];
t = [1 1 1];

global k

% UIWAIT makes IC wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.edit_e,'enable','off');
set(handles.edit_i,'enable','off');
set(handles.edit_w,'enable','off');
set(handles.edit_omega,'enable','off');
set(handles.simular_K,'enable','off','String','Simular');
    
% --- Outputs from this function are returned to the command line.
function varargout = IC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1
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
    
function edit_e_Callback(hObject, eventdata, handles)
% hObject    handle to edit_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_e as text
%        eval(get(hObject,'String')) returns contents of edit_e as a double
global e
e = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_e_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_i_Callback(hObject, eventdata, handles)
% hObject    handle to edit_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_i as text
%        eval(get(hObject,'String')) returns contents of edit_i as a double
global i
i = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_w_Callback(hObject, eventdata, handles)
% hObject    handle to edit_omega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_omega as text
%        eval(get(hObject,'String')) returns contents of edit_omega as a double
global w
w = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_w_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_omega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_omega_Callback(hObject, eventdata, handles)
% hObject    handle to edit_omega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_omega as text
%        eval(get(hObject,'String')) returns contents of edit_omega as a double
global omega
omega = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_omega_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_omega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in simular_K.
function simular_K_Callback(hObject, eventdata, handles)

% global ut 
% ut = 3.9860040*(10^5);
% global Rt 
% Rt = 6378.16;
% 
% a = 1.5*Rt;

global e
global i
global w
global omega
global X0
global Y0
global Z0
global Vx0
global Vy0
global Vz0
global D
global d
global T
global t

[Pos_0,V_0] = ValoresIniciais(e,i,omega,w,T,D,t,d);
Delta_JD = juliano(T, D, t, d);

if Delta_JD == 0
    errordlg('O instante inicial não pode ser idêntico ao instante final.','Erro');
    return
end

if Delta_JD < 0 
    errordlg('O instante inicial não pode ser maior que o instante final.','Erro');
    return
end
%Aplicação do ode45 para os valores iniciais calculados
options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6); 

[tempo,valores_saida] = ode45( @odefun,[0 Delta_JD],[Pos_0(1) Pos_0(2) Pos_0(3) V_0(1) V_0(2) V_0(3)],options);
[lon,lat,r] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));
%set(handles.plot3D,'Xlim',[-2.5 2.5]*10^4,'Ylim',[-2.5 2.5]*10^4,'Zlim',[-1.5 1.5]*10^4);
global p
delete(p);
p = plot3(handles.plot3D,valores_saida(:,1),valores_saida(:,2),valores_saida(:,3), 'r');
%hold on;
%Earth(handles);
axis equal;
grid on;
rotate3d on;
%Seta os valores em suas respectivas caixas de texto
X0 = Pos_0(1);
Y0 = Pos_0(2);
Z0 = Pos_0(3);

Vx0 = V_0(1);
Vy0 = V_0(2);
Vz0 = V_0(3);

set(handles.edit_X0,'String',X0);
set(handles.edit_Y0,'String',Y0);
set(handles.edit_Z0,'String',Z0);

set(handles.edit_Vx0,'String',Vx0);
set(handles.edit_Vy0,'String',Vy0);
set(handles.edit_Vz0,'String',Vz0);


function edit_X0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_X0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_X0 as text
%        eval(get(hObject,'String')) returns contents of edit_X0 as a double
global X0
X0 = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_X0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_X0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Y0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Y0 as text
%        eval(get(hObject,'String')) returns contents of edit_Y0 as a double
global Y0
Y0 = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_Y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Z0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Z0 as text
%        eval(get(hObject,'String')) returns contents of edit_Z0 as a double
global Z0
Z0 = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_Z0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Z0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vx0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Vx0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Vx0 as text
%        eval(get(hObject,'String')) returns contents of edit_Vx0 as a double
global Vx0
Vx0 = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_Vx0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Vx0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vy0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Vy0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Vy0 as text
%        eval(get(hObject,'String')) returns contents of edit_Vy0 as a double
global Vy0
Vy0 = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_Vy0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Vy0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Vz0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Vz0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Vz0 as text
%        eval(get(hObject,'String')) returns contents of edit_Vz0 as a double
global Vz0
Vz0 = eval(get(gcbo,'string'));

% --- Executes during object creation, after setting all properties.
function edit_Vz0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Vz0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in simular_C.
function simular_C_Callback(hObject, eventdata, handles)
% hObject    handle to simular_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global X0
global Y0
global Z0
global Vx0
global Vy0
global Vz0
global i
global omega
global e 
global w

R= [X0 Y0 Z0];
V =[Vx0 Vy0 Vz0];

disp(R);
disp(V);

r = norm(R);
v = norm(V);
vr= dot(V,R)/r;

H = cross(R,V);
h = norm(H);

Hz = H(3);

i = acos(Hz/h);
i = radtodeg(i);

set(handles.edit_i,'String',i);


N = cross([0 0 1], H);
n = norm(N);
Nx = N(1);
Ny = N(2);
if(Ny >= 0)
    omega = acos(Nx/n);
else
    omega = 2*pi - acos(Ny/n);
end
omega = radtodeg(omega);
set(handles.edit_omega,'String',omega);


u = 3.9860040*(10^5);
E = (1/u) * ((v^2 - (u/r))*R- (r*vr*V));
e = norm(E);
set(handles.edit_e,'String',e);

Ez = E(3);
if Ez >= 0
    w = acos((dot(N,E))/(n*e));
else
    w = 2*pi - acos(dot(N,E)/(n*e));
end

w = radtodeg(w);

set(handles.edit_w,'String',w);

global T
global D
global t
global d

Delta_JD = juliano(T, D, t, d);
%Aplicação do ode45 para os valores iniciais calculados
options = odeset('Abstol', [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6], 'Reltol', 1e-6); 

[tempo,valores_saida] = ode45( @odefun,[0 Delta_JD],[X0 Y0 Z0 Vx0 Vy0 Vz0],options);
[lon,lat,r] = cart2sph(valores_saida(:,1),valores_saida(:,2),valores_saida(:,3));


plot3(handles.plot3D,valores_saida(:,1),valores_saida(:,2),valores_saida(:,3), 'r');
set(handles.plot3D,'Xlim',[-2.5 2.5]*10^4,'Ylim',[-2.5 2.5]*10^4,'Zlim',[-1.5 1.5]*10^4);
grid on;
set(handles.plot3D,'GridLineStyle','--');
global cor_eixo
set(handles.plot3D,'color','k','Xcolor','w','Ycolor','w','Zcolor','w');

% if vr >=0
%     O = acos((sum(E.*R))/(e*r));
% else
%     O = 2*pi - acos((sum(E.*R))/(e*r));
% end
% O = radtodeg(O);
% disp('O(graus) = ')
% disp(O);
% 
% rp = (h^2/u)*(1/(1+e*cos(0)));
% disp('rp(km) = ')
% disp(rp);
% 
% ra = (h^2/u)*(1/(1+e*cos(pi)));
% disp('ra(km) = ')
% disp(ra);
% 
% a = (1/2)*(rp+ra);
% disp('a(km) = ')
% disp(a);
% 
% T= (((2*pi)/(sqrt(u)))* (a^(3/2)))/3600;
% disp('T(hr) = ')
% disp(T);
set(handles.edit_w,'String',w);

function Rot = rotation(omega,i,w)
    Rz_omega = [cosd(-omega) sind(-omega) 0;-sind(-omega) cosd(-omega) 0;0 0 1];
    Rx_i = [1 0 0;0 cosd(-i) sind(-i);0 -sind(-i) cosd(-i)];
    Rz_w = [cosd(-w) sind(-w) 0;-sind(-w) cosd(-w) 0;0 0 1];

    Rot = Rz_omega*Rx_i;
    Rot = Rot*Rz_w;
    
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


% --- Executes on selection change in Min_0.
function Min_0_Callback(hObject, eventdata, handles)
% hObject    handle to Min_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global T
T(2) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Min_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Min_0


% --- Executes during object creation, after setting all properties.
function Min_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Seg_0.
function Seg_0_Callback(hObject, eventdata, handles)
% hObject    handle to Seg_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global T
T(3) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Seg_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Seg_0


% --- Executes during object creation, after setting all properties.
function Seg_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Seg_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Hora_0.
function Hora_0_Callback(hObject, eventdata, handles)
% hObject    handle to Hora_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global T
T(1) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Hora_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Hora_0


% --- Executes during object creation, after setting all properties.
function Hora_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hora_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Hora_F.
function Hora_F_Callback(hObject, eventdata, handles)
% hObject    handle to Hora_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
t(1) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Hora_F contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Hora_F


% --- Executes during object creation, after setting all properties.
function Hora_F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hora_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Min_F.
function Min_F_Callback(hObject, eventdata, handles)
% hObject    handle to Min_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
t(2) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Min_F contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Min_F


% --- Executes during object creation, after setting all properties.
function Min_F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Seg_F.
function Seg_F_Callback(hObject, eventdata, handles)
% hObject    handle to Seg_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
t(3) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Seg_F contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Seg_F


% --- Executes during object creation, after setting all properties.
function Seg_F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Seg_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ano_F.
function Ano_F_Callback(hObject, eventdata, handles)
% hObject    handle to Ano_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
d(3) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Ano_F contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ano_F


% --- Executes during object creation, after setting all properties.
function Ano_F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ano_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Mes_F.
function Mes_F_Callback(hObject, eventdata, handles)
% hObject    handle to Mes_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
d(2) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Mes_F contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Mes_F


% --- Executes during object creation, after setting all properties.
function Mes_F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mes_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Dia_F.
function Dia_F_Callback(hObject, eventdata, handles)
% hObject    handle to Dia_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global d
d(1) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Dia_F contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Dia_F


% --- Executes during object creation, after setting all properties.
function Dia_F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dia_F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Dia_0.
function Dia_0_Callback(hObject, eventdata, handles)
% hObject    handle to Dia_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D
D(1) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Dia_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Dia_0


% --- Executes during object creation, after setting all properties.
function Dia_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dia_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ano_0.
function Ano_0_Callback(hObject, eventdata, handles)
% hObject    handle to Ano_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D
D(3) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Ano_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ano_0


% --- Executes during object creation, after setting all properties.
function Ano_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ano_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Mes_0.
function Mes_0_Callback(hObject, eventdata, handles)
% hObject    handle to Mes_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D
D(2) = get(hObject,'Value');
% Hints: contents = cellstr(get(hObject,'String')) returns Mes_0 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Mes_0

function Mes_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mes_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Earth(handles)
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible

image_file = 'https://www.evl.uic.edu/pape/data/Earth/2048/BigEarth.jpg';
erad    = 6718; % equatorial radius (km)
prad    = 6718; % polar radius (km)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
globe = surf(handles.plot3D,x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata = imread(image_file);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

% --- Executes on button press in Save_button.
function Save_button_Callback(hObject, eventdata, handles)
    % hObject    handle to Save_button (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    global i
    global e
    global w
    global omega
    global T
    global t
    global D
    global d
    %Organiza os dados em uma struct
    [Pos_0,V_0] = ValoresIniciais(e,i,omega,w,T,D,t,d);
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
    %Pega o horário e a data para gerar o arquivo
    c = clock;
    nome = strcat('SetKepler_(',num2str(c(3)),'_',num2str(c(2)),'_',num2str(c(1)),')_(',num2str(c(4)),'_',num2str(c(5)),'_',num2str(c(6)),').mat');
    save(nome,'Configuracao');

function [Pos_0,V_0] = ValoresIniciais(e,i,omega,w,T,D,t,d)

global ut 
ut = 3.9860040*(10^5);
global Rt 
Rt = 6378.16;

a = 1.5*Rt;
n = sqrt(ut/a^3);

Delta_JD = juliano(T, D, t, d);
M = n*Delta_JD;

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

r = a*(1 - e*cosd(E));
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

%Correção dos valores para o modelo inercial
global wt 
wt = 7.29*10^-5;
Pos_0(1) = Pos_0(1)+wt*V_0(2);
Pos_0(2) = Pos_0(2)-wt*V_0(1);

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
    
 
