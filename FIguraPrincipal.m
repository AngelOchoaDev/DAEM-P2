function varargout = FIguraPrincipal(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FIguraPrincipal_OpeningFcn, ...
                   'gui_OutputFcn',  @FIguraPrincipal_OutputFcn, ...
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


% --- Executes just before FIguraPrincipal is made visible.
function FIguraPrincipal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose default command line output for FIguraPrincipal
handles.output = hObject;
%%%%%%%%%%%%%%aqui es el codigo que se ejecuta al inicio%%%%%%%%%%%%%%%%%%
filas={'N_dtes' 'PasoDiam' 'F' 'AngHelice' ...
'Potencia' 'Vel_Ang' 'DiamEje' 'LongEje'};
Datos={14 54 'd=' 0 false;
10 10 'v=' 0 false;
1.0 1.0 'w=' 0 false;
30 30 'wt' 0 false;
10 10 'wr' 0 true;
900 0 'wa' 0 true;
0.75 1.75 'nr' 5 false;
4 4 'nd' 9 false};%
set(handles.uitable1,'data',Datos,'RowName',filas,...
'ColumnName',{'ValoresP' 'ValoresE' 'Variables' ...
'Valores' 'Carga(+)Sign'},'ColumnEditable',true,...
'ForegroundColor',[0 0 1],'FontSize',10,...
'ColumnWidth',{65 75 55 75 85});
%ahora se setea la info en la tabla 2 y 3
Posiones={0 0 0 true;
0 0 -5/2 false;
-0.7 0 -5/4 false};
set(handles.uitable2,'FontSize',10,'data',Posiones,...
'RowName',{'PuntoA' 'PuntoB' 'PosCarga'},...
'ColumnName',{'x/Or' 'y/Or' 'z/Or' 'x/A' 'y/A' 'z/A' 'x/B' 'y/B' 'z/B'},...
'ColumnEditable',true,...
'ForegroundColor',[0 0 1],...
'ColumnWidth',{75});
Fuerzas={0 0 0;0 0 0;0 0 0;0 0 0};
set(handles.uitable3,'FontSize',10,'data',Fuerzas,...
'RowName',{'PuntoA' 'PuntoB' 'PosCarga'},...
'RowName',{'ReaccA' 'ReaccB' 'Carga'},...
'ColumnName',{'Datax' 'Datay' 'Dataz'},...
'ColumnEditable',true,...
'ForegroundColor',[0 0 1],...
'ColumnWidth',{75});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FIguraPrincipal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FIguraPrincipal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
%uuusuuhuusu
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%se setea la informacion en la tabla 1
Z=get(handles.uitable1,'data'); %datos de engrane
P=get(handles.uitable2,'data'); %datos de posicion
F=get(handles.uitable3,'data'); %datos de fuerzas
if get(handles.radiobutton1,'value')==1
    Z{6,2}=Z{1,1}*Z{6,1}/Z{1,2};
    Z{1,4}=Z{1,1}/Z{2,1}; %d=N/P
    Z{2,4}=pi*Z{1,4}*Z{6,1}/2; %vel tangencial (pies/min)
    Z{3,4}=(33000*Z{5,1}/Z{2,4})/(cosd(20)*cosd(Z{4,1})); %W=(33000*H/V)/(cos(fi_n)*cos(psi))
    Z{4,4}=33000*Z{5,1}/Z{2,4}; %Wt=33000*H/V
    Z{5,4}=(33000*Z{5,1}/Z{2,4})*tand(20)/cosd(Z{4,1}); %Wr
    Z{6,4}=(33000*Z{5,1}/Z{2,4})*tand(Z{4,1}); %Wa
    F{3,1}=Z{5,4};F{3,2}=Z{4,4};F{3,3}=Z{4,1};
    %correccion de signo para los valores de las tablas
    %Para x de la carga
    if Z{5,5}==1
        F{3,1}=Z{5,4};
    else
        F{3,1}=-Z{5,4};
    end
    %Para y de la carga
    if Z{4,5}==1
        F{3,2}=Z{4,4};
    else
        F{3,2}=-Z{4,4};
    end
    %Para z de la carga
    if Z{6,5}==1
        F{3,3}=Z{6,4};
    else
        F{3,3}=-Z{6,4};
    end
    %Calculo de posiciones relativas Respecto de A
    P{1,4}=P{1,1}-P{1,1};P{1,5}=P{1,2}-P{1,2};P{1,6}=P{1,3}-P{1,3};
    P{2,4}=P{2,1}-P{1,1};P{2,5}=P{2,2}-P{1,2};P{2,6}=P{2,3}-P{1,3};
    P{3,4}=P{3,1}-P{1,1};P{3,5}=P{3,2}-P{1,2};P{3,6}=P{3,3}-P{1,3};
    %Calculo de posiciones relativas Respecto de B
    P{1,7}=P{1,1}-P{2,1};P{1,8}=P{1,2}-P{2,2};P{1,9}=P{1,3}-P{2,3};
    P{2,7}=P{2,1}-P{2,1};P{2,8}=P{2,2}-P{2,2};P{2,9}=P{2,3}-P{2,3};
    P{3,7}=P{3,1}-P{2,1};P{3,8}=P{3,2}-P{2,2};P{3,9}=P{3,3}-P{2,3};
    MB=cross([P{3,7:9}],[F{3,:}]); %MB=PF1/BxF1;
    F{4,1}=MB(1,1); F{4,2}=MB(1,2); F{4,3}=MB(1,3);
    F{1,1}=-MB(1,2)/P{1,9}; F{1,2}=MB(1,1)/P{1,9};
    MA=cross([P{3,4:6}],[F{3,1:3}]); %MA = (r1/A)xF1
    F{5,1}=MA(1,1); F{5,2}=MA(1,2); F{5,3}=MA(1,3);
    F{2,1}=-MA(1,2)/P{2,6}; F{2,2}=MA(1,1)/P{2,6};
    %Se grafican los vectores de fuerza
    Vector3D202202({F{3,1}/max([F{3,1} F{3,2} F{3,3}]) 0.5 'x'},...
        [P{3,1} P{3,2} P{3,3}]);
end
set(handles.uitable1,'data',Z);
set(handles.uitable2,'data',P);
set(handles.uitable3,'data',F);
%%%%%%%-------------DATOS DE LOS ENGRANES--------------%%%%%%%%%
%The parameters of gear
% n=10;			%The number of teeth
% pd=6;			%The diametral pitch
n=Z{1,1};%str2double(get(handles.edit1,'string')); %dientes
pd=Z{2,1};%str2double(get(handles.edit2,'string'));	%paso diametral
phi_d=20;	%The pressure angle in degrees 
% ----------------------------------------------------------------------------
r_fillet=0.05;		%The radius of fillet
% ----------------------------------------------------------------------------
%To declare variables
xp=zeros(10,1);yp=zeros(10,1);
xo=zeros(5,1);yo=zeros(5,1);
xr=zeros(3,1);yr=zeros(3,1);
xro=zeros(5,1);yro=zeros(5,1);
xf=zeros(5,1);yf=zeros(5,1);
theta=zeros(10,1);
f=zeros(2,28);
M=[];c=[];e=[];g=[];h=[];
% ----------------------------------------------------------------------------
%To calculate the basic parameters of a gear
d=n/pd;				%pitch diamter
phi=phi_d*pi/180;	%pressure angle in radians
db=d*cos(phi);		%diameter of base circle
do=d+2/pd;			%addendum (outside) diameter
tt=pi/(2*pd);		%tooth thickness at the pitch circle
dr=d-2*1.25/pd;	%dedendum (root) diameter
% ----------------------------------------------------------------------------
%To calculate the coordinates of the involute profile
n1=10;
tp=pi*d/(2*n);
for i=1:n1;
	r=do/2-(do-db)*(i-1)/(2*(n1-1));
	pha=acos(db/(2*r));
	t=2*r*(tp/d+(tan(phi)-phi)-(tan(pha)-pha));	%tooth tickness at any angle phi
																%involute equation - 
																%refer to Shigley's book
	theta(i)=t/(2*r);
	xp(i)=r*sin(theta(i));		%change from polar coordinates to cartesian coordinates
	yp(i)=r*cos(theta(i));
end
xp=xp';yp=yp';
% ----------------------------------------------------------------------------
%To calculate the addendum circle
n2=5;
for i=1:n2;
	theta_o=theta(1)*(i-1)/(n2-1);
	xo(i)=(do/2)*sin(theta_o);
	yo(i)=(do/2)*cos(theta_o);
end
xo=xo';yo=yo';
% ----------------------------------------------------------------------------
%To calculate the non-involute portion of the curve- between the base circle and 
% dedendum circle - in this case, a straight line parallel to the y axis and connects 
% to the fillet arc
for i=1:3;
	theta0=asin((xp(1,n1)+r_fillet)/(dr/2)); 
%to find the angle between the central line (y-axis) and the line from the center 
%to the last point of the involute curve.
	xr(i)=xp(1,10);
	yr(i)=yp(1,10)-(yp(1,10)-r_fillet-(dr/2)*cos(theta0))*i/3;
end
xr=xr';yr=yr';
% ----------------------------------------------------------------------------
%To calculate the dedendum circle
n3=5;
for i=1:n3;
   thetar=theta0+(pi/n-theta0)*(i-1)/(n3-1);	
   %(pi/n-theta0) angle subtended for dededem arc
   xro(i)=dr*sin(thetar)/2;	%xro(1) is the beginning point
	yro(i)=dr*cos(thetar)/2;
end
xro=xro';yro=yro';
% ----------------------------------------------------------------------------
%To calculate fillet
% to draw the quarter of a circle from the last point of the non-involute part to 
% the tangent of the dedenum circle.
n4=5;
for i=1:n4;
   xf(i)=xro(1)-r_fillet*cos((i-1)*pi/(2*n4-2));
   yf(i)=yro(1)+r_fillet*(1-sin((i-1)*pi/(2*n4-2)));	%yf(5)=yro(1)-r_fillet*sin(4*pi/8)
end
xf=xf';yf=yf';
% ----------------------------------------------------------------------------
%To append each piece of curve to generate one-half of a tooth profile
c=[c,xo,xp,xr,xf,xro];
e=[e,yo,yp,yr,yf,yro];
g=[c',e'];
g=g';						%the one-half  tooth profile
% ----------------------------------------------------------------------------
%To reflecte the involute curve about y axis to get the whole tooth
ff=[-1 0;0 1]*g;		%reflection 
n5=n1+n2+n3+n4+3
for i=1:n5;				%n4 points =n1(involute)+n2(addendum)+n4(fillet)
   						%			+3(noninvolute)+n3(dedendum)
   f(1,i)=ff(1,n5+1-i);	%reverse the order of the points, easy for plotting 
	f(2,i)=ff(2,n5+1-i);
end
h=[h,f,g];				%the whole tooth profile
% ----------------------------------------------------------------------------
%To rotate and append the tooth to generate the gear
for i=1:n;
	kk=[cos(2*pi*(i-1)/n) sin(2*pi*(i-1)/n);-sin(2*pi*(i-1)/n) cos(2*pi*(i-1)/n)];
  												 		%rotation matrix
	mm=kk*h;		%rotate
	M=[M,mm];	%append
end
M=[M,h(:,1)]; %add the first point, so the curve returns to the original point
Q=atan2d(M(:,2),M(:,1));
x=Z{7,1}*cosd(Q)/2;
y=Z{7,1}*sind(Q)/2;
% ----------------------------------------------------------------------------
cla; %Para limpiar la figura, hasta este punto solo se tenia el diente rojo ploteado
hold on
%-----------------------------------------------------------------------------
%plot (g(1,:),g(2,:))	%plot one-half tooth
%plot (h(1,:),h(2,:))	%plot one tooth
%plot (M(1,:),M(2,:))	%plot the whole gear - the first row (x) and second row (y)
%%%plot (g(1,:),g(2,:),'-.b',xo,yo,'-r', 'linewidth',4)  %plot one-half tooth, the addendum part is red
%plot(M(1,:),M(2,:),'.b', h(1,:),h(2,:),'-r', 'linewidth',2) %plot whole gear, the original tooth is red
%plot(M(1,:),M(2,:),'.b','linewidth',2) %plot whole gear, without the red tooth
%-------------------------------------------------------------------------------
%ploteando para ver las dos caras del engrane
%Primera cara
%plot3(M(1,:),M(2,:),zeros(size(M(2,:)))+Z{3,1},'.b', 'linewidth',2)

%se intenta plotear todo en dos plots nada mas
%cara exterior

%plot3([M(1,:) M(1,:)],[M(2,:) M(2,:)],[zeros(size(M(2,:)))+Z{3,1}/2 zeros(size(M(2,:)))-Z{3,1}/2], '.b','linewidth',2)
%axis('equal')
%hold on

surf([M(1,:);M(1,:)],...
     [M(2,:);M(2,:)],...
     [zeros(size(M(2,:)))+Z{3,1}/2;...
      zeros(size(M(2,:)))-Z{3,1}/2],...
      'FaceAlpha',0.75)
shading interp
axis('equal')
hold on

%Cara interior
M=M';
%plot3([x x],[y y],[zeros(size(y))+Z{3,1}/2 zeros(size(y))-Z{3,1}/2],'.b','linewidth',0.5)
surf([x x x],[y y y],[zeros(size(y))+Z{3,1}/2 zeros(size(y))+Z{3,1}/2 zeros(size(y))-Z{3,1}/2],'FaceAlpha',0.75)
shading interp
h = rotate3d;set(h,'RotateStyle','box','Enable','on');hold on;
%----------------------------------------------------
%transpose the matrix to get only two columns,
