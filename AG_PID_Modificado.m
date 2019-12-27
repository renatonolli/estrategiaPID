
function varargout = AG_PID_Modificado(varargin)
% AG_PID_Modificado MATLAB code for AG_PID_Modificado.fig
%      AG_PID_Modificado, by itself, creates a new AG_PID_Modificado or raises the existing
%      singleton*.
%
%      H = AG_PID_Modificado returns the handle to a new AG_PID_Modificado or the handle to
%      the existing singleton*.
%
%      AG_PID_Modificado('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AG_PID_Modificado.M with the given input arguments.
%
%      AG_PID_Modificado('Property','Value',...) creates a new AG_PID_Modificado or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AG_PID_Modificado_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AG_PID_Modificado_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AG_PID_Modificado

% Last Modified by GUIDE v2.5 09-Apr-2016 16:04:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AG_PID_Modificado_OpeningFcn, ...
                   'gui_OutputFcn',  @AG_PID_Modificado_OutputFcn, ...
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


% --- Executes just before AG_PID_Modificado is made visible.
function AG_PID_Modificado_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AG_PID_Modificado (see VARARGIN)

% Choose default command line output for AG_PID_Modificado
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AG_PID_Modificado wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AG_PID_Modificado_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edtBits_Callback(hObject, eventdata, handles)
% hObject    handle to edtBits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBits as text
%        str2double(get(hObject,'String')) returns contents of edtBits as a double


% --- Executes during object creation, after setting all properties.
function edtBits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3



function edtPopulacao_Callback(hObject, eventdata, handles)
% hObject    handle to edtPopulacao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPopulacao as text
%        str2double(get(hObject,'String')) returns contents of edtPopulacao as a double


% --- Executes during object creation, after setting all properties.
function edtPopulacao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPopulacao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtGeracoes_Callback(hObject, eventdata, handles)
% hObject    handle to edtGeracoes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtGeracoes as text
%        str2double(get(hObject,'String')) returns contents of edtGeracoes as a double


% --- Executes during object creation, after setting all properties.
function edtGeracoes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtGeracoes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtPcrate_Callback(hObject, eventdata, handles)
% hObject    handle to edtPcrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPcrate as text
%        str2double(get(hObject,'String')) returns contents of edtPcrate as a double


% --- Executes during object creation, after setting all properties.
function edtPcrate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPcrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtPmrate_Callback(hObject, eventdata, handles)
% hObject    handle to edtPmrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPmrate as text
%        str2double(get(hObject,'String')) returns contents of edtPmrate as a double


% --- Executes during object creation, after setting all properties.
function edtPmrate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPmrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pbtIniciar2.
function pbtIniciar2_Callback(hObject, eventdata, handles)
% hObject    handle to pbtIniciar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
bits=3;
warning('off','all')

tampop=str2num(get(handles.edtPopulacao2, 'String'));
geracoes=str2num(get(handles.edtGeracoes2, 'String'));
pcrate=str2num(get(handles.edtPcrate2, 'String'));
pmrate=str2num(get(handles.edtPmrate2, 'String'));
csttorneio=str2num(get(handles.edtCsttorneio2, 'String'));
RCMin=str2num(get(handles.edtRCMin, 'String'));
RCMax=str2num(get(handles.edtRCMax, 'String'));
beta=str2num(get(handles.edtBeta2, 'String'));
flagelit=get(handles.rbtElitismo2, 'Value');
flagtorneio=get(handles.rbtTorneio2, 'Value');
flagroleta=get(handles.rbtRoleta2, 'Value');
flagrank=get(handles.rbtRank2, 'Value');
flagradcliff=get(handles.rbtRadcliff2, 'Value');
flagwright=get(handles.rbtWright2, 'Value');

global dados

PMin=RCMin;
PMax=RCMax;
IMin=RCMin;
IMax=RCMax;
DMin=RCMin;
DMax=RCMax;

geratual(1)=1;
popiniccro=gerapopinicial2(tampop,RCMin,RCMax);

tempo=dados(:,3);
amostra=dados(:,2);
plot(dados(:,3),dados(:,2),'g')

fitnesspop=testepop2(popiniccro,tampop,bits,dados);
mediafitness(1)=mean(fitnesspop);

[d f]=sort(fitnesspop,'ascend');
for i=1:tampop
    g(i,:)=popiniccro(f(i),:);
end
popiniccro=g;
fitnesspop=d;

melhorindividuo=popiniccro(1,:);

melhorZ(1)=fitnesspop(1);

%-------------------Parte Gráfica------------------------------------------
axes(handles.axsGrafico1);
plot(geratual,melhorZ,'r.',geratual,mediafitness,'c.')
set(handles.txErro2, 'string', num2str(melhorZ(1)))
set(handles.edtAcalc, 'string',popiniccro(1,1).')
set(handles.edtBcalc, 'string',popiniccro(1,2).')
set(handles.edtCcalc, 'string',popiniccro(1,3).')
set(handles.txtMediafxy2, 'string', num2str(mediafitness(1)));
axes(handles.axes5);
G=tf(popiniccro(1,1),[1 popiniccro(1,2) popiniccro(1,3)]);
ft=lsim(G,dados(:,1),dados(:,3),'zoh');
hold on
graf=plot(dados(:,3),ft);
hold off
%--------------------------------------------------------------------------

fitnesspopinv=-1*fitnesspop+fitnesspop(tampop);

for j=2:geracoes
    geratual(j)=j;
    disp(popiniccro);
%------------------------------Elistismo-----------------------------------
    if flagelit==1    
        filhos(1,:)=popiniccro(1,:);
        filhos(2,:)=popiniccro(2,:);
        kelitismo=3;
    else
        kelitismo=1;
    end    
    %--------------------------------------------------------------------------

        for i=kelitismo:2:tampop
    %-------------------Torneio------------------------------------------------        

            if flagtorneio==1
                comp1=seletorneio(fitnesspopinv,tampop,csttorneio);
                comp2=seletorneio(fitnesspopinv,tampop,csttorneio);
            end
    %--------------------------------------------------------------------------         

    %----------------------------Roleta----------------------------------------
            if flagroleta==1
                comp1=seleroleta(fitnesspopinv,tampop);
                comp2=seleroleta(fitnesspopinv,tampop);
            end
    %--------------------------------------------------------------------------

    %----------------------------Rank------------------------------------------
            if flagrank==1
                comp1=selerank(fitnesspopinv,tampop);
                comp2=selerank(fitnesspopinv,tampop);

            end
    %--------------------------------------------------------------------------

    %----------------------------Crossover-------------------------------------    
            pai(1,:)=popiniccro(comp1,:);
            mae(1,:)=popiniccro(comp2,:) ;  
            
            if flagradcliff==1
                filho=crossoverradcliff(pai,mae,pcrate,bits,beta);
            end
            
            
            if flagwright==1
                filho=crossoverwright(pai,mae,pcrate,bits,PMin,PMax,IMin,IMax,DMin,DMax);
            end
            
            filho1=filho(1,:);
            filho2=filho(2,:);
          
    %--------------------------------------------------------------------------

    %----------------------------Mutação---------------------------------------
            filho1=mutacao(filho1,pmrate,bits,PMax,PMin,IMax,IMin,DMax,DMin);
            filho2=mutacao(filho2,pmrate,bits,PMax,PMin,IMax,IMin,DMax,DMin);       

            filhos(i,:)=filho1;
            filhos(i+1,:)=filho2;
    %--------------------------------------------------------------------------        
        end

        popiniccro=filhos;
        fitnesspop=testepop2(popiniccro,tampop,bits,dados);
        mediafitness(j)=mean(fitnesspop);
        [d f]=sort(fitnesspop,'ascend');
        for i=1:tampop
            g(i,:)=popiniccro(f(i),:);
        end
        popiniccro=g;
        fitnesspop=d;
        
        
        melhorindividuo=popiniccro(1,:);
        melhorZ(j)=fitnesspop(1);
        fitnesspopinv=-1*fitnesspop+fitnesspop(tampop);
        axes(handles.axsGrafico1);
        plot(geratual,melhorZ,'r.',geratual,mediafitness,'c.')
        
        set(handles.txErro2, 'string', num2str(melhorZ(j)))
        set(handles.edtAcalc, 'string',popiniccro(1,1).')
        set(handles.edtBcalc, 'string',popiniccro(1,2).')
        set(handles.edtCcalc, 'string',popiniccro(1,3).')
        set(handles.txtMediafxy2, 'string', num2str(mediafitness(j)));
        
        axes(handles.axes5);
        G=tf(popiniccro(1,1),[1 popiniccro(1,2) popiniccro(1,3)]);
        ft=lsim(G,dados(:,1),dados(:,3),'zoh');
        delete(graf)
        hold on
        graf=plot(dados(:,3),ft);
        hold off  
end
        warning('on','all')
        axes(handles.axsGrafico1);
        legend('Melhor Indivíduo','Média da População');
        axes(handles.axes5);
        legend('Sistema Original','Sistema Identificado','Location','SouthEast');

        
function  popinic2=gerapopinicial2(a,b,c)
%a=tampop, b=Pmin, c=Pmax
popinic2(:,1)=rand(a,1)*(c-b)+b;
popinic2(:,2)=rand(a,1)*(c-b)+b;
popinic2(:,3)=rand(a,1)*(c-b)+b;

function  fx2=testepop2(a,b,c,d)
%a=população, b=tampop, c=bits, d=dados
for i=1:b
    if sum(a(i,:))==0
        a(i,:)=a(i,:)+0.1;
    end
    xo=d(1,2);
    G=tf(a(i,1),[1 a(i,2) a(i,3)]);
    ft=lsim(G,d(:,1),d(:,3),xo,'zoh');
    erro=ft-d(:,2);
    erroquad=erro.*erro;
    testeNAN = isfinite(erroquad);
    for cont=1:length(erroquad)
        if testeNAN(cont)==0
           erroquad(cont) = 1e300;
        end 
    end
    fx2(i)=sum(erroquad)/length(d);
    fx2(i)=sqrt(fx2(i));
end

% --- Executes on button press in pbtIniciar.
function pbtIniciar_Callback(hObject, eventdata, handles)
% hObject    handle to pbtIniciar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc;
bits=3;
tampop=str2num(get(handles.edtPopulacao, 'String'));
geracoes=str2num(get(handles.edtGeracoes, 'String'));
pcrate=str2num(get(handles.edtPcrate, 'String'));
pmrate=str2num(get(handles.edtPmrate, 'String'));
csttorneio=str2num(get(handles.edtCsttorneio, 'String'));
PMin=str2num(get(handles.edtPMin, 'String'));
PMax=str2num(get(handles.edtPMax, 'String'));
IMin=str2num(get(handles.edtIMin, 'String'));
IMax=str2num(get(handles.edtIMax, 'String'));
DMin=str2num(get(handles.edtDMin, 'String'));
DMax=str2num(get(handles.edtDMax, 'String'));
beta=str2num(get(handles.edtBeta, 'String'));
flagelit=get(handles.rbtElitismo, 'Value');
flagtorneio=get(handles.rbtTorneio, 'Value');
flagroleta=get(handles.rbtRoleta, 'Value');
flagrank=get(handles.rbtRank, 'Value');
flagradcliff=get(handles.rbtRadcliff, 'Value');
flagwright=get(handles.rbtWright, 'Value');
over=str2num(get(handles.edtOver, 'String'));
AA=str2num(get(handles.edtAcalc, 'String'));
BB=str2num(get(handles.edtBcalc, 'String'));
CC=str2num(get(handles.edtCcalc, 'String'));
Kganho=str2num(get(handles.edtKPganho, 'String'));
tal=str2num(get(handles.edtTal, 'String'));
numamostra=str2num(get(handles.edtAmostra, 'String'));
tempamostra=str2num(get(handles.edtTempo, 'String'));

min=PMin;
max=PMax;

if PMin==0  
PMin=0.0001;
end

geratual(1)=1;
popiniccro=gerapopinicial(tampop,PMin,PMax,IMin,IMax,DMin,DMax);

H=tf(AA,[1 BB CC]);

T=linspace(0,tempamostra,numamostra);
GPad=tf(Kganho,[tal 1]);
%GPad=feedback(GPad,1);
Gamostra=step(GPad,T);

axes(handles.axes4);
plot(T,Gamostra,'g')
axis([0 tempamostra 0 over]);

fitnesspop=testepop(popiniccro,tampop,bits,H,T,over,Gamostra);
mediafitness(1)=mean(fitnesspop);

[d f]=sort(fitnesspop,'ascend');
for i=1:tampop
    g(i,:)=popiniccro(f(i),:);
end
popiniccro=g;
fitnesspop=d;

melhorindividuo=popiniccro(1,:);

melhorZ(1)=fitnesspop(1);

%-------------------Parte Gráfica------------------------------------------
axes(handles.axsGrafico1);
plot(geratual,melhorZ,'r.',geratual,mediafitness,'c.')
set(handles.txErro, 'string', num2str(melhorZ(1)))
set(handles.edtPcalc, 'string',popiniccro(1,1).')
set(handles.edtIcalc, 'string',popiniccro(1,2).')
set(handles.edtDcalc, 'string',popiniccro(1,3).')
set(handles.txtMediafxy, 'string', num2str(mediafitness(1)));

axes(handles.axes4);
GPID=tf([popiniccro(1,3) popiniccro(1,1) popiniccro(1,2)],[1 0]);
HLC2=feedback(GPID*H,1); %ganho em malha fechada com realimentação unitária
[tempy]=step(HLC2,T);
hold on
graf=plot(T,tempy);
hold off
%--------------------------------------------------------------------------

fitnesspopinv=-1*fitnesspop+fitnesspop(tampop);

for j=2:geracoes
    geratual(j)=j;
  
%------------------------------Elistismo-----------------------------------
    if flagelit==1    
        filhos(1,:)=popiniccro(1,:);
        filhos(2,:)=popiniccro(2,:);
        kelitismo=3;
    else
        kelitismo=1;
    end    
    %--------------------------------------------------------------------------

        for i=kelitismo:2:tampop
    %-------------------Torneio------------------------------------------------        

            if flagtorneio==1
                comp1=seletorneio(fitnesspopinv,tampop,csttorneio);
                comp2=seletorneio(fitnesspopinv,tampop,csttorneio);
            end
    %--------------------------------------------------------------------------         

    %----------------------------Roleta----------------------------------------
            if flagroleta==1
                comp1=seleroleta(fitnesspopinv,tampop);
                comp2=seleroleta(fitnesspopinv,tampop);
            end
    %--------------------------------------------------------------------------

    %----------------------------Rank------------------------------------------
            if flagrank==1
                comp1=selerank(fitnesspopinv,tampop);
                comp2=selerank(fitnesspopinv,tampop);

            end
    %--------------------------------------------------------------------------

    %----------------------------Crossover-------------------------------------    
            pai(1,:)=popiniccro(comp1,:);
            mae(1,:)=popiniccro(comp2,:) ;  
            
            if flagradcliff==1
                filho=crossoverradcliff(pai,mae,pcrate,bits,beta);
            end
            
            
            if flagwright==1
                filho=crossoverwright(pai,mae,pcrate,bits,PMin,PMax,IMin,IMax,DMin,DMax);
            end
            
            filho1=filho(1,:);
            filho2=filho(2,:);
    %--------------------------------------------------------------------------

    %----------------------------Mutação---------------------------------------
            filho1=mutacao(filho1,pmrate,bits,PMax,PMin,IMax,IMin,DMax,DMin);
            filho2=mutacao(filho2,pmrate,bits,PMax,PMin,IMax,IMin,DMax,DMin);       

            filhos(i,:)=filho1;
            filhos(i+1,:)=filho2;
    %--------------------------------------------------------------------------        
        end

        popiniccro=filhos;
        fitnesspop=testepop(popiniccro,tampop,bits,H,T,over,Gamostra);
        mediafitness(j)=mean(fitnesspop);
        [d f]=sort(fitnesspop,'ascend');
        for i=1:tampop
            g(i,:)=popiniccro(f(i),:);
        end
        popiniccro=g;
        fitnesspop=d;
                
        melhorindividuo=popiniccro(1,:);
        melhorZ(j)=fitnesspop(1);
        fitnesspopinv=-1*fitnesspop+fitnesspop(tampop);
        axes(handles.axsGrafico1);
        plot(geratual,melhorZ,'r.',geratual,mediafitness,'c.')
        
        set(handles.txErro, 'string', num2str(melhorZ(j)))
        set(handles.txtMediafxy, 'string', num2str(mediafitness(j)))
        set(handles.edtPcalc, 'string',popiniccro(1,1).')
        set(handles.edtIcalc, 'string',popiniccro(1,2).')
        set(handles.edtDcalc, 'string',popiniccro(1,3).')
               
        axes(handles.axes4);
        GPID=tf([popiniccro(1,3) popiniccro(1,1) popiniccro(1,2)],[1 0]);
        HLC2=feedback(GPID*H,1);%ganho em malha fechada com realimentação unitária
        delete(graf);
        [tempy]=step(HLC2,T);
        hold on
        graf=plot(T,tempy);
        hold off
    
end
        axes(handles.axsGrafico1);
        legend('Melhor Indivíduo','Média da População');
        axes(handles.axes4);
        legend('Resposta Desejada','Resposta Obtida','Location','SouthEast');
        
function  popinic=gerapopinicial(a,b,c,d,e,f,g,h)
%a=tampop, b=Pmin, c=Pmax, d=Imin, e=Imax, f=Dmin, g=Dmax
popinic(:,1)=rand(a,1)*(c-b)+b;
popinic(:,2)=rand(a,1)*(e-d)+d;
popinic(:,3)=rand(a,1)*(g-f)+f;

function  fx=testepop(a,b,c,d,e,f,g)
%a=população, b=tampop, c=bits, d=Planta, e=vetor tempo, f=overshut,
%g=amostra

for i=1:b 
    GPID=tf([a(i,3) a(i,1) a(i,2)],[1 0]);
    HLC1=feedback(GPID*d,1); %malha fechada com realimentação unitária
    ft=step(HLC1,e);
    if max(ft)>f
        ft=ft+2; %punição individuo overshoot
    end
    erro=ft-g;
    erroquad=erro.*erro;%erro
    fx(i)=sum(erroquad)/length(e);   
    fx(i)=sqrt(fx(i));
end

function m=melhorindiv(a,b)
%a=fitnesspop
%b=população
pmelhor=0;
smelhor=0;
for i=1:b    
    if(a(i)>=pmelhor)
        pmelhor=a(i);
        m(1)=i;
    end
end
for i=1:b
    if(a(i)>=smelhor)
        if(a(i)<pmelhor)
            smelhor=a(i);
            m(2)=i;
        end
    end
end

function media=mediafit(a,b)
%a=fitnesspop
%b=população
soma=0;
for i=1:b    
    soma=soma+a(i);
end
media=soma/b;

function t=seletorneio(a,b,d)
%a=fitnesspop
%b=população
%d=individuos do torneio
ind(1)=randi([1,b],1);
for i=1:d-1 
    c=0; 
    while c==0
        sorteado=randi([1,b],1);
        flag=0;
        for j=1:i
            if sorteado == ind(j)
                flag=1;
            end
        end
        if flag == 1
            c=0;
        else
            c=1;
            ind(i+1)=sorteado;
        end 
    end
end
compet=0;
for i=1:d
    if a(ind(i))>=compet
        compet=a(ind(i));
        t=ind(i);
    end
end

function muta=mutacao(a,b,c,d,e,f,g,h,ii)
%a=filho1,b=pmrate,c=bits,d=pmax,e=pmin,f=Imax,g=Imin,h=Dmax,ii=Dmin
pm=rand;
if pm<=b
    jj=randi(c);
    if jj==1
        a(1)=rand()*(d-e)+e; 
    end
    if jj==2
        a(2)=rand()*(f-g)+g; 
    end
    if jj==3
        a(3)=rand()*(h-ii)+ii; 
    end
    if jj>3
        a(jj)=rand()*(d-e)+e; 
    end
    
end
muta=a;

function indroleta=seleroleta(fit,tam)
p_i=(1/sum(fit))*fit;
ci(1)=p_i(1);
for k=2:tam
    ci(k)=ci(k-1)+p_i(k);   
end
apont=rand;
for k=1:tam
    if apont<ci(k)
        indroleta=k;
        break;       
    end
end

function indrank=selerank(fit,tam)
[xcrescente,indzdec]=sort(fit,'descend');
min=0.5;
max=1.5;
for k=1:tam
    p_i(k)= min+(max-min)*(tam-k)/(tam-1);
end
p_i=p_i/sum(p_i);
ci(1)=p_i(1);
for k=2:tam
    ci(k)=ci(k-1)+p_i(k);   
end
apont=rand;
for k=1:tam
    if apont<ci(k)
        indice=k;
        break;       
    end
end
indrank=indzdec(indice);

function filhos=crossoverradcliff(a,b,c,d,e)
%a=pai, b=mae, c=pcrate, d=bits, e=beta
for i=1:d
    filhos(1,i)=e*a(i)+(1-e)*b(i);
end
for i=1:d
    filhos(2,i)=(1-e)*a(i)+e*b(i);
end

function filhos=crossoverwright(pai,mae,taxa,a,pmi,pma,imi,ima,dmi,dma)
pc=rand;
if pc<=taxa
    for k=1:a
        f1(1,k)=0.5*pai(k)+0.5*mae(k);
        f1(2,k)=1.5*pai(k)-0.5*mae(k);
        f1(3,k)=1.5*mae(k)-0.5*pai(k);
        
        if k==1
            for l=1:3
            if f1(l,k)>pma
                f1(l,k)=pma;
            end
            if f1(l,k)<pmi
                f1(l,k)=pmi;
            end
            end
        end
        if k==2
            for l=1:3
            if f1(l,k)>ima
                f1(l,k)=ima;
            end
            if f1(l,k)<imi
                f1(l,k)=imi;
            end
            end
        end
        if k==3
            for l=1:3
            if f1(l,k)>dma
                f1(l,k)=dma;
            end
            if f1(l,k)<dmi
                f1(l,k)=dmi;
            end
            end
        end
        if k>3
            for l=1:3
            if f1(l,k)>pma
                f1(l,k)=pma;
            end
            if f1(l,k)<pmi
                f1(l,k)=pmi;
            end
            end
        end 
              
    end
    corte=randi([1,3],1);
    g=1;
    for h=1:3
        if h==corte
        else 
            filhos(g,:)=f1(h,:);
            g=g+1;
        end
    end
else
    filhos(1,:)=pai;
    filhos(2,:)=mae;
end


function edtCsttorneio_Callback(hObject, eventdata, handles)
% hObject    handle to edtCsttorneio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtCsttorneio as text
%        str2double(get(hObject,'String')) returns contents of edtCsttorneio as a double


% --- Executes during object creation, after setting all properties.
function edtCsttorneio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtCsttorneio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtTemptorneio_Callback(hObject, eventdata, handles)
% hObject    handle to edtTemptorneio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTemptorneio as text
%        str2double(get(hObject,'String')) returns contents of edtTemptorneio as a double


% --- Executes during object creation, after setting all properties.
function edtTemptorneio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTemptorneio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtDelay_Callback(hObject, eventdata, handles)
% hObject    handle to edtDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtDelay as text
%        str2double(get(hObject,'String')) returns contents of edtDelay as a double


% --- Executes during object creation, after setting all properties.
function edtDelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbtElitismo.
function rbtElitismo_Callback(hObject, eventdata, handles)
% hObject    handle to rbtElitismo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtElitismo


% --- Executes on button press in rbtTorneio.
function rbtTorneio_Callback(hObject, eventdata, handles)
% hObject    handle to rbtTorneio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtTorneio
set(handles.rbtTorneio,'value',1);
set(handles.rbtRank,'value',0); 
set(handles.rbtRoleta,'value',0); 
set(handles.edtCsttorneio,'visible','on');



% --- Executes on button press in rbtRank.
function rbtRank_Callback(hObject, eventdata, handles)
% hObject    handle to rbtRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtRank
set(handles.rbtTorneio,'value',0);
set(handles.rbtRank,'value',1); 
set(handles.rbtRoleta,'value',0);
set(handles.edtCsttorneio,'visible','off');
% --- Executes on button press in rbtRoleta.
function rbtRoleta_Callback(hObject, eventdata, handles)
% hObject    handle to rbtRoleta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtRoleta
set(handles.rbtTorneio,'value',0);
set(handles.rbtRank,'value',0); 
set(handles.rbtRoleta,'value',1);
set(handles.edtCsttorneio,'visible','off');



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtMiny_Callback(hObject, eventdata, handles)
% hObject    handle to edtMiny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtMiny as text
%        str2double(get(hObject,'String')) returns contents of edtMiny as a double


% --- Executes during object creation, after setting all properties.
function edtMiny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtMiny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtPMax_Callback(hObject, eventdata, handles)
% hObject    handle to edtPMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPMax as text
%        str2double(get(hObject,'String')) returns contents of edtPMax as a double


% --- Executes during object creation, after setting all properties.
function edtPMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtMaxy_Callback(hObject, eventdata, handles)
% hObject    handle to edtMaxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtMaxy as text
%        str2double(get(hObject,'String')) returns contents of edtMaxy as a double


% --- Executes during object creation, after setting all properties.
function edtMaxy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtMaxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtPMin_Callback(hObject, eventdata, handles)
% hObject    handle to edtPMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPMin as text
%        str2double(get(hObject,'String')) returns contents of edtPMin as a double


% --- Executes during object creation, after setting all properties.
function edtPMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbtRadcliff.
function rbtRadcliff_Callback(hObject, eventdata, handles)
% hObject    handle to rbtRadcliff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtRadcliff
set(handles.rbtWright,'value',0);
set(handles.rbtRadcliff,'value',1); 

% --- Executes on button press in rbtWright.
function rbtWright_Callback(hObject, eventdata, handles)
% hObject    handle to rbtWright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtWright
set(handles.rbtWright,'value',1);
set(handles.rbtRadcliff,'value',0);

% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton12



function edtBeta_Callback(hObject, eventdata, handles)
% hObject    handle to edtBeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBeta as text
%        str2double(get(hObject,'String')) returns contents of edtBeta as a double


% --- Executes during object creation, after setting all properties.
function edtBeta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ltbX.
function ltbX_Callback(hObject, eventdata, handles)
% hObject    handle to ltbX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ltbX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ltbX


% --- Executes during object creation, after setting all properties.
function ltbX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ltbX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtIMax_Callback(hObject, eventdata, handles)
% hObject    handle to edtIMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtIMax as text
%        str2double(get(hObject,'String')) returns contents of edtIMax as a double


% --- Executes during object creation, after setting all properties.
function edtIMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtIMin_Callback(hObject, eventdata, handles)
% hObject    handle to edtIMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtIMin as text
%        str2double(get(hObject,'String')) returns contents of edtIMin as a double


% --- Executes during object creation, after setting all properties.
function edtIMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtDMax_Callback(hObject, eventdata, handles)
% hObject    handle to edtDMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtDMax as text
%        str2double(get(hObject,'String')) returns contents of edtDMax as a double


% --- Executes during object creation, after setting all properties.
function edtDMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtDMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtDMin_Callback(hObject, eventdata, handles)
% hObject    handle to edtDMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtDMin as text
%        str2double(get(hObject,'String')) returns contents of edtDMin as a double


% --- Executes during object creation, after setting all properties.
function edtDMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtDMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtOver_Callback(hObject, eventdata, handles)
% hObject    handle to edtOver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtOver as text
%        str2double(get(hObject,'String')) returns contents of edtOver as a double


% --- Executes during object creation, after setting all properties.
function edtOver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtOver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbtRoleta2.
function rbtRoleta2_Callback(hObject, eventdata, handles)
% hObject    handle to rbtRoleta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtRoleta2
set(handles.rbtTorneio2,'value',0);
set(handles.rbtRank2,'value',0); 
set(handles.rbtRoleta2,'value',1);
set(handles.edtCsttorneio2,'visible','off');

% --- Executes on button press in radiobutton16.
function radiobutton16_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton16


% --- Executes on button press in radiobutton17.
function radiobutton17_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton17



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton18.
function radiobutton18_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton18


% --- Executes on button press in radiobutton19.
function radiobutton19_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton19



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtPcalc_Callback(hObject, eventdata, handles)
% hObject    handle to edtPcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPcalc as text
%        str2double(get(hObject,'String')) returns contents of edtPcalc as a double


% --- Executes during object creation, after setting all properties.
function edtPcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtDcalc_Callback(hObject, eventdata, handles)
% hObject    handle to edtDcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtDcalc as text
%        str2double(get(hObject,'String')) returns contents of edtDcalc as a double


% --- Executes during object creation, after setting all properties.
function edtDcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtDcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtB_Callback(hObject, eventdata, handles)
% hObject    handle to edtB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtB as text
%        str2double(get(hObject,'String')) returns contents of edtB as a double


% --- Executes during object creation, after setting all properties.
function edtB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtBeta2_Callback(hObject, eventdata, handles)
% hObject    handle to edtBeta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBeta2 as text
%        str2double(get(hObject,'String')) returns contents of edtBeta2 as a double


% --- Executes during object creation, after setting all properties.
function edtBeta2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBeta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton14.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton14



function edtPopulacao2_Callback(hObject, eventdata, handles)
% hObject    handle to edtPopulacao2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPopulacao2 as text
%        str2double(get(hObject,'String')) returns contents of edtPopulacao2 as a double


% --- Executes on button press in rbtTorneio2.
function rbtTorneio2_Callback(hObject, eventdata, handles)
% hObject    handle to rbtTorneio2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtTorneio2
set(handles.rbtTorneio2,'value',1);
set(handles.rbtRank2,'value',0); 
set(handles.rbtRoleta2,'value',0); 
set(handles.edtCsttorneio2,'visible','on');

% --- Executes on button press in rbtRank2.
function rbtRank2_Callback(hObject, eventdata, handles)
% hObject    handle to rbtRank2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtRank2
set(handles.rbtTorneio2,'value',0);
set(handles.rbtRank2,'value',1); 
set(handles.rbtRoleta2,'value',0);
set(handles.edtCsttorneio2,'visible','off');


% --- Executes during object creation, after setting all properties.
function edtRCMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRCMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtRCMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRCMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtPopulacao2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPopulacao2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtGeracoes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtGeracoes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtPcrate2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPcrate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtPmrate2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPmrate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtCsttorneio2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtCsttorneio2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtBcalc_Callback(hObject, eventdata, handles)
% hObject    handle to edtBcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBcalc as text
%        str2double(get(hObject,'String')) returns contents of edtBcalc as a double


% --- Executes during object creation, after setting all properties.
function edtBcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtC_Callback(hObject, eventdata, handles)
% hObject    handle to edtC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtC as text
%        str2double(get(hObject,'String')) returns contents of edtC as a double


% --- Executes during object creation, after setting all properties.
function edtC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtA_Callback(hObject, eventdata, handles)
% hObject    handle to edtA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtA as text
%        str2double(get(hObject,'String')) returns contents of edtA as a double


% --- Executes during object creation, after setting all properties.
function edtA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtAcalc_Callback(hObject, eventdata, handles)
% hObject    handle to edtAcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtAcalc as text
%        str2double(get(hObject,'String')) returns contents of edtAcalc as a double


% --- Executes during object creation, after setting all properties.
function edtAcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtAcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtCcalc_Callback(hObject, eventdata, handles)
% hObject    handle to edtCcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtCcalc as text
%        str2double(get(hObject,'String')) returns contents of edtCcalc as a double


% --- Executes during object creation, after setting all properties.
function edtCcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtCcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edtIcalc_Callback(hObject, eventdata, handles)
% hObject    handle to edtIcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtIcalc as text
%        str2double(get(hObject,'String')) returns contents of edtIcalc as a double


% --- Executes during object creation, after setting all properties.
function edtIcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbtRadcliff2.
function rbtRadcliff2_Callback(hObject, eventdata, handles)
% hObject    handle to rbtRadcliff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtRadcliff2
set(handles.rbtWright2,'value',0);
set(handles.rbtRadcliff2,'value',1); 


% --- Executes on button press in rbtWright2.
function rbtWright2_Callback(hObject, eventdata, handles)
% hObject    handle to rbtWright2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtWright2
set(handles.rbtWright2,'value',1);
set(handles.rbtRadcliff2,'value',0); 


function edtDDcalc_Callback(hObject, eventdata, handles)
% hObject    handle to edtDDcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtDDcalc as text
%        str2double(get(hObject,'String')) returns contents of edtDDcalc as a double


% --- Executes during object creation, after setting all properties.
function edtDDcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtDDcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtKPganho_Callback(hObject, eventdata, handles)
% hObject    handle to edtKPganho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtKPganho as text
%        str2double(get(hObject,'String')) returns contents of edtKPganho as a double


% --- Executes during object creation, after setting all properties.
function edtKPganho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtKPganho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtTal_Callback(hObject, eventdata, handles)
% hObject    handle to edtTal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTal as text
%        str2double(get(hObject,'String')) returns contents of edtTal as a double


% --- Executes during object creation, after setting all properties.
function edtTal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbtImportarDados.
function pbtImportarDados_Callback(hObject, eventdata, handles)
% hObject    handle to pbtImportarDados (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dados
arquivo = uigetfile('*.csv','Select the CSV code file');
dados=csvread(arquivo); 
% dados=csvread('dados.csv');
axes(handles.axes5);
plot(dados(:,3),dados(:,2),'g')



function edtAmostra_Callback(hObject, eventdata, handles)
% hObject    handle to edtAmostra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtAmostra as text
%        str2double(get(hObject,'String')) returns contents of edtAmostra as a double


% --- Executes during object creation, after setting all properties.
function edtAmostra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtAmostra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtTempo_Callback(hObject, eventdata, handles)
% hObject    handle to edtTempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTempo as text
%        str2double(get(hObject,'String')) returns contents of edtTempo as a double


% --- Executes during object creation, after setting all properties.
function edtTempo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbtResposta.
function pbtResposta_Callback(hObject, eventdata, handles)
% hObject    handle to pbtResposta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Kganho=str2num(get(handles.edtKPganho, 'String'));
tal=str2num(get(handles.edtTal, 'String'));
numamostra=str2num(get(handles.edtAmostra, 'String'));
tempamostra=str2num(get(handles.edtTempo, 'String'));
over=str2num(get(handles.edtOver, 'String'));

T=linspace(0,tempamostra,numamostra);
GPad=tf(Kganho,[tal 1]);
Gamostra=step(GPad,T);
axes(handles.axes4);
plot(T,Gamostra,'g')
axis([0 tempamostra 0 over]);
