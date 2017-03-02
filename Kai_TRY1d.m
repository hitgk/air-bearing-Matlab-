function varargout = Kai_TRY1d(varargin)
% KAI_TRY1D MATLAB code for Kai_TRY1d.fig
%      KAI_TRY1D, by itself, creates a new KAI_TRY1D or raises the existing
%      singleton*.
%
%      H = KAI_TRY1D returns the handle to a new KAI_TRY1D or the handle to
%      the existing singleton*.
%
%      KAI_TRY1D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KAI_TRY1D.M with the given input arguments.
%
%      KAI_TRY1D('Property','Value',...) creates a new KAI_TRY1D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Kai_TRY1d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Kai_TRY1d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Kai_TRY1d

% Last Modified by GUIDE v2.5 05-May-2013 10:07:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Kai_TRY1d_OpeningFcn, ...
                   'gui_OutputFcn',  @Kai_TRY1d_OutputFcn, ...
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


% --- Executes just before Kai_TRY1d is made visible.
function Kai_TRY1d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Kai_TRY1d (see VARARGIN)

% Choose default command line output for Kai_TRY1d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Kai_TRY1d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Kai_TRY1d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function i_rb_Callback(hObject, eventdata, handles)
% hObject    handle to i_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_rb as text
%        str2double(get(hObject,'String')) returns contents of i_rb as a double


% --- Executes during object creation, after setting all properties.
function i_rb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_rb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function i_ra_Callback(hObject, eventdata, handles)
% hObject    handle to i_ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_ra as text
%        str2double(get(hObject,'String')) returns contents of i_ra as a double


% --- Executes during object creation, after setting all properties.
function i_ra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function i_r_Callback(hObject, eventdata, handles)
% hObject    handle to i_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_r as text
%        str2double(get(hObject,'String')) returns contents of i_r as a double


% --- Executes during object creation, after setting all properties.
function i_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function i_n_Callback(hObject, eventdata, handles)
% hObject    handle to i_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_n as text
%        str2double(get(hObject,'String')) returns contents of i_n as a double


% --- Executes during object creation, after setting all properties.
function i_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in r_1.
function r_1_Callback(hObject, eventdata, handles)
% hObject    handle to r_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%��������Ԫ������ֹ����еľ�̬����
%%%���ա���ѹ�����󻬡����ӣ����������ئȷ��򻮷�Ϊ8�ȷ֣���r������ra��r0���仮��Ϊ8�ȷ֣���r0��rb���仮��Ϊ8�ȷ�
n=str2num(get(handles.i_n,'String'));%����ֹ����н����׸���
Pa=1.01325e5;   %һ������ѹPa
P0=str2num(get(handles.i_P,'String'));   %����ѹ��6������ѹ
fa=(Pa/P0)^2;   %1����ѹѹ��

h0=0.010;%��Ĥ���0.02mm�����ڴ�ֵ��������Ĥ��ȣ�Ҳ�ǵ���ƫ����
ra=str2num(get(handles.i_ra,'String'));%�ھ�7mmi_ra
rb=str2num(get(handles.i_rb,'String'));% i_rb->;%�⾶14mm
r0=str2num(get(handles.i_r0,'String'));%����Բ�뾶10mm
ras=ra/rb;%ra����ΰ뾶��rbΪ�ο�����
rbs=rb/rb;%rb����ΰ뾶

hs=1;%�������Ĥ���   hs=h0/h0

deltatheta=pi/(n*8);%����
deltars=(rbs-ras)/16;%ra��rb����εȷ�ֵ

rs=zeros(1,15);
for i=1:15
    rs(i)=ras+i*deltars;%����r1��r15����ΰ뾶
end

%%%%%%%%%�任��=lnr
xisa=log(ras);%ras�任
xisb=log(rbs);%rbs�任
xis=zeros(1,15);
for i=1:15
    xis(i)=log(rs(i));%����r1��r15�任
end

%%%%%%%���㦤��(i)s
deltaxis=zeros(1,16);
deltaxis(1)=xis(1)-xisa;%����1
for i=2:15
    deltaxis(i)=xis(i)-xis(i-1);%����2������15
end
deltaxis(16)=xisb-xis(15);%����16

%%%����G(i��
G=zeros(1,16);
for i=1:16
    G(i)=(hs^3)/(2*deltatheta*deltaxis(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����նȾ���K
A=zeros(1,16);
B=zeros(1,16);
C=zeros(1,16);
for i=1:16
    A(i)=G(i)*((deltaxis(i))^2+deltatheta^2);
    B(i)=G(i)*(deltaxis(i))^2;
    C(i)=G(i)*deltatheta^2;
end

aa=zeros(15,15);%�Ӿ����ʼֵ
k=cell(9);%9x9�ֿ����

for i=1:9
    for j=1:9
        k{i,j}=aa;%�ֿ���󸳳�ֵ
    end
end

for i=1:15           %k11��k99
    k{1,1}(i,i)=A(i)+A(i+1);
    k{9,9}(i,i)=A(i)+A(i+1);
    if i<15
        k{1,1}(i,i+1)=-C(i+1);
        k{1,1}(i+1,i)=-C(i+1);
        
        k{9,9}(i,i+1)=-C(i+1);
        k{9,9}(i+1,i)=-C(i+1);
    else
        break
    end
end
        
 for j=2:8       %k22,k33,...,k88
     for i=1:15
         k{j,j}(i,i)=2*(A(i)+A(i+1));
         if i<15
             k{j,j}(i,i+1)=-2*C(i+1);
             k{j,j}(i+1,i)=-2*C(i+1);
         else
             break
         end
     end
 end
             
 for j=1:8     %k12,k21,k23,k32,...,k89,k98
     for i=1:15
         k{j,j+1}(i,i)=-B(i)-B(i+1);
         k{j+1,j}(i,i)=-B(i)-B(i+1);
     end
 end

K=cell2mat(k);%�ϲ���һ������

%%%%�������T
T=zeros(135,1);

T(1)=G(1)*(deltatheta^2)*fa;
T(121)=T(1);
T(15)=G(16)*(deltatheta^2)*fa;
T(135)=T(15);
T(16)=2*G(1)*(deltatheta^2)*fa;
T(31)=T(16);
T(46)=T(16);
T(61)=T(16);
T(76)=T(16);
T(91)=T(16);
T(106)=T(16);
T(30)=2*G(16)*(deltatheta^2)*fa;
T(45)=T(30);
T(60)=T(30);
T(75)=T(30);
T(90)=T(30);
T(105)=T(30);
T(120)=T(30);

%%%����T(8)

rhoa=1.205;%20���϶ȿ����ܶ�kg/m^3
eta=1.79e-5;%20���϶ȿ�������ճ��Pa.s
d=str2num(get(handles.i_r,'String'));%������ֱ��0.2mm
phi=0.8;%����ϵ��
coT8=3*eta*pi*((d/1000)^2)*phi/((h0/1000)^3*P0)*sqrt(2*Pa/rhoa);%t8�ļ���ʽϵ��t8=coT8*��8

ka=1.401;%�������峣��
Bk=(2/(ka+1))^(ka/(ka-1));%��k

psi1=sqrt(ka/2*((2/(ka+1))^((ka+1)/(ka-1))));%��8��1��

%%%%%%%%����F��ʼֵ
%%%��ڵ�1��15�ĳ�ʼ����ѹ�������Էֲ���������׺��ѹ������γ�ʼֵP8=0.7��
%�ذ뾶����������ڵ��ѹ�����Ե��½��������߽�ڵ㴦��Pa=1/6.������뾶�ϵĶ�Ӧ�ڵ�����ʼѹ����ڵ�1��15��ͬ��
pas=Pa/P0;
ps=zeros(1,15);
for i=1:15
    if i<=8
        ps(i)=pas+i*(0.7-pas)/8;
    else
        ps(i)=pas+(16-i)*(0.7-pas)/8;
    end
end

fs=zeros(1,15);
for i=1:15
    fs(i)=ps(i)^2;
end
    
F=zeros(135,1);%F��ʼֵ
for i=1:135
    if i<=15
        F(i)=fs(i);
    elseif i<=30
        F(i)=fs(i-15);
    elseif i<=45
        F(i)=fs(i-30);
    elseif i<=60
        F(i)=fs(i-45);
    elseif i<=75
        F(i)=fs(i-60);
    elseif i<=90
        F(i)=fs(i-75);
    elseif i<=105
        F(i)=fs(i-90);
    elseif i<=120
        F(i)=fs(i-105);        
    else
        F(i)=fs(i-120);
    end
end


%%%%%%%%%ѭ������

ERROR=1;
while ERROR>=1e-6
    
    F8=F(8); %F8���F(8)���ɳڵ�����ʼֵ
    FB=F;    %FB���F���ɳڵ�����ʼֵ
    TOL=1;   %ǰ�����ε������    
    while TOL>=1e-6      %���ɳڵ���SOR
        if F(8)^(1/2)<=Bk       %%%%%%%%����T(8)
            T(8)=coT8*psi1;
        else
            psi2=sqrt(ka/(ka-1)*((F(8)^(1/2))^(2/ka)-(F(8)^(1/2))^((ka+1)/ka)));
            T(8)=coT8*psi2;
        end
        omega=0.1;%�ɳ�����
        for i=1:135   
            M=0;
            N=0;
            if i<=1
                M=0;
            else
                for j=1:(i-1)
                    M=M+K(i,j)*F(j);
                end
            end
            for j=i:135
                N=N+K(i,j)*F(j);
            end
            F(i)=F(i)+omega/K(i,i)*(T(i)-M-N);
        end
        TOL=norm(F-FB);%ǰ������������
        FB=F;        
    end
    
    ERROR=abs((F(8)-F8)/F8);%��Ծ���
  
    if F(8)^(1/2)<=Bk       %%%%%%%%����T(8)
        T(8)=coT8*psi1;
    else
        psi2=sqrt(ka/(ka-1)*((F(8)^(1/2))^(2/ka)-(F(8)^(1/2))^((ka+1)/ka)));
        T(8)=coT8*psi2;
    end
    
    PRC=1+T(8)/(2*ka*K(8,8)*F(8))*abs(1-(ka-1)/2/(F(8)^(-(ka-1)/(2*ka))-1));   %�����ָ�����
    alpha=1.3;%����
    F(8)=1/(PRC*alpha)*(F(8)-F8)+F8;
end
%%%%�����������ս��Fֵ����FB��

trueps=zeros(135,1);
for i=1:135
    trueps(i)=sqrt(FB(i));  %������ڵ������ѹ��
end

%%%%�����������ϵ��
    %%%%%%�������Ԫ���
elementarea=zeros(1,16);%��Ԫ�������16��
for i=1:16
    elementarea(i)=1/2*deltatheta*deltaxis(i);
end

 %%%%%%��Ԫ1��33��65��97��129��161��193��225��32��64��96��128��160��192��224��256���ڵ�����ͬѹ�������൥Ԫ�ڵ�ѹ������ͬ
CW=zeros(1,256);%����256����Ԫ

for i=0:7
    CW(1+32*i)=8*elementarea(1)/15*(3/2*fa^(5/2)-5/2*FB(1+15*i)*fa^(3/2)+FB(1+15*i)^(5/2))/(fa-FB(1+15*i))^2;%��Ԫ1��33��65��97��129��161��193��225
    CW(2+32*i)=8*elementarea(1)/15*((FB(1+15*i)^(5/2))/((FB(1+15*i)-FB(16+15*i))*(FB(1+15*i)-fa))+(FB(16+15*i)^(5/2))/((FB(16+15*i)-FB(1+15*i))*(FB(16+15*i)-fa))+(fa^(5/2))/((fa-FB(1+15*i))*(fa-FB(16+15*i))));%��Ԫ2,34,66,98,130,162,194,226
end

for i=1:8
    CW(32*i)=8*elementarea(16)/15*(3/2*fa^(5/2)-5/2*FB(15+15*i)*fa^(3/2)+FB(15+15*i)^(5/2))/(fa-FB(15+15*i))^2;%��Ԫ32��64��96��128��160��192��224��256
    CW(32*i-1)=8*elementarea(16)/15*((FB(15*i)^(5/2))/((FB(15*i)-FB(15*i+15))*(FB(15*i)-fa))+(FB(15*i+15)^(5/2))/((FB(15*i+15)-FB(15*i))*(FB(15*i+15)-fa))+(fa^(5/2))/((fa-FB(15*i))*(fa-FB(15*i+15))));%��Ԫ31,63,95,127,159,191,223,255
end

for i=0:7
    CW(3+32*i)=8*elementarea(2)/15*((FB(1+15*i)^(5/2))/((FB(1+15*i)-FB(2+15*i))*(FB(1+15*i)-FB(16+15*i)))+(FB(2+15*i)^(5/2))/((FB(2+15*i)-FB(1+15*i))*(FB(2+15*i)-FB(16+15*i)))+(FB(16+15*i)^(5/2))/((FB(16+15*i)-FB(1+15*i))*(FB(16+15*i)-FB(2+15*i))));
    CW(4+32*i)=8*elementarea(2)/15*((FB(17+15*i)^(5/2))/((FB(17+15*i)-FB(2+15*i))*(FB(17+15*i)-FB(16+15*i)))+(FB(2+15*i)^(5/2))/((FB(2+15*i)-FB(17+15*i))*(FB(2+15*i)-FB(16+15*i)))+(FB(16+15*i)^(5/2))/((FB(16+15*i)-FB(17+15*i))*(FB(16+15*i)-FB(2+15*i))));
    
    CW(5+32*i)=8*elementarea(3)/15*((FB(17+15*i)^(5/2))/((FB(17+15*i)-FB(2+15*i))*(FB(17+15*i)-FB(3+15*i)))+(FB(2+15*i)^(5/2))/((FB(2+15*i)-FB(17+15*i))*(FB(2+15*i)-FB(3+15*i)))+(FB(3+15*i)^(5/2))/((FB(3+15*i)-FB(17+15*i))*(FB(3+15*i)-FB(2+15*i))));
    CW(6+32*i)=8*elementarea(3)/15*((FB(17+15*i)^(5/2))/((FB(17+15*i)-FB(18+15*i))*(FB(17+15*i)-FB(3+15*i)))+(FB(18+15*i)^(5/2))/((FB(18+15*i)-FB(17+15*i))*(FB(18+15*i)-FB(3+15*i)))+(FB(3+15*i)^(5/2))/((FB(3+15*i)-FB(17+15*i))*(FB(3+15*i)-FB(18+15*i))));
    
    CW(7+32*i)=8*elementarea(4)/15*((FB(4+15*i)^(5/2))/((FB(4+15*i)-FB(18+15*i))*(FB(4+15*i)-FB(3+15*i)))+(FB(18+15*i)^(5/2))/((FB(18+15*i)-FB(4+15*i))*(FB(18+15*i)-FB(3+15*i)))+(FB(3+15*i)^(5/2))/((FB(3+15*i)-FB(4+15*i))*(FB(3+15*i)-FB(18+15*i))));
    CW(8+32*i)=8*elementarea(4)/15*((FB(4+15*i)^(5/2))/((FB(4+15*i)-FB(18+15*i))*(FB(4+15*i)-FB(19+15*i)))+(FB(18+15*i)^(5/2))/((FB(18+15*i)-FB(4+15*i))*(FB(18+15*i)-FB(19+15*i)))+(FB(19+15*i)^(5/2))/((FB(19+15*i)-FB(4+15*i))*(FB(19+15*i)-FB(18+15*i))));
    
    CW(9+32*i)=8*elementarea(5)/15*((FB(4+15*i)^(5/2))/((FB(4+15*i)-FB(5+15*i))*(FB(4+15*i)-FB(19+15*i)))+(FB(5+15*i)^(5/2))/((FB(5+15*i)-FB(4+15*i))*(FB(5+15*i)-FB(19+15*i)))+(FB(19+15*i)^(5/2))/((FB(19+15*i)-FB(4+15*i))*(FB(19+15*i)-FB(5+15*i))));
    CW(10+32*i)=8*elementarea(5)/15*((FB(20+15*i)^(5/2))/((FB(20+15*i)-FB(5+15*i))*(FB(20+15*i)-FB(19+15*i)))+(FB(5+15*i)^(5/2))/((FB(5+15*i)-FB(20+15*i))*(FB(5+15*i)-FB(19+15*i)))+(FB(19+15*i)^(5/2))/((FB(19+15*i)-FB(20+15*i))*(FB(19+15*i)-FB(5+15*i))));
    
    CW(11+32*i)=8*elementarea(6)/15*((FB(20+15*i)^(5/2))/((FB(20+15*i)-FB(5+15*i))*(FB(20+15*i)-FB(6+15*i)))+(FB(5+15*i)^(5/2))/((FB(5+15*i)-FB(20+15*i))*(FB(5+15*i)-FB(6+15*i)))+(FB(6+15*i)^(5/2))/((FB(6+15*i)-FB(20+15*i))*(FB(6+15*i)-FB(5+15*i))));
    CW(12+32*i)=8*elementarea(6)/15*((FB(20+15*i)^(5/2))/((FB(20+15*i)-FB(21+15*i))*(FB(20+15*i)-FB(6+15*i)))+(FB(21+15*i)^(5/2))/((FB(21+15*i)-FB(20+15*i))*(FB(21+15*i)-FB(6+15*i)))+(FB(6+15*i)^(5/2))/((FB(6+15*i)-FB(20+15*i))*(FB(6+15*i)-FB(21+15*i))));
    
    CW(13+32*i)=8*elementarea(7)/15*((FB(7+15*i)^(5/2))/((FB(7+15*i)-FB(21+15*i))*(FB(7+15*i)-FB(6+15*i)))+(FB(21+15*i)^(5/2))/((FB(21+15*i)-FB(7+15*i))*(FB(21+15*i)-FB(6+15*i)))+(FB(6+15*i)^(5/2))/((FB(6+15*i)-FB(7+15*i))*(FB(6+15*i)-FB(21+15*i))));
    CW(14+32*i)=8*elementarea(7)/15*((FB(7+15*i)^(5/2))/((FB(7+15*i)-FB(21+15*i))*(FB(7+15*i)-FB(22+15*i)))+(FB(21+15*i)^(5/2))/((FB(21+15*i)-FB(7+15*i))*(FB(21+15*i)-FB(22+15*i)))+(FB(22+15*i)^(5/2))/((FB(22+15*i)-FB(7+15*i))*(FB(22+15*i)-FB(21+15*i))));
    
    CW(15+32*i)=8*elementarea(8)/15*((FB(7+15*i)^(5/2))/((FB(7+15*i)-FB(8+15*i))*(FB(7+15*i)-FB(22+15*i)))+(FB(8+15*i)^(5/2))/((FB(8+15*i)-FB(7+15*i))*(FB(8+15*i)-FB(22+15*i)))+(FB(22+15*i)^(5/2))/((FB(22+15*i)-FB(7+15*i))*(FB(22+15*i)-FB(8+15*i))));
    CW(16+32*i)=8*elementarea(8)/15*((FB(23+15*i)^(5/2))/((FB(23+15*i)-FB(8+15*i))*(FB(23+15*i)-FB(22+15*i)))+(FB(8+15*i)^(5/2))/((FB(8+15*i)-FB(23+15*i))*(FB(8+15*i)-FB(22+15*i)))+(FB(22+15*i)^(5/2))/((FB(22+15*i)-FB(23+15*i))*(FB(22+15*i)-FB(8+15*i))));
    
    CW(17+32*i)=8*elementarea(9)/15*((FB(23+15*i)^(5/2))/((FB(23+15*i)-FB(8+15*i))*(FB(23+15*i)-FB(9+15*i)))+(FB(8+15*i)^(5/2))/((FB(8+15*i)-FB(23+15*i))*(FB(8+15*i)-FB(9+15*i)))+(FB(9+15*i)^(5/2))/((FB(9+15*i)-FB(23+15*i))*(FB(9+15*i)-FB(8+15*i))));
    CW(18+32*i)=8*elementarea(9)/15*((FB(23+15*i)^(5/2))/((FB(23+15*i)-FB(24+15*i))*(FB(23+15*i)-FB(9+15*i)))+(FB(24+15*i)^(5/2))/((FB(24+15*i)-FB(23+15*i))*(FB(24+15*i)-FB(9+15*i)))+(FB(9+15*i)^(5/2))/((FB(9+15*i)-FB(23+15*i))*(FB(9+15*i)-FB(24+15*i))));
    
    CW(19+32*i)=8*elementarea(10)/15*((FB(10+15*i)^(5/2))/((FB(10+15*i)-FB(24+15*i))*(FB(10+15*i)-FB(9+15*i)))+(FB(24+15*i)^(5/2))/((FB(24+15*i)-FB(10+15*i))*(FB(24+15*i)-FB(9+15*i)))+(FB(9+15*i)^(5/2))/((FB(9+15*i)-FB(10+15*i))*(FB(9+15*i)-FB(24+15*i))));
    CW(20+32*i)=8*elementarea(10)/15*((FB(10+15*i)^(5/2))/((FB(10+15*i)-FB(24+15*i))*(FB(10+15*i)-FB(25+15*i)))+(FB(24+15*i)^(5/2))/((FB(24+15*i)-FB(10+15*i))*(FB(24+15*i)-FB(25+15*i)))+(FB(25+15*i)^(5/2))/((FB(25+15*i)-FB(10+15*i))*(FB(25+15*i)-FB(24+15*i))));
    
    CW(21+32*i)=8*elementarea(11)/15*((FB(10+15*i)^(5/2))/((FB(10+15*i)-FB(11+15*i))*(FB(10+15*i)-FB(25+15*i)))+(FB(11+15*i)^(5/2))/((FB(11+15*i)-FB(10+15*i))*(FB(11+15*i)-FB(25+15*i)))+(FB(25+15*i)^(5/2))/((FB(25+15*i)-FB(10+15*i))*(FB(25+15*i)-FB(11+15*i))));
    CW(22+32*i)=8*elementarea(11)/15*((FB(26+15*i)^(5/2))/((FB(26+15*i)-FB(11+15*i))*(FB(26+15*i)-FB(25+15*i)))+(FB(11+15*i)^(5/2))/((FB(11+15*i)-FB(26+15*i))*(FB(11+15*i)-FB(25+15*i)))+(FB(25+15*i)^(5/2))/((FB(25+15*i)-FB(26+15*i))*(FB(25+15*i)-FB(11+15*i))));
    
    CW(23+32*i)=8*elementarea(12)/15*((FB(26+15*i)^(5/2))/((FB(26+15*i)-FB(11+15*i))*(FB(26+15*i)-FB(12+15*i)))+(FB(11+15*i)^(5/2))/((FB(11+15*i)-FB(26+15*i))*(FB(11+15*i)-FB(12+15*i)))+(FB(12+15*i)^(5/2))/((FB(12+15*i)-FB(26+15*i))*(FB(12+15*i)-FB(11+15*i))));
    CW(24+32*i)=8*elementarea(12)/15*((FB(26+15*i)^(5/2))/((FB(26+15*i)-FB(27+15*i))*(FB(26+15*i)-FB(12+15*i)))+(FB(27+15*i)^(5/2))/((FB(27+15*i)-FB(26+15*i))*(FB(27+15*i)-FB(12+15*i)))+(FB(12+15*i)^(5/2))/((FB(12+15*i)-FB(26+15*i))*(FB(12+15*i)-FB(27+15*i))));
    
    CW(25+32*i)=8*elementarea(13)/15*((FB(13+15*i)^(5/2))/((FB(13+15*i)-FB(27+15*i))*(FB(13+15*i)-FB(12+15*i)))+(FB(27+15*i)^(5/2))/((FB(27+15*i)-FB(13+15*i))*(FB(27+15*i)-FB(12+15*i)))+(FB(12+15*i)^(5/2))/((FB(12+15*i)-FB(13+15*i))*(FB(12+15*i)-FB(27+15*i))));
    CW(26+32*i)=8*elementarea(13)/15*((FB(13+15*i)^(5/2))/((FB(13+15*i)-FB(27+15*i))*(FB(13+15*i)-FB(28+15*i)))+(FB(27+15*i)^(5/2))/((FB(27+15*i)-FB(13+15*i))*(FB(27+15*i)-FB(28+15*i)))+(FB(28+15*i)^(5/2))/((FB(28+15*i)-FB(13+15*i))*(FB(28+15*i)-FB(27+15*i))));
    
    CW(27+32*i)=8*elementarea(14)/15*((FB(13+15*i)^(5/2))/((FB(13+15*i)-FB(14+15*i))*(FB(13+15*i)-FB(28+15*i)))+(FB(14+15*i)^(5/2))/((FB(14+15*i)-FB(13+15*i))*(FB(14+15*i)-FB(28+15*i)))+(FB(28+15*i)^(5/2))/((FB(28+15*i)-FB(13+15*i))*(FB(28+15*i)-FB(14+15*i))));
    CW(28+32*i)=8*elementarea(14)/15*((FB(29+15*i)^(5/2))/((FB(29+15*i)-FB(14+15*i))*(FB(29+15*i)-FB(28+15*i)))+(FB(14+15*i)^(5/2))/((FB(14+15*i)-FB(29+15*i))*(FB(14+15*i)-FB(28+15*i)))+(FB(28+15*i)^(5/2))/((FB(28+15*i)-FB(29+15*i))*(FB(28+15*i)-FB(14+15*i))));
    
    CW(29+32*i)=8*elementarea(15)/15*((FB(29+15*i)^(5/2))/((FB(29+15*i)-FB(14+15*i))*(FB(29+15*i)-FB(15+15*i)))+(FB(14+15*i)^(5/2))/((FB(14+15*i)-FB(29+15*i))*(FB(14+15*i)-FB(15+15*i)))+(FB(15+15*i)^(5/2))/((FB(15+15*i)-FB(29+15*i))*(FB(15+15*i)-FB(14+15*i))));
    CW(30+32*i)=8*elementarea(15)/15*((FB(29+15*i)^(5/2))/((FB(29+15*i)-FB(30+15*i))*(FB(29+15*i)-FB(15+15*i)))+(FB(30+15*i)^(5/2))/((FB(30+15*i)-FB(29+15*i))*(FB(30+15*i)-FB(15+15*i)))+(FB(15+15*i)^(5/2))/((FB(15+15*i)-FB(29+15*i))*(FB(15+15*i)-FB(30+15*i))));
    
end
 
CWTOTAL=sum(CW);%���е�Ԫ����ϵ��֮��

%%%%�����������
W=20*CWTOTAL*pi*(rb^2-ra^2)/20*P0*10^(-6);%��λN


%%%%�����������������
com=phi*pi/4*(d/1000)^2*P0*sqrt(2*rhoa/Pa);%%����������������ϵ��

if trueps(8)<=Bk
    psi8=psi1;
else
    psi8=sqrt(ka/(ka-1)*((trueps(8))^(2/ka)-(trueps(8))^((ka+1)/ka)));
end
m8=com*psi8; %�ڵ�8��������
mtotal=n*2*m8;%����ֹ����������ܵ���������g/s


%%%��ͼ

theta=0:pi/(8*n):20*pi/10;

r=zeros(1,15);
for i=1:15
    r(i)=rb*rs(i);%����r1��r15
end

R=[ra r(1) r(2) r(3) r(4) r(5) r(6) r(7) r(8) r(9) r(10) r(11) r(12) r(13) r(14) r(15) rb];

[RHO,THETA]=meshgrid(R,theta);


Z=zeros(2*8*n+1,17);%16ΪԲ���������ָ�����160Ϊ�������ֲ�����
for i=1:2*8*n+1
    Z(i,1)=pas;
    Z(i,17)=pas;
end



j=0;
    for i=2:16
        Z(j+1,i)=trueps(15*(j)+i-1);
    end


for k_i=0:n-1

for j=k_i*16+1:k_i*16+8
    for i=2:16
        Z(j+1,i)=trueps(15*(j-k_i*16)+i-1);
    end
end

for j=k_i*16+8+1:(k_i+1)*16
    for i=2:16
        Z(j+1,i)=trueps(15*((k_i+1)*16-j)+i-1);
    end
end
end

[XX,YY,ZZ]=pol2cart(THETA,RHO,Z);
surf(XX,YY,ZZ);

xlabel('X');
ylabel('Y');
zlabel('P/P0');













% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over r_1.
function r_1_ButtonDownFcn(hObject, eventdata, handles)
{
    
OK_ZT.m
    }
% hObject    handle to r_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function i_r0_Callback(hObject, eventdata, handles)
% hObject    handle to i_r0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_r0 as text
%        str2double(get(hObject,'String')) returns contents of i_r0 as a double


% --- Executes during object creation, after setting all properties.
function i_r0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_r0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2



function i_P_Callback(hObject, eventdata, handles)
% hObject    handle to i_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_P as text
%        str2double(get(hObject,'String')) returns contents of i_P as a double


% --- Executes during object creation, after setting all properties.
function i_P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function i_T_Callback(hObject, eventdata, handles)
% hObject    handle to i_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_T as text
%        str2double(get(hObject,'String')) returns contents of i_T as a double


% --- Executes during object creation, after setting all properties.
function i_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function i_P0_Callback(hObject, eventdata, handles)
% hObject    handle to i_P0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_P0 as text
%        str2double(get(hObject,'String')) returns contents of i_P0 as a double


% --- Executes during object creation, after setting all properties.
function i_P0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i_P0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes during object creation, after setting all properties.
function r_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
A_P=imread('ZG.png');
imshow(A_P);

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

javaFrame = get(hObject, 'JavaFrame');
javaFrame.setFigureIcon(javax.swing.ImageIcon('favicon.jpg'));

% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
