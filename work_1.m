function varargout = work_1(varargin)
% WORK_1 MATLAB code for work_1.fig
%      WORK_1, by itself, creates a new WORK_1 or raises the existing
%      singleton*.
%
%      H = WORK_1 returns the handle to a new WORK_1 or the handle to
%      the existing singleton*.
%
%      WORK_1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORK_1.M with the given input arguments.
%
%      WORK_1('Property','Value',...) creates a new WORK_1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before work_1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to work_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help work_1

% Last Modified by GUIDE v2.5 25-Nov-2018 16:46:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @work_1_OpeningFcn, ...
                   'gui_OutputFcn',  @work_1_OutputFcn, ...
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


% --- Executes just before work_1 is made visible.
function work_1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to work_1 (see VARARGIN)

% Choose default command line output for work_1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes work_1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


backgroundlmage1 = importdata('jiemian1.jpg');
axes(handles.axes3);
image(backgroundlmage1);
axis off
% 设置显示窗口图案
backgroundlmage2 = importdata('1.jpg');
axes(handles.axes1);
image(backgroundlmage2);
axis off

backgroundlmage3 = importdata('2.jpg');
axes(handles.axes2);
image(backgroundlmage3);
axis off


set(handles.slider1,'enable','off');



% --- Outputs from this function are returned to the command line.
function varargout = work_1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 导入图片

global filename pathname I
% 选择图片
[filename, pathname] = uigetfile({'*.png'; '*.jpg'; '*.bmp'; '*.gif'; '*.fts'}, '选择图片');
% 合成路径+文件名
name = [pathname filename];
I = imread(name);
if ndims(I) == 3  % 如果是彩图，转化为灰度图
    I = rgb2gray(I);
end
axes(handles.axes1);
imshow(I);
axis equal
axis off
axes(handles.axes2);
imhist(I);
% axis([0,255,min(imhist(I)),max(imhist(I))]);
% axis 'auto x' % x轴的范围自动调整
% axis equal
axis off % 去掉坐标轴

set(handles.slider1,'enable','off');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,filespec] = uiputfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.png','(*.png)';...
                                 '*.*',  'All Files (*.*)'},...
                                 'Save Picture','Untitled');
h = getframe(handles.axes2);
imwrite(h.cdata,fullfile(PathName,FileName));
hs = msgbox('图片保存成功！','图片保存','help','modal');
ht = findobj(hs,'Type','text');
set(ht,'FontSize',20,'Unit','normal');
% 改变对话框大小
set(hs,'Resize','on');% 自动改变

set(handles.slider1,'enable','off');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf);
close(figure(2));


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 阈值分割类型
global filename pathname I
str = get(hObject,'string');

switch str
    case '原图'
        set(handles.slider1,'enable','off');
        rotate3d off;
        close(figure(2));
        
        cla reset
        axes(handles.axes1);
        imshow(I);
        
        axes(handles.axes2);
        cla reset  %清空axes2
        axes(handles.axes2);
        imhist(I);
        % axis([0,255,min(imhist(I)),max(imhist(I))]);
        % axis 'auto x' % x轴的范围自动调整
        % axis equal
        axis off % 去掉坐标轴
    case '人工选择'
        rotate3d off;
        close(figure(2));
        set(handles.slider1,'enable','on');
        
        while 1
            prompt = {'人工选择阈值为(0<= x <= 255)：'};%设置提示字符串
            title = '确定人工选择阈值';%设置标题
            % numlines = 1;%指定输入数据的行数
            defAns = {'128'};%设定默认值
            Resize = 'on';%设定对话框尺寸可调节
            answer = inputdlg(prompt,title,[1 45],defAns,'on');%创建输入对话框
            arr = answer(1,1);
            arr = str2num(char(arr));  % 手动输入的阈值
            if arr >= 0 && arr <= 255
                break;
            end
            uiwait(errordlg('请输入合理的阈值！','阈值错误','model'));
            % strh=findobj(h,'tag','MessageBox');
            % set(strh,'FontSize',15,'FontName','楷体');
        end
        [row,col] = size(I);
        for i = 1:row
            for j = 1:col
                if (I(i,j) < arr)
                    num1(i,j) = 0;
                else
                    num1(i,j) = 1;
                end
            end
        end
        cla reset
        axes(handles.axes1);
        imshow(I);
        axis off
        
        axes(handles.axes2);
        cla reset  %清空axes2
        axes(handles.axes2);
        imshow(num1);
        axis off
        set(handles.text4,'string',[num2str(arr)]);
        
        set(handles.slider1,'value',arr);
        
        figure(2)
        imhist(I);
        axis([0,255,min(imhist(I)),max(imhist(I))]);
        hold on
        plot([arr,arr],[0,max(imhist(I))],'r','linewidth',2)
        text(arr+10,round(max(imhist(I))*0.9),['阈值为',num2str(arr)],'horiz','left','color','r','fontsize',15)
        xlabel('灰度级');ylabel('次数（次）')

    case '迭代法'
        set(handles.slider1,'enable','off');
        rotate3d off;
        close(figure(2));
        
        f = double(I);
        T = (min(f(:)) + max(f(:))) / 2;
        done = false;
        i = 0;
        while ~done
            r1 = find(f <= T);
            r2 = find(f > T);
            Tnew = (mean(f(r1)) + mean(f(r2))) / 2;
            done = abs(Tnew-T) < 0.1;
            T = Tnew;  % 阈值为T
            i = i + 1;
        end
        f(r1) = 0;
        f(r2) = 1;
        
        cla reset
        axes(handles.axes1);
        imshow(I);
        axis off
        
        axes(handles.axes2);
        cla reset  %清空axes2
        axes(handles.axes2);
        imshow(f);
        axis off
        set(handles.text4,'string',[num2str(T)]);
        
        figure(2)
        imhist(I)
        axis([0,255,min(imhist(I)),max(imhist(I))]);
        hold on
        plot([T,T],[0,max(imhist(I))],'r','linewidth',2)
        text(T+10,round(max(imhist(I))*0.9),['阈值为',num2str(T)],'horiz','left','color','r','fontsize',15)
        xlabel('灰度级');ylabel('次数（次）');

    case 'Otsu算法'
        set(handles.slider1,'enable','off');
        rotate3d off;
        close(figure(2));
        
        T2 = graythresh(I);
        num2 = im2bw(I,T2);
        
        cla reset
        axes(handles.axes1);
        imshow(I);
        axis off
        
        axes(handles.axes2);
        cla reset  %清空axes2
        axes(handles.axes2);
        imshow(num2);
        axis off
        set(handles.text4,'string',[num2str(T2*255)]);
        
        ans_gray_j = T2 * 255;
        figure(2)
        imhist(I)
        axis([0,255,min(imhist(I)),max(imhist(I))]);
        hold on
        plot([ans_gray_j,ans_gray_j],[0,max(imhist(I))],'r','linewidth',2)
        text(ans_gray_j+10,round(max(imhist(I))*0.9),['阈值为',num2str(ans_gray_j)],'horiz','left','color','r','fontsize',15)
        xlabel('灰度级');ylabel('次数（次）');
    case '二维Ostu算法'
        set(handles.slider1,'enable','off');
        close(figure(2))
        
        [row,col]=size(I);
        I0=double(I);
        hh=1;
        I1=zeros(row,col);
        for i=1:row
            for j=1:col
                for k=-hh:hh
                    for w=-hh:hh
                        p=i+k;
                        q=j+w;
                        if (p<=0) || (p>row)
                            p = i;
                        end
                        if (q <= 0) || (q>col)
                            q = j;
                        end
                        I1(i,j) = I0(p,q)+I1(i,j);
                    end
                end
                a2(i,j) = uint8(1/9*I1(i,j));
            end
        end
        fxy=zeros(256,256);
        for i=1:row
            for j=1:col
                c=I0(i,j);
                d=double(a2(i,j));
                fxy(c+1,d+1)=fxy(c+1,d+1)+1;
            end
        end
%         figure
%         mesh(fxy)
        pxy=fxy/row/col;
        w0=zeros(256,256);
        ui=zeros(256,256);
        uj=zeros(256,256);
        w0(1,1)=pxy(1,1);
        for i=2:256
            w0(i,1)=w0(i-1,1)+pxy(i,1);
        end
        for i=2:256
            w0(1,i)=w0(1,i-1)+pxy(1,i);
        end
        for i=2:256
            for j=2:256
                w0(i,j)=w0(i-1,j)+w0(i,j-1)-w0(i-1,j-1)+pxy(i,j);
            end
        end
        w1=ones(256,256)-w0;
        ui(1,1)=0;
        for i=2:256
            ui(1,i)=ui(1,i-1)+(1-1)*pxy(1,i);
        end
        for i=2:256
            ui(i,1)=ui(i-1,1)+(i-1)*pxy(i,1);
        end
        for i=2:256
            for j=2:256
                ui(i,j)=ui(i-1,j)+ui(i,j-1)-ui(i-1,j-1)+(i-1)*pxy(i,j);
            end
        end
        uj(1,1)=0;
        for i=2:256
            uj(1,i)=uj(1,i-1)+(i-1)*pxy(1,i);
        end
        for i=2:256
            uj(i,1)=uj(i-1,1)+(1-1)*pxy(i,1);
        end
        for i=2:256
            for j=2:256
                uj(i,j)=uj(i-1,j)+uj(i,j-1)-uj(i-1,j-1)+(j-1)*pxy(i,j);
            end
        end
        uti=0;
        utj=0;
        for i=1:256
            for j=1:256
                uti=uti+(i-1)*pxy(i,j);
                utj=utj+(j-1)*pxy(i,j);
            end
        end
        for i=1:256
            for j=1:256
                if w0(i,j)~=0 && w1(i,j)~=0
                    hh(i,j)=((uti*w0(i,j)-ui(i,j))^2+(utj*w0(i,j)-uj(i,j))^2)/(w0(i,j)*w1(i,j));
                else
                    hh(i,j)=0;
                end
            end
        end
        hmax=max(hh(:));
        for i=1:256
            for j=1:256
                if hh(i,j)==hmax
                    s = i-1;
                    k = j-1;
                    continue;
                end
            end
        end
        z=ones(row,col);
        for i=1:row
            for j=1:col
                if I(i,j)<=s || a2(i,j)<=k
                    z(i,j)=0;
                end
            end
        end
        cla reset
        
        set(handles.text4,'string',['  领域均值: ',num2str(s),'   灰度值: ',num2str(k)]);
        axes(handles.axes1);
        mesh(fxy);
        grid on
        zlabel('像素值');
        colormap jet
        rotate3d on;  % 可以旋转
        
        axes(handles.axes2);
        cla reset  %清空axes2
        axes(handles.axes2);
        imhist(I)
        axis([0,255,min(imhist(I)),max(imhist(I))]);
        hold on
        plot([k,k],[0,max(imhist(I))],'r','linewidth',2)
        text(k+10,round(max(imhist(I))*0.9),['阈值为',num2str(k)],'horiz','left','color','r','fontsize',15);
        xlabel('灰度级');ylabel('次数（次）');
        
        figure(2);
        subplot(1,2,1),imshow(I);
        subplot(1,2,2),imshow(z);
        
%         axes(handles.axes2);
%         imshow(z);
%         axis off
%         set(handles.text4,'string',['x: ',num2str(s),'  y: ',num2str(k)]);
%         figure(2)
%         mesh(fxy);
%         figure(3)
%         imhist(I)
%         axis([0,255,min(imhist(I)),max(imhist(I))]);
%         hold on
%         plot([k,k],[0,max(imhist(I))],'r','linewidth',2)
%         text(k+10,round(max(imhist(I))*0.9),['阈值为',num2str(k)],'horiz','left','color','r','fontsize',15);
%         xlabel('灰度级');ylabel('次数（次）');

    case '模拟退火'
        set(handles.slider1,'enable','off');
        rotate3d off;
        close(figure(2));
        
        [row,col] = size(I);
        A = 256;
        N = row * col;   %图像的总像素数
        J = imhist(I,A);
        ANS = 0;
        ans_gray_i = 0;
        ans_gray_j = 0;
        ut = 0;   %平均像素
        tic
        for i = 1:1:A
             ut = ut + i * (J(i,:)/N);
        end
        
        for i = 1:1:A-1
            for j = i:1:A
                J1 = J(1 : i , :);
                w0 = sum(J1) ./ N;
                J2 = J(i+1 : j,:);
                w1 = sum(J2) ./ N;
                J3 = J(j+1 : A,:);
                w2 = sum(J3) ./ N;
                u0 = 0;
                u1 = 0;
                u2 = 0;
                for k = 1 : 1 : i
                    M = J(k);
                    P = M ./ N;
                    u0 = u0 + k * P;
                end
                for k = i + 1 : 1 : j
                    M = J(k);
                    P = M ./ N;
                    u1 = u1 + k * P;
                end
                for k = j + 1 : 1 : A
                    M = J(k);
                    P = M ./ N;
                    u2 = u2 + k * P;
                end
                u0 = u0 / w0;
                u1 = u1 / w1;
                u2 = u2 / w2;
                jiao = w0 * (u0 - ut)^2 + w1 * (u1 - ut)^2 + w2 * (u2 - ut)^2;
                if ANS < jiao
                    ANS = jiao;
                    ans_gray_i = i;
                    ans_gray_j = j;
                end
            end
        end
        
        im2_ans = ones(row,col);
        for i = 1:row
            for j = 1:col
                if I(i,j) <= ans_gray_j
                    im2_ans(i,j) = 0;
                end
            end
        end
        toc
        cla reset
        
        axes(handles.axes1);
        imshow(I);
        axis off
        
        axes(handles.axes2);
        cla reset  %清空axes2
        
        axes(handles.axes2);
        imshow(im2_ans);
        axis off
        set(handles.text4,'string',[num2str(ans_gray_j)]);
        
        figure(2)
        imhist(I)
        axis([0,255,min(imhist(I)),max(imhist(I))]);
        hold on
        plot([ans_gray_j,ans_gray_j],[0,max(imhist(I))],'r','linewidth',2)
        text(ans_gray_j+10,round(max(imhist(I))*0.9),['阈值为',num2str(ans_gray_j)],'horiz','left','color','r','fontsize',15)
        xlabel('灰度级');ylabel('次数（次）');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global I

close(figure(2));

k = round(get(hObject,'value'));
set(handles.text4,'string',num2str(k));

[row,col] = size(I);
        for i = 1:row
            for j = 1:col
                if (I(i,j) < k)
                    num1(i,j) = 0;
                else
                    num1(i,j) = 1;
                end
            end
        end
        
        axes(handles.axes1);
        cla reset
        axes(handles.axes1);
        imshow(I);
        axis off
        
        axes(handles.axes2);
        cla reset  %清空axes2
        axes(handles.axes2);
        imshow(num1);
        axis off

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
