clear all;
% Stuff to edit for different cores 
% 1. finalcorename
% 2. filenames
% 3. dst array 
% 4. elements you want to plot (at the bottom) 
% 5. can also play around with "n" value (low-pass filter) 

%combines XRF files

finalcorename='JPC05';

% % Slicker to load in files, but not necessary 
% corename{1}='VBM1_0to100';dst(1)=0;
% corename{2}='VBM1_100to185';dst(2)=nan;

% Easier way to load in files 
filename{1} = 'JPC05_1of12_xrf.mat'; dst(1)=0;
filename{2} = 'JPC05_2of12_xrf.mat'; dst(2)=nan;
filename{3} = 'JPC05_3of12_xrf.mat'; dst(3)=nan;
filename{4} = 'JPC05_4of12_xrf.mat'; dst(4)=nan;
filename{5} = 'JPC05_5of12_xrf.mat'; dst(5)=nan;
filename{6} = 'JPC05_6of12_xrf.mat'; dst(6)=nan;
filename{7} = 'JPC05_7of12_xrf.mat'; dst(7)=nan;
filename{8} = 'JPC05_8of12_xrf.mat'; dst(8)=nan;
filename{9} = 'JPC05_9of12_xrf.mat'; dst(9)=nan;
filename{10} = 'JPC05_10of12_xrf.mat'; dst(10)=nan;
filename{11} = 'JPC05_11of12_xrf.mat'; dst(11)=nan;
filename{12} = 'JPC05_12of12_xrf.mat'; dst(12)=nan;

% Array of start depths
% Write "nan" if there is no gap between drives
% Fill in a value if there's a gap, and you want the drive to start at a specific depth 
% dst = [0 nan]; 

for i=1:length(filename);
    
    C=load(filename{i}); 
%     C=load(['./' corename{i} '/' corename{i} '_xrf.mat']); %use this line if you go with the slick file loading method 
    
    if i==1; %for first drive 
        vname=C.vname;
%         Xrad=C.Xrad;
        Xopt=C.Xopt;
        data=C.data;
        d=C.d;
        xres=C.xres;
    else %for subsequent drives that don't start at 0
        
        if isnan(dst(i));
            d=[d C.d+max(d)+C.xres/10000];
        else
            d=[d C.d+dst(i)+C.xres/10000];
        end
        

        
        
%         Xrad=[Xrad; C.Xrad];
        Xopt=[Xopt; C.Xopt];
        data=[data C.data];


        if xres~=C.xres;
            display(['WARNING: resolution of scans are different, x-ray images not plotted at accurate depths'])
        end
    end
end

%         %remove suspect data
%         ind=96<d & d<103; %d is depth
%         data(:,[ind])=nan; %zn is in row 24 in data array
% 
%         ind=195<d &d<199;
%         data(:,[ind])=nan;


%Defining Elements to plot, can change these to look at other elements
res=10;


% var{1}='Sr';n(1)=res;%var is the name of the element in vname, n is the low pass filter (n=1 is for no filter)
% var{1}='Ca';n(1)=res;
% var{2}='Ti';n(1)=res;
% var{3}='Sr';n(1)=res;
% var{4}='Fe';n(1)=res;
% var{5}='Zn';n(1)=res;
% var{6}='Pb';n(1)=res;
% var{7}='Br';n(1)=res;
% var{8}='Ba';n(1)=res;
% var{6}='Zr';n(1)=res;
var{1}='TS';n(1)=res;
var{2}='TC';n(1)=res;
var{3}='TZ';n(1)=res;


%Get ratios
% [Ti]=getvar(var{2},data,vname,n); 
% [Ca]=getvar(var{3},data,vname,n); 
% [Sr]=getvar(var{1},data,vname,n); 
% 
% TiSr=Ti./Sr;
% data=vertcat(data,TiSr);
% TiCa=Ti./Ca;
% data=vertcat(data,TiSr);
% 
% var{6}='Ti/Sr';
% vname{48}='Ti/Sr';
% var{7}='Ti/Ca';
% vname{49}='Ti/Ca';

figure
plt_xrf_SKQ(finalcorename,Xopt,vname,data,d,var,n)
ax=gca
grid on
ax.FontSize = 16;
yticks(0:100:1500)
grid on
hold off


eval(['save ' finalcorename ' Xopt vname var data d finalcorename n']);

Ca_bottom = nanmean(data(16,[11000:15097]))
Ca_top = nanmean(data(16,[1:11000]))

TiSr_bottom = nanmean(data(50,[11000:15097]))
TiSr_top = nanmean(data(50,[1:11000]))

TiCa_bottom = nanmean(data(51,[11000:15097]))
TiCa_top = nanmean(data(51,[1:11000]))

TiZr_bottom = nanmean(data(52,[11000:15097]))
TiZr_top = nanmean(data(52,[1:11000]))
    