function []=plt_xrf_SKQ(CoreName,Xopt,vname,data,d,var,n)

%defining number of subplots based on number of variables to plot
np=length(var)+3;

%getting max and min depths to define y extent for optical and rad images
dst=min(d);
dend=max(d);

%% adjusting title (taking out '_' and replacing with ' '
a=find(CoreName=='_');
if sum(a)>0;
CoreName(a)=' ';
end

%% optical image
subplot(1,np,1);
image([0 1],[dst dend],Xopt);%ploting image and defining the dstart and end, arbitrary position for x extent
% set(gca,'DataAspectRatio',[1 10 1]);% sets aspect ratio so image is not stretched too much

%below is a trick to replace the '_' with a space since matlab used _ as an
%indication to supscript the text after it
CoreTitle=CoreName;
a=findstr(CoreTitle,'_');
CoreTitle(a)=' ';

ttl{1}=CoreTitle;%title for subplot is core name
ttl{2}='photo';%photo for second line in title
title(ttl);%putting title on subplot
ylabel('Depth (cm)');%ylabel
set(gca,'xtick',[]);%makes it so there are not ticks and numbers on x axis
axis([0 1 dst dend]);%setting x axis limits to 0-1 and y axis limits dst and dend
set(gca,'fontweight','bold');%making font weight for axis bold
set(gca, 'FontSize', 16);
yticks(0:100:1500)
grid on

%% plotting xradiograph
% subplot(1,np,2);
% image([0 1],[dst dend],Xrad);% same deal as optical
% set(gca,'DataAspectRatio',[1 3 1]);% adjusting axis for xradiograph so it is skinnier
% set(gca,'xtick',[],'fontweight','bold');%taking xticks off and making font weight for subplot bold
% title('xray');%second line for title to plot

%% plotting xrf


col='bymgcrkGwbymgcrkGwbymgcrkGw';%defining color for each variable
col='rbgmrbgmrbgmrbgmrbgmrbgmrbgm';

%loop to plot all 4 variables
for ii=1:length(var);
    subplot(1,np,ii+2);
    
    %getvar command belwo is another mfile that gets the data for a specific element from "data" based on order defined in vname
    %n is the lowpass filter length, V is the output data for the particular element
    [V]=getvar(var{ii},data,vname,n);
    
    ind=find(isfinite(V)==1);%only plotting validity=1
    hd=plot(V(ind),d(ind)); %plotting V data with color for i =  1 to 4 (i.e. red=1, blue=2, green=3, magenta=4)
    set(hd,'linewidth',3);%thicking line to 1
    set(hd,'color',cjw_v2(col(ii))); 
    set(gca,'ydir','reverse'); %flipping yaxis so 0 is on the top and depth increases down
    xlabel(var{ii});%labeling axis with element name
    set(gca,'ylim',[dst dend]);%setting ylimits to dst and dend
    grid on
    set(gca, 'FontSize', 16);
    yticks(dst:100:dend) %ticks every 100cm


    %COMENT THIS OUT IF NOT DOING SIKULIAQ CORES
    % if ii==5;
    %    set(gca,'xlim',[0 150])
    % end
           
    set(gca,'xaxislocation','top','fontweight','bold'); %making text for subplot bold and moving axis location to top
end




