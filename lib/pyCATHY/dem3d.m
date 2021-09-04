% DEM3D creates a 3D representation from a Grass DEM file
% Default = 0 or -9999 ?
clear all
close all
delete dem3d.m~;
load cape
clear X
clear caption

% DEM 1
dem = fopen('/Users/campo/Work/papers/panola/deep/prepro_bedstraight/dem','r');
DEM=[];
DEm=[];
a=0;
b=0;
Dem=0;
highest=-100;
lowest=10000;
length=1;
width=1;
north=0;
south=0;
east=0;
west=0;
cols=0;
rows=0;

% Read the Header
%for a=1:2
%    fgets(dem);
%end
north = fscanf(dem,'%*s %g',1);
south = fscanf(dem,'%*s %g',1);
east = fscanf(dem,'%*s %g',1);
west = fscanf(dem,'%*s %g',1);
rows = fscanf(dem,'%*s %g',1);
cols = fscanf(dem,'%*s %g',1);

% Find out Highest and Lowest
for a=1:rows
    for b=1:cols
        Dem = fscanf(dem,'%g',1);
        if Dem ~=0
            if Dem > highest
                highest = Dem;
            end
            if Dem < lowest
                lowest = Dem;
            end
        end
    end
end
frewind(dem);

% Read the DEM
for a=1:6
    fgets(dem);
end
for a=1:rows
    for b=1:cols
        Dem = fscanf(dem,'%g',1);
        if Dem==0
            %Dem = -9999;
            Dem = NaN;
        end
        DEm = [DEm;Dem];
    end
    DEm = DEm';
    DEM = [DEM;DEm];
    DEm = [];
end
% DEM=DEM';
% Multiplies by scale factor
DEM=DEM.*1.0E+00;
DEM=flipud(DEM);
DEM1=DEM;
fclose(dem);

% % Determine grid size
% length = (north-south)/rows;
% width = (east-west)/cols;

for a=1:cols
    x(a)=west+length*a;
end
for a=1:rows
    y(a)=south+width*a;
end
x=x-width/2;
y=y-length/2;
% y=fliplr(y);
% Create the plot
%imagesc(x,y,DEM)

% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1])
% set(gcf,'Renderer','Zbuffer');

% subplot(1,3,1);
figure(1)
set(gcf,'Renderer','Zbuffer');
surf(x,y,DEM1);
% caxis([12.3 59.5]);
view(2)
shading flat
colorbar
% line([3 240],[123 123],[100 100],'LineStyle','--','Color','k','LineWidth',4);
% line([6 6],[3 240],[100 100],'LineStyle','--','Color','k','LineWidth',4);
axis image
colormap(map)
% axis image
% set(gca,'YDir','normal')
% set(gca,'FontSize',18);
xlabel('Easting (m)')
ylabel('Northing (m)')
% print('-djpeg','-r300','/Users/campo/Work/papers/journals/cathy-LCC/review/VCatchment/dem1');
% 
% % print('-djpeg','-zbuffer','-r300','/Users/campo/Work/papers/journals/cathy-LCC/review/VCatchment/dem');