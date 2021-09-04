% MESH3D creates a 3D representation of the 3D grid
close all
clear all
load cape

fgrid = fopen('/Users/campo/Work/papers/transcathy-num-diff/cathy/bea/output/grid3d','r');
TETRA=[];
NODES=[];
NNOD=0;
N=0;
NT=0;

A = fscanf(fgrid,'%u %u %u',3);
NNOD = A(1);
N = A(2);
NT = A(3);
TETRA = fscanf(fgrid,'%u',[5,NT]); % Read data
TETRA = TETRA';
TETRA = TETRA(:,1:4);
NODES = fscanf(fgrid,'%g',[3,N]);
NODES = NODES';
% tetramesh(TETRA,NODES);%,'CData',NODES(:,3));
% % colormap(map)
% axis image;
% set(gca,'FontSize',14)
% xlabel('Easting (m)')
% ylabel('Northing (m)')
% zlabel('Elevation (m)')
clear S

S.Vertices = NODES;
S.Faces = TETRA;
S.FaceVertexCData = NODES(:,3);
S.FaceColor = 'interp';
% S.FaceAlpha = 0.5;
% S.LineStyle = ':';
S.LineWidth = 0.25;
% S.EdgeColor = 'red';
% S.LineWidth = 2;
figure
patch(S)
axis image
view(3)
set(gca,'FontSize',14);
xlabel('Easting (m)');
ylabel('Northing (m)');
zlabel('Elevation (m)');
% colormap(map);
colorbar;
