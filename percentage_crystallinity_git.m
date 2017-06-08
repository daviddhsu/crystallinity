% Calculates percentage crystallinity in boxed finite regions of a
% cross-linked polyethalene shape memory polymer.
% author: Chesterton Schuchardt, date: 05/11/2017

clear all, close all 

%% Import Data from File
% % Read xyz data from DCD file
% chains=20;
% beads_chain=500;
% total_bead=chains*beads_chain;
% 
% addpath 'F:\Box Sync\Research\Shape Memory Polymer 2017\Prior Research Codes\readDCD'
% name = sprintf('PE_20_500_muthukumar_lite.dcd');
% atom = readdcd(name,1:total_bead);
    
% Read from excel file
[num] = xlsread('xyz.xlsx');

%% Initialize data and box conditions
xpos = num(1:10000,1); % gets x values from excel
ypos = num(1:10000,2); % gets y values from excel
zpos = num(1:10000,3); % gets z values from excel

% xpos = atom(1450,1:3:end); % gets x values from dcd
% ypos = atom(1450,2:3:end); % gets y values from dcd
% zpos = atom(1450,3:3:end); % gets z values from dcd

xmin = min(xpos); xmax = max(xpos); 
ymin = min(ypos); ymax = max(ypos);     % finds bounds of SMP
zmin = min(zpos); zmax = max(zpos);

%Calculate total box dimensions
xlength = xmax-xmin;
ylength = ymax-ymin;
zlength = zmax-zmin;

nx = 20; % number of x slices
ny = 5; % number of y slices
nz = 5; % number of z slices

xgrid = linspace(xmin,xmax,nx);
ygrid = linspace(ymin,ymax,ny);     % creates grid (boxes)
zgrid = linspace(zmin,zmax,nz);

nchains = 20;        % ns = number of chains
nbeadspchain = 500;     % nbps = number of beads per chain
nbeadstotal = nchains*nbeadspchain;  % nbt = number of beads in total

%% Plot simulation box

%Plot chains w/ different colors, no wrap
c=jet(nchains); % Set up color matrix
figure;
for i=1:nchains
    for j=1:nbeadspchain-1
        length=norm([xpos((i-1)*nbeadspchain+j+1)-xpos((i-1)*nbeadspchain+j),ypos((i-1)*nbeadspchain+j+1)-ypos((i-1)*nbeadspchain+j),zpos((i-1)*nbeadspchain+j+1)-zpos((i-1)*nbeadspchain+j)]);
        if abs(length) <=5
        line(xpos((i-1)*nbeadspchain+j:(i-1)*nbeadspchain+j+1),ypos((i-1)*nbeadspchain+j:(i-1)*nbeadspchain+j+1),zpos((i-1)*nbeadspchain+j:(i-1)*nbeadspchain+j+1),'color',c(i,:))
        hold on
        else
        end
    end
end
axis equal




%% Get chain direction vectors use midpoint method

%Create new position array that follows the centers of every other bond
%Loop over chains
midbondxyz = zeros((nbeadspchain-1)*nchains,3);
for i=1:nchains
    for j=1:nbeadspchain-1
        point1=[xpos((i-1)*nbeadspchain+j),ypos((i-1)*nbeadspchain+j),zpos((i-1)*nbeadspchain+j)];
%         point2=[xpos((i-1)*nbeadspchain+j+1),ypos((i-1)*nbeadspchain+j+1),zpos((i-1)*nbeadspchain+j+1)];
        %use minimum image convention to prevent locating avg bonds in the
        %middle of the box
        umid=[(xpos(((i-1)*nbeadspchain)+j+1)+xlength)-xpos(((i-1)*nbeadspchain)+j), xpos(((i-1)*nbeadspchain)+j+1)-xpos(((i-1)*nbeadspchain)+j), (xpos(((i-1)*nbeadspchain)+j+1)-xlength)-xpos(((i-1)*nbeadspchain)+j)];
        vmid=[(ypos(((i-1)*nbeadspchain)+j+1)+ylength)-ypos(((i-1)*nbeadspchain)+j), ypos(((i-1)*nbeadspchain)+j+1)-ypos(((i-1)*nbeadspchain)+j), (ypos(((i-1)*nbeadspchain)+j+1)-ylength)-ypos(((i-1)*nbeadspchain)+j)];
        wmid=[(zpos(((i-1)*nbeadspchain)+j+1)+zlength)-zpos(((i-1)*nbeadspchain)+j), zpos(((i-1)*nbeadspchain)+j+1)-zpos(((i-1)*nbeadspchain)+j), (zpos(((i-1)*nbeadspchain)+j+1)-zlength)-zpos(((i-1)*nbeadspchain)+j)];
        point2vector=[umid(find(abs(umid) == min(abs(umid)))),vmid(find(abs(vmid) == min(abs(vmid)))),wmid(find(abs(wmid) == min(abs(wmid))))];
        point2=point1+point2vector;
        midbondxyz((i-1)*(nbeadspchain-1)+j,:)=mean([point1;point2]);
    end
end

plot3(midbondxyz(:,1),midbondxyz(:,2),midbondxyz(:,3),'o')

%Create vectors that follow centers of every bond
midvectors = zeros((nbeadspchain-2)*nchains,6);
for i=1:nchains
    for j=1:nbeadspchain-2
        midvectors(((i-1)*(nbeadspchain-2))+j,1:3) = midbondxyz(((i-1)*(nbeadspchain-1))+j,:);
        
        u=[(midbondxyz(((i-1)*(nbeadspchain-1))+j+1,1)+xlength)-midbondxyz(((i-1)*(nbeadspchain-1))+j,1), midbondxyz(((i-1)*(nbeadspchain-1))+j+1,1)-midbondxyz(((i-1)*(nbeadspchain-1))+j,1), (midbondxyz(((i-1)*(nbeadspchain-1))+j+1,1)-xlength)-midbondxyz(((i-1)*(nbeadspchain-1))+j,1)];
        midvectors((i-1)*((nbeadspchain-2))+j,4)=u(find(abs(u) == min(abs(u))));
        
        v=[(midbondxyz(((i-1)*(nbeadspchain-1))+j+1,2)+ylength)-midbondxyz(((i-1)*(nbeadspchain-1))+j,2), midbondxyz(((i-1)*(nbeadspchain-1))+j+1,2)-midbondxyz(((i-1)*(nbeadspchain-1))+j,2), (midbondxyz(((i-1)*(nbeadspchain-1))+j+1,2)-ylength)-midbondxyz(((i-1)*(nbeadspchain-1))+j,2)];
        midvectors((i-1)*((nbeadspchain-2))+j,5)=v(find(abs(v) == min(abs(v))));
        
        w=[(midbondxyz(((i-1)*(nbeadspchain-1))+j+1,3)+zlength)-midbondxyz(((i-1)*(nbeadspchain-1))+j,3), midbondxyz(((i-1)*(nbeadspchain-1))+j+1,3)-midbondxyz(((i-1)*(nbeadspchain-1))+j,3), (midbondxyz(((i-1)*(nbeadspchain-1))+j+1,3)-zlength)-midbondxyz(((i-1)*(nbeadspchain-1))+j,3)];
        midvectors((i-1)*((nbeadspchain-2))+j,6)=w(find(abs(w) == min(abs(w))));
    end
end

for i=1:size(midvectors,1)
    lengths(i,:)=norm(midvectors(i,[4 5 6]));
end

% figure;
c=jet(nchains);
for i=1:nchains
    for j=1:nbeadspchain-2
        quiver3(midvectors((i-1)*(nbeadspchain-2)+j,1),midvectors((i-1)*(nbeadspchain-2)+j,2),midvectors((i-1)*(nbeadspchain-2)+j,3),midvectors((i-1)*(nbeadspchain-2)+j,4),midvectors((i-1)*(nbeadspchain-2)+j,5),midvectors((i-1)*(nbeadspchain-2)+j,6),'color',c(i,:))
        hold on
    end
end
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%% Calculate P2 order parameter using neighbor sphere algorithm

% add stuff in here




%% Calculate volume of crystals 
% Sort vector IDs into bins along X-direction
xbins=10;





% try initializing structure instead of building it 

% Sorts vector IDs into bins in 3D structure
for x = 1:nx-1
    for y = 1:ny-1
        for z = 1:nz-1
            grids(x,y,z).data = find(Mvectors(:,1) >= xgrid(x) & Mvectors(:,1) < xgrid(x+1) & Mvectors(:,2) >= ygrid(y) & Mvectors(:,2) < ygrid(y+1) & Mvectors(:,3) >= zgrid(z) & Mvectors(:,3) < zgrid(z+1));
        end 
    end 
end 

%B = zeros(10000,10000);

%% Sort vectors into grid and calculate P2 order parameter

for x = 1:nx-1      % cycles through x coord of bins
    for y = 1:ny-1      % cycles through y coord of bins
        for z = 1:nz-1      % cycles through z coord of bins
            bindotprod = [];
            
            for row  = 1:length(grids(x,y,z).data) % cycles through indices in an array at a given x,y,z coordinate
                vectorsinbin = grids(x,y,z).data; % takes the array of indices from a given bin
                u = [Mvectors(vectorsinbin(row),4) Mvectors(vectorsinbin(row),5) Mvectors(vectorsinbin(row),6)]; % creates a vector row(index) of M, column 4,5,6
                unitu = u./norm(u); % takes unit vector
                
                for r = row:length(grids(x,y,z).data)-1  % cycles through the rest of the indeces
                    v = [Mvectors(vectorsinbin(r+1),4) Mvectors(vectorsinbin(r+1),5) Mvectors(vectorsinbin(r+1),6)];   % generates vectors to take dot prod with
                    unitv = v./norm(v);  % takes unit vector
                    bindotprod(row,r) = dot(unitv,unitu); % takes dot product and stores it in matrix
                    
                end 
            end 
             P2(y,x,z) = 3/2.*(mean(bindotprod(bindotprod~=0))).^2-1/2;
        end 
    end 
end 
%x = -2:0.2:2;
%y = -2:0.25:2;
%z = -2:0.16:2;
%[X,Y,Z] = meshgrid(x,y,z);

%Use X, Y, and Z to define V as a matrix of volume data.
%V = X.*exp(-X.^2-Y.^2-Z.^2);

%Return matrices Xi, Yi, and Zi from the sphere function.
%[Xi,Yi,Zi] = sphere;

%Draw contours through the volume V along the surface defined by Xi, Yi, and Zi. Change the plot view to a 3-D view.
%contourslice(X,Y,Z,V,Xi,Yi,Zi)
%view(3)

%% Visualize order parameter data

[X,Y,Z] = meshgrid((xgrid(1:end-1)+0.5.*(xmax-xmin)./nx),(ygrid(1:end-1)+0.5.*(ymax-ymin)./ny),(zgrid(1:end-1)+0.5.*(zmax-zmin)./nz));
% visualization method #1 
xbin = [1:nx-1];
ybin = [1:ny-1];
zbin = [1:nz-1];
Parameter = P2(ybin,xbin,zbin);
Xi = xgrid;
Yi = ygrid;
Zi = [];
figure
contourslice(X,Y,Z,Parameter,Xi,Yi,Zi)
view(3)

% visualization method #2
figure
contourf(X(:,:,1),Y(:,:,1),Parameter(:,:,1),1000,'edgecolor','none')
[x,y,z,v] = flow;

% visualization method #3
figure;
Xi = [-90 -10 80];
Yi = 10;
Zi = -10;
slice(X,Y,Z,Parameter,Xi,Yi,Zi);
view(3);
axis on;
grid on;
light;
lighting phong;
camlight('left');
shading interp;
