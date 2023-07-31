function [A,Gr] = Parse_LaMEM_bin(Parallel_partition,Paraview_output,LaMEM_Parallel_output,npart)
%==========================================================================
% OUTPUT OPTIONS
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO

Parallel_partition     = 'ProcessorPartitioning_8cpu_4.1.2.bin';
RandomNoise             =   logical(0);
Is64BIT                 =   logical(0);
%==========================================================================
% LOAD MESH GRID FROM LaMEM PARTITIONING FILE
%==========================================================================
npart_x =   npart(1);
npart_y =   npart(2);
npart_z =   npart(3);

% Load grid from parallel partitioning file
[X,Y,Z,x,y,z, Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition, RandomNoise, Is64BIT);
Gr.x_g = [min(x),max(x)];
Gr.z_g =[min(z),max(z)];
Gr.y_g = [min(y),max(y)];
% Update variables (size of grid in [x,y,z] direction
nump_x  =   size(X,2);
nump_y  =   size(X,1);
nump_z  =   size(X,3);
W       =   max(x)-min(x);
mW      =   abs(min(x));
L       =   max(y)-min(y);
mL      =   abs(min(y));
H       =   max(z)-min(z);
mH      =   abs(min(z));
Xvec    =   squeeze(X(1,:,1));
Yvec    =   squeeze(Y(:,1,1));
Zvec    =   squeeze(Z(1,1,:));
% Temporary save memory
clear Z Xpart Ypart Zpart
% Prepare data for visualization/output
A           =   struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);
% Linear vectors containing coords
[A.Xpart,A.Ypart,A.Zpart] =meshgrid(single(Xvec),single(Yvec),single(Zvec));
end