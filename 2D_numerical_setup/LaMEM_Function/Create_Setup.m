function Create_Setup(Terranes,ph,TI,A,npart,Gr,Parallel_partition)

RandomNoise             =   logical(0);
Is64BIT                 =   logical(0);
Phase = 0.0.*((A.Xpart));
Temp  = 0.0.*((A.Xpart));

% Set Mantle Phase
Phase(:,:,:)  = nan;
Temp(:,:,:)   = TI.TP;

% Function in each classes of terranes (i.e., trench, passive margin, and
% so forth) has the same name, but different internal structure. It is
% easier to call automatically in the future.
%=====================================================================%
[Phase,Temp]=Terranes.Continent1.fill_terranes(A,Phase,Temp);
[Phase,Temp]=Terranes.Continent2.fill_terranes(A,Phase,Temp);
[Phase,Temp]=Terranes.T.fill_terranes(A,Phase,Temp); 


terranes_list = fieldnames(Terranes);
for it =1:length(terranes_list)
    t = Terranes.(terranes_list{it});
    disp(['Filling the ',terranes_list{it} ,' terranes .... '])

    a = cputime;
    [Phase,Temp] =  Set_Phase_Temperature(A,Phase,Temp,t,Gen);
    b = cputime;
    time = b-a;
    disp(['took ', num2str(round(time)), ' s'])
    disp('=================================')


end
%===========================================================================%
% Set Air Phase
ind = Phase == 0 & A.Zpart<0.0;
Phase(ind)  = Gen.Ph_UM;
Temp(ind)   = Gen.T_P;
ind = Phase == 0 & A.Zpart>0.0;
Temp(ind)   = Gen.T_S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A.Xpart  =  permute(A.Xpart,[2 1 3]);
A.Ypart  =  permute(A.Ypart,[2 1 3]);
A.Zpart  =  permute(A.Zpart,[2 1 3]);


% We can still manually change all this & include, say, a lower mantle
A.nump_x = npart(1);
A.nump_y = npart(2);
A.nump_z = npart(3);
A.Phase  = double(Phase); clear Phase
A.Temp   = double(Temp);  clear Temp
A.Phase  = permute(A.Phase,[2 1 3]);
A.Temp   = permute(A.Temp, [2 1 3]);

x = squeeze(A.Xpart(:,1,1));
y = squeeze(A.Ypart(1,:,1));
z = squeeze(A.Zpart(1,1,:));
A.x      =  double(x(:));
A.y      =  double(y(:));
A.z      =  double(z(:));

[A,surf] = displace_phase_isostasy(ph,A,Gr,Gen);

plot_initial_setup2D(A,surf);

A.RandomNoise = logical(0);

clear Temp Phase


% PARAVIEW VISUALIZATION
FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option

% SAVE PARALLEL DATA (parallel)
FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition, Is64BIT);

end

