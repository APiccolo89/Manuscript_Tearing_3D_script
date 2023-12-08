function Create_Setup(Terranes,ph,TI,A,npart,Gr,Parallel_partition,Save_Test_Benchmark)
%=========================================================================
% Benchmark {Since waiting for 30 minutes everytime I fucked up something,
% I will save A,ph,Gr, in benchmark3D.m}
%==========================================================================
if nargin ~= 0 
RandomNoise             =   logical(0);
Is64BIT                 =   logical(0);
%A.Phase = 0.0.*((A.Xpart));
%A.Temp  = 0.0.*((A.Xpart));
% Set Mantle Phase
A.Phase(:,:,:)  = nan;
A.Temp(:,:,:)   = TI.TP;
% 
terranes_list = fieldnames(Terranes);
TA=cputime;
for it =1:length(terranes_list)
    t = Terranes.(terranes_list{it});
    disp('=================================')
    [A] =  t.fill_terranes(A);
    disp('=================================')
end
TB = cputime;
disp('=================================TERRANES ARE FINISHED==============')
disp(['The full operation took ', num2str(round((TB-TA)./60)), ' minutes'])
disp('====================================================================')
save_h5file_data(Terranes,A.Phase,A.Temp,A.Zpart)

% Final correction 
ind =   abs(A.Temp-TI.TP)<0.1 & A.Zpart<0.0; 
A.Phase(ind)  = TI.Ph_Ast;
A.Temp(ind)   = TI.TP;

ind = A.Zpart>=0.0;
A.Temp(ind)   = TI.TS;
A.Phase(ind)  = TI.Ph_Air;


A.Xpart  =  permute(A.Xpart,[2 1 3]);
A.Ypart  =  permute(A.Ypart,[2 1 3]);
A.Zpart  =  permute(A.Zpart,[2 1 3]);

A.nump_x = npart(1);
A.nump_y = npart(2);
A.nump_z = npart(3);
A.Phase  = double(A.Phase); 
A.Temp   = double(A.Temp);  %clear Temp
A.Phase  = permute(A.Phase,[2 1 3]);
A.Temp   = permute(A.Temp, [2 1 3]);

x = squeeze(A.Xpart(:,1,1));
y = squeeze(A.Ypart(1,:,1));
z = squeeze(A.Zpart(1,1,:));
A.x      =  double(x(:));
A.y      =  double(y(:));
A.z      =  double(z(:));
if Save_Test_Benchmark == 1
    save('Benchmark3D.mat',"A","Gr","ph","TI")
end
else
    clear all;
    close all;
    A_time = cputime;
    load('Benchmark3D.mat');
    B_time = cputime; 
    disp(['Benchmark3D took',num2str(B_time-A_time,3) ,' s, to load']); 
end
[A] = displace_phase_isostasy(ph,A,Gr,TI);
A.Phase(A.Phase==0) = 12;


A.RandomNoise = logical(0);



% PARAVIEW VISUALIZATION
FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option

% SAVE PARALLEL DATA (parallel)
FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition, Is64BIT);

end

