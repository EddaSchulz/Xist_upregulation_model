% ODEs for modeling Xist upregulation at the onset of XCI
% Nelly Kanata and Edda Schulz
% OWL Schulz, Max Planck Institute for Molecular Genetics
% Created: 20.09.2023
% Modified: 05.01.2024


function dx = feedback_model(t,x,k) %t: time, x:%of cells, k=array with k_upx2, k_up, k_silx2, k_res parameters

% initialize dx
dx = [0; 0; 0; 0; 0; 0; 0];

%parameters
k_upx2=k(1); % Xist upregulation rate when there are two Xs to upregulate Xist
k_up=k_upx2/2; % Xist upregulation rate
k_silx2=k(2); % x-chromosome silencing rate when there are two Xs to be silenced 
k_sil=k_silx2/2; % x-chromosome silencing rate
k_res=k(3); % biallelic resolution rate
k_diff=k(4); % differentiation rate

%equations
dx(1)=k_diff*x(7)-k_upx2*x(1); %no Xist
dx(2)=k_upx2*x(1)-k_up*x(2)-k_sil*x(2); %monoallelic not silenced
dx(3)=k_res*x(6)+k_sil*x(2); %monoallelic silenced
dx(4)=k_up*x(2)-k_silx2*x(4); %biallelic not silenced
dx(5)=k_silx2*x(4)-k_sil*x(5); %biallelic one silenced
dx(6)=k_sil*x(5)-k_res*x(6); %biallelic both silenced
dx(7)=-k_diff*x(7); %undifferentiated

end
