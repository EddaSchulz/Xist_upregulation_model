% distance function
% Nelly Kanata and Edda Schulz
% OWL Schulz, Max Planck Institute for Molecular Genetics
% Created: 20.09.2023
% Modified: 12.10.2023



function dist= model_fit(parameters, x_data,init_cond,k,pp,tspan) % init_cond: initial conditions

k(pp)=parameters;

tspan = [0 tspan(2:end)-k(5)]; %fit delay of experimental data vs start of phenomenon


[t,x]= ode45(@feedback_model, tspan, init_cond, [], k);

% sum up the subgroups (negative, mono, bi) of cells (object 'x') that are given by the
% model with the selected parameters

x_summed=x;
x_summed(:,1) = x_summed(:,1)+x_summed(:,7);
%sum up monoallelic silenced and not silenced
x_summed(:,2)=x_summed(:,2)+x_summed(:,3);
%sum up biallelic silenced and not silenced
x_summed(:,3)=x_summed(:,4)+x_summed(:,5)+x_summed(:,6);
%remove old columns
x_summed=x_summed(:,1:3);


% distance of experimental data (x_data) and modelled data (x_summed)
dist_unflattened=(x_summed'-x_data);

dist=dist_unflattened(:); %flatten

end