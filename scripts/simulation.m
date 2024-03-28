% Script for modeling Xist upregulation at the onset of XCI
% Nelly Kanata and Edda Schulz
% OWL Schulz, Max Planck Institute for Molecular Genetics
% Created: 20.09.2023
% Modified: 27.03.2024

% read data

x_data_summed=readtable('../data/shiura_abe_2019.txt', 'Delimiter', '\t', 'ReadVariableNames', true);
columnNames= x_data_summed.Properties.VariableNames;
t_data =str2double(strrep(columnNames, 'x', ''));
transposed_x_data_summed=x_data_summed{:,:}';


tspan=[0 48]+5;
x0=[0 0 0 0 0 0 100]; % order of cell groups: xist negative differentiated, 
% monoallelic not silenced, monoallelic silenced, biallelic not silenced, 
% biallelic one allele silenced, biallelic both alleles silenced, xist negative undifferentiated

% solve with preset parameters
k=[0.2 0.3 1 0.1];

% get percentages of cells ("x") that belong to each subgroup across time ("t")
% for the selected parameters ("k") in the timespan "tspan", based on the ODEs of
% "feedback_model.m"

[t,x] = ode45(@feedback_model, tspan, x0,[], k);

% sum up predicted data (all negative, monoallelic and biallelic
% subgroups together)
x_summed=x;
x_summed(:,1) = x_summed(:,1)+x_summed(:,7);
%sum up monoallelic silenced and not silenced
x_summed(:,2)=x_summed(:,2)+x_summed(:,3);
%sum up biallelic silenced and not silenced
x_summed(:,3)=x_summed(:,4)+x_summed(:,5)+x_summed(:,6);
%remove old columns
x_summed=x_summed(:,1:3);

mycolors=[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];

plot(t_data,transposed_x_data_summed,'.'); % plot data points

hold on; % to keep the same plot and overlay

plot(t,x_summed);
legend('no Xist', 'monoallelic','biallelic', 'sim no xist', 'sim monoallelic', 'sim biallelic'); 
ylabel('% of cells');
xlabel('Time (hours)');
ax = gca;
ax.ColorOrder = mycolors;
hold off;

saveas(gcf, '../output/model_random_par.pdf')


%% Fit data

options = optimset('display','iter','MaxIter',1000,'TolX',1e-7, 'TolFun',2e-8);

options.MaxFunEvals = 2000; % to prevent solver from stopping prematurely


tspan=[0 48];
x0=[0 0 0 0 0 0 100]; %set x0 as the experimental data on t1

k=[10.^(rand(1,4)-1) rand(1,1)*6]; % use random starting values for the parameter

lb=zeros(length(k),1); %lower boundary (we want non negative values)
up = 100*ones(length(k),1); %upper boundary (arbitrary?)
up(5) = 6; %maxium 6 hours "delay"

pp=[1:5]; %number of parameters

fp=k(pp);
up_sel = up(pp);
lb_sel = lb(pp);

%use lsqnonlin to minimize the distance function model_fit.m
[k_neu,resnorm,residual,exitflag] =lsqnonlin(@model_fit,fp,lb_sel,up_sel,options,...
    x_data_summed{:,:},x0,k,pp,t_data);

k(pp)=k_neu;


%model with estimated k and x0 initial conditions
[t,x]=ode45(@feedback_model,tspan,x0,[],k);


%plot and overlay data (sum up silenced and not silenced fractions)

% sum up predicted data
x_summed=x;
x_summed(:,1) = x_summed(:,1)+x_summed(:,7);
%sum up monoallelic silenced and not silenced
x_summed(:,2)=x_summed(:,2)+x_summed(:,3);
%sum up biallelic silenced and not silenced
x_summed(:,3)=x_summed(:,4)+x_summed(:,5)+x_summed(:,6);
%remove old columns
x_summed=x_summed(:,1:3);


plot(t_data-k(5),transposed_x_data_summed,'.'); % plot data points
hold on; % to keep the same plot and overlay curves
plot(t,x_summed); % plot model curves
legend('no Xist', 'monoallelic','biallelic', 'sim no xist', 'sim monoallelic', 'sim biallelic'); 
ylabel('% of cells');
xlabel('Time (hours)');
ax = gca;
ax.ColorOrder = mycolors;
hold off;

%display estimated k
k

% save plot as pdf
saveas(gcf, '../output/fitted_model.pdf')

% save parameters
filename='../output/fitted_parameters.txt';
writematrix(round(k, 2, "significant"), filename);

%% Determine how robust the fitted parameters are

parameter_index = [1:1:5]; % for each one of the 5 parameters
for param=parameter_index
    
% define testing range

    if param == 3 %k3 only has a lower bound
        step=0.05;
        spread=1.1;
        
    else
        step=0.05;
        spread=0.5;
    end
    
width_of_test= [-spread:step:spread];

    
    
k_iter=[10.^(rand(1,4)-1) rand(1,1)*6];

resnorm_iter=[];

for i=width_of_test
    k1i=k(param)*(1+i)
    
    pp=[1:5]; %number of parameters
    pp=pp(pp~=param);
    fp=k(pp);
    up_sel = up(pp);
    lb_sel = lb(pp);
    k_iter(param)=k1i;
    
    %use lsqnonlin to minimize the distance function model_fit.m
    [k_neu,resnorm,residual,exitflag] =lsqnonlin(@model_fit,fp,lb_sel,up_sel,options,...
    x_data_summed{:,:},x0,k_iter,pp,t_data);


k_iter(pp)=k_neu;


%model with estimated k and x0 initial conditions
[t,x]=ode45(@feedback_model,tspan,x0,[],k_iter);

%plot and overlay data (sum up silenced and not silenced fractions)

% sum up predicted data
x_summed=x;
x_summed(:,1) = x_summed(:,1)+x_summed(:,7);
%sum up monoallelic silenced and not silenced
x_summed(:,2)=x_summed(:,2)+x_summed(:,3);
%sum up biallelic silenced and not silenced
x_summed(:,3)=x_summed(:,4)+x_summed(:,5)+x_summed(:,6);
% remove old columns
x_summed=x_summed(:,1:3);


k_iter
resnorm_iter=[resnorm_iter; [i, resnorm]];
end

resnorm_iter

% save residuals for each tested value within the test range
filename=sprintf('../robustness_data/resnorm_k%.0f_%.2f-%.2f.xlsx', param, spread, step);
writematrix(resnorm_iter, filename);
end