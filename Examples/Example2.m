% t_nonlinearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = f(x(k),u(k)) + w(k)
% The approach is based on [1]. The nonlinear system is found in [2]
%
%
%
% Syntax:
%    t_nonlinearDT
%
% Inputs:
%    no
%
% Outputs:
%    no
%
% Example:
%
% References:

%
% Author:       Alireza Naderi Akhormeh
% Written:      07-May-2025
% Last update:
% Last revision:---


%------------- BEGIN CODE --------------

clear
close all
rand('seed',1);
params.dt =0.015;
NN=4;
params.tFinal = params.dt*NN;
params.tStart = 0;
%input set
params.U = zonotope([[1.1;-1.3],diag([0.1;0.2])]);

%initial set
params.R0 = zonotope([[-1.9;-20],diag([0.005;0.3])]);

% dimension of x
options.dim_x=2;

%Number of trajectories
initpoints=35;
%Number of time steps
Win_Step = 5;
steps=15;

%Totoal number of samples
totalsamples = steps*initpoints;

%noise zonotope
sigma_v = 0.0001;
options.W = zonotope(zeros(options.dim_x,1),sigma_v*diag(ones(options.dim_x,1))); % disturbance

dim_x = options.dim_x;
GW = cell(dim_x * totalsamples, 1);
index = 1;
for row = 1:dim_x
    for col = 1:totalsamples
        G = zeros(dim_x, totalsamples);
        G(row, col) = sigma_v;
        GW{index} = G;
        index = index + 1;
    end
end
options.Wmatzono = matZonotope(zeros(dim_x, totalsamples), GW);

% Reachability Settings

options.zonotopeOrder = 100;
options.tensorOrder = 2;
options.errorOrder = 5;


% System Dynamics
fun = @(x,u) cstrDiscr(x,u,params.dt);

%input random sample points
for i=1:totalsamples
    u(:,i) = randPointExtreme(params.U);
end

%get state trajectories
index=1;
for j=1:options.dim_x:initpoints*options.dim_x
 
    x(j:j+options.dim_x-1,1) = randPoint(params.R0);

        x_free(j:j+options.dim_x-1,1) = x(j:j+options.dim_x-1,1);
    for i=1:steps
        x_free(j:j+options.dim_x-1,i+1) = fun(x(j:j+options.dim_x-1,i),u(:,index));
        x(j:j+options.dim_x-1,i+1) = fun(x(j:j+options.dim_x-1,i),u(:,index)) +randPoint(options.W);
        index=index+1;
    end
end


%combine trajectories
index_0 =1;
index_1 =1;
for j=1:options.dim_x:initpoints*options.dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x(j:j+options.dim_x-1,i);
        x_free_vec_1(:,index_1) = x_free(j:j+options.dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_free_vec_0(:,index_0) = x_free(j:j+options.dim_x-1,i);
        x_meas_vec_0(:,index_0) = x(j:j+options.dim_x-1,i);
        index_0 = index_0 +1;
    end
end

stepsLip=1;
initpointsLip=1000;
[gamma,L]= compLipConst(fun,params.U,params.R0,stepsLip,initpointsLip,options.dim_x);
eps(1)= L(1) .* gamma/2;
eps(2)= L(2) .* gamma/2;
options.Zeps = zonotope([zeros(2,1),diag(eps)]);
Zeps=options.Zeps;

% X_+ is X_1T
% X_- is X_0T
options.U_full = u(:,1:totalsamples);
options.X_0T = x_meas_vec_0(:,1:totalsamples);
options.X_1T = x_meas_vec_1(:,1:totalsamples);

% define system
sysDisc = nonlinearSysDT('stirredTankReactor',fun,params.dt,2,2);

%% Reachability Analysis ---------------------------------------------------
% compute model based reachability (R) and data driven one (R_data)
tic
steps = steps-4;
params.R0 = zonotope([x(1:2,steps),diag([0.005;0.3])]);

% model based Reachability Analysis
R_model= reach(sysDisc,params,options);
% data-driven(LS) based Reachability Analysis
R_data = reach_LS(params,options);

tComp = toc;
disp("Computation time: " + tComp);


%% ZRLS
n = 5;
m = 2;
c = zeros(n, m);
initial_uncertainty_radius_per_element = 1.5;

Gcell = {};
for i = 1:n
    for j = 1:m
        G_v = zeros(n, m);
        G_v(i,j) = initial_uncertainty_radius_per_element/sqrt(n*m);
        Gcell{end+1} = G_v;
    end
end


G = Gcell ;
totalsteps = 3;
P = eye(n) * 1e7;
Q_v = (sigma_v)^2; % Scalar variance for single measurement noise
q_v = [sqrt(Q_v);sqrt(Q_v)];
lamda = 0.96;

for i = 1:steps
    phi = [ones(1,1);options.X_0T(:,i);options.U_full(:,i)]';
    y = options.X_1T(:,i)';

    [c_next, G_next, P_next, theta_next, kappa] = ZRLS(c, G, P, phi, y, Q_v,lamda,sigma_v,m);
    c = c_next;
    G = G_next;
    theta = theta_next;
    P = P_next;

end

IAB = theta;
oneMat = repmat([1],1,size(options.U_full(:,steps-Win_Step:steps),2));

V =  options.X_1T(:,steps-Win_Step:steps) + -1*(IAB.center*[ oneMat;options.X_0T(:,steps-Win_Step:steps);options.U_full(:,steps-Win_Step:steps) ]);


VInt = intervalMatrix(V);
leftLimit = infimum(VInt);
rightLimit = supremum(VInt);

V_one= zonotope(interval(min(leftLimit')',max(rightLimit')'));
R_data_RLS{1} = params.R0;

for k=1:NN

    R_data_RLS{k} = reduce(R_data_RLS{k},'girard',100);

    R_data_RLS{k + 1} = (IAB)*cartProd([ 1],cartProd(R_data_RLS{k},params.U)) +V_one+options.W+Zeps ;
end
params.tStart = 0;
t = params.tStart:params.dt:params.tFinal;

timePoint_data.set = R_data_RLS(2:end)';
timePoint_data.time = num2cell(t(2:end)');

R_data_RLS1 = reachSet(timePoint_data);


%% Visualization -----------------------------------------------------------
    figure('Renderer', 'painters', 'Position', [10 10 700 900]);hold on; box on;


% plot initial set
handleX0=plot(params.R0,[1,2],'k-','LineWidth',2);

% plot model based reachable set
handleModel=plot(R_model,[1 2],'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b','LineWidth',0.7);

% plot data driven reachable set
 handleData =  plot(R_data,[1 2], ...
     'EdgeColor','r', ...
     'LineWidth',0.7, ...
     'FaceColor','none');   % outline only, no fill

 handleDataRLS =  plot(R_data_RLS1,[1 2], ...
     'EdgeColor','g', ...
     'LineWidth',0.7, ...
     'FaceColor','none');   % outline only, no fill

% After all plots are done
   ylim([-11.35 -9.75]); % Example: from -2 to 5 in x_2


% formatting
xlabel('$x_1$','Interpreter', 'latex', 'FontSize', 30);
ylabel('$x_2$','Interpreter', 'latex', 'FontSize', 20);

% skip warning for extra legend entries
warOrig = warning; warning('off','all');
legend([handleX0,handleModel,handleData,handleDataRLS],...
    'Initial set ($\mathcal{X}_0$)','Reachable Sets ($\mathcal{R}_k$) from model','Reachable Sets ($\bar{\mathcal{R}}_k$) from data using LS','Reachable Sets ($\mathcal{R}^\prime_k$) from data using EF-ZRLS','Location','northwest','Interpreter','latex');
warning(warOrig);
legend boxoff

ax = gca;
ax.FontSize = 16;
%set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+0.04 ax_width ax_height-0.05];
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% example completed
completed = 1;


%------------- END OF CODE --------------