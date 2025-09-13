% t_linearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = Ax(k) + Bu(k) + w(k)
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author:       Alireza Naderi Akhormeh
% Written:      September-2025
% Last update:
% Last revision:---

%------------- BEGIN CODE --------------

rand('seed',1);

clear all close all
clc
%% system dynamics
dim_x = 5;
dim_u = 1;
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(5,1);
C = [1,0,0,0,0];
D = 0;
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.1;
sys_d = c2d(sys_c,samplingtime);
c0 = [sys_d.A,sys_d.B];
%Number of trajectories
initpoints =1;
%Number of time steps
steps = 20;
totalsamples = initpoints*steps;
%% initial set and input
X0 = zonotope(ones(dim_x,1),0.1*diag(ones(dim_x,1)));
U = zonotope(10,2.25);

%noise zontope W
sigma_v = 0.0005;
W = zonotope(zeros(dim_x,1),sigma_v*diag(ones(dim_x,1)));
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
Wmatzono = matZonotope(zeros(dim_x, totalsamples), GW);

% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(i) = randPoint(U);
end


%simulate the system to get the data
x0 = X0.center;
x(:,1) = x0;
index=1;
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    for i=1:steps
        utraj(j,i) = u(index);
        x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*u(index) + randPoint(W);
        index=index+1;
    end

end


% concatenate the data trajectories
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

% X_+ is X_1T
% X_- is X_0T
U_full = u_mean_vec_0(:,1:totalsamples); %same as u
X_0T = x_meas_vec_0(:,1:totalsamples);
X_1T = x_meas_vec_1(:,1:totalsamples);

%get state trajectories
index=1;
for j=1:dim_x:initpoints*dim_x
    x_free(j:j+dim_x-1,1) = randPoint(X0);
    for i=1:steps
        x_free(j:j+dim_x-1,i+1) = sys_d.A*x_free(j:j+dim_x-1,i) + sys_d.B*u(:,index);
        index=index+1;
    end
end

%combine trajectories
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_free_vec_1(:,index_1) = x_free(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_free_vec_0(:,index_0) = x_free(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end
x_free_vec_1= normalize(x_free_vec_1);
x_free_vec_0 = normalize(x_free_vec_0);

%% Compute data-driven model using LS method

X1W_cen =  X_1T - Wmatzono.center;
X1W = matZonotope(X1W_cen,Wmatzono.generator);

% set of A and B
AB = X1W  *pinv([X_0T;U_full]);

% Compute Lip Constant
steps1 = 1;
initpoints = 1000;
[gamma, L] = compLipConstLinear(sys_d,AB,U,X0,steps1,initpoints,dim_x);

for i = 1:dim_x
    eps(i)= L(i) .* gamma/2;
end
Zeps_LS = zonotope([zeros(dim_x,1),diag(eps)]);

%% ZRLS
n = dim_x + dim_u;
m = dim_x;

% Initial center
c = zeros(n, m);
initial_uncertainty_radius_per_element = 0.9;

Gcell = {};
for i = 1:n
    for j = 1:m
        G0 = zeros(n, m);
        G0(i,j) = initial_uncertainty_radius_per_element/sqrt(n*m);
        Gcell{end+1} = G0;
    end
end


G = Gcell ;
totalsteps = 2;
% Initial Covarience
P = eye(n) * 1e7;
Q_v = (sigma_v)^2; % Scalar variance for single measurement noise
q_v = [sqrt(Q_v);sqrt(Q_v);sqrt(Q_v);sqrt(Q_v);sqrt(Q_v)];
% Forgetting factor
lamda = 1;
l = 1;

for i = 1:length(G)

    G_NEXT{i} = G{i}';
end

AB_ZRLS = matZonotope(c', G_NEXT);

%% compute next step sets from model / data
for i = 1:totalsamples
    % Compute Lip Constant
    steps1 = 1;
    initpoints = 1000;
    [gamma, L] = compLipConstLinear(sys_d,AB_ZRLS,U,X0,steps1,initpoints,dim_x);

    for j = 1:dim_x
        eps(j)= L(j) .* gamma/2;
    end

    Zeps_ZRLS = zonotope([zeros(dim_x,1),diag(eps)]);

    % set number of steps in analysis
    X_model = cell(totalsteps+1,1);
    X_data = cell(totalsteps+1,1);
    % init sets for loop
    X_model{1} = X0; X_data{1} = X0;X_data_zrls{1} = X0;

    for j=1:totalsteps

        % 1) model-based computation
        X_model{j,1}=reduce(X_model{j,1},'girard',400);
        X_model{j+1,1} = sys_d.A * X_model{j} + sys_d.B * U+W;

        % 2) Data Driven approach LS
        X_data{j,1}=reduce(X_data{j,1},'girard',400);
        X_data{j+1,1} = AB * (cartProd(X_data{j},U)) + Zeps_LS +W;

        % 3)  Data Driven approach ZRLS
        X_data_zrls{j,1}=reduce(X_data_zrls{j,1},'girard',400);
        X_data_zrls{j+1,1} = AB_ZRLS * (cartProd(X_data_zrls{j},U)) +Zeps_ZRLS +W;


    end

    phi = [X_0T(:,i);U_full(i)]';
    y = X_1T(:,i)';

    [c_next, G_NEXT, P_next, theta_next, kappa] = ZRLS(c, G, P, phi, y, Q_v,lamda,sigma_v,m);
    c = c_next;
    G = G_NEXT;
    theta = theta_next;
    P = P_next;

    AB_ZRLS = theta;


    %% Visualization

    projectedDims = 3;


    % Plot initial set
    int_ = interval(X_model{totalsteps});
    lower = int_.inf(projectedDims);
    upper = int_.sup(projectedDims);
    model_lower(l) = lower;
    model_upper(l) = upper;


    int_ = interval(X_data{totalsteps});
    lower = int_.inf(projectedDims);
    upper = int_.sup(projectedDims);
    LS_lower(l) = lower;
    LS_upper(l) = upper;




    int_ = interval(X_data_zrls{totalsteps});
    lower = int_.inf(projectedDims);
    upper = int_.sup(projectedDims);
    ZRLS_lower(l) = lower;
    ZRLS_upper(l) = upper;


    l = l +1

end

% Create a figure and explicitly number it as 1
h_fig1 = figure(1);

% Now you can set its properties using the handle
set(h_fig1, 'Renderer', 'painters', 'Position', [10 10 700 900]);
% Set axis
model_lower = [model_lower model_lower(end)];
model_upper = [model_upper model_upper(end)];
LS_lower = [LS_lower LS_lower(end)];
LS_upper = [LS_upper LS_upper(end)];
ZRLS_lower = [ZRLS_lower ZRLS_lower(end)];
ZRLS_upper = [ZRLS_upper ZRLS_upper(end)];
time = 0:totalsamples;
handleModelArray = plot(time,model_lower, '-k', 'LineWidth', 1.2);hold on
handleModelArray = plot(time,model_upper, '-k', 'LineWidth', 1.2);hold on
handleLSArray = plot(time,LS_lower, '- r', 'LineWidth', 1.2);hold on
handleLSArray = plot(time,LS_upper, '- r', 'LineWidth', 1.2);hold on
handleZRLSArray = plot(time,ZRLS_lower, '-g', 'LineWidth', 1.2);hold on
handleZRLSArray = plot(time,ZRLS_upper, '-g', 'LineWidth', 1.2);hold on
% Label plot
xlabel('$i$', 'Interpreter', 'latex', 'FontSize', 20);

ylabel(sprintf('$x_{%d}$', projectedDims), ...
    'Interpreter', 'latex', 'FontSize', 20);

% Skip warning for extra legend entries
warOrig = warning;
warning('off', 'all');

% Use representative handles for the legend
legend([handleModelArray, handleLSArray, handleZRLSArray], ...
    'Reachable Sets ($\mathcal{R}_k$) from model','Reachable Sets ($\bar{\mathcal{R}}_k$) from data using LS ','Reachable Sets ($\hat{\mathcal{R}}_k$) from data using ZRLS ','Location','northwest','Interpreter','latex');
ylim([-3 4.1]); % Example: from -2 to 5 in x_2

% warning(warOrig);
legend boxoff
% Formatting plot
ax = gca;
ax.FontSize = 22;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3) - 0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom+0.03 ax_width ax_height-0.03];

