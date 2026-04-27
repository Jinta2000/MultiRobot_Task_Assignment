%% Distributed Multirobot Task Assignment via Consensus ADMM %%
% Persistant Survillance % 
clc; clear; close all;

% parameters
N = 12;            % robots
m_s = 4;           % surveillance stations
m_c = 8;           % charging stations
T = 100;            % time steps

rho = 0.01;
max_iter = 10000;
step_size = 0.9;
b_max = 100;
b_min = 20;
X = zeros(N,m_s+m_c);

% Initialization
robot_pos = rand(N,2)*20;
battery = (rand(N,1))*100;

theta = linspace(0,2*pi,m_s+1); theta(end)=[];
radius = 4;
center = [10 10];

charging_pos = [...
    0 0; 0 20; 20 0; 20 20; ...
    10 0; 0 10; 20 10; 10 20];

Adj = ones(N) - eye(N);
Neighbors = cell(N,1);
for i=1:N
    Neighbors{i} = find(Adj(i,:) == 1);
end

battery_history = zeros(N,T);
min_battery = zeros(1,T);
avg_battery = zeros(1,T);

figure;

for t = 1:T
    
    surv_pos = zeros(m_s,2);
    for j=1:m_s
        surv_pos(j,:) = center + ...
            radius * [cos(theta(j)+0.05*t), ...
                      sin(theta(j)+0.05*t)];
    end
    
    task_pos = [surv_pos; charging_pos];
    m = size(task_pos,1);
    
    % Cost function update
    C = zeros(N,m);
    
    for i=1:N
        for j=1:m
            dist = norm(robot_pos(i,:) - task_pos(j,:));
            
            if j <= m_s
                % Surveillance → prefer closer + higher battery
                C(i,j) = dist^2 - 2*battery(i);
            else
                [~, idx] = max(X(:,j));
                % Charging → prefer low battery + closer
                C(i,j) = -(battery(idx) - battery(i))^2;
            end
        end
    end

    % MURID-TAP call
    X_relaxed = MURID_TAP(C, Neighbors, rho, max_iter);
    X = zeros(N,m);

    %----------- STEP 1: Assign Surveillance (1 robot per station) ----------
    assigned_robots = [];
    
    for j = 1:m_s
        [~, idx] = max(X_relaxed(:,j));   % best robot for station j
        
        while ismember(idx, assigned_robots)
            X_relaxed(idx,j) = -inf;      % remove already used robot
            [~, idx] = max(X_relaxed(:,j));
        end
        
        X(idx,j) = 1;
        assigned_robots = [assigned_robots idx];
    end
    
    %----------- STEP 2: Assign Remaining Robots to Charging ----------
    remaining = setdiff(1:N, assigned_robots);
    charging_taken = zeros(m_c,1);
    
    for i = remaining
        
        costs = X_relaxed(i, m_s+1:end);   % preference for charging
        [~, sorted_idx] = sort(costs, 'descend');
        
        for k = 1:m_c
            j = sorted_idx(k);
            
            if charging_taken(j) == 0
                X(i, m_s + j) = 1;
                charging_taken(j) = 1;
                break;
            end
        end
    end

    for i=1:N
        
        [~, task] = max(X(i,:));
        
        direction = task_pos(task,:) - robot_pos(i,:);
        robot_pos(i,:) = robot_pos(i,:) + step_size * direction;
       
        if task <= m_s
            battery(i) = max(b_min, battery(i) - 3);
        else
            battery(i) = min(b_max, battery(i) + 5);
        end
    end
    
    battery_history(:,t) = battery;
    min_battery(:,t) = min(battery);
    avg_battery(:,t) = mean(battery(:));
    
    % Visualization
    clf; hold on;
    
    scatter(surv_pos(:,1), surv_pos(:,2), 150, 'r', 'filled');
    scatter(charging_pos(:,1), charging_pos(:,2), 120, 'g', 'filled');
    scatter(robot_pos(:,1), robot_pos(:,2), 100, 'b', 'filled');
    
    title(['Persistent Surveillance (MURID) | t = ', num2str(t)]);
    legend('Surveillance','Charging','Robots');
    
    xlim([0 20]); ylim([0 20]);
    grid on;
    
    drawnow;
    pause(0.05);
end

% plot 
figure; hold on;

for i = 1:N
    plot(1:T, battery_history(i,:), 'LineWidth', 1.5);
end
plot(1:T, avg_battery(1,:), 'k','LineWidth', 2);
plot(1:T, min_battery(1,:), 'b','LineWidth', 2);

yline(b_min, '--r', 'Recharge Threshold');

title('Battery Levels of Robots');
xlabel('Time Step');
ylabel('Battery Level');
ylim([0 100]);
grid on;


%% MURID-TAP function %%
function X = MURID_TAP(C, Neighbors, rho, max_iter)

[N, m] = size(C);
beta = 0.05;

G = cell(N, 1);
for i = 1:N
    G{i} = zeros(N, m);
    G{i}(i, :) = ones(1, m);
end

x      = rand(m, N);   
y      = zeros(m, N);  
lambda = zeros(N, N);   
eta    = zeros(m, N);   
psi    = zeros(N, N); 
k = 0;
comp_time   = 0;

for k = 1:max_iter
    x_prev = x;
    y_prev = y;
    lam_prev = lambda;

    t_start = tic;

    x_new      = zeros(m, N);   
    y_new      = zeros(m, N);  
    lambda_new = zeros(N, N);

    for i = 1:N
        deg_i = length(Neighbors{i});
        ci    = C(i,:)';

        if deg_i == 0
            nu_i    = zeros(m, 1);
            kappa_i = zeros(N, 1);
        else
            sum_y_pairs   = zeros(m, 1);
            sum_lam_pairs = zeros(N, 1);
            for j = Neighbors{i}
                sum_y_pairs   = sum_y_pairs   + (y(:,i) + y(:,j));
                sum_lam_pairs = sum_lam_pairs + (lambda(:,i) + lambda(:,j));
            end
            coeff = 1 / (2 * rho * deg_i);
            nu_i    = coeff * ((1/N) * ones(m,1) - x(:,i) - eta(:,i) + rho * sum_y_pairs);
            kappa_i = coeff * (G{i} * x(:,i) - (1/N) * ones(N,1) - psi(:,i) + rho * sum_lam_pairs);
        end

        grad_fi = ci;
        kappa_term = G{i}' * kappa_i;
        kappa_nu_term = - nu_i + kappa_term;
        x_hat = x(:,i) - beta * (grad_fi + kappa_nu_term);

        x_new(:,i) = min(1, max(0, x_hat));

        if deg_i == 0
            nu_new_i    = zeros(m, 1);
            kappa_new_i = zeros(N, 1);
        else
            coeff = 1 / (2 * rho * deg_i);
            nu_new_i    = coeff * ((1/N) * ones(m,1) - x_new(:,i) - eta(:,i) + rho * sum_y_pairs);
            kappa_new_i = coeff * (G{i} * x_new(:,i) - (1/N) * ones(N,1) - psi(:,i) + rho * sum_lam_pairs);
        end

        y_new(:,i)      = max(0, nu_new_i);         % y_i = max(0, nu_i(x_new))
        lambda_new(:,i) = kappa_new_i;              % lambda_i = kappa_i(x_new)
    end
    eta_new = zeros(m, N);
    psi_new = zeros(N, N);
    for i = 1:N
        eta_new(:,i) = eta(:,i);
        psi_new(:,i) = psi(:,i);
        for j = Neighbors{i}
            eta_new(:,i) = eta_new(:,i) + rho * (y_new(:,i) - y_new(:,j));
            psi_new(:,i) = psi_new(:,i) + rho * (lambda_new(:,i) - lambda_new(:,j));
        end
    end

    comp_time   = comp_time + toc(t_start);

    x      = x_new;
    y      = y_new;
    lambda = lambda_new;
    eta    = eta_new;
    psi    = psi_new;
    res = 0;
    for i = 1:N
        dx = norm(x(:,i) - x_prev(:,i));
        res = max(res, dx / (norm(x_prev(:,i)) + 1e-12));
        x_plot(i,:) = x(:,i)';
    end

    if res < 10^-5
        break;
    end
end
X = zeros(N, m);
for i = 1:N
    X(i,:) = x(:,i)';
end

end