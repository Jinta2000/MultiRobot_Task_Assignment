%% Distributed Multirobot Task Assignment via Consensus ADMM %%
% Linear Multirobot Task Assignment Problem %
clear; clc; close all;

Repeat = 30;               % repeating number
N = 10;                    % No. of robots
m = N;                    % No. of tasks
rho = 0.15;               % step size
max_iter = 10000;          % for limiting No. of iterations
tol = 1e-4;               % tolerance for converging
Xedges = 1;               % No. of edges extra required for graph
beta = 0.05*ones(N,1);    % beta value for MURID-TAP

time = zeros(3,Repeat);
iter_no = zeros(3,Repeat);

fprintf('--- N=%d robots, m=%d tasks ---\n', N, m);
for rep = 1:Repeat
    fprintf('--------------- Iteration -%d -------------------\n',rep);
    C = rand(m,N);
    IN = eye(N);
    
    %% Connected graph creation randomly %%
    while true 
        edges = [];
        for i = 2:N
            j = randi(i-1);
            edges = [edges; i j];
        end
        
        k = 0;
        while k < Xedges
            u = randi(N);
            v = randi(N);
            
            if u ~= v   
                edges = [edges; u v];
                k = k + 1;
            end
        end
        
        G = graph(edges(:,1), edges(:,2));
        Adj = full(adjacency(G));
        if max(conncomp(G)) == 1
            %plot(G)
            break;
        end
    end
    
    %  Neighbourhood sets
    Ni = cell(N, 1);
    nE = size(edges,1);
    for i = 1:N, Ni{i} = []; end
    for e = 1:nE
        i = edges(e,1); j = edges(e,2);
        Ni{i} = [Ni{i}, j];
        Ni{j} = [Ni{j}, i];
    end
    for i = 1:N
        Ni{i} = unique(Ni{i}); 
    end

    %% Iterations for MUR-TAP %%
    fprintf('----------MUR-TAP----------');

    % Initialization
    c_all = zeros(N*N,N);
    x = ones(N*N,N)/N;
    q = zeros(N*N,N);
    x_mem_mur = cell(max_iter, 1);
    
    for i = 1:N
        ci_all = zeros(N*N, 1);
        ci_all((i-1)*N+1 : i*N) = C(:,i);
        c_all(:,i) = ci_all;
    end
    mur = tic;
    for k = 1:max_iter
        x_mem_mur{k} = x(:,N);
        x_prev = x;
    
        % X updation
        for i = 1:N

            cardi = numel(Ni{i});
            sum_mid = zeros(N*N,1);
            for j = Ni{i}
                sum_mid = sum_mid + (x_prev(:,i) + x_prev(:,j));
            end
            H = rho * 2 * cardi * eye(N*N);
            f = c_all(:,i) + q(:,i) - rho*sum_mid;
            A = [];
            Q = [];
            for l = 1:N
                A = [A -IN];
                if l == i
                   Q = [Q ones(N,1)']; 
                else
                   Q = [Q zeros(N,1)']; 
                end
            end
            A =[A; Q; -Q; -eye(N*N);eye(N*N)];
            b = [-ones(N,1);1;-1;zeros(N*N,1);ones(N*N, 1)];
            options = optimoptions('quadprog','Display','off');
            [x_val, fval, exitflag, output] = quadprog(H, f, A, b, [], [], [], [], x_prev(:,i), options);
            if exitflag <= 0
                x_val = zeros(N*N, 1);
                disp("Solution is not existing")
            end
            x(:,i) = x_val;
            error_mur(i) = norm(x_prev(:,i) - x(:,i))/(norm(x(:,i))+ 1e-12);
        end
    
        % q updation
        for i = 1:N
            diff = ones(N*N, 1);
            for j = Ni{i}
                diff = diff + (x(:,i) - x(:,j));
            end
            q(:,i) = q(:,i) + rho * diff;
        end

        if (error_mur < (tol*ones(N,1))) 
            time_mur = toc(mur);
            x_opt_mur = [];
            for i = 1:N
                x_opt_mur = [x_opt_mur x(((i-1)*N+1 : i*N),1)];
            end
            x_opt_mur
            fprintf("No. of iterations: %d\n", k)
            fprintf("Duration in ms : %.2f\n",time_mur*10^3)

            time(1, rep) = time_mur*10^3;
            iter_no(1, rep) = k;
            break;
        end
    end

    x_it = 1: k;
    y_val = [];
    for i = 1:k
        x_diff = x(:,N) - x_mem_mur{i};
        y_val = [y_val norm(x_diff)/ norm(x(:,N))];
    end

    figure;
    plot(x_it,y_val,'-','Color',"r","LineWidth",2);
    hold on

    %% Iterations for MURD-TAP %%
    fprintf('----------MURD-TAP----------');
    % Initialization
    y = zeros(N, m);          
    lambda = zeros(N, 1);     
    r = zeros(N, m);      
    x = zeros(N, m);
    time_murd = 0;
    murd = tic;

    for k = 1:max_iter
        x_mem_murd_2{k} = round(x,1);
        
        y_old = y;
        for i = 1:N
            
            deg = length(Ni{i});
           
            dim = m + 1;
            
            H = zeros(dim);
            H(1:m,1:m) = 2*rho*deg*eye(m);
            
            f = zeros(dim,1);
            
            f(1:m) = r(i,:)' - (1/N)*ones(m,1);
            f(end) = 1;   

            avg_sum = zeros(m,1);
            for j = Ni{i}
                avg_sum = avg_sum + (y_old(i,:)' + y_old(j,:)')/2;
            end
            
            f(1:m) = f(1:m) - 2*rho*avg_sum;
           
            A1 = [-eye(m), zeros(m,1)];
            b1 = zeros(m,1);
            
            A2 = [eye(m), -ones(m,1)];
            b2 = C(i,:)';
            
            A = [A1; A2];
            b = [b1; b2];
           
            options = optimoptions('quadprog','Display','off');
            z = quadprog(H, f, A, b, [], [], [], [], [], options);
            
            y(i,:) = z(1:m)';
            lambda(i) = z(end);
            
        end
        
        for i = 1:N
            
            for j = Ni{i}
                r(i,:) = r(i,:) + rho*(y(i,:) - y(j,:));
            end
        end
        
        err = norm(y - y_old, 'fro');
        x_prev = x; 

        for i = 1:N
            
            deg = length(Ni{i});
            
            sum_y = zeros(m,1);
            for j = Ni{i}
                sum_y = sum_y + (y(i,:)' + y(j,:)');
            end
            
            const = (1/N)*ones(m,1) - r(i,:)' + rho * sum_y;
            
            A_nu = -(1/(2*rho*deg)) * eye(m);
            b_nu = (1/(2*rho*deg)) * const;
            
            H = rho*deg*(A_nu'*A_nu);
            f = C(:,i) + rho*deg*(A_nu'*b_nu);
            
            Aeq = ones(1,m);
            beq = 1;
            
            lb = zeros(m,1);
            ub = ones(m,1);
            
            options = optimoptions('quadprog','Display','off');
            x(i,:) = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options)';
            
        end

        if err < tol
            x_mem_murd_2{k} = round(x,0);
            break;
        end
        
    end
    time_murd = time_murd + toc(murd)* 10^3;

    round(x)

    x_it = 1: k;
    y_val = [];
    for i = 1:k
        x_diff = x_opt_mur - x_mem_murd_2{i};
        y_val = [y_val norm(x_diff)/ norm(x_opt_mur)];
    end
    y_val = [y_val];
    plot(x_it,y_val,'-','Color',"g","LineWidth",2);
    
    fprintf ("Duration in ms : %.2f\n",time_murd)
    fprintf ("No. of iterations : %d\n", k)

    time(2, rep) = time_murd;
    iter_no(2, rep) = k;
    
    %% Iterations for MURID-TAP %%
    fprintf('----------MURID-TAP----------');
    
    G = cell(N, 1);
    for i = 1:N
        G{i} = zeros(N, m);
        G{i}(i, :) = ones(1, m);
    end

    x      = ones(m, N) /m;   
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
            deg_i = length(Ni{i});
            ci    = C(:,i);

            if deg_i == 0
                nu_i    = zeros(m, 1);
                kappa_i = zeros(N, 1);
            else
                sum_y_pairs   = zeros(m, 1);
                sum_lam_pairs = zeros(N, 1);
                for j = Ni{i}
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
            x_hat = x(:,i) - beta(i) * (grad_fi + kappa_nu_term);

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
            for j = Ni{i}
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
        x_mem_murid{k} = x_plot';

        if res < tol
            break;
        end
    end
    x_opt = zeros(N, m);
    for i = 1:N
        x_opt(i,:) = x(:,i)';
    end
    fprintf("Time taken in ms : %.2f\n",comp_time* 10^3)
    fprintf("NO. of Iterations : %d\n",k)
    x_opt'

    time(3, rep) = comp_time* 10^3;
    iter_no(3, rep) = k;

    x_it = 1: k-1;
    y_val = [];
    for i = 1:k-1
        x_diff = x_opt_mur - x_mem_murid{i};
        y_val = [y_val norm(x_diff)/ norm(x_opt_mur)];
    end
    y_val = [y_val];
    plot(x_it,y_val,'-','Color',"y","LineWidth",2);
    title("C-ADMM Error Analysis - N = m = 5");
    xlabel("No. of iterations");
    ylabel("Relative error");
    legend('MUR-TAP', 'MURD-TAP', 'MURID-TAP')
    
end