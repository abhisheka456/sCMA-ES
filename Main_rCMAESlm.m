%% ============ A Reference Vector-Based Simplified Covariance ============
%% =============== Matrix Adaptation Evolution Strategy for ===============
%% ==================== Constrained Global Optimization ===================
% Should you have any queries, please contact
% Dr. Abhishek Kumar
% email-id: abhishek.kumar.eee13@iitbhu.ac.in
%%=========================================================================
%% Main routine for running the CMA-ES variants on the 2017 CEC Benchmarks
clear all;
clc;   


        global  initial_flag
        initial_flag = 0;
        global opt_f
        global func_num
        % input -- initial (fixed) strategy parameter setting
        input.budget            = 5e5;
        input.maxIter           = 5e4;
        input.delta             = 10^-4;                    % error margin for equality constraints
        input.runs              = 25;                       % number of repetitions
       


   

        problem.constr_fun_name = 'cec2006';                    % objective function class

        strategy='rCMAESlm';                                    % algorithm variant
        D = [13 20 10 5 4 2 10 2 7 8 2 3 5 10 3 5 6 9 15 24 7 22 9 2];

        disp(['ES variant ' strategy ' --- 25 independent runs!'])

        for k=1:24                                              % on each CEC 2018 test function do
            func_num=k;                                         % test function number
            initial_flag = 0;
            [L,U,opt_f,ineq,eq] = get_Fun_info(func_num,D(func_num));
            problem.gn = ineq;
            problem.hn = eq;
            input.dim = D(func_num);
            input.lb = L;
            input.ub = U;
            problem.opt_f = opt_f;

            problem.upper_bounds    = input.ub';
            problem.lower_bounds    = input.lb';
      
            % CMA specific paprameters
            input.lambda            = 2*D(func_num);                          % population size
            input.sigma             = 1;                            % initial mutation strength
            input.mu                = floor(input.lambda/3);        % parental ppopulation size
            input.weights = log(input.mu+1/2)-log(1:input.mu)';     % muXone array for weighted recombination
            input.weights = input.weights./sum(input.weights);      % normalize recombination weights array
            input.mueff=1/sum(input.weights.^2);                    % variance-effectiveness of sum w_i x_i
            input.cc = (4+input.mueff/D(func_num))/(D(func_num)+4+2*input.mueff/D(func_num));
            input.cs = (input.mueff+2) / (D(func_num)+input.mueff+5);         % t-const for cumulation for sigma control
            input.c1 = 2 / ((D(func_num)+1.3)^2+input.mueff);                 % learning rate for rank-one update of M
            input.cmu = min(1-input.c1, 2 * (input.mueff-2+1/input.mueff) / ((D(func_num)+2)^2+input.mueff));     % and for rank-mu update of M
            input.damps = 1 + 2*max(0, sqrt((input.mueff-1)/(D(func_num)+1))-1) + input.cs;  
            SP = zeros(input.runs,1);
            for j=1:input.runs                                  % perform multiple runs on each test function
                eval(['[tab,SPj]=' strategy '(problem,input,func_num);']); % run epsMAgES
                FitT(j,:)=[tab(1,:) tab(2,:) tab(3,:)];
                SP(j) = SPj;
                disp([func_num j tab(1,:) tab(2,:) tab(3,:) SPj]);
            end
            S(k) = mean(SP);
            Tab=build_stats(FitT,input);                        % build statistics according to specification of the CEC2018 competition 
            Stats(k,:)=[func_num Tab(1,:) Tab(2,:) Tab(3,:)];
        end
              
        filename = strcat(['rCMAESlm_on_CEC2017_D' num2str(input.dim) '.txt']);
        fileID = fopen(filename,'w');             
        fprintf(fileID,'%14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\n',Stats);
        fclose(fileID)
        
