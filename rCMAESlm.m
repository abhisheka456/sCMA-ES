%% ============ A Reference Vector-Based Simplified Covariance ============
%% =============== Matrix Adaptation Evolution Strategy for ===============
%% ==================== Constrained Global Optimization ===================
% Should you have any queries, please contact
% Dr. Abhishek Kumar
% email-id: abhishek.kumar.eee13@iitbhu.ac.in
%%=========================================================================
%% Main file of algorithm
function [out,global_best] = rCMAESlm(problem,input,CEC_fun_no)

    dim         = length(problem.lower_bounds);
    sigma       = input.sigma;
    mu          = input.mu;
    lambda      = input.lambda;
    newpop.y    = zeros(dim,lambda);    % initialize new population matrix (n times NP)
	newpop.f    = 0;
	newpop.conv = 0;
	evals.fun   = 0;

    g           = 0;
    termination = 0;
    
    %% Initialize dynamic (internal) strategy parameters and constants
    ps      = zeros(dim,1);                                 % evolution paths for sigma
    pc      = zeros(dim,1);
    MM      = eye(dim);
    oMM     = eye(dim);
    sqrt_s  = sqrt(input.cs*(2-input.cs)*input.mueff);      % factor in path update
    sqrt_c  = sqrt(input.cc*(2-input.cc)*input.mueff);
    
      
    constraint_number = problem.gn(CEC_fun_no) + problem.hn(CEC_fun_no);
    %% Initialize random population of lambda candidate solutions
    for k=1:lambda
        newpop.y(:,k)   = problem.lower_bounds...
                            +(problem.upper_bounds-problem.lower_bounds).*rand(dim,1); 
        [fval, gv, hv]  = feval(problem.constr_fun_name,newpop.y(:,k)',CEC_fun_no);
        newpop.f(k)     = fval;                             % fitness vector 
        newpop.conv(k)  = sum([sum(gv.*(gv>0)), sum(abs(hv).*(abs(hv)>input.delta))])./constraint_number;  % mean constraint violations)
		newpop.g(:,k) = gv;
        newpop.h(:,k) = hv;
        evals.fun       = evals.fun + 1;                    % count objective function evaluations
    end

    %% Initial parameter for epsilon level ordering
    [TC,Eg,Eh,CPg,CPh] = ETAintialization(newpop);
    Eg0 = Eg;
    Eh0 = Eh;
             
    
    %% Rank initial population 
    [ranking]           = eta_sort_ref(newpop,Eg,Eh);
    ParentPop           = newpop.y(:,ranking(1:mu));
    yParent             = sum(ParentPop,2)./mu;
    
    %% Best individual of current population
    best_ind            = ranking(1);
    best_val            = newpop.f(best_ind);       % best fitness value of current population
    best_y              = newpop.y(:,best_ind);
    best_conv           = newpop.conv(:,best_ind);
          
    %% Best solution found so far
    global_best.y       = best_y; 				% best solution found so far
	global_best.val     = best_val;
    global_best.conv    = best_conv;
    %% local intialization 
    local_best.y       = best_y; 				% best solution found so far
	local_best.val     = best_val;
    local_best.conv    = best_conv;
    tolY               = 0;
    tolF               = inf;
    %% Upper mutation strength bound
    sigmaMAX            = 100;
    flag1               = 0;
    flag2               = 0;
    restart             = 0;
    while ~termination
        if restart == 1
            [sigma, g, ps, pc, MM, oMM, FEs, Eg0, Eh0, yParent, TC, Eg, Eh, CPg, CPh, input ] = restarts(problem,input,CEC_fun_no);
            evals.fun       = evals.fun + FEs;
            restart         = 0;
            sqrt_s          = sqrt(input.cs*(2-input.cs)*input.mueff);      % factor in path update
            sqrt_c          = sqrt(input.cc*(2-input.cc)*input.mueff);
            lambda          = input.lambda;
            mu              = input.mu;
        end
        %% Compute pseudo inverse of transformation matrix MM
        piM       = pinv(MM,1e-12);
        repi      = zeros(1,lambda);   
        %% Sample lambda offspring distributed around yParent
        newpop.z  = randn(dim,lambda);
        newpop.d  = MM*newpop.z;
        newpopy   = repmat(yParent,1,lambda) + sigma.*newpop.d;
        
        for k=1:lambda 
                                        
            %% ensure box constraint satisfaction 
            newpop.y(:,k)   = keep_range(newpopy(:,k),problem.lower_bounds,problem.upper_bounds);
            % log whether an offspring was repaired
            repi(k) = (sum(newpop.y(:,k)~=newpopy(:,k)) > 0);
            
            %% evaluation  
            [fval, gv, hv]  = feval(problem.constr_fun_name,newpop.y(:,k)',CEC_fun_no);
            evals.fun       = evals.fun + 1;
            fitval          = fval;                                                                 % fitness vector 
            convio          = sum([sum(gv.*(gv>0)), sum(abs(hv).*(abs(hv)>input.delta))])./(problem.gn(CEC_fun_no) + problem.hn(CEC_fun_no));                  % vector of the corresponding constraint violations)
            
            %% Levenberg-Marquadt repair step
            if mod(g,dim)==0 && rand(1) <=0.2 && convio > 0
               [new_mutant,FE] = Broyden_LM(problem, newpop.y(:,k),gv,hv,CEC_fun_no);
               [fval, gv, hv]  = feval(problem.constr_fun_name,new_mutant',CEC_fun_no);
               fitval          = fval;                                                                 % fitness vector 
               convio          = sum([sum(gv.*(gv>0)), sum(abs(hv).*(abs(hv)>input.delta))])./(problem.gn(CEC_fun_no) + problem.hn(CEC_fun_no));
               evals.fun = evals.fun + FE;
               
               newpop.y(:,k) = new_mutant;
               repi(k)= repi(k)+1;                     % log repaired offspring individuals
            end
            newpop.f(k)     = fitval;                                                                 % fitness vector 
            newpop.conv(k)  = convio;
            newpop.g(:,k) = gv;
            newpop.h(:,k) = hv;
            if repi(k) > 0
                newpop.d(:,k) = (newpop.y(:,k)    - yParent)./sigma;
                newpop.z(:,k) = piM*newpop.d(:,k);
            end
        end

        
        %% Rank current population w.r.t. recent epsilon level        
        [ran]   = eps_sort(newpop.f,newpop.conv,0);
        [ranking]  = eta_sort_ref(newpop,Eg,Eh);
        %% Best individual of current population
        best_ind    = ran(1);
        tolF        = abs(best_val-newpop.f(best_ind));
        best_val    = newpop.f(best_ind);            
        best_y      = newpop.y(:,best_ind);          
        best_conv   = newpop.conv(:,best_ind);
        
        %% Recombination of mutation vectors and centroid update
        parent_z = newpop.z(:,ranking(1:mu)) * input.weights;  
        parent_d = newpop.d(:,ranking(1:mu)) * input.weights;
        yParent  = yParent + sigma * parent_d; 
        tolX     = max(abs(sigma * parent_d));
        
                       
        %% Update evolution path and transformation matrix
        ps = (1-input.cs) * ps + sqrt_s * parent_z; 
        hsig = sum(ps.^2)/(1-(1-input.cs)^(2*evals.fun/lambda))/dim < 2+4/(dim+1);
        hsig2 = piM*oMM;
        pc = (1-input.cc) * hsig2 * pc + hsig * sqrt_c * parent_z;
        CC = (1-input.c1-input.cmu)*eye(dim)+input.c1*(pc*pc'+(1-hsig)*input.cc*(2-input.cc)*eye(dim))+input.cmu...
            *newpop.z(:,ranking(1:mu))*diag(input.weights)*(newpop.z(:,ranking(1:mu)))';


        oMM = MM;
        MM = 0.5*(oMM+oMM*CC);


        liMM = MM > 1e+12;
        siMM = MM < -1e+12;
        if sum(sum(liMM))>1 || sum(sum(siMM))>1
            MM = eye(input.dim);
            oMM = eye(input.dim);
        end
        %% Adapt the mutation strength 
        sigma = min(sigma  * exp((input.cs/2)*(norm(ps)^2/input.dim - 1)),sigmaMAX); 
                   
       
        %% update best solution found so far                    
        if (best_conv==0 && global_best.conv==0 && best_val <= global_best.val) ||...
                (best_conv==global_best.conv && best_val <= global_best.val) || best_conv<global_best.conv
            global_best.y   = best_y; 				
            global_best.val = best_val;
            global_best.conv = best_conv;
            global_best.fevs = evals.fun;
        end

        if evals.fun>=input.budget                      % check termination criterion
            termination = 1;
        end
        %% update local best solution found so far
        if (best_conv==0 && local_best.conv==0 && best_val <= local_best.val) ||...
                (best_conv==local_best.conv && best_val <= local_best.val) || best_conv<local_best.conv
            local_best.y   = best_y; 
            local_best.val = best_val;
            local_best.conv = best_conv;
            local_best.fevs = evals.fun;
            tolY            = 0;
        else
            tolY            = tolY+1;
        end
        
        
        if (tolX < 1e-10 || tolY > 300 || tolF < 1e-25) && g > TC
            restart = 1;
            local_best.y   = []; 				
            local_best.val = inf;
            local_best.conv = inf;
        end
        %% Update epsilon value 
        g=g+1;  
        if(g>1 && g<TC)
          Eg=Eg0.*((1-g./TC).^CPg);
          Eh=Eh0.*(1-g./TC).^CPh;
        elseif(g+1>=TC)
          Eg = zeros(size(Eg));
          Eh = zeros(size(Eh));
        end 
        
        %% log global best after having used 10%, and 50% of the evaluation budget
        if evals.fun>=input.budget*10/100 && flag1==0
            fit10=global_best.val;
            con10=global_best.conv;
            [ff,gg,hh]=feval(problem.constr_fun_name,global_best.y',CEC_fun_no);
            c10_1    = sum(gg>1)                  + sum(abs(hh)>1);
            c10_2    = sum((gg>0.01) & (gg<1))    + sum(abs(hh)>0.01 & abs(hh)<1);
            c10_3    = sum((gg>0.0001)&(gg<0.01)) + sum(abs(hh)>0.0001 & abs(hh)<0.01);  
            flag1=1;
        elseif evals.fun>=input.budget*50/100 && flag2==0
            fit50=global_best.val;
            con50=global_best.conv;
            [ff,gg,hh]=feval(problem.constr_fun_name,global_best.y',CEC_fun_no);
            c50_1    = sum(gg>1)                  + sum(abs(hh)>1);
            c50_2    = sum((gg>0.01)&(gg<1))      + sum(abs(hh)>0.01 & abs(hh)<1);
            c50_3    = sum((gg>0.0001)&(gg<0.01)) + sum(abs(hh)>0.0001 & abs(hh)<0.01)  ;
            flag2=1;
        end

    end
    
    %% log final global best solution
    fit100=global_best.val;
    con100=global_best.conv;
    [ff,gg,hh]=feval(problem.constr_fun_name,global_best.y',CEC_fun_no);
    c100_1    = sum(gg>1)                   + sum(abs(hh)>1);
    c100_2    = sum((gg>0.01)&(gg<1))       + sum(abs(hh)>0.01 & abs(hh)<1);
    c100_3    = sum((gg>0.0001)&(gg<0.01))  + sum(abs(hh)>0.0001 &abs(hh)<0.01);

    out = [fit10 con10 c10_1 c10_2 c10_3;
             fit50 con50 c50_1 c50_2 c50_3;
             fit100 con100 c100_1 c100_2 c100_3];
    
end
function [sigma, g, ps, pc, MM, oMM, FEs, Eg0, Eh0, yParent, TC, Eg, Eh, CPg, CPh, input ] = restarts(problem,input,CEC_fun_no)
    dim         = length(problem.lower_bounds);
    sigma       = input.sigma;
    lambda      = 2*input.lambda;
    input.lambda = lambda;
    mu          = floor(lambda/3);
    input.mu    = mu;
   
    input.weights   = log(mu+1/2)-log(1:mu)';     % muXone array for weighted recombination
    input.weights   = input.weights./sum(input.weights);      % normalize recombination weights array
    input.mueff     = 1/sum(input.weights.^2);                    % variance-effectiveness of sum w_i x_i
    input.cc        = (4+input.mueff/dim)/(dim+4+2*input.mueff/dim);
    input.cs        = (input.mueff+2) / (dim+input.mueff+5);         % t-const for cumulation for sigma control
    input.c1        = 2 / ((dim+1.3)^2+input.mueff);                 % learning rate for rank-one update of M
    input.cmu       = min(1-input.c1, 2 * (input.mueff-2+1/input.mueff) / ((dim+2)^2+input.mueff));     % and for rank-mu update of M
    input.damps     = 1 + 2*max(0, sqrt((input.mueff-1)/(dim+1))-1) + input.cs;
            
    newpop.y    = zeros(dim,lambda);    % initialize new population matrix (n times NP)
	newpop.f    = 0;
	newpop.conv = 0;
	FEs         = 0;

    g           = 0;

    
    %% Initialize dynamic (internal) strategy parameters and constants
    ps      = zeros(dim,1);                                 % evolution paths for sigma
    pc      = zeros(dim,1);
    MM      = eye(dim);
    oMM     = eye(dim);   
    constraint_number = problem.gn(CEC_fun_no) + problem.hn(CEC_fun_no);
    %% Initialize random population of lambda candidate solutions
    for k=1:lambda
        newpop.y(:,k)   = problem.lower_bounds...
                            +(problem.upper_bounds-problem.lower_bounds).*rand(dim,1); 
        [fval, gv, hv]  = feval(problem.constr_fun_name,newpop.y(:,k)',CEC_fun_no);
        newpop.f(k)     = fval;                             % fitness vector 
        newpop.conv(k)  = sum([sum(gv.*(gv>0)), sum(abs(hv).*(abs(hv)>input.delta))])./constraint_number;  % mean constraint violations)
		newpop.g(:,k) = gv;
        newpop.h(:,k) = hv;
        FEs           = FEs + 1;                    % count objective function evaluations
    end

    %% Initial parameter for epsilon level ordering
    [TC,Eg,Eh,CPg,CPh] = ETAintialization(newpop);
    Eg0 = Eg;
    Eh0 = Eh;

    [ranking]           = eta_sort_ref(newpop,Eg,Eh);
    ParentPop           = newpop.y(:,ranking(1:mu));
    yParent             = sum(ParentPop,2)./mu;
end