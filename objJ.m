%% ============ A Reference Vector-Based Simplified Covariance ============
%% =============== Matrix Adaptation Evolution Strategy for ===============
%% ==================== Constrained Global Optimization ===================
% Should you have any queries, please contact
% Dr. Abhishek Kumar
% email-id: abhishek.kumar.eee13@iitbhu.ac.in
%%=========================================================================
function f = objJ(x)
global func_num
        gn=[2 2 0 0 0 0 1 1 0 0 0 1 3 3 3 2 2 1];
        hn=[0 1 1 4 2 2 0 0 1 1 1 1 0 0 0 2 1 1];


[~,g,h] = CEC2010(x,func_num);
%g = max(g,0);
if gn(func_num) == 0
f = h;
elseif hn(func_num) == 0
    f = g;
else
    f = [g;h];
end
