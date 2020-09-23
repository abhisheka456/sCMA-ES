%% ============ A Reference Vector-Based Simplified Covariance ============
%% =============== Matrix Adaptation Evolution Strategy for ===============
%% ==================== Constrained Global Optimization ===================
% Should you have any queries, please contact
% Dr. Abhishek Kumar
% email-id: abhishek.kumar.eee13@iitbhu.ac.in
%%=========================================================================
%% ranking based on reference vector
function [ranking]= eta_sort_ref(pop,Eg,Eh)
f = pop.f;
g = pop.g;
h = pop.h;
G = sum(max(0,g-Eg),1);
H = sum(max(0,abs(h)-(Eh+0.0001)),1);
conv = G+H;
ranking = eps_sort(f,conv,0);
f1 = f(:);
c1 = conv(:);
f2 = f1';
c2 = c1';
ff = f1(:,ones(1,length(f)))-f2(ones(length(f),1),:);
cc = c1(:,ones(1,length(f)))-c2(ones(length(f),1),:);

tt = -ff./cc;
tt(tt<0|tt==inf|tt==-inf) = 0;
al = max(tt,[],2);
% al(find(al==inf)) = 0;
% a1(find(al==-inf)) = 0;
alpha = max(al);
[~,ranking2] = sort(f+alpha*conv);
% if alpha < 0
%     alpha = 0;
% end
tt = (-ff)./(cc-ff);
tt(cc==0) = 0;
tt(tt>1|tt<=0) = 0;
a2 = max(tt,[],2);
aa = [a2(ranking(1)), max(a2)];
nn = floor(rand*2)+1;
alpha2 = aa(nn);
alpha2(alpha2==0) = 0;
if alpha2 < 0
    alpha2 = 0;
elseif alpha2 > 1
    alpha2 = 0;
end
% [f' conv']
[~,ranking] = sort((1-alpha2)*f+(alpha2)*conv);%alpha
% [ranking', ranking2', ranking1']
% alpha2*f+(1-alpha2)*conv
% if rand < 0.5
%     ranking = ranking1;
% end

end
