%% ============ A Reference Vector-Based Simplified Covariance ============
%% =============== Matrix Adaptation Evolution Strategy for ===============
%% ==================== Constrained Global Optimization ===================
% Should you have any queries, please contact
% Dr. Abhishek Kumar
% email-id: abhishek.kumar.eee13@iitbhu.ac.in
%%=========================================================================
function [L,U,opt_f,ineq,eq] = get_Fun_info(fun,D)
% Get the lower-, upper-bound and the optimal value for the function 'fun'

% Get the global optimal value of fun
switch(fun) % fun from 1 to 25 are CEC'05 functions
    case {1}, % g01
        LB = zeros(1,D);UB = zeros(1,D);
        LB(1:9) = 0; UB(1:9) = 1;
        LB(10:12) = 0; UB(10:12) = 100;
        LB(13)=0; UB(13)=1;
        opt_f = -15;
        ineq = 9; eq =0;
    case {2}, % g02
        LB =1e-2; UB = 10;
        opt_f = -0.80361910412559;
         ineq = 2; eq =0;
    case {3}, % g03
        LB = 0; UB =1;
        opt_f = -1.00050010001000;
        ineq = 0; eq =1;
    case {4}, % g04
       LB = zeros(1,D);UB = zeros(1,D);
       LB(1)=78;UB(1)=102;
       LB(2)=33;UB(2)=45;
       LB(3:5)=27;UB(3:5)=45;
       opt_f = -3.066553867178332e04;
        ineq = 6; eq =0;
    case {5}, % g05
        LB = zeros(1,D); UB = zeros(1,D);
        LB(1:2)=0; UB(1:2) = 1200;
        LB(3:4)=-0.55;UB(3:4)=0.55;
        opt_f = 5126.4967140071;
          ineq = 2; eq =3;
    case {6}, % g06
        LB = zeros(1,D); UB = zeros(1,D);
        LB(1) = 13; UB(1) = 100;
        LB(2) = 0; UB(2) = 100;
        opt_f = -6961.81387558015;
         ineq = 2; eq =0;
    case {7}, %g07
        LB = -10;UB=10;
        opt_f = 24.30620906818;
        ineq = 8; eq =0;
    case {8}, % g08
        LB = zeros(1,D); UB = zeros(1,D);
        LB(1:2) = 0 ; UB(1:2) = 10;
        opt_f = -0.09582504141;
         ineq = 2; eq =0;
    case {9}, % g09
        LB = -10; UB = 10;
        opt_f = 680.630057374402;
        ineq = 4; eq =0;
    case {10},%g10
         LB = zeros(1,D); UB = zeros(1,D);
         LB(1)=100;UB(1)=10000;
         LB(2:3)=1000;UB(2:3)=10000;
         LB(4:8)=10;UB(4:8)=1000;
         opt_f = 7049.24802052867;
          ineq = 6; eq =0;
    case {11}, % g11
        LB = -1; UB = 1;
        opt_f = 0.7499;
         ineq = 0; eq =1;
    case {12}, % g12
        LB =0; UB = 10;
        opt_f = -1;
          ineq = 9^3; eq =0;
    case {13} % g13
         LB = zeros(1,D); UB = zeros(1,D);
         LB(1:2)=-2.3;UB(1:2)=2.3;
         LB(3:5)=-3.2;UB(3:5)=3.2;
         opt_f = 0.053941514041898;
          ineq = 0; eq =3;
    case {14}, %g14
        LB = 1e-5;UB = 10;
        opt_f = -47.7648884594915;
          ineq = 0; eq =3;
    case {15}, % g15
        LB = 0; UB = 10;
        opt_f = 961.715022289961;
           ineq = 0; eq =2;
    case {16}, % g16
        LB = zeros(1,D); UB = zeros(1,D);
        LB(1)=704.4148;UB(1)=906.3855;
        LB(2)=68.6;UB(2)=288.88;
        LB(3)=0;UB(3)=134.75;
        LB(4)=193;UB(4)=287.0966;
        LB(5)=25;UB(5)=84.1988;
        opt_f = -1.90515525853479;
           ineq = 38; eq =0;
    case {17}, % g17
        LB=zeros(1,D);UB = zeros(1,D);
        LB(1)=0;UB(1) = 400;
        LB(2)=0;UB(2) =1000;
        LB(3)=340;UB(3)=420;
        LB(4)=340;UB(4)=420;
        LB(5)=-1000;UB(5)=1000;
        LB(6)=0;UB(6)=0.5236;
        opt_f =8853.53967480648;
%         opt_f =8853.533874806484;
        ineq = 0; eq =4;
    case{18}, % g18
           LB=zeros(1,D);UB = zeros(1,D);
           LB(1:8)=-10;UB(1:8)=10;
           LB(9)=0;UB(9)=20;
           opt_f = -0.866025403784439;
            ineq = 13; eq =0;
    case {19}, % g19
         LB=0;UB =10;
         opt_f = 32.6555929502463;
         ineq = 5;eq =0;
    case {20}, % g20
        LB = 0; UB = 10;
        opt_f = 0.159268357;
           ineq = 6; eq =14;
    case {21}, %g21
        LB=zeros(1,D);UB = zeros(1,D);
        LB(1)=100;UB(1)=200;
        LB(2)=1e-28;LB(3)=0.1;UB(2:3)=40;
        LB(4)=100;UB(4)=300;
        LB(5)=6.3;UB(5)=6.7;
        LB(6)=5.9;UB(6)=6.4;
        LB(7)=4.5;UB(7)=6.25;
        opt_f = 193.724510070035;
        ineq = 1; eq =5;
    case {22}, %g22
         LB=zeros(1,D);UB = zeros(1,D);
         LB(1)=0;UB(1)=20000;
         LB(2:4)=0;UB(2:4)=1e6;
         LB(5:7)=0;UB(5:7)=4e7;
         LB(8)=100;UB(8)=299.99;
         LB(9)=100;UB(9)=399.99;
         LB(10)=100.01;UB(10)=300;
         LB(11)=100;UB(11)=400;
         LB(12)=100;UB(12)=600;
         LB(13:15)=0;UB(13:15)=500;
         LB(16)=0.01;UB(16)=300;
         LB(17)=0.01;UB(17)=400;
         LB(18:22)=-4.7;UB(18:22)=6.25;
         opt_f = 236.430975504001; ineq = 1; eq =19;
    case {23}, %g23
          LB = zeros(1,D); UB = zeros(1,D);
          LB(1:8)=0;UB(1:2)=300;UB(3)=100;UB(4)=200;UB(5)=100;UB(6)=300;UB(7)=100;UB(8)=200;
          LB(9)=0.01;UB(9)=0.03;
          opt_f = -400.055099999999584;
          ineq = 2; eq =4;
    case {24}, % g24
         LB = zeros(1,D); UB = zeros(1,D);
         LB(1)=0;UB(1)=3;
        LB(2) = 0; UB(2) = 4;
        opt_f = -5.50801327159536;
        ineq = 2; eq =0;
    
end

% If LB and UB are not vectors make them vectors
sl = size(LB);

if (sl(1)*sl(2) == 1) % LB and UB are scalers
    L = LB*ones(1,D);
    U = UB*ones(1,D);
else
    L = LB;
    U = UB;
end
end