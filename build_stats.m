%% ============ A Reference Vector-Based Simplified Covariance ============
%% =============== Matrix Adaptation Evolution Strategy for ===============
%% ==================== Constrained Global Optimization ===================
% Should you have any queries, please contact
% Dr. Abhishek Kumar
% email-id: abhishek.kumar.eee13@iitbhu.ac.in
%%=========================================================================
%% Main routine for saving the results
function [out] = build_stats(FitT,input)
            r10       = eps_sort(FitT(:,1),FitT(:,2),0);
            r50       = eps_sort(FitT(:,6),FitT(:,7),0);
            r100      = eps_sort(FitT(:,11),FitT(:,12),0);
            
            rFit10    = FitT(r10,1);
            rCon10    = FitT(r10,2);
            rFit50    = FitT(r50,6);
            rCon50    = FitT(r50,7);
            rFit100   = FitT(r100,11);
            rCon100   = FitT(r100,12);
            
            bestS10   = rFit10(1);
            bestS50   = rFit50(1);
            bestS100  = rFit100(1);
            
            mediS10   = rFit10((input.runs+1)/2);
            mediS50   = rFit50((input.runs+1)/2);
            mediS100  = rFit100((input.runs+1)/2);
            
            c10_1     = FitT(r10((input.runs+1)/2),3);
            c10_2     = FitT(r10((input.runs+1)/2),4);
            c10_3     = FitT(r10((input.runs+1)/2),5);
            c50_1     = FitT(r50((input.runs+1)/2),8);
            c50_2     = FitT(r50((input.runs+1)/2),9);
            c50_3     = FitT(r50((input.runs+1)/2),10);
            c100_1    = FitT(r100((input.runs+1)/2),13);
            c100_2    = FitT(r100((input.runs+1)/2),14);
            c100_3    = FitT(r100((input.runs+1)/2),15);
            
            mediC10   = rCon10((input.runs+1)/2);
            mediC50   = rCon50((input.runs+1)/2);
            mediC100  = rCon100((input.runs+1)/2);
            
            meanS10   = mean(FitT(:,1));
            stdS10    = std(FitT(:,1));
            meanS50   = mean(FitT(:,6));
            stdS50    = std(FitT(:,6));
            meanS100  = mean(FitT(:,11));
            stdS100   = std(FitT(:,11));
            
            worstS10  = rFit10(input.runs);
            worstS50  = rFit50(input.runs);
            worstS100 = rFit100(input.runs);
            
            num_fea10 = sum(rCon10==0);
            FR10      = num_fea10./input.runs;
            num_fea50 = sum(rCon50==0);
            FR50      = num_fea50./input.runs;
            num_fea100 = sum(rCon100==0);
            FR100      = num_fea100./input.runs;
            
            vio10      = mean(rCon10);
            vio50      = mean(rCon50);
            vio100     = mean(rCon100);
            
            out=[bestS10 mediS10 c10_1 c10_2 c10_3 mediC10 meanS10 stdS10 worstS10 FR10 vio10;
                 bestS50 mediS50 c50_1 c50_2 c50_3 mediC50 meanS50 stdS50 worstS50 FR50 vio50;
                 bestS100 mediS100 c100_1 c100_2 c100_3 mediC100 meanS100 stdS100 worstS100 FR100 vio100];
end