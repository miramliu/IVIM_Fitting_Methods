%% algorithm 5 is similar to algorithm 2, and carries over f and D. only fits D* to the entire set of b-values
%algorithm 5, segmented fit with f and D held, D* fit in the second run
%but also forced D* upper bound (unlike algo5)
function Output = Algorithm6(bvalues, signal,bval_cutoff_idx)


    if ~exist('bval_cutoff_idx','var')
         % third parameter does not exist, so default it to something
          bval_cutoff_idx = 6;%the threshold of bvalues where b>bvalue(blim) is slow only.
     end
    vec = log(signal/signal(1)); % ln(S/S0) normalize all signal
    vec=double(vec(:)); %force to be column vector
  
    bvalues = double(bvalues(:));
    N_bvalues=length(bvalues);


    %generate D map
    [mono_fitresult, ~] = fit(bvalues(bval_cutoff_idx:N_bvalues),vec(bval_cutoff_idx:N_bvalues),'poly1');
    
    D_fit = -mono_fitresult.p1;
    f_fit = 1-exp(mono_fitresult.p2); % as  ln(Ae^-bx) = ln(a) - bx.



    %generate D* and f map
    ft_bi = fittype('(1-f)*exp(-x*D)+f*exp(-x*(Dstar))','dependent',{'y'},'independent',{'x'},'problem',{'D', 'f'},'coefficients',{'Dstar'});
    fo_bi = 0.005; %startpoints, in alphabetical order, so D, Dstar, f. Check Algorithm 1 startpoint order!

    try
        [fitmod_bi,good_bi,~]=fit(bvalues,signal(:)/signal(1),ft_bi,'startpoint',fo_bi, 'upper',0.1, 'problem', {D_fit, f_fit}); %no longer in log space
        Output.D=D_fit;
        Output.Dstar=fitmod_bi.Dstar;
        Output.f = f_fit;
        Output.SSE = good_bi.sse;
        Output.rsq = good_bi.rsquare;
        Output.adj_rsq = good_bi.adjrsquare;
    catch
        Output.D=0;
        Output.Dstar=0;
        Output.f = 0;
        Output.SSE = NaN;
        Output.rsq = NaN;
        Output.adj_rsq = NaN;
    end