%% algorithm 5 is similar to algorithm 2, and carries over f and D. only fits D* to the entire set of b-values
% algorithm 5, segmented fit with f and D held, D* fit in the second run
% and not in log space
% but also forced D* upper bound (unlike algo5)
function Output = Algorithm7(bvalues, signal)

%{
for slice_ct=11 %1:25
        for i=1:max
            for j=1:max
                for c=1:5
                    vec(c)=log(slice(i,j,c+5,slice_ct)/slice(i,j,b_0,slice_ct)); %ln(S/S0) from b=250-900
                    if isnan(vec(c))||vec(c)==Inf||vec(c)==-Inf
                        vec(c)=0;
                    end
                end
%}
    vec = signal/signal(1); % ln(S/S0) normalize all signal
    vec=double(vec(:)); %force to be column vector
  
    bvalues = double(bvalues(:));
    N_bvalues=length(bvalues);
    b_split = 6; %the threshold of bvalues where b>bvalue(blim) is slow only.


    %generate D map
    %[mono_fitresult, gof] = fit(b_D',vec',ft_mono);
    fit_seg1 = fittype('(1-f)*exp(-x*D)','dependent',{'y'},'independent',{'x'});
    [seg1fit,~]=fit(bvalues(b_split:N_bvalues),vec(b_split:N_bvalues), fit_seg1,'startpoint',[0.005,0.2],'lower',[0,0],'upper',[0.1,1]);
    
    D_fit = seg1fit.D;
    f_fit = seg1fit.f; % as  ln(Ae^-bx) = ln(a) - bx.


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