%% really messy code copy-pasted from different functions to provide to salman
% So sorry Salman - Mira Liu 1/27/2025
% runs in Matlab R2022b
function Output = Algorithm1(bvalues,Signal,bval_cutoff_idx)
%{
%% assumes given a volume of nx by ny by b-value 
for i=1:nx
    for j=1:ny
    %-----------------------------------------------------------------%
    % precondition                                                    %
    %-----------------------------------------------------------------%
    % Add signal conditioning- catch noise spikes.                    %
    %-----------------------------------------------------------------%
        if( ImageStack(1,i,j) > 100 ) 
            Signal = double(ImageStack(1:N_Bvalues,i,j)); 
            %scatter(Bvalues, Signal/Signal(1))
%}

    if ~exist('bval_cutoff_idx','var')
             % third parameter does not exist, so default it to something
        bval_cutoff_idx = 6;%the threshold of bvalues where b>bvalue(blim) is slow only.
    end
    Signal = double(Signal(:));
    bvalues = double(bvalues(:));
    [f_meas, Dstar_meas, D_meas, SSE, rsq, adj_rsq]  = IVIM_UC_LM_Fit_2step(Signal,bvalues,bval_cutoff_idx);

    Output.D=D_meas;
    Output.Dstar=Dstar_meas;
    Output.f = f_meas;
    Output.SSE = SSE;
    Output.rsq = rsq;
    Output.adj_rsq=adj_rsq;
    %{
    if( (D_meas > 0 ) && (D_meas < .01 ) )
        Output.D=D_meas;
        Output.Dstar=Dstar_meas;
        Output.f = f_meas;
        Output.RSSE = Residual;
        %{
       %fprintf(string(i)+', ' + string(j)+'\n')
       D(i,j)     = D_meas;
       Dstar(i,j) = Dstar_meas;
       f(i,j)     = f_meas;
       RSSE(i,j)  = Residual;
        %}
        
    else
        Output.D=0;
        Output.Dstar=0;
        Output.f = 0;
        Output.RSSE = Residual;
        %{
       D(i,j)     = 0.;
       Dstar(i,j) = 0.;
       f(i,j)     = 0.;
       RSSE(i,j)  = Residual;
        %}
    end
    %}
end



function [f_meas, Dstar_meas, D_meas, sse, rsq, adj_rsq] = IVIM_UC_LM_Fit_2step( IVIM, Bvalues,bval_cutoff_idx)
%-------------------------------------------------------------------------%
% Fit a single Diffusion vs B-values curve optimized for IVIM modelling   %
% %TWO STEP FIT %
%-------------------------------------------------------------------------%

f_best     = 0.05;
Dstar_best = 0.009;
D_best     = 0.0009;

N_Bvalues = length(Bvalues);

Bvalues(bval_cutoff_idx:N_Bvalues);

if( IVIM(1) == 0 )
    f_meas     = 0.;
    Dstar_meas = 0.;
    D_meas     = 0.;
    Residual   = 0.;
    IFLAG      = 0.;
    return
end

[f0, G0 ] = fit(Bvalues(bval_cutoff_idx:N_Bvalues), (IVIM(bval_cutoff_idx:N_Bvalues)./IVIM(1))  , ...   
                          '(1-f)*exp(-x*D)'                                    , ... % (1-f)e^-bD (diffusion, high b values)
                          'Startpoint', [f_best, D_best]                     , ... % (1-f) and D best fits
                          'Lower'     , [0.0 0.0]                          , ... 
                          'MaxIter'   , 100                                , ...
                          'Upper'     , [0.04 1]                         , ...
                          'TolFun'    , 10e-30                             );
   
IFLAG = 0;
   if( G0.sse > 10e-08) %goodness of fit, sse is sum of squares due to error

%       fprintf('SSE %f D input %f D_fit %f \n', G0.sse, D_best, f0.D);
       IFLAG = 1;
       [f0, G0 ] = fit(Bvalues(bval_cutoff_idx:N_Bvalues), (IVIM(bval_cutoff_idx:N_Bvalues)./IVIM(1)) , ...
                          '(1-f)*exp(-x*D)'                                    , ... %(1-f)e^-bD
                          'Startpoint', [f0.f f0.D]                    , ... % (1-f) and D
                          'Lower'     , [0.0 0.0]                          , ...
                          'MaxIter'   , 1000                                , ...
                          'Upper'     , [0.04 1]                         , ...
                          'TolFun'    , 10e-30                             );
   end
   
%----------------------------------------------------------------------%
%  Now subtract and fit (two step!).   %
%----------------------------------------------------------------------%                  
IVIM2 = IVIM./IVIM(1)-((1-f0.f)*exp(-Bvalues.*f0.D));
[f1, G1] = fit(Bvalues(1:bval_cutoff_idx-1), (IVIM2(1:bval_cutoff_idx-1))    , ...
                          'f*exp(-x*Dstar)'                                        , ...
                          'Startpoint', [f_best Dstar_best]               , ...
                          'Lower', [0  0 ]                                , ...
                          'MaxIter',100000                                , ...
                          'Upper', [.1 1]                              , ...
                          'TolFun', 10e-30                                ); ...

%-------------------------------------------------------------------------%
% Get the Residuals                                                        %
%-------------------------------------------------------------------------%
     f_meas     = f1.f;
     Dstar_meas = f1.Dstar;
     D_meas     = f0.D; 
     sse        = G1.sse;
     adj_rsq = G1.adjrsquare;
     rsq = G1.rsquare;
     
     
end
