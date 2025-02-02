
function Output = Algorithm2(bvalues, signal)

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
    vec = log(signal/signal(1)); % ln(S/S0) normalize all signal
    vec=vec(:); %force to be column vector
  
    bvalues = bvalues(:);
    N_bvalues=length(bvalues);
    b_split = 6; %the threshold of bvalues where b>bvalue(blim) is slow only.


    %generate D map
    %[mono_fitresult, gof] = fit(b_D',vec',ft_mono);
    [mono_fitresult, ~] = fit(bvalues(b_split:N_bvalues),vec(b_split:N_bvalues),'poly1');
    D_fit = -mono_fitresult.p1


    %generate D* and f map
    ft_bi = fittype('(1-f)*exp(-x*D)+f*exp(-x*(Dstar))','dependent',{'y'},'independent',{'x'},'problem',{'D'},'coefficients',{'Dstar','f'});
    fo_bi = [0.005, 0.07]; %startpoints, in alphabetical order, so D, Dstar, f. Check Algorithm 1 startpoint order!

    try
        [fitmod_bi,good_bi,~]=fit(bvalues,signal(:)/signal(1),ft_bi,'startpoint',fo_bi, 'problem', D_fit); %no longer in log space
        Output.D=D_fit;
        Output.Dstar=fitmod_bi.Dstar;
        Output.f = fitmod_bi.f;
        Output.RSSE = good_bi.adjrsquare;
    catch
        Output.D=0;
        Output.Dstar=0;
        Output.f = 0;
        Output.RSSE = NaN;
    end
                    %{
%                     if good_bi.adjrsquare < 0.80
%                         Dstar_mat(i,j,slice_ct) = 0;
%                         f_mat(i,j,slice_ct) = 0;
%                     else
                        Dstar_mat(i,j,slice_ct) = fitmod_bi.Dstar;
                        f_mat(i,j,slice_ct) = fitmod_bi.f;
                        good_bi_fit(i,j,slice_ct) = good_bi.adjrsquare;
                        
%                     end
                catch
                    error = error + 1;
                end
                
                raw_counter=raw_counter+1;
                raw_normvec(raw_counter,:)=norm_vec;
                for d=1:10
                    bifit_raw_normvec(raw_counter,d)=(1-fitmod_bi.f)*exp(-b_val(d)*fitmod_bi.D)+fitmod_bi.f*exp(-b_val(d)*fitmod_bi.Dstar);
                end

%            else
%                 Dstar_mat(i,j,slice_ct) = 0;
%                 f_mat(i,j,slice_ct) = 0;
%             end
        end
        end
                %}