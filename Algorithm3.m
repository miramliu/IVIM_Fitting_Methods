    for slice_ct=11 %1:25
        for i=1:max
            for j=1:max
                for c=1:5
                    vec(c)=log(blur_slice(i,j,c+5,slice_ct)/blur_slice(i,j,b_0,slice_ct)); %ln(S/S0) from b=250-900
                    if isnan(vec(c))||vec(c)==Inf||vec(c)==-Inf
                        vec(c)=0;
                    end
                end

                %generate D map
                [mono_fitresult, gof] = fit(b_D',vec',ft_mono);
                D_mat(i,j,slice_ct) = -mono_fitresult.p1;
                f0_mat(i,j,slice_ct) = 1-exp(mono_fitresult.p2);
                for d=1:10
                    norm_vec(d)=blur_slice(i,j,d,slice_ct)/blur_slice(i,j,b_0,slice_ct); %S/S0 from b=0-900
                    if isnan(norm_vec(d))||norm_vec(d)==Inf||norm_vec(d)==-Inf
                        norm_vec(d)=0;
                    end
                end
                

    %            [maxi,loc]=max(norm_vec);
    %             if maxi==1 && loc==1
                    curr_D=D_mat(i,j,slice_ct);
                    curr_f0 = f0_mat(i,j,slice_ct);
                    %generate D* and f map
                    ft_bi = fittype('(1-f0)*exp(-x*D)+f*exp(-x*(Dstar))','dependent',{'y'},'independent',{'x'},'problem',{'D','f0'},'coefficients',{'Dstar','f'});
                    try
                        [fitmod_bi,good_bi,out_bi]=fit(b_val',norm_vec',ft_bi,fo_bi, 'problem', {curr_D curr_f0});
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
                        bifit_raw_normvec(raw_counter,d)=(1-fitmod_bi.f0)*exp(-b_val(d)*fitmod_bi.D)+fitmod_bi.f*exp(-b_val(d)*fitmod_bi.Dstar);
                    end
    %            else
    %                 Dstar_mat(i,j,slice_ct) = 0;
    %                 f_mat(i,j,slice_ct) = 0;
    %             end
            end
        end