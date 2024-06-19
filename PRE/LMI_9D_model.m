function controller = LMI_9D_model(flag,param,ind_Re)

%%obtain the region of attraction (ROA) from the linear matrix inequality
%%formulation.

%%add the flag of setting lambda and kappa to zero to see what is happening
%%if we do not have enough null space...

try 
    load([flag.path_data,flag.result_name,'_Re',num2str(ind_Re),'.mat'],'controller');
    disp(['Successfully load the results at the ',num2str(ind_Re),'th Reynolds number']);
catch 
    disp(['Compute the results at the ',num2str(ind_Re),'th Reynolds number']);

    %%initialize the 9D model...
    %%This is disturbance independent...
    
    %%get the formulation from RHS() function
    param.Re=param.Re_list(ind_Re); %%This is the critical Reynolds number for the 
    RHS_out=RHS(flag,param);
    A=double(RHS_out.A_linear);
    B=RHS_out.B; %Obtain A, B, C and J_norm_weight and n, null space of nonlinear term
    C=RHS_out.C;
    n=RHS_out.n;
   J_norm_weight=RHS_out.J_norm_weight;
   
    epsilon=param.epsilon; %set up the SDP option
    sdp_options = sdpsettings('verbose',flag.verbose,'cachesolvers',1,...
        'solver',flag.lmi_solver,'savesolverinput',flag.savesolverinput, 'savesolveroutput',flag.savesolveroutput);
    
    %the option for Mosek
    if strcmp(flag.lmi_solver,'mosek')
        sdp_options.sedumi.eps=flag.sedumi_epsilon;
    end

    for ind_delta=1:length(param.delta_list)
        yalmip('clear');%clear one SDP problem
        delta2=param.delta_list(ind_delta)^2; %fix the delta2 value
        
        %%input output property of lambda and kappa
            lambda_af=zeros(size(C,1),size(C,1)); %%modify the size of lambda_af
            kappa_ff=zeros(size(B,2),size(B,2)); %%modify the size of kappa_ff as the size of the input %%Update 2019/12/25
            if ~isempty(n)  %%we have a non-zero null space of nonlinear term.
                lambda=sdpvar(length(A),1);
                kappa=sdpvar(length(A),1);
                n=double(n);
                %implementing the inequality constraint like the Lagrange
                %multiplier
                for ind_e=1:length(A)
                      e=zeros(length(A),1);
                      e(ind_e)=1;
                      lambda_af=lambda_af+lambda(ind_e,1)*e*n';
                      kappa_ff=kappa_ff+kappa(ind_e,1)*(e*n'+n*e');
                end
            else
                lambda=sdpvar(1,1);
                lambda_af=zeros(size(A));
                kappa=sdpvar(1,1);
                kappa_ff=zeros(size(A));
            end
        
            if flag.IO_ff==0
                kappa_ff=zeros(size(B,2),size(B,2));
            end

            %%set up the input output property a^T M f=0
            lambda_energy=sdpvar(1,1);
            if flag.IO_af==1
                lambda_af=lambda_af+lambda_energy*RHS_out.af_energy;
            elseif flag.IO_af=='energy'
                %%implement then energy conserving property
                    lambaf_af=lambda_energy*RHS_out.af_energy;
            else
                 lambda_af=zeros(size(C,1),size(C,1));
                %%size(C,1) should equal to size(B,2)
            end

            
            if flag.IO_ff_bound %%set up the sdp on the bounds of ff.
                s=sdpvar(1);
                bound=delta2*double(J_norm_weight);
            else 
                s=0;
            end
                 
            
       if flag.lmi_roa_IO==1
        %This is the main input-output based LMI formulation
        %------start SDP formulation
        P=sdpvar(length(A),length(A)); %Set up the P, PSD cone
        I_A=eye(size(A));


         s=sdpvar(length(RHS_out.F_square),1);
         s_bound=zeros(size(A));
         for m_ind=1:length(RHS_out.F_square)
                %Update 2019/11/29, use the square of the F matrix, which
                %gives a sharper bound...
                %s(1)=0; s(3)=0; s(4)=0; s(5)=0;

                s_bound=s_bound+s(m_ind)*...
                    delta2*double(RHS_out.F_square{m_ind});
         end
         
         if length(s)<=length(A)
            diag_s=diag(s);
         else
            diag_s=diag(s(1:length(A)))+eye(size(A))*s(length(A)+1);
         end
        
         dV_ineq=[P*A+A'*P+s_bound+flag.strict*epsilon*I_A, P*B+C'*lambda_af;
            B'*P+lambda_af'*C,zeros(size(B,2),size(B,2))-diag_s+kappa_ff ];

        V_ineq=P-epsilon*I_A;
        F=[V_ineq>=0,...%epsilon*gamma2*I_A-P>=0,... %%This gamma2 restrict the condition number of P, but seems not necessary
            dV_ineq<=0,...
            s>=0];
        result_yalmip=optimize(F,[],sdp_options);
        
        P=value(P);
        try 
            [V,D]=eig(P);%compute the eigenvalue of P to obtain the delta_min
            V=fliplr(V);
            D=flipud(diag(D)); %%sort the eigenvalue and eigenvector from the max to min...
            delta_min=param.delta_list(ind_delta)/sqrt(D(1)/D(end));
        catch 
            V=zeros(size(P));
            D=zeros(length(P));
            delta_min=param.delta_list(ind_delta);
        end
        
        %output results...
        controller{ind_delta}.lmi_delta_max=param.delta_list(ind_delta);
        controller{ind_delta}.lmi_delta_min=delta_min;
        controller{ind_delta}.lmi_P_eig_val=D;
        controller{ind_delta}.lmi_P_eig_vec=V;
        controller{ind_delta}.lmi_P=value(P);
        controller{ind_delta}.lmi_dV_ineq_eig=eig(value(dV_ineq));
        controller{ind_delta}.lmi_time_yalmip=result_yalmip.yalmiptime;
        controller{ind_delta}.lmi_time_solver=result_yalmip.solvertime;
        controller{ind_delta}.lmi_lambda=value(lambda);
        controller{ind_delta}.lmi_kappa=value(kappa);
        controller{ind_delta}.lmi_lambda_energy=value(lambda_energy);
        controller{ind_delta}.lmi_s=value(s);
        controller{ind_delta}.lmi_solverinput=result_yalmip.solverinput;
        controller{ind_delta}.lmi_solveroutput=result_yalmip.solveroutput;
        [controller{ind_delta}.lmi_primalfeas,controller{ind_delta}.lmi_dualfeas] = ...
            check(F);
       end
       
    end
    
    %save after all computation.
    save([flag.path_data,flag.result_name,'_Re',num2str(ind_Re),'.mat'],'controller');
end

end


