function  RHS_out  = RHS(flag, param)

%%This function return the right hand side of the governing equation
%%Author: Chang Liu
%%Date: 2019/4/27

%%Output include: 
%%the total right hand side f 
%%The right hand side represented in the state dependent linear system
%%The port-Hamiltonian description..
%%The output are default to shift the mean to the zero... otherwise, it has
%%a ori(gin) suffix.

%%2D,3D, 4D and 6D shear flow model are based on:
%%Lebovitz, N., & Mariotti, G. (2013). Edges in models of shear flow. Journal of Fluid Mechanics, 721, 386-402.



switch flag.model
    case {'9D_GalerkinNS_4pi_2pi','9D_GalerkinNS_1p75pi_1p2pi'}
        %%Read the 9D model coefficients from GalerkinNS_basis_9mode_Couette.mat
        addpath('./GalerkinNS/');
        try 
            %disp(flag.model);
            load(['./GalerkinNS/',flag.model,'.mat']);
            disp('Successfully load the coefficients of 9D model');
        catch
            disp('Do not find the coefficients of the 9D model');
            if strcmp(flag.model,'9D_GalerkinNS_4pi_2pi')
                Lx=4*pi; height=2; Lz=2*pi;
            elseif strcmp(flag.model,'9D_GalerkinNS_1p75pi_1p2pi')
                Lx=1.75*pi; height=2; Lz=1.2*pi;
            end
            GalerkinNS_mode_gen(Lx,height,Lz);
            load(['./GalerkinNS/',flag.model,'.mat']);
        end
        
        Re=param.Re;%%Reynolds number

        %%generate the J and R terms in port-Hamiltonian description
        [RHS_J_NL,RHS_J_convec,RHS_J_mean_shear,RHS_R_mean_shear,RHS_R_viscous] =...
            GalerkinNS_RHS_J_R(viscous_coeff,J_NL,convec_coeff,mean_shear_skew_symmetric,mean_shear_symmetric,inner_product,Re,a);
 
        %%These are original J and R matrix, just incase that we may need
        %%to use... typically, we move the equilibrium to zero
        RHS_out.J_ori=RHS_J_NL;
        RHS_out.R_ori=RHS_R_viscous;
        RHS_out.ulaw_mean=RHS_out.R_ori*param.mean-subs(RHS_out.J_ori,a,param.mean)*param.mean;

        %%For the original dynamics, we may need a force to balance this
        %%mean....
        
        RHS_out.H=1/2*transpose(a)*a;
        RHS_out.f=(RHS_J_NL+RHS_J_convec+RHS_J_mean_shear...
            -RHS_R_mean_shear-RHS_R_viscous)*a;
       
        %%output the state dependent representation.
        RHS_out.A=RHS_J_NL+RHS_J_convec+RHS_J_mean_shear...
         -RHS_R_mean_shear-RHS_R_viscous;
        RHS_out.Z=a;
        
        RHS_out.A_linear=RHS_J_convec+RHS_J_mean_shear-RHS_R_mean_shear-RHS_R_viscous;
        %%Also output the mean shear term, which is important for the shear
        %%production.
        RHS_out.J_MS=RHS_J_mean_shear;
        RHS_out.R_MS=RHS_R_mean_shear;
        
        %%Build the symbolic table of the state, input and output
        try
            param.C*RHS_out.R_ori*param.B;
            RHS_out.B=param.B;
            RHS_out.C=param.C;
        catch
            
            error('The size of A, B, and C are not consistent, set up as the identity matrices');
            RHS_out.B=eye(size(RHS_out.A_linear));
            RHS_out.C=eye(size(RHS_out.A_linear));
        end
        
        RHS_out.x=a; %%original state table
        [~,num_input]=size(param.B); %%They are states, u is input, y is output
        RHS_out.u=sym('u',[num_input,1]); %%define the decision variable
        RHS_out.y=param.C*RHS_out.x;
        
        
    RHS_out.af_energy=eye(9,9); %energy conserving of the nonlinearity
    RHS_out.n=[1;0;0;0;0;0;0;0;-1];
    %RHS_out.n=null(RHS_out.J_ori); %%The null space of J_NL^T
    
    
    %Test another norm, where we search the minimum distance
    switch flag.J_norm_weight
         
        case 'explicit'
            nonlinear=RHS_out.J_ori*RHS_out.x;
            for n_ind=1:length(nonlinear)
%                  RHS_out.F{n_ind}=zeros(size(RHS_out.A_linear));
                 for x_ind=1:length(RHS_out.x)
                    for y_ind=1:length(RHS_out.x)
                        RHS_out.F{n_ind}(x_ind,y_ind)=1/2*diff(diff(nonlinear(n_ind),RHS_out.x(x_ind)),RHS_out.x(y_ind));
                    end
                 end
                 [V,D]=eig(double(RHS_out.F{n_ind}));
                 RHS_out.F_eig_vec{n_ind}=V;
                 RHS_out.F_eig{n_ind}=D;
                 RHS_out.F_plus_minus{n_ind}=V*abs(D)*inv(V);
                 RHS_out.F_square{n_ind}=V*D^2*inv(V);
            end
            
            bound=zeros(9,9);
            for n_ind=1:9
                bound=bound+max(diag(RHS_out.F_eig{n_ind}))*RHS_out.F_plus_minus{n_ind};
            end
            
            J_norm=norm(RHS_out.J_ori,'fro')^2;
            % J_norm_fun=matlabFunction(J_norm,'Vars',{RHS.x});
            assume(RHS_out.x,'real');
            J_norm=simplify(J_norm);
            for x_ind=1:length(RHS_out.x)
                for y_ind=1:length(RHS_out.x)
                    J_norm_weight(x_ind,y_ind)=1/2*diff(diff(J_norm,RHS_out.x(x_ind)),RHS_out.x(y_ind));
                end
            end

            RHS_out.J_norm_weight=J_norm_weight;    
            RHS_out.F_square{length(nonlinear)+1}=J_norm_weight;
        otherwise
            RHS_out.J_norm_weight=zeros(size(RHS_out.A_linear));
    end
       
    case 'TTRD'
        %(TTRD) Trefethen, L. N., Trefethen, A. E., Reddy, S. C., & Driscoll, T. A. (1993). Hydrodynamic stability without eigenvalues.燬cience,�261(5121), 578-584.

        Re=param.Re;%%Reynolds number
        a=sym('a',[2,1]);        
        RHS_R_viscous=[1/Re, 0 ;
        0 , 1/Re];
        mean_shear=[0, 1;
                    0, 0];
%         B=eye(size(A));
%         C=eye(size(A));
        RHS_J_NL=norm(a)*[0, -1;
                          1, 0];
        RHS_out.af_energy=eye(2,2); %The input-output property that it conserve energy
        RHS_out.J_norm_weight=([0 -1; 1,0])'*[0 -1; 1,0];
        RHS_out.n=[]; %%The null space of J_NL^T
        
        %%Update: 2019.11.29, bound the nonlinear term one by one
        if strcmp(flag.J_norm_weight,'explicit')
            RHS_out.F_plus_minus{1}=([0, -1]'*[0, -1]);
            RHS_out.F_plus_minus{2}=([1, 0]'*[1,0]);
            for n_ind=1:2
                RHS_out.F_eig{n_ind}=eye(2,2);
            end
            RHS_out.F_square{1}=([0, -1]'*[0, -1]);
            RHS_out.F_square{2}=([1, 0]'*[1,0]);
        end
        
    case 'TTRDp'
        Re=param.Re;%%Reynolds number
        a=sym('a',[2,1]);        
        RHS_R_viscous=[1/Re, 0 ;
        0 , 1/Re];
        mean_shear=[0, 1;
                    0, 0];
        RHS_J_NL=[0, -a(1);
                  a(1), 0];
        RHS_out.af_energy=eye(2,2);
        RHS_out.J_norm_weight=[2, 0;
                                0, 0]; %%This is the correct frobenius norm bound
        RHS_out.n=[]; %%The null space of J_NL^T
        
        
        
    case 'TTRDpp'
        Re=param.Re;%%Reynolds number
        a=sym('a',[2,1]);        
        RHS_R_viscous=[1/Re, 0 ;
        0 , 1/Re];
        mean_shear=[0, 1;
                    0, 0];
        RHS_J_NL=[0, -a(2);
                  a(2), 0];
        RHS_out.af_energy=eye(2,2);
        RHS_out.J_norm_weight=[0, 0;
                                0, 2];%This is the correct Frobenius bound.
        RHS_out.n=[]; %%The null space of J_NL^T
        
        %%Update: 2019.11.29, bound the nonlinear term one by one
        if strcmp(flag.J_norm_weight,'explicit')
            nonlinear=RHS_J_NL*a;
            for n_ind=1:length(nonlinear)
%                  RHS_out.F{n_ind}=zeros(size(RHS_out.A_linear));
                 for x_ind=1:length(a)
                    for y_ind=1:length(a)
                        RHS_out.F{n_ind}(x_ind,y_ind)=1/2*diff(diff(nonlinear(n_ind),a(x_ind)),a(y_ind));
                    end
                 end
                 [V,D]=eig(double(RHS_out.F{n_ind}));
                 RHS_out.F_eig_vec{n_ind}=V;
                 RHS_out.F_eig{n_ind}=D;
                 RHS_out.F_plus_minus{n_ind}=V*abs(D)*inv(V);
                 RHS_out.F_square{n_ind}=V*D^2*inv(V);

            end
             RHS_out.F_plus_minus{3}=[0,0;0,1];
             RHS_out.F_square{3}=[0,0;0,1];
             RHS_out.F_eig{3}=eye(2,2);
        end
    case 'BDT'
        Re=param.Re;
        RHS_R_viscous=[1/Re, 0, 0;
            0, 1/Re, 0;
            0, 0, 1/Re];
        a=sym('a',[3,1]);        

        mean_shear=[0, Re^(-1/2), 0;
            0, 0, Re^(-1/2);
            0, 0, 0];
        RHS_J_NL=norm(a)*[0,-1,1;
                        1,0,1;
                        -1, -1,0];
        RHS_out.J_norm_weight=([0,-1,1;
                        1,0,1;
                        -1, -1,0])'*[0,-1,1;
                                    1,0,1;
                                    -1, -1,0];
        RHS_out.n=null(([0,-1,1;
                1,0,1;
               -1, -1,0])');
        RHS_out.af_energy=eye(3,3);

        %incorporate the energy conserving property.
        if strcmp(flag.J_norm_weight,'explicit')
            RHS_out.F_plus_minus{1}=([0,-1,1]'*[0,-1,1]);
            RHS_out.F_plus_minus{2}=([1,0,1]'*[1,0,1]);
            RHS_out.F_plus_minus{3}=([-1, -1,0]'*[-1, -1,0]);
            for n_ind=1:3
                RHS_out.F_eig{n_ind}=eye(3,3);
            end
            RHS_out.F_square{1}=([0,-1,1]'*[0,-1,1]);
            RHS_out.F_square{2}=([1,0,1]'*[1,0,1]);
            RHS_out.F_square{3}=([-1, -1,0]'*[-1, -1,0]);
        end
        RHS_out.F_square{4}=RHS_out.J_norm_weight;

    case 'W'
        Re=param.Re;
        RHS_R_viscous=[1/Re, 0, 0, 0;
                        0, 1/Re, 0, 0;
                        0, 0, 1/Re, 0;
                        0, 0, 0, 1/Re];
        mean_shear=[0, 1, 0, 0;
                    0, 0, 0, 0;
                    0, 0, 0, 0;
                    0, 0, 0, 0];
        a=sym('a',[4,1]);    
        RHS_J_NL=[0, 0, -a(3), -a(2);
                0, 0, a(3), 0 ;
                a(3), -a(3), 0, 0;
                a(2), 0, 0, 0];
        RHS_out.n=null(RHS_J_NL');
        RHS_out.af_energy=eye(4,4);
        
    case 'Wp'
        Re=param.Re;
        RHS_R_viscous=[1/Re, 0, 0;
            0, 1/Re, 0;
            0, 0, 1/Re];
        mean_shear=[0, 1, 0;
                    0, 0, 0;
                    0, 0, 0];
        a=sym('a',[3,1]);    
        RHS_J_NL=[0, 0, -a(3);
            0, 0, a(3);
            a(3), -a(3), 0];
        RHS_out.n=null(RHS_J_NL'); %%null space of the nonlinearity
        RHS_out.af_energy=eye(3,3);
        
end
switch flag.model

    case {'TTRDp','W','Wp','WKH','KLH','KLHp'}

        %%Update: 2019.11.29, bound the nonlinear term one by one
        if strcmp(flag.J_norm_weight,'explicit')
            nonlinear=RHS_J_NL*a;
            for n_ind=1:length(nonlinear)
    %                  RHS_out.F{n_ind}=zeros(size(RHS_out.A_linear));
                 for x_ind=1:length(a)
                    for y_ind=1:length(a)
                        RHS_out.F{n_ind}(x_ind,y_ind)=1/2*diff(diff(nonlinear(n_ind),a(x_ind)),a(y_ind));
                    end
                 end
                 [V,D]=eig(double(RHS_out.F{n_ind}));
                 RHS_out.F_eig_vec{n_ind}=V;
                 RHS_out.F_eig{n_ind}=D;
                 RHS_out.F_plus_minus{n_ind}=V*abs(D)*inv(V);
                 RHS_out.F_square{n_ind}=V*D^2*inv(V);
    
            end
        end

end

switch flag.model
    case {'TTRD','TTRDp','TTRDpp','BDT','W','Wp',}
        RHS_J_mean_shear=(mean_shear-transpose(mean_shear))/2;
        RHS_R_mean_shear=-(mean_shear+transpose(mean_shear))/2;
        RHS_out.f=(RHS_J_NL+RHS_J_mean_shear...
            -RHS_R_mean_shear-RHS_R_viscous)*a;
        RHS_out.J_ori=RHS_J_NL;
        %%output the state dependent representation.
        RHS_out.A=RHS_J_NL+RHS_J_mean_shear...
         -RHS_R_mean_shear-RHS_R_viscous;
        RHS_out.A_linear=RHS_J_mean_shear...
         -RHS_R_mean_shear-RHS_R_viscous;
        RHS_out.Z=a;
        RHS_out.H=1/2*transpose(a)*a;
        RHS_out.x=a;
        %%Also output the mean shear term, which is important for the shear
        %%production.
        RHS_out.J_MS=RHS_J_mean_shear;
        RHS_out.R_MS=RHS_R_mean_shear;
        RHS_out.x=a; %%original state table
        [~,num_input]=size(param.B); %%They are states, u is input, y is output
        RHS_out.u=sym('u',[num_input,1]); %%define the decision variable
        try 
            RHS_out.y=param.C*RHS_out.x;
        catch 
            RHS_out.y=RHS_out.x;
        end
        RHS_out.B=eye(size(RHS_out.A_linear));
        RHS_out.C=eye(size(RHS_out.A_linear));
     
        
end
