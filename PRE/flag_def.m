classdef flag_def

    %%This file define the class of flag and set up the default and all
    %%possible value of the flag.
    %%Author: Chang Liu
    %%Date: 2019/04/27
    
    properties
    model={'9D_GalerkinNS_4pi_2pi','9D_GalerkinNS_1p75pi_1p2pi','Lorenz','CRTBP','Duffing',...
            'TTRD','TTRDp','TTRDpp','BDT','W','Wp','KLH','KLHp','JRB','Landau'};
    %%9D_GalerkinNS_4pi2pi: The 9D galerkin model for free shear flow,
    %%obtained for domain 4pi and 2pi
    %%9D_GalerkinNS_175pi12pi: the 9D galerkin model for free shear flow,
    %%obtained for domain 1.75 pi and 1.2pi 
    %%Lorenz: the lorenz system
    %%CRTBP: The circular restricted three body problem
    %%Duffing: The duffing system
    %%Landau equation: Stuart–Landau equation following Chapter 5.3 of
    %%Schmid Henningson (2001)
    %%Use this to study the nonlinear numerical stability 
    
    result_name='9D_GalerkinNS_4pi_2pi';
   
    %%The name that you want to save data, figure, and parameters.    
        
    %%The path to save data.
    controller={'shear','port_hamiltonian','LQR','sos','lmi','no','finished'};
    simulation='finished';
    
    path_fig='C:/Figure/feedback_control_clean/'; %%my local lattop path
    path_data='C:/Data/feedback_control_clean/';%This is the path to read results. As marcc cannot plot figures....

    %%The path to save figures.
    
    par_num=0; %%parallization core number for parfor.
    post_shear_2D={'ROA','control','gain'};
    post_shear_4D={'ROA','control','gain'};
    post_9D_GalerkinNS={'ROA'};
    print=0;
    
    %%default simulation domain of 9D model. used in the integration..
    Lx=4*pi;
    height=2;
    Lz=2*pi; 
    
    IO_ff=1; %%include the input output of ff in the lmi formulation
    IO_af=1; %%include the input output of af in the lmi formulation
    IO_ff_bound=1; %%include the bound on ff in the lmi formulation
    J_norm_weight={'sos','fro','orth'}; %%if 1, try another weight to bound the nonlinear term.
    lmi_solver={'sedumi','mosek'};
    lmi_roa_sos=0;
    lmi_roa_IO=0;
    lmi_roa_IO_fro=0;
    lmi_roa_explicit=0;
    lmi_roa_sos_mix=0;
    lmi_roa_sostools=0;
    lmi_roa_sostools_fast=0;
    lmi_roa_KSH=0;
    
    %%These are flag updated in 2020/08/19
    lmi_roa_mu=0; %%compute delta_max using mu analysis
    lmi_roa_IO_bisection=0; %%compute using lmi with bisection
    lmi_roa_sos_bisection=0; %%compute using sos with bisection
    lmi_roa_sostools_bisection=0; %%compute using sostools with bisection 
    
    post_lmi_roa_sos_com=0;
    post_lmi_roa_com_KSH=0;
    
    sedumi_epsilon=10^(-13);
    timing=0; %track the time of the eigenvalue, lmi solution...
    verbose=0; %the flag for display/printing of the solver information.s
    savesolverinput=0;
    savesolveroutput=0;
    %%The flag of post processing... isolate this flag from the model..
    %%Update 2019/10/29
    post={'lmi_roa','finished','9D_GalerkinNS_4pi_2pi','9D_GalerkinNS_1p75pi_1p2pi','Lorenz','CRTBP','Duffing'};
    strict=0; %%The strict inequality in dV/dt
    
    %%Update 2020/07/06 Here add the flag that are used for numerical
    %%stability testing
    
    debug='master'; %%flag for debug.
    
    %%Update 2020/10/05, add flag that is used for the mu_rotation branch
    %%post processing
    post_mu_rotation_angle=0;
   
    end
    
end

