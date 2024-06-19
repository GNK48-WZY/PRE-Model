clear all;
close all;
clc;

%%This version of code clean up many previous work and make it clear as
%%different function
%%Also use the objective oriented programming to make the structure clear
%%Author: Chang Liu
%%Date: 2019/09/30, clean up
%%update: try to using the property that the nonlinear forcing is
%%orthogonal to the velocity and vorticity to add more flexiblity for the
%%Lyapunov function...

%%include the system analysis based on dissipation inequality and Lur'e
%%decomposition.
%%include the input output analysis for the 9D model.
flag=flag_def;
flag.model='W';
flag.J_norm_weight='explicit'; 
flag.lmi_solver='mosek';
flag.controller='lmi';
flag.result_name=['ROA_',flag.model,'_bound_separate_',flag.controller];%,'_bound_sos'];

flag.simulation='finished';
flag.strict=1;
%------My desktop path
flag.path_data='C:/Data/feedback_control_clean/';
flag.path_fig='C:/Figure/feedback_control_clean/'; %%my local lattop path
flag.print=1;
flag.sedumi_epsilon=10^(-9);
flag.savesolveroutput=1;
flag.savesolverinput=1;
flag.lmi_roa_sos=0; %%Use the version of Sum of Square, most automatic, but most computational expensive
flag.lmi_roa_IO=1;  %%Use the implicit IO formulation to describe the nonlinear term.
flag.lmi_roa_KSH=0;
flag.par_num=8;
flag.IO_ff=1;
flag.IO_af=1;
flag.IO_ff_bound=1;
flag.verbose=0;
flag.post_9D_GalerkinNS='ROA';

flag.post='lmi_roa';

param=param_def;

param.Re_list=logspace(0,3.3,40);%logspace(0,3.3,40)
param.epsilon=10^(-2);
param.delta_list=logspace(-6,0,400);
%transverse over the Reynolds number and delta, perturbation distance...

%design the controller
RHS=RHS(flag,param);
if ~strcmp(flag.controller,'finished')
    controller = control_law(flag,param,RHS);
else
    disp('User has finished designing a controller. Load parameters and go to the next step of post-processing.');
end

if ~strcmp(flag.post,'finished')
    addpath('./post_processing');
    post_out=post_main(flag,param,RHS,0,0);
else
    disp('User has finished post-processing.')
end