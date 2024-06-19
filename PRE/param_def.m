classdef param_def
   %%This class define the parameters that we will use for further
   %%computation.
   %%
   
   
    properties
        Re=400;
        %%Reynolds number used for 9D Glaerkin model
        
        %%The input and output matrix.
        B=eye(9,9); 
        C=eye(9,9);
                
        %%These two parameters are used for the bounds bounds on
        %%perturbation.
        Re_list=logspace(1,5,40);
        delta_list=logspace(-5,-1,20);
        epsilon=0.01; %%small positive number to gurantee the strict positive definiteness.
 
    end
end