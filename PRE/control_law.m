function controller = control_law(flag,param,RHS)
%%Design the controller according to user choose law..


switch flag.controller
    case 'lmi'
      %RHS from the input is not used for this case. We need to set up the RHS for each Reynolds number
      %This code currently is limited to the 9D model
      %They can be parallelized or not... I need a flag to do so...
      if flag.par_num>0 %%
          delete(gcp('nocreate'));
          parpool(flag.par_num); %%specify parfor worker num
          parfor ind_Re=1:length(param.Re_list)
            controller=lmi_roa(flag,param,ind_Re);
          end
          %%parallel computation make results not valid..
          controller=flag.par_num;
      else 
          for ind_Re=1:length(param.Re_list)
            controller=lmi_roa(flag,param,ind_Re);
          end
      end
    
end


end

