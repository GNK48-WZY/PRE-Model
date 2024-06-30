
%y0 = [0,0,0,0,0,0,0,0,0];  % 9d model %
num_vectors = 10000; %should be large enough
Re_R = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000];  
min_distances = zeros(size(Re_R));
laminar_state=[1,0,0,0,0,0,0,0,0];
delta_list=logspace(-5,-1,128); %the last value of logspace means spacing should also be large enough
delete(gcp('nocreate'));
parpool(64); 
for i = 1:length(Re_R)
    Re = Re_R(i);
    laminar_delta_ind=zeros(size(delta_list));
    parfor delta_ind=1:length(delta_list)
        distances = zeros(1, num_vectors);
        delta=delta_list(delta_ind);
        for j = 1:num_vectors  
            perturbation=randn(1,9);
            y0=laminar_state+delta*perturbation/norm(perturbation);
            [t, a] = ode45(@(t, a) dynamics(t, a, Re), [0, 5000], y0, []);  
            % compute the distance 
            distances(delta_ind,j) = compute_distance(a);
        end  
        if all(distances(delta_ind,:)<10^(-3))
            laminar_delta_ind(delta_ind)=1; 
        end
    end
    % min distance 
    min_distances(i) = max(laminar_delta_ind.*delta_list);
end
    
save('ODE45_9D_model.mat');
loglog(Re_R, min_distances, '-o')
xlabel('Re')
ylabel('Minimum boundary distance')
title('Minimum boundary distance vs Reynolds number')
grid on
saveas(gcf,'figure.fig');
function adot = dynamics(~,a,Re) 
 Lx = 1.75*pi;
 Lz = 1.2*pi;
 alpha = (2*pi)/Lx;
 Beta = pi/2;
 Gamma = 2*pi/Lz;
 KBG = sqrt(Beta^2+Gamma^2);
 KAG = sqrt(alpha^2+Gamma^2);
 KABG = (alpha^2+Beta^2+Gamma^2)^(1/2);
 
 adot = zeros(9,1);
 adot(1) = Beta^2/Re - Beta^2 * a(1)/Re - sqrt(3/2)*Beta*Gamma*a(6)*a(8)/KABG + sqrt(3/2)*Beta*Gamma*a(2)*a(3)/KBG;
 adot(2) = -(4*Beta^2/3 + Gamma^2)*a(2)/Re + (5/3)*sqrt(2/3)*Gamma^2*a(4)*a(6)/KAG - Gamma^2*a(5)*a(7)/(sqrt(6)*KAG) ...
           - alpha*Beta*Gamma*a(5)*a(8)/(sqrt(6)*KAG*KABG) - sqrt(3/2)*Beta*Gamma*a(1)*a(3)/KBG - sqrt(3/2)*Beta*Gamma*a(3)*a(9)/KBG;
 adot(3) = -(Beta^2+Gamma^2)*a(3)/Re + 2*alpha*Beta*Gamma*(a(4)*a(7)+a(5)*a(6))/(sqrt(6)*KAG*KBG) + (Beta^2*(3*alpha^2 + Gamma^2) - 3*Gamma^2*(alpha^2 + Gamma^2))*a(4)*a(8)/(sqrt(6)*KAG*KBG*KABG);
 adot(4) = -(3*alpha^2+4*Beta^2)*a(4)/(3*Re) - alpha*a(1)*a(5)/sqrt(6) - 10*alpha^2*a(2)*a(6)/(3*sqrt(6)*KAG)  ...
           - sqrt(3/2)*alpha*Beta*Gamma*a(3)*a(7)/KAG*KBG - sqrt(3/2)*alpha^2*Beta^2*a(3)*a(8)/KAG*KBG*KABG - alpha*a(5)*a(9)/sqrt(6);
 adot(5) = -(alpha^2+Beta^2)*a(5)/Re + alpha*a(1)*a(4)/sqrt(6) + alpha^2*a(2)*a(7)/(sqrt(6)*KAG) - alpha*Beta*Gamma*a(2)*a(8)/(sqrt(6)*KAG*KABG) + alpha*a(4)*a(9)/sqrt(6) + 2*alpha*Beta*Gamma*a(3)*a(6)/(sqrt(6)*KAG*KBG);
 
 adot(6) = -(3*alpha^2+4*Beta^2+3*Gamma^2)*a(6)/(3*Re) + alpha*a(1)*a(7)/sqrt(6) + sqrt(3/2)*Beta*Gamma*a(1)*a(8)/KABG  ...
            + 10*(alpha^2-Gamma^2)*a(2)*a(4)/(KAG*3*sqrt(6)) - 2*sqrt(2/3)*a(3)*a(5)*alpha*Beta*Gamma/(KAG*KBG) + alpha*a(7)*a(9)/sqrt(6) + sqrt(3/2)*Beta*Gamma*a(8)*a(9)/KABG ;
 
 adot(7) = -(alpha^2+Beta^2+Gamma^2)*a(7)/Re - alpha*(a(1)*a(6)+ a(6)*a(9))/sqrt(6) + (Gamma^2- alpha^2)*a(2)*a(5)/(sqrt(6)*KAG) + alpha*Beta*Gamma*a(3)*a(4)/(sqrt(6)*KAG*KBG);
 
 adot(8) = -(alpha^2+Beta^2+Gamma^2)*a(8)/Re + 2*alpha*Beta*Gamma*a(2)*a(5)/(sqrt(6)*KAG*KABG) + Gamma^2*(3*alpha^2-Beta^2+3*Gamma^2)*a(3)*a(4)/(sqrt(6)*KAG*KBG*KABG);
 adot(9) = -9*Beta^2*a(9)/Re +sqrt(3/2)*Beta*Gamma*a(2)*a(3)/KBG -sqrt(3/2)*Beta*Gamma*a(6)*a(8)/KABG;
 
end
function d = compute_distance(a)
    laminar_state=[1,0,0,0,0,0,0,0,0];
    for i=1:20
        d_last(1) = norm(a(end-i+1, :)-laminar_state); 
    end
    d=mean(d_last);
end
