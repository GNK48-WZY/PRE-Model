clear all;
close all;
clc;

syms a [9,1];
laminar_state=[1,0,0,0,0,0,0,0,0];
Re_R = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000];  


for Re_ind=1:length(Re_R)
     Re=Re_R(Re_ind);
     Lx = 1.75*pi;
     Lz = 1.2*pi;
     alpha = (2*pi)/Lx;
     Beta = pi/2;
     Gamma = 2*pi/Lz;
     KBG = sqrt(Beta^2+Gamma^2);
     KAG = sqrt(alpha^2+Gamma^2);
     KABG = (alpha^2+Beta^2+Gamma^2)^(1/2);
     
     adot(1,1) = Beta^2/Re - Beta^2 * a(1)/Re - sqrt(3/2)*Beta*Gamma*a(6)*a(8)/KABG + sqrt(3/2)*Beta*Gamma*a(2)*a(3)/KBG;
     adot(2,1) = -(4*Beta^2/3 + Gamma^2)*a(2)/Re + (5/3)*sqrt(2/3)*Gamma^2*a(4)*a(6)/KAG - Gamma^2*a(5)*a(7)/(sqrt(6)*KAG) ...
               - alpha*Beta*Gamma*a(5)*a(8)/(sqrt(6)*KAG*KABG) - sqrt(3/2)*Beta*Gamma*a(1)*a(3)/KBG - sqrt(3/2)*Beta*Gamma*a(3)*a(9)/KBG;
     adot(3,1) = -(Beta^2+Gamma^2)*a(3)/Re + 2*alpha*Beta*Gamma*(a(4)*a(7)+a(5)*a(6))/(sqrt(6)*KAG*KBG) + (Beta^2*(3*alpha^2 + Gamma^2) - 3*Gamma^2*(alpha^2 + Gamma^2))*a(4)*a(8)/(sqrt(6)*KAG*KBG*KABG);
     adot(4,1) = -(3*alpha^2+4*Beta^2)*a(4)/(3*Re) - alpha*a(1)*a(5)/sqrt(6) - 10*alpha^2*a(2)*a(6)/(3*sqrt(6)*KAG)  ...
               - sqrt(3/2)*alpha*Beta*Gamma*a(3)*a(7)/KAG*KBG - sqrt(3/2)*alpha^2*Beta^2*a(3)*a(8)/KAG*KBG*KABG - alpha*a(5)*a(9)/sqrt(6);
    
     adot(5,1) = -(alpha^2+Beta^2)*a(5)/Re + alpha*a(1)*a(4)/sqrt(6) + alpha^2*a(2)*a(7)/(sqrt(6)*KAG) - alpha*Beta*Gamma*a(2)*a(8)/(sqrt(6)*KAG*KABG) + alpha*a(4)*a(9)/sqrt(6) + 2*alpha*Beta*Gamma*a(3)*a(6)/(sqrt(6)*KAG*KBG);
     
     adot(6,1) = -(3*alpha^2+4*Beta^2+3*Gamma^2)*a(6)/(3*Re) + alpha*a(1)*a(7)/sqrt(6) + sqrt(3/2)*Beta*Gamma*a(1)*a(8)/KABG  ...
                + 10*(alpha^2-Gamma^2)*a(2)*a(4)/(KAG*3*sqrt(6)) - 2*sqrt(2/3)*a(3)*a(5)*alpha*Beta*Gamma/(KAG*KBG) + alpha*a(7)*a(9)/sqrt(6) + sqrt(3/2)*Beta*Gamma*a(8)*a(9)/KABG ;
     
     adot(7,1) = -(alpha^2+Beta^2+Gamma^2)*a(7)/Re - alpha*(a(1)*a(6)+ a(6)*a(9))/sqrt(6) + (Gamma^2- alpha^2)*a(2)*a(5)/(sqrt(6)*KAG) + alpha*Beta*Gamma*a(3)*a(4)/(sqrt(6)*KAG*KBG);
     
     adot(8,1) = -(alpha^2+Beta^2+Gamma^2)*a(8)/Re + 2*alpha*Beta*Gamma*a(2)*a(5)/(sqrt(6)*KAG*KABG) + Gamma^2*(3*alpha^2-Beta^2+3*Gamma^2)*a(3)*a(4)/(sqrt(6)*KAG*KBG*KABG);
    
     adot(9,1) = -9*Beta^2*a(9)/Re +sqrt(3/2)*Beta*Gamma*a(2)*a(3)/KBG -sqrt(3/2)*Beta*Gamma*a(6)*a(8)/KABG;


end
adot_jacobian=jacobian(adot);
% adot_hessian=hessian(adot,a);

