function [Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1)
%This function calculate the location and the fluid properties at the dowstream
%intersection of a characteristic with the axis X.
%Irrotational Flow
%Axisymmetric flow
%Assume the ambient pressure to be zero

%Geometrical Data of points

%Tolerance
Tol = 0.0001;     %Tolerance for the geometrical coordinates
Tol2 = 0.1;       %Tolerance for the velocity

%Known features of point 4
y_4 = 0;
v_4 = 0;
theta_4 = 0;

%PREDICTOR ALGORITHM

%Field features of Negative characteristic
Vg_1 = round(sqrt(u_1^2 + v_1^2),2);
tetha_1 = round(atand(v_1/u_1),3);
a_1 = round((a_0^2 - (gamma - 1)/2*Vg_1^2)^(1/2),2);
alpha_1 = round(asind(a_1/Vg_1),3);
lambda_1 = tand(tetha_1 - alpha_1);

%Coefficients in finite difference equations
Q_1 = u_1^2 - a_1^2;
R_1 = 2*u_1*v_1 - Q_1*lambda_1;
S_1 = delta*(a_1^2*v_1)/y_1;
    
%Geometrical Data of point 4th for the predictor
x_4 = -(y_1 - lambda_1*x_1 - y_4)/lambda_1;

%Coefficient T in finite difference equations
T_1 = S_1*(x_4 - x_1) + Q_1*u_1 + R_1*v_1;

%Flowfield properties of point 4th for the predictor
u_4  = (T_1 - R_1*v_4)/Q_1;
    
%Coefficients in finite difference equations
Q_2 = NaN;
R_2 = NaN;
S_2 = NaN;
T_2 = NaN;
lambda_2 = NaN;
      
%Predictor = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
    
%CORRECTOR ALGORITHM WITH THREE ITERATIONS

iter = 0;

%Calculo valores promedio
u_12 = (u_1 + u_4)/2;
v_12 = (v_1 + v_4)/2;
y_12 = (y_1 + y_4)/2;

%Corrector variables to save results and compare
u_4c1 = 1;
x_4c1 = 1;
u_4c2 = 0;
x_4c2 = 0;
    
while(abs(x_4c2 - x_4c1) > Tol && abs(u_4c2 - u_4c1) > Tol2)
    
    for i=1:2
        
        %Field features of Negative characteristic
        Vg_1 = round(sqrt(u_12^2 + v_12^2),3);
        tetha_1 = round(atand(v_12/u_12),3);
        a_1 = round((a_0^2 - (gamma - 1)/2*Vg_1^2)^(1/2),2);
        alpha_1 = round(asind(a_1/Vg_1),3);
        lambda_1 = tand(tetha_1 - alpha_1);
        
        %Coefficients in finite difference equations
        Q_1 = u_12^2 - a_1^2;
        R_1 = 2*u_12*v_12 - Q_1*lambda_1;
        S_1 = delta*(a_1^2*v_12)/y_12;
        
        %Geometrical Data of point 4th for the predictor
        x_4 = -(y_1 - lambda_1*x_1 - y_4)/lambda_1;
        
        %Coefficient T in finite difference equations
        T_1 = S_1*(x_4 - x_1) + Q_1*u_1 + R_1*v_1;
        
        %Flowfield properties of point 4th for the predictor
        u_4  = (T_1 - R_1*v_4)/Q_1;
        
        %Coefficients in finite difference equations
        Q_2 = NaN;
        R_2 = NaN;
        S_2 = NaN;
        T_2 = NaN;
        lambda_2 = NaN;
        
        if i == 1
            u_4c1 = u_4;
            x_4c1 = x_4;
            
            u_12 = (u_1 + u_4)/2;
            v_12 = (v_1 + v_4)/2;
            y_12 = (y_1 + y_4)/2;
            
            %Corrector_1 = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
            
        else
            u_4c2 = u_4;
            x_4c2 = x_4;
            
            u_12 = (u_1 + u_4)/2;
            v_12 = (v_1 + v_4)/2;
            y_12 = (y_1 + y_4)/2;
            
            %Corrector_2 = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
            
        end
    end
    
    iter = iter + 1;
    
end

%Corrector = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];

Vg_4 = round(sqrt(u_4^2 + v_4^2),1);

[Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);

Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];

%name = {'lambda 2'; 'lambda 1'; 'Q 2'; 'R 2'; 'S 2'; 'T 2'; 'Q 1'; 'R 1'; 'S 1'; 'T 1'; 'x 4'; 'y 4'; 'u 4'; 'v 4'};
%Point4_Data = table(Predictor, Corrector_1, Corrector_2, Corrector, 'RowNames', name)
end