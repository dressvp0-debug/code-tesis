function [Array_2, Array_4] = InverseWallPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_3, y_3, u_3, v_3, x_4, y_4, theta_4)

%Employing the unit process for an inverse wall 

%Tolerance
Tol = 0.00001;     %Tolerance for the geometrical coordinates
Tol2 = 0.00001;       %Tolerance for the velocity

%PREDICTOR ALGORITHM

%FLOW PROPERTIES AT POINT 2 MAKEN PROPERTIES 2 = PROPERTIES 3

u_2 = u_3;
v_2 = v_3;

Vg_2 = (u_2^2 + v_2^2)^(1/2);                %Magnitude of the velocity
a_2 = (a_0^2 - (gamma - 1)/2*Vg_2^2)^(1/2);  %Local speed of sound
tetha_2 = atand(v_2/u_2);                    %Flow or stramline angle
alpha_2 = asind(a_2/Vg_2);                   %Mach angle
lambda_2 = tand(tetha_2 + alpha_2);          %Slope of the characteristics

%POSITION OF POINT 2
clear x y
syms x y
eqn1 = y == ((y_3 - y_1)/(x_3 - x_1))*(x - x_1) + y_1; %chord point 3 - 1
eqn2 = y - y_4 == lambda_2*(x - x_4);
sol = solve([eqn1, eqn2], [x,y]);
x_2 = double(sol.x);
y_2 = double(sol.y);

%INTERPOLATING BETWEEN POINT 3 AND 1
u_22 = u_1+((x_2 - x_1)/(x_3 - x_1))*(u_3 - u_1);
v_22 = v_1+((y_2 - y_1)/(y_3 - y_1))*(v_3 - v_1);

iter_0 = 0;

while(abs(u_2 - u_22) > Tol2 && abs(v_2 - v_22) > Tol2)
    
    iter_0 = iter_0 +1;
    
    %Repeating all the procedure again
    Vg_2 = (u_22^2 + v_22^2)^(1/2);                %Magnitude of the velocity
    a_2 = (a_0^2 - (gamma - 1)/2*Vg_2^2)^(1/2);  %Local speed of sound
    tetha_2 = atand(v_22/u_22);                    %Flow or stramline angle
    alpha_2 = asind(a_2/Vg_2);                   %Mach angle
    lambda_2 = tand(tetha_2 + alpha_2);          %Slope of the characteristics
    
    %POSITION OF POINT 2
    clear x y
    syms x y
    eqn1 = y == ((y_3 - y_1)/(x_3 - x_1))*(x - x_1) + y_1; %chord point 3 - 1
    eqn2 = y - y_4 == lambda_2*(x - x_4);
    sol = solve([eqn1, eqn2], [x,y]);
    x_2 = double(sol.x);
    y_2 = double(sol.y);
    
    u_2 = u_22;
    v_2 = v_22;
    
    u_22 = u_1+((x_2 - x_1)/(x_3 - x_1))*(u_3 - u_1);
    v_22 = v_1+((y_2 - y_1)/(y_3 - y_1))*(v_3 - v_1);
end

u_2 = u_22;
v_2 = v_22;

%COEFFICIENTS FOR THE PREDICTOR

%Flow properties
Vg_2 = (u_2^2 + v_2^2)^(1/2);                %Magnitude of the velocity
a_2 = (a_0^2 - (gamma - 1)/2*Vg_2^2)^(1/2);  %Local speed of sound
tetha_2 = atand(v_2/u_2);                    %Flow or stramline angle
alpha_2 = asind(a_2/Vg_2);                   %Mach angle
lambda_2 = tand(tetha_2 + alpha_2);          %Slope of the characteristics

%Coefficients in finite difference equations
Q_2 = u_2^2 - a_2^2;
R_2 = 2*u_2*v_2 - Q_2*lambda_2;
S_2 = delta*(a_2^2*v_2)/y_2;
T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
 
clear x y
syms x y
eqn1 = x*tand(theta_4) == y;
eqn2 = Q_2*x + R_2*y == T_2;
sol = solve([eqn1, eqn2], [x, y]);
u_4 = double(sol.x);
v_4 = double(sol.y);

%Predictor = [lambda_2; Q_2; R_2; S_2; T_2; x_2; y_2; u_2; v_2; u_4; v_4];

%Termodynamics properties at point 2

%CORRECTOR ALGORITHM WITH THREE ITERATIONS

iter = 0;

%Calculo valores promedio
u_22 = (u_2 + u_4)/2;
v_22 = (v_2 + v_4)/2;
y_22 = (y_2 + y_4)/2;

%Corrector values
u_4c1 = 1;
v_4c1 = 1;
y_2c1 = 1;
x_2c1 = 1;
u_4c2 = 0;
v_4c2 = 0;
y_2c2 = 0;
x_2c2 = 0;
    
while(abs(x_2c2 - x_2c1) > Tol && abs(y_2c2 - y_2c1) > Tol && abs(u_4c2 - u_4c1) > Tol2 && abs(v_4c2 - v_4c1) > Tol2)
    
    for i=1:2
        
        %Field features of Positive characteristic
        Vg_2 = round(sqrt(u_22^2 + v_22^2),1);
        tetha_2 = round(atand(v_22/u_22),3);
        a_2 = round((a_0^2 - (gamma - 1)/2*Vg_2^2)^(1/2),2);
        alpha_2 = round(asind(a_2/Vg_2),3);
        lambda_2 = tand(tetha_2 + alpha_2);
        
        %POSITION OF POINT 2
        clear x y
        syms x y
        eqn1 = y == ((y_3 - y_1)/(x_3 - x_1))*(x - x_1) + y_1; %chord point 3 - 1
        eqn2 = y - y_4 == lambda_2*(x - x_4);
        sol = solve([eqn1, eqn2], [x,y]);
        x_2 = double(sol.x);
        y_2 = double(sol.y);
        
        %INTERPOLATING BETWEEN POINT 3 AND 1
        u_2 = u_1+((x_2 - x_1)/(x_3 - x_1))*(u_3 - u_1);
        v_2 = v_1+((y_2 - y_1)/(y_3 - y_1))*(v_3 - v_1);
        
        %Coefficients in finite difference equations
        Q_2 = u_22^2 - a_2^2;
        R_2 = 2*u_22*v_22 - Q_2*lambda_2;
        S_2 = delta*(a_2^2*v_22)/y_22;
        
        %POSITION OF POINT 2
        clear x y
        syms x y
        eqn1 = y == ((y_3 - y_1)/(x_3 - x_1))*(x - x_1) + y_1; %chord point 3 - 1
        eqn2 = y - y_4 == lambda_2*(x - x_4);
        sol = solve([eqn1, eqn2], [x,y]);
        x_2 = double(sol.x);
        y_2 = double(sol.y);
        
        %Coefficient T in finite difference equations
        T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
        
        %Flowfield properties of point 4th for the predictor
        clear x y
        syms x y
        eqn1 = x*tand(theta_4) == y;
        eqn2 = Q_2*x + R_2*y == T_2;
        sol = solve([eqn1, eqn2], [x, y]);
        u_4 = double(sol.x);
        v_4 = double(sol.y);
        
        if i == 1
            u_4c1 = u_4;
            v_4c1 = v_4;
            y_2c1 = y_2;
            x_2c1 = x_2;
            
            u_22 = (u_2 + u_4)/2;
            v_22 = (v_2 + v_4)/2;
            y_22 = (y_2 + y_4)/2;
            
            %Corrector_1 = [lambda_2; Q_2; R_2; S_2; T_2; x_2; y_2; u_2; v_2; u_4; v_4];
            
        else
            u_4c2 = u_4;
            v_4c2 = v_4;
            y_2c2 = y_2;
            x_2c2 = x_2;
            
            u_22 = (u_2 + u_4)/2;
            v_22 = (v_2 + v_4)/2;
            y_22 = (y_2 + y_4)/2;
            
            %Corrector_2 = [lambda_2; Q_2; R_2; S_2; T_2; x_2; y_2; u_2; v_2; u_4; v_4];
            
        end
    end
    
    iter = iter + 1;
    
end

%Corrector = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];

Vg_4 = Vg_2;
[Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
Array_2 = [y_2, x_2, u_2, v_2, Ma t, rho, p];

Vg_4 = round(sqrt(u_4^2 + v_4^2),1);
[Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
Array_4 = [y_4, x_4, u_4, v_4, Ma t, rho, p];

%name = {'lambda 2'; 'lambda 1'; 'Q 2'; 'R 2'; 'S 2'; 'T 2'; 'Q 1'; 'R 1'; 'S 1'; 'T 1'; 'x 4'; 'y 4'; 'u 4'; 'v 4'};
%Point4_Data = table(Predictor, Corrector_1, Corrector_2, Corrector, 'RowNames', name)

end