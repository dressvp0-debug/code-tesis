function [Array] = WallPointXX(gamma, R, delta, p_st, T_st, a_0, x_2, y_2, u_2, v_2, DataPoints, InMa, Nmf, counter)

%Tolerance
Tol = 0.0001;     %Tolerance for the geometrical coordinates
Tol2 = 0.1;       %Tolerance for the velocity

%Field features of Positive characteristic
Vg_2 = round(sqrt(u_2^2 + v_2^2),2);     %Magnitude of the velocity
tetha_2 = round(atand(v_2/u_2),3);       %Flow or stramline angle
a_2 = round((a_0^2 - (gamma - 1)/2*Vg_2^2)^(1/2),2);
alpha_2 = round(asind(a_2/Vg_2),3);      %Mach angle
lambda_2 = tand(tetha_2 + alpha_2);      %Slope of the characteristics

%Coefficients in finite difference equations
Q_2 = u_2^2 - a_2^2;
R_2 = 2*u_2*v_2 - Q_2*lambda_2;
S_2 = delta*a_2^2;
if y_2 == 0
    S_2 = S_2*v_1/y_1;
else
    S_2 = S_2*v_2/y_2;
end

%Geometrical Data of point 4th for the predictor
syms x y
eqn1 = y - lambda_2*x == y_2 - lambda_2*x_2;
eqn2 = y == InMa(Nmf,1) + (DataPoints(counter,1) - InMa(Nmf,1))/(DataPoints(counter,2) - InMa(Nmf,2))*(x - InMa(Nmf,2));
sol = solve([eqn1, eqn2], [x, y]);
x_4 = double(sol.x);
y_4 = double(sol.y);


T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;

%Flowfield properties of point 4th for the predictor
clear x y
syms x y
eqn1 = Q_2*x + R_2*y == T_2;
eqn2 = y == (DataPoints(counter,1) - InMa(Nmf,1))/(DataPoints(counter,2) - InMa(Nmf,2));
sol = solve([eqn1, eqn2], [x, y]);
u_4 = double(sol.x);
v_4 = double(sol.y);

iter = 0;
    
%Calculo valores promedio
u_22 = (u_2 + u_4)/2;
v_22 = (v_2 + v_4)/2;
y_22 = (y_2 + y_4)/2;
    
%Corrector values
u_4c1 = 1;
v_4c1 = 1;
y_4c1 = 1;
x_4c1 = 1;
u_4c2 = 0;
v_4c2 = 0;
y_4c2 = 0;
x_4c2 = 0;
    
while(abs(x_4c2 - x_4c1) > Tol && abs(y_4c2 - y_4c1) > Tol && abs(u_4c2 - u_4c1) > Tol2 && abs(v_4c2 - v_4c1) > Tol2)
    
    for i=1:2
        
        %Field features of Positive characteristic
        Vg_2 = round(sqrt(u_22^2 + v_22^2),1);
        tetha_2 = round(atand(v_22/u_22),3);
        a_2 = round((a_0^2 - (gamma - 1)/2*Vg_2^2)^(1/2),2);
        alpha_2 = round(asind(a_2/Vg_2),3);
        lambda_2 = tand(tetha_2 + alpha_2);
        
        %Coefficients in finite difference equations
        Q_2 = u_22^2 - a_2^2;
        R_2 = 2*u_22*v_22 - Q_2*lambda_2;
        S_2 = delta*(a_2^2*v_22)/y_22;
        

        %Geometrical Data of point 4th for the predictor
        syms x y
        eqn1 = y - lambda_2*x == y_2 - lambda_2*x_2;
        eqn2 = y == InMa(Nmf,1) + (DataPoints(counter,1) - InMa(Nmf,1))/(DataPoints(counter,2) - InMa(Nmf,2))*(x - InMa(Nmf,2));
        sol = solve([eqn1, eqn2], [x, y]);
        x_4 = double(sol.x);
        y_4 = double(sol.y);
        
        T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
        
        %Flowfield properties of point 4th for the predictor
        clear x y
        syms x y
        eqn1 = Q_2*x + R_2*y == T_2;
        eqn2 = y == (DataPoints(counter,1) - InMa(Nmf,1))/(DataPoints(counter,2) - InMa(Nmf,2));
        sol = solve([eqn1, eqn2], [x, y]);
        u_4 = double(sol.x);
        v_4 = double(sol.y);

        
        if i == 1
            u_4c1 = u_4;
            v_4c1 = v_4;
            y_4c1 = y_4;
            x_4c1 = x_4;
            
            u_22 = (u_2 + u_4)/2;
            v_22 = (v_2 + v_4)/2;
            y_22 = (y_2 + y_4)/2;
                        
        else
            u_4c2 = u_4;
            v_4c2 = v_4;
            y_4c2 = y_4;
            x_4c2 = x_4;
            
            u_22 = (u_2 + u_4)/2;
            v_22 = (v_2 + v_4)/2;
            y_22 = (y_2 + y_4)/2;
                       
        end
    end
    
    iter = iter + 1;
    
end

Vg_4 = round(sqrt(u_4^2 + v_4^2),1);

[Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);

Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];

end