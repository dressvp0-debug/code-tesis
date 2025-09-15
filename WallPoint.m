function [Array] = WallPoint(gamma, R, delta, p_st, T_st, a_0, x_2, y_2, u_2, v_2, p)

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
eqn2 = y == p(1,1)*x^15 + p(1,2)*x^14 + p(1,3)*x^13 + p(1,4)*x^12 + p(1,5)*x^11 + p(1,6)*x^10 + p(1,7)*x^9 + p(1,8)*x^8 + p(1,9)*x^7 + p(1,10)*x^6 + p(1,11)*x^5 + p(1,12)*x^4 + p(1,13)*x^3 + p(1,14)*x^2 + p(1,15)*x + p(1,16);
sol = solve([eqn1, eqn2], [x, y]);
x_44 = double(sol.x);
for ii = 1 : numel(x_44)
    z = x_44(ii,1);
    if imag(z) == 0
        x_4 = z;
        break
    end
end
y_44 = double(sol.y);
y_4 = y_44(ii,1);

T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;

%Flowfield properties of point 4th for the predictor
clear x y
syms x y
eqn1 = Q_2*x + R_2*y == T_2;
eqn2 = (15*p(1,1)*x_4^14 + 14*p(1,2)*x_4^13 + 13*p(1,3)*x_4^12 + 12*p(1,4)*x_4^11 + 11*p(1,5)*x_4^10 + 10*p(1,6)*x_4^9 + 9*p(1,7)*x_4^8 + 8*p(1,8)*x_4^7 + 7*p(1,9)*x_4^6 + 6*p(1,10)*x_4^5 + 5*p(1,11)*x_4^4 + 4*p(1,12)*x_4^3 + 3*p(1,13)*x_4^2 + 2*p(1,14)*x_4 + p(1,15))*x - y == 0;
sol = solve([eqn1, eqn2], [x, y]);
u_44 = double(sol.x);
for iii = 1 : numel(u_44)
    zz = u_44(iii,1);
    if imag(zz) == 0
        u_4 = zz;
        break
    end
end
v_44 = double(sol.y);
v_4 = v_44(iii,1);

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
        eqn2 = y == p(1,1)*x^15 + p(1,2)*x^14 + p(1,3)*x^13 + p(1,4)*x^12 + p(1,5)*x^11 + p(1,6)*x^10 + p(1,7)*x^9 + p(1,8)*x^8 + p(1,9)*x^7 + p(1,10)*x^6 + p(1,11)*x^5 + p(1,12)*x^4 + p(1,13)*x^3 + p(1,14)*x^2 + p(1,15)*x + p(1,16);
        sol = solve([eqn1, eqn2], [x, y]);
        x_44 = double(sol.x);
        for ii = 1 : numel(x_44)
            z = x_44(ii,1);
            if imag(z) == 0
                x_4 = z;
                break
            end
        end
        y_44 = double(sol.y);
        y_4 = y_44(ii,1);
        
        T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
        
        %Flowfield properties of point 4th for the predictor
        clear x y
        syms x y
        eqn1 = Q_2*x + R_2*y == T_2;
        eqn2 = (15*p(1,1)*x_4^14 + 14*p(1,2)*x_4^13 + 13*p(1,3)*x_4^12 + 12*p(1,4)*x_4^11 + 11*p(1,5)*x_4^10 + 10*p(1,6)*x_4^9 + 9*p(1,7)*x_4^8 + 8*p(1,8)*x_4^7 + 7*p(1,9)*x_4^6 + 6*p(1,10)*x_4^5 + 5*p(1,11)*x_4^4 + 4*p(1,12)*x_4^3 + 3*p(1,13)*x_4^2 + 2*p(1,14)*x_4 + p(1,15))*x - y == 0;
        sol = solve([eqn1, eqn2], [x, y]);
        u_44 = double(sol.x);
        u_4 = u_44(iii,1);
        v_44 = double(sol.y);
        v_4 = v_44(iii,1);
        
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