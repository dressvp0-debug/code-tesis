function [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2)
%This function calculate the location and the fluid properties at the dowstream
%intersection of the characteristics through points 1 and 2 and result in a
%interior point
%Irrotational Flow
%Axisymmetric flow
%Assume the ambient pressure to be zero

%Geometrical Data of points

%Tolerance
Tol = 0.0001;     %Tolerance for the geometrical coordinates
Tol2 = 0.1;       %Tolerance for the velocity

%PREDICTOR ALGORITHM

%Field features of Positive characteristic
Vg_2 = round(sqrt(u_2^2 + v_2^2),2);     %Magnitude of the velocity
tetha_2 = round(atand(v_2/u_2),3);       %Flow or stramline angle
a_2 = round((a_0^2 - (gamma - 1)/2*Vg_2^2)^(1/2),2);
alpha_2 = round(asind(a_2/Vg_2),3);      %Mach angle
lambda_2 = tand(tetha_2 + alpha_2);      %Slope of the characteristics

%Field features of Negative characteristic
Vg_1 = round(sqrt(u_1^2 + v_1^2),2);
tetha_1 = round(atand(v_1/u_1),3);
a_1 = round((a_0^2 - (gamma - 1)/2*Vg_1^2)^(1/2),2);
alpha_1 = round(asind(a_1/Vg_1),3);
lambda_1 = tand(tetha_1 - alpha_1);

if lambda_2 == Inf
    
    %Aproximation
    x_4 = x_2;
    u_4 = u_2;
    
    %Coefficients in finite difference equations
    Q_1 = u_1^2 - a_1^2;
    R_1 = 2*u_1*v_1 - Q_1*lambda_1;
    S_1 = delta*(a_1^2*v_1)/y_1;
    
    %Geometrical Data of point 4th for the predictor
    y_4 = y_1 + lambda_1*(x_4 - x_1);
    
    %Coefficient T in finite difference equations
    T_1 = S_1*(x_4 - x_1) + Q_1*u_1 + R_1*v_1;  
    
    %Flowfield properties of point 4th for the predictor
    v_4 = (T_1 - Q_1*u_4)/R_1;
    
    %Coefficients in finite difference equations
    Q_2 = NaN;
    R_2 = NaN; 
    S_2 = NaN;
    T_2 = NaN;
      
    %Predictor = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
    
    %CORRECTOR ALGORITHM WITH THREE ITERATIONS

    iter = 0;
    
    %Calculo valores promedio
    u_12 = (u_1 + u_4)/2;
    v_12 = (v_1 + v_4)/2;
    y_12 = (y_1 + y_4)/2;
    
    %Corrector variables to save results and compare
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
                        
            %Field features of Negative characteristic
            Vg_1 = round(sqrt(u_12^2 + v_12^2),1);
            tetha_1 = round(atand(v_12/u_12),3);
            a_1 = round((a_0^2 - (gamma - 1)/2*Vg_1^2)^(1/2),2);
            alpha_1 = round(asind(a_1/Vg_1),3);
            lambda_1 = tand(tetha_1 - alpha_1);
            
            %Coefficients in finite difference equations
            Q_1 = u_12^2 - a_1^2;
            R_1 = 2*u_12*v_12 - Q_1*lambda_1;
            S_1 = delta*(a_1^2*v_12)/y_12;
            
            %Geometrical Data of point 4th for the predictor
            y_4 = y_1 + lambda_1*(x_4 - x_1);
    
            %Coefficient T in finite difference equations
            T_1 = S_1*(x_4 - x_1) + Q_1*u_1 + R_1*v_1;
    
            %Flowfield properties of point 4th for the predictor
            v_4 = (T_1 - Q_1*u_4)/R_1;
            
            if i == 1
                u_4c1 = u_4;
                v_4c1 = v_4;
                y_4c1 = y_4;
                x_4c1 = x_4;
                
                u_12 = (u_1 + u_4)/2;
                v_12 = (v_1 + v_4)/2;
                y_12 = (y_1 + y_4)/2;
                
                %Corrector_1 = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
                
            else
                u_4c2 = u_4;
                v_4c2 = v_4;
                y_4c2 = y_4;
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

elseif lambda_1 == Inf
        
    %Aproximation
    x_4 = x_1;
    u_4 = u_1;
    
    %Coefficients in finite difference equations
    Q_2 = u_2^2 - a_2^2;
    R_2 = 2*u_2*v_2 - Q_2*lambda_2; 
    S_2 = delta*(a_2^2*v_2)/y_2;
    
    %Geometrical Data of point 4th for the predictor
    y_4 = y_2 + lambda_2*(x_4 - x_2);
    
    %Coefficient T in finite difference equations
    T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
    
    %Flowfield properties of point 4th for the predictor
    v_4 = (T_2 - Q_2*u_4)/R_2;
    
    %Coefficients in finite difference equations
    Q_1 = NaN;
    R_1 = NaN; 
    S_1 = NaN;
    T_1 = NaN;
    
    %Predictor = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
    
    %CORRECTOR ALGORITHM WITH THREE ITERATIONS

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
            y_4 = y_2 + lambda_2*(x_4 - x_2);
            
            %Coefficient T in finite difference equations
            T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
               
            %Flowfield properties of point 4th for the predictor
            v_4 = (T_2 - Q_2*u_4)/R_2;

            
            if i == 1
                u_4c1 = u_4;
                v_4c1 = v_4;
                y_4c1 = y_4;
                x_4c1 = x_4;
                
                u_22 = (u_2 + u_4)/2;
                v_22 = (v_2 + v_4)/2;
                y_22 = (y_2 + y_4)/2;
                
                %Corrector_1 = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
                
            else
                u_4c2 = u_4;
                v_4c2 = v_4;
                y_4c2 = y_4;
                x_4c2 = x_4;
                
                u_22 = (u_2 + u_4)/2;
                v_22 = (v_2 + v_4)/2;
                y_22 = (y_2 + y_4)/2;
                
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
    
else
    %Coefficients in finite difference equations
    Q_2 = u_2^2 - a_2^2;
    R_2 = 2*u_2*v_2 - Q_2*lambda_2; 
    S_2 = delta*a_2^2;
    if y_2 == 0
        S_2 = S_2*v_1/y_1;
    else
        S_2 = S_2*v_2/y_2;
    end
    
    %Coefficients in finite difference equations
    Q_1 = u_1^2 - a_1^2;
    R_1 = 2*u_1*v_1 - Q_1*lambda_1;
    S_1 = delta*(a_1^2*v_1)/y_1;
    
    %Geometrical Data of point 4th for the predictor
    syms x y
    eqn1 = y - lambda_2*x == y_2 - lambda_2*x_2;
    eqn2 = y - lambda_1*x == y_1 - lambda_1*x_1;
    sol = solve([eqn1, eqn2], [x, y]);
    x_4 = double(sol.x);
    y_4 = double(sol.y);
    
    %Coefficient T in finite difference equations
    T_1 = S_1*(x_4 - x_1) + Q_1*u_1 + R_1*v_1;
    T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
    
    %Flowfield properties of point 4th for the predictor
    clear x y
    syms x y
    eqn1 = Q_2*x + R_2*y == T_2;
    eqn2 = Q_1*x + R_1*y == T_1;
    sol = solve([eqn1, eqn2], [x, y]);
    u_4 = double(sol.x);
    v_4 = double(sol.y);
    
    %Predictor = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
    
    %CORRECTOR ALGORITHM WITH THREE ITERATIONS

    iter = 0;
    
    %Calculo valores promedio
    u_22 = (u_2 + u_4)/2;
    v_22 = (v_2 + v_4)/2;
    y_22 = (y_2 + y_4)/2;
    u_12 = (u_1 + u_4)/2;
    v_12 = (v_1 + v_4)/2;
    y_12 = (y_1 + y_4)/2;
    
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
            
            %Field features of Negative characteristic
            Vg_1 = round(sqrt(u_12^2 + v_12^2),1);
            tetha_1 = round(atand(v_12/u_12),3);
            a_1 = round((a_0^2 - (gamma - 1)/2*Vg_1^2)^(1/2),2);
            alpha_1 = round(asind(a_1/Vg_1),3);
            lambda_1 = tand(tetha_1 - alpha_1);
            
            %Coefficients in finite difference equations
            Q_1 = u_12^2 - a_1^2;
            R_1 = 2*u_12*v_12 - Q_1*lambda_1;
            S_1 = delta*(a_1^2*v_12)/y_12;
            
            %Geometrical Data of point 4th for the predictor
            syms x y
            eqn1 = y - lambda_2*x == y_2 - lambda_2*x_2;
            eqn2 = y - lambda_1*x == y_1 - lambda_1*x_1;
            sol = solve([eqn1, eqn2], [x, y]);
            x_4 = double(sol.x);
            y_4 = double(sol.y);
            
            %Coefficient T in finite difference equations
            T_1 = S_1*(x_4 - x_1) + Q_1*u_1 + R_1*v_1;
            T_2 = S_2*(x_4 - x_2) + Q_2*u_2 + R_2*v_2;
            
            %Flowfield properties of point 4th for the predictor
            clear x y
            syms x y
            eqn1 = Q_2*x + R_2*y == T_2;
            eqn2 = Q_1*x + R_1*y == T_1;
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
                u_12 = (u_1 + u_4)/2;
                v_12 = (v_1 + v_4)/2;
                y_12 = (y_1 + y_4)/2;
                
                %Corrector_1 = [lambda_2; lambda_1; Q_2; R_2; S_2; T_2; Q_1; R_1; S_1; T_1; x_4; y_4; u_4; v_4];
                
            else
                u_4c2 = u_4;
                v_4c2 = v_4;
                y_4c2 = y_4;
                x_4c2 = x_4;
                
                u_22 = (u_2 + u_4)/2;
                v_22 = (v_2 + v_4)/2;
                y_22 = (y_2 + y_4)/2;
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