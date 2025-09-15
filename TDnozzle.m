%%
% OmarFnts - Student of the Universidad Industrial de Santander

%DATA
clear;
clc;

%Thermodynamics Data
gamma = 1.2;                        %Specific heat ratio
R = 287.04;                         %Gas constant
T_st = 3000;                        %Stagnation Temperature
p_st = 70*10^5;                     %Stagnation Pressure
rho_st = p_st/(R*T_st);   

M_D = 2.5;                          %Mach numbrer design'
M_e = M_D;                          %Exit Mach number
theta_e = 0;                        %Exit angle
delta = 1;
y_cham = 2;                         %Combustion chamber radius - arbitrarily chosen
x_sl = 1.637;                       %Length of the subsonic nozzle, calculated in EES
                                    %to make the point of the sonic line
                                    %coincide with the Vitoshinsky formula

%Geometrical Data
y_t = 1;                            %Throat radius
rho_tup = 2.0;                      %Circular arc throat upstream
rho_td = 2.0;                       %Circular arc throat downstram
Npoints = 11;                       %Equally spaced points along the v~ = 0 line

a_0 = (gamma*R*T_st)^(1/2);         %Stagnation or total acustic speed

f1 = figure('Name','Nozzle contour design','NumberTitle','off');
f2 = figure('Name','Total expansion of the Nozzle','NumberTitle','off');
f3 = figure('Name','Analysis of the flow field - One-dimensional','NumberTitle','off');
f4 = figure('Name','Analysis of the flow field - Two-dimensional - Ideal contour','NumberTitle','off');
f5 = figure('Name','Analysis of the flow field - Two-dimensional - 16th Polynomial contour','NumberTitle','off');

%%
for q=1:2
    %%
    figure(q);
    
    if q == 2
        
        for qq=1:2
            
            if qq==1
                
                SP21 = subplot(2,1,qq);
                
            else
                
                SP22 = subplot(2,1,qq);
                
            end
            
            WallPoints = zeros(1,8);
            NWallPoints = 0;
            
            AxisPoints = zeros(1,8);
            NAxisPoints = 0;
            
            %%TRANSONIC DESIGN
            
            %Initial-value line
            [A_t, Delta_y, CriPro, epsilon, alpha, TableVariables, MassFlow_IVL, Thrust_IVL, MassFlow_1D, C_D, Thrust_1D, lambda] = SonicLine(gamma, R, delta, p_st, T_st, y_t, rho_tup, Npoints, a_0);
            
            %Boundary value - IV Line extend
            
            IVarray = table2array(TableVariables);
            
            %Internsection of arc
            syms x y
            eqn1 = x == -(((gamma + 1)*alpha)/(2*(1 + delta)))*y^2 - epsilon;
            eqn2 = x^2 + (y-3)^2 == 2^2;
            sol = solve([eqn1, eqn2], [x, y]);
            x_inter_1 = double(sol.x);
            x_inter = x_inter_1(1,1);
            y_inter_1 = double(sol.y);
            y_inter = y_inter_1(1,1);
            
            x_ss = [];
            y_ss = [];
            M_ss = [];
            t_ss = [];
            rho_ss = [];
            p_ss = [];
            
            for x = -1.637:0.01637:0
                
                y = 1/(1-(1-(1/y_cham)^2)*(1-(x+x_sl).^2/x_sl.^2).^2/(1+(x+x_sl).^2/(3.*x_sl.^2)).^3)^(1/2);
                A = pi*y^2;
                
                syms z
                eqn = A/A_t == 1/z*((2/(gamma + 1))*(1 + ((gamma - 1)/2)*z^2))^((gamma + 1)/(2*(gamma - 1)));
                sol = vpasolve(eqn, z, [0 1]);
                M = double(sol);
                
                t = T_st/(1 + (gamma - 1)/2*M^2);                 % Static temperature
                p = p_st/(1 + (gamma - 1)/2*M^2)^(gamma/(gamma - 1));
                rho = p/(R*t);
                
                x_ss = [x_ss;x];
                y_ss = [y_ss;y];
                M_ss = [M_ss;M];
                t_ss = [t_ss;t];
                rho_ss = [rho_ss;rho];
                p_ss = [p_ss;p];
                
            end
            subsonic = [y_ss,x_ss,M_ss,t_ss,rho_ss,p_ss];
            
            plot(subsonic(1:numel(subsonic(:,1)),2),subsonic(1:numel(subsonic(:,1)),1),'k','LineWidth',2)
            
            Array = [0 y_inter 0 x_inter 0 0 0 0 0 0 0 0 0];
            IVarray = [IVarray; Array];
            
            IVextend = zeros(Npoints^2,8);
            IVextend(1:Npoints, 1) = IVarray(1:Npoints, 2);    %Fill the zero matrix with the IV Line matrix for y position
            IVextend(1:Npoints, 2) = IVarray(1:Npoints, 6);    %Fill the zero matrix with the IV Line matrix for x position
            IVextend(1:Npoints, 3) = IVarray(1:Npoints, 7);    %Fill the zero matrix with the IV Line matrix for velocity in x (u)
            IVextend(1:Npoints, 4) = IVarray(1:Npoints, 8);    %Fill the zero matrix with the IV Line matrix for velocity in y (v)
            IVextend(1:Npoints, 5) = IVarray(1:Npoints, 10);   %Fill the zero matrix with the IV Line matrix for the Mach number M
            IVextend(1:Npoints, 6) = IVarray(1:Npoints, 11);   %Fill the zero matrix with the IV Line matrix for the Temperature T
            IVextend(1:Npoints, 7) = IVarray(1:Npoints, 12);   %Fill the zero matrix with the IV Line matrix for the Density (rho)
            IVextend(1:Npoints, 8) = IVarray(1:Npoints, 13);   %Fill the zero matrix with the IV Line matrix for the Presure P
            
            counter = 11;
            
            WallPoints(1,1:8) = IVextend(11,1:8);
            NWallPoints = NWallPoints + 1;
            
            AxisPoints(1,1:8) = IVextend(1,1:8);
            NAxisPoints = NAxisPoints + 1;
            
            hold on
            plot(IVarray(1:numel(IVarray(:,1)),4),IVarray(1:numel(IVarray(:,1)),2),'k:','LineWidth',2)
            plot(IVextend(1:Npoints,2),IVextend(1:Npoints,1),'k--','LineWidth',2)

            for i=1:Npoints-1
                
                J_1 = i;
                J_2 = 3*i - 1;
                
                for j=J_1:J_2
                    
                    if j==J_1
                        
                        %Ponit 1 - Negative charactristic
                        y_1 = IVarray(i+1,2);
                        x_1 = IVarray(i+1,6);
                        u_1 = IVarray(i+1,7);
                        v_1 = IVarray(i+1,8);
                        
                        %Ponit 2 - Positive characteristic
                        y_2 = IVarray(i,2);
                        x_2 = IVarray(i,6);
                        u_2 = IVarray(i,7);
                        v_2 = IVarray(i,8);
                        
                        [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                        
                        counter = counter +1;
                        IVextend(counter,1:8) = Array(1,1:8);
                                               
                        plot([IVextend(i,2),IVextend(counter,2)],[IVextend(i,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                        plot([IVextend(i+1,2),IVextend(counter,2)],[IVextend(i+1,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                        
                    elseif j==J_2
                        
                        %Ponit 1 - Negative charactristic
                        y_1 = IVextend(counter,1);
                        x_1 = IVextend(counter,2);
                        u_1 = IVextend(counter,3);
                        v_1 = IVextend(counter,4);
                        
                        [Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1);
                        
                        counter = counter +1;
                        IVextend(counter,1:8) = Array(1,1:8);
                        
                        NAxisPoints = NAxisPoints + 1;
                        AxisPoints(NAxisPoints,1:8) = Array(1,1:8);
                        
                        plot([IVextend(counter - 1,2),IVextend(counter,2)],[IVextend(counter - 1,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                        
                    else
                        
                        %Ponit 1 - Negative charactristic
                        y_1 = IVextend(counter,1);
                        x_1 = IVextend(counter,2);
                        u_1 = IVextend(counter,3);
                        v_1 = IVextend(counter,4);
                        
                        %Ponit 2 - Positive characteristic
                        y_2 = IVextend(counter - (2*i - 2),1);
                        x_2 = IVextend(counter - (2*i - 2),2);
                        u_2 = IVextend(counter - (2*i - 2),3);
                        v_2 = IVextend(counter - (2*i - 2),4);
                        
                        [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                        
                        counter = counter +1;
                        
                        IVextend(counter,1:8) = Array(1,1:8);
                        
                        plot([IVextend(counter - 1,2),IVextend(counter,2)],[IVextend(counter - 1,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                        plot([IVextend(counter - (2*i - 1),2),IVextend(counter,2)],[IVextend(counter - (2*i - 1),1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                        
                    end
                end
            end
            
            %%
            %KERNELL DESIGN
            
            Betta = 1;
            counter_2 = 0;
            DataPoints = IVextend;
            
            InMa = IVextend(Npoints,:);
            InMa = [InMa; DataPoints(counter - ((3*(Npoints - 1) - 1) - (Npoints - 1)):counter,:)];
            
            counter7 = 1;
            
            %Initial value points
            
            %Ponit 1
            y_3 = InMa(1,1);
            x_3 = InMa(1,2);
            u_3 = InMa(1,3);
            v_3 = InMa(1,4);
            
            %Ponit 3
            y_1 = InMa(2,1);
            x_1 = InMa(2,2);
            u_1 = InMa(2,3);
            v_1 = InMa(2,4);
            
            M_axis = 0;
            counter_7 = 0;
            
            while M_axis < M_D
                
                %Isosceles triangle of the arc and properties at point 4
                b = sqrt(2*rho_td^2*(1 - cosd(Betta)));
                alpha_i = (180 - Betta)/2;
                theta_i = 90 - alpha_i;
                x_4 = b*cosd(theta_i);
                y_4 = 1 + b*sind(theta_i);
                theta_4 = Betta;
                
                [Array_2, Array_4] = InverseWallPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_3, y_3, u_3, v_3, x_4, y_4, theta_4);
                
                counter = counter +1;
                
                DataPoints(counter,1:8) = Array_2(1,1:8);
                
                counter = counter +1;
                
                DataPoints(counter,1:8) = Array_4(1,1:8);
                
                NWallPoints = NWallPoints + 1;
                WallPoints(NWallPoints,1:8) = Array_4(1,1:8);
                
                ReMa = Array_4;
                
                plot([x_3,DataPoints(counter,2)],[y_3,DataPoints(counter,1)],'k','LineWidth',2)
                plot([DataPoints(counter - 1,2),DataPoints(counter,2)],[DataPoints(counter - 1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3])
                
                J_1 = 1;
                J_2 = numel(InMa(:,1));
                
                for j = J_1:J_2-1
                    
                    if Array_2(1,1) < InMa(2,1)
                        
                        %Ponit 1 - Negative charactristic
                        y_1 = ReMa(counter7,1);
                        x_1 = ReMa(counter7,2);
                        u_1 = ReMa(counter7,3);
                        v_1 = ReMa(counter7,4);
                        
                        %Ponit 2 - Positive characteristic
                        y_2 = InMa(counter7 + 2,1);
                        x_2 = InMa(counter7 + 2,2);
                        u_2 = InMa(counter7 + 2,3);
                        v_2 = InMa(counter7 + 2,4);
                        
                        [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                        
                        counter = counter + 1;
                        counter7 = counter7 + 1;
                        
                        ReMa(counter7,1:8) = Array(1,1:8);
                        DataPoints(counter,1:8) = Array(1,1:8);
                        
                        plot([ReMa(counter7,2),ReMa(counter7-1,2)],[ReMa(counter7,1),ReMa(counter7-1,1)],'color',[0.3 0.3 0.3])
                        plot([ReMa(counter7,2),InMa(counter7 + 1,2)],[ReMa(counter7,1),InMa(counter7 + 1,1)],'color',[0.3 0.3 0.3])
                        
                        if j==J_2-2
                            counter_7 = counter_7 + 1;
                            break
                        end
                        
                    else
                        
                        %Ponit 1 - Negative charactristic
                        y_1 = ReMa(counter7,1);
                        x_1 = ReMa(counter7,2);
                        u_1 = ReMa(counter7,3);
                        v_1 = ReMa(counter7,4);
                        
                        %Ponit 2 - Positive characteristic
                        y_2 = InMa(counter7 + 1,1);
                        x_2 = InMa(counter7 + 1,2);
                        u_2 = InMa(counter7 + 1,3);
                        v_2 = InMa(counter7 + 1,4);
                        
                        [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                        
                        counter = counter + 1;
                        counter7 = counter7 + 1;
                        
                        ReMa(counter7,1:8) = Array(1,1:8);
                        DataPoints(counter,1:8) = Array(1,1:8);
                        
                        plot([ReMa(counter7,2),ReMa(counter7-1,2)],[ReMa(counter7,1),ReMa(counter7-1,1)],'color',[0.3 0.3 0.3])
                        plot([ReMa(counter7,2),InMa(counter7,2)],[ReMa(counter7,1),InMa(counter7,1)],'color',[0.3 0.3 0.3])
                        
                    end
                end
                
                %Ponit 1 - Negative charactristic
                y_1 = ReMa(counter7,1);
                x_1 = ReMa(counter7,2);
                u_1 = ReMa(counter7,3);
                v_1 = ReMa(counter7,4);
                
                [Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1);
                
                counter = counter + 1;
                counter7 = counter7 + 1;
                
                ReMa(counter7,1:8) = Array(1,1:8);
                DataPoints(counter,1:8) = Array(1,1:8);
                
                NAxisPoints = NAxisPoints + 1;
                AxisPoints(NAxisPoints,1:8) = Array(1,1:8);
                
                M_axis = Array(1,5);
                
                plot([ReMa(counter7,2),ReMa(counter7-1,2)],[ReMa(counter7,1),ReMa(counter7-1,1)],'color',[0.3 0.3 0.3])
                
                Betta = Betta + 1;
                
                InMa = ReMa;
                counter7 = 1;
                
                %Ponit 1
                y_3 = InMa(1,1);
                x_3 = InMa(1,2);
                u_3 = InMa(1,3);
                v_3 = InMa(1,4);
                
                %Ponit 3
                y_1 = InMa(2,1);
                x_1 = InMa(2,2);
                u_1 = InMa(2,3);
                v_1 = InMa(2,4);
                
                counter_2 = counter_2 + 1;
                
            end
            
            NKernell = counter;
            
            NArc = NWallPoints;
            
            hold off
            
        end
        
    else
        
        
        WallPoints = zeros(1,8);
        NWallPoints = 0;
        AxisPoints = zeros(1,8);
        NAxisPoints = 0;
        
        %%TRANSONIC DESIGN
        
        %Initial-value line
        [A_t, Delta_y, CriPro, epsilon, alpha, TableVariables, MassFlow_IVL, Thrust_IVL, MassFlow_1D, C_D, Thrust_1D, lambda] = SonicLine(gamma, R, delta, p_st, T_st, y_t, rho_tup, Npoints, a_0);
        
        %Boundary value - IV Line extend
        
        IVarray = table2array(TableVariables);
        
        %Internsection of arc
        syms x y
        eqn1 = x == -(((gamma + 1)*alpha)/(2*(1 + delta)))*y^2 - epsilon;
        eqn2 = x^2 + (y-3)^2 == 2^2;
        sol = solve([eqn1, eqn2], [x, y]);
        x_inter_1 = double(sol.x);
        x_inter = x_inter_1(1,1);
        y_inter_1 = double(sol.y);
        y_inter = y_inter_1(1,1);
        
        syms x
        eqn1 = y_inter == 1/(1-(1-(1/y_cham)^2)*(1-(x_inter+x)^2/x^2)^2/(1+(x_inter+x)^2/(3*x^2))^3)^(1/2);
        sol = vpasolve(eqn1, x);
        
        x_ss = [];
        y_ss = [];
        M_ss = [];
        t_ss = [];
        rho_ss = [];
        p_ss = [];
        
        for x = -1.637:0.01637:0
            
            y = 1/(1-(1-(1/y_cham)^2)*(1-(x+x_sl).^2/x_sl.^2).^2/(1+(x+x_sl).^2/(3.*x_sl.^2)).^3)^(1/2);
            A = pi*y^2;
            
            syms z
            eqn = A/A_t == 1/z*((2/(gamma + 1))*(1 + ((gamma - 1)/2)*z^2))^((gamma + 1)/(2*(gamma - 1)));
            sol = vpasolve(eqn, z, [0 1]);
            M = double(sol);
            
            t = T_st/(1 + (gamma - 1)/2*M^2);                 % Static temperature
            p = p_st/(1 + (gamma - 1)/2*M^2)^(gamma/(gamma - 1));
            rho = p/(R*t);
            
            x_ss = [x_ss;x];
            y_ss = [y_ss;y];
            M_ss = [M_ss;M];
            t_ss = [t_ss;t];
            rho_ss = [rho_ss;rho];
            p_ss = [p_ss;p];
            
        end
        subsonic = [y_ss,x_ss,M_ss,t_ss,rho_ss,p_ss];
        
        plot(subsonic(1:numel(subsonic(:,1)),2),subsonic(1:numel(subsonic(:,1)),1),'k','LineWidth',2)
        
        Array = [0 y_inter 0 x_inter 0 0 0 0 0 0 0 0 0];
        
        IVarray = [IVarray; Array];
        
        IVextend = zeros(Npoints,8);
        IVextend(1:Npoints, 1) = IVarray(1:Npoints, 2);    %Fill the zero matrix with the IV Line matrix for y position
        IVextend(1:Npoints, 2) = IVarray(1:Npoints, 6);    %Fill the zero matrix with the IV Line matrix for x position
        IVextend(1:Npoints, 3) = IVarray(1:Npoints, 7);    %Fill the zero matrix with the IV Line matrix for velocity in x (u)
        IVextend(1:Npoints, 4) = IVarray(1:Npoints, 8);    %Fill the zero matrix with the IV Line matrix for velocity in y (v)
        IVextend(1:Npoints, 5) = IVarray(1:Npoints, 10);   %Fill the zero matrix with the IV Line matrix for the Mach number M
        IVextend(1:Npoints, 6) = IVarray(1:Npoints, 11);   %Fill the zero matrix with the IV Line matrix for the Temperature T
        IVextend(1:Npoints, 7) = IVarray(1:Npoints, 12);   %Fill the zero matrix with the IV Line matrix for the Density (rho)
        IVextend(1:Npoints, 8) = IVarray(1:Npoints, 13);   %Fill the zero matrix with the IV Line matrix for the Presure P
        
        counter = 11;
        
        WallPoints(1,1:8) = IVextend(11,1:8);
        NWallPoints = NWallPoints + 1;
        
        AxisPoints(1,1:8) = IVextend(1,1:8);
        NAxisPoints = NAxisPoints + 1;
        
        hold on
        plot(IVextend(1:Npoints,2),IVextend(1:Npoints,1),'k--','LineWidth',2)
        plot(IVarray(1:numel(IVarray(:,1)),4),IVarray(1:numel(IVarray(:,1)),2),'k:','LineWidth',2)
    
        for i=1:Npoints-1
            
            J_1 = i;
            J_2 = 3*i - 1;
            
            for j=J_1:J_2
                
                if j==J_1
                    
                    %Ponit 1 - Negative charactristic
                    y_1 = IVarray(i+1,2);
                    x_1 = IVarray(i+1,6);
                    u_1 = IVarray(i+1,7);
                    v_1 = IVarray(i+1,8);
                    
                    %Ponit 2 - Positive characteristic
                    y_2 = IVarray(i,2);
                    x_2 = IVarray(i,6);
                    u_2 = IVarray(i,7);
                    v_2 = IVarray(i,8);
                    
                    [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                    
                    counter = counter +1;
                    
                    IVextend(counter,1:8) = Array(1,1:8);
                    
                    plot([IVextend(i,2),IVextend(counter,2)],[IVextend(i,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                    plot([IVextend(i+1,2),IVextend(counter,2)],[IVextend(i+1,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                    
                elseif j==J_2
                    
                    %Ponit 1 - Negative charactristic
                    y_1 = IVextend(counter,1);
                    x_1 = IVextend(counter,2);
                    u_1 = IVextend(counter,3);
                    v_1 = IVextend(counter,4);
                    
                    [Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1);
                    
                    counter = counter +1;
                    IVextend(counter,1:8) = Array(1,1:8);
                    
                    NAxisPoints = NAxisPoints + 1;
                    AxisPoints(NAxisPoints,1:8) = Array(1,1:8);
                                        
                    plot([IVextend(counter - 1,2),IVextend(counter,2)],[IVextend(counter - 1,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                    
                else
                    
                    %Ponit 1 - Negative charactristic
                    y_1 = IVextend(counter,1);
                    x_1 = IVextend(counter,2);
                    u_1 = IVextend(counter,3);
                    v_1 = IVextend(counter,4);
                    
                    %Ponit 2 - Positive characteristic
                    y_2 = IVextend(counter - (2*i - 2),1);
                    x_2 = IVextend(counter - (2*i - 2),2);
                    u_2 = IVextend(counter - (2*i - 2),3);
                    v_2 = IVextend(counter - (2*i - 2),4);
                    
                    [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                    
                    counter = counter +1;
                    
                    IVextend(counter,1:8) = Array(1,1:8);
                    
                    plot([IVextend(counter - 1,2),IVextend(counter,2)],[IVextend(counter - 1,1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                    plot([IVextend(counter - (2*i - 1),2),IVextend(counter,2)],[IVextend(counter - (2*i - 1),1),IVextend(counter,1)],'color',[0.3 0.3 0.3])
                    
                end
            end
        end
        
        %%
        %KERNELL DESIGN
        
        Betta = 1;
        counter_2 = 0;
        DataPoints = IVextend;
        
        InMa = IVextend(Npoints,:);
        InMa = [InMa; DataPoints(counter - ((3*(Npoints - 1) - 1) - (Npoints - 1)):counter,:)];
        
        counter7 = 1;
        
        %Initial value points
        
        %Ponit 1
        y_3 = InMa(1,1);
        x_3 = InMa(1,2);
        u_3 = InMa(1,3);
        v_3 = InMa(1,4);
        
        %Ponit 3
        y_1 = InMa(2,1);
        x_1 = InMa(2,2);
        u_1 = InMa(2,3);
        v_1 = InMa(2,4);
        
        M_axis = 0;
        counter_7 = 0;
        
        while M_axis < M_D
            
            %Isosceles triangle of the arc and properties at point 4
            b = sqrt(2*rho_td^2*(1 - cosd(Betta)));
            alpha_i = (180 - Betta)/2;
            theta_i = 90 - alpha_i;
            x_4 = b*cosd(theta_i);
            y_4 = 1 + b*sind(theta_i);
            theta_4 = Betta;
            
            [Array_2, Array_4] = InverseWallPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_3, y_3, u_3, v_3, x_4, y_4, theta_4);
            
            counter = counter +1;
            
            DataPoints(counter,1:8) = Array_2(1,1:8);
            
            counter = counter +1;
            
            DataPoints(counter,1:8) = Array_4(1,1:8);
            
            NWallPoints = NWallPoints + 1;
            WallPoints(NWallPoints,1:8) = Array_4(1,1:8);
            
            ReMa = Array_4;
            
            plot([x_3,DataPoints(counter,2)],[y_3,DataPoints(counter,1)],'k','LineWidth',2)
            plot([DataPoints(counter - 1,2),DataPoints(counter,2)],[DataPoints(counter - 1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3])
            
            J_1 = 1;
            J_2 = numel(InMa(:,1));
            
            for j = J_1:J_2-1
                
                if Array_2(1,1) < InMa(2,1)
                    
                    %Ponit 1 - Negative charactristic
                    y_1 = ReMa(counter7,1);
                    x_1 = ReMa(counter7,2);
                    u_1 = ReMa(counter7,3);
                    v_1 = ReMa(counter7,4);
                    
                    %Ponit 2 - Positive characteristic
                    y_2 = InMa(counter7 + 2,1);
                    x_2 = InMa(counter7 + 2,2);
                    u_2 = InMa(counter7 + 2,3);
                    v_2 = InMa(counter7 + 2,4);
                    
                    [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                    
                    counter = counter + 1;
                    counter7 = counter7 + 1;
                    
                    ReMa(counter7,1:8) = Array(1,1:8);
                    DataPoints(counter,1:8) = Array(1,1:8);
                    
                    plot([ReMa(counter7,2),ReMa(counter7-1,2)],[ReMa(counter7,1),ReMa(counter7-1,1)],'color',[0.3 0.3 0.3])
                    plot([ReMa(counter7,2),InMa(counter7 + 1,2)],[ReMa(counter7,1),InMa(counter7 + 1,1)],'color',[0.3 0.3 0.3])
                    
                    if j==J_2-2
                        counter_7 = counter_7 + 1;
                        break
                    end
                    
                else
                    
                    %Ponit 1 - Negative charactristic
                    y_1 = ReMa(counter7,1);
                    x_1 = ReMa(counter7,2);
                    u_1 = ReMa(counter7,3);
                    v_1 = ReMa(counter7,4);
                    
                    %Ponit 2 - Positive characteristic
                    y_2 = InMa(counter7 + 1,1);
                    x_2 = InMa(counter7 + 1,2);
                    u_2 = InMa(counter7 + 1,3);
                    v_2 = InMa(counter7 + 1,4);
                    
                    [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
                    
                    counter = counter + 1;
                    counter7 = counter7 + 1;
                    
                    ReMa(counter7,1:8) = Array(1,1:8);
                    DataPoints(counter,1:8) = Array(1,1:8);
                    
                    plot([ReMa(counter7,2),ReMa(counter7-1,2)],[ReMa(counter7,1),ReMa(counter7-1,1)],'color',[0.3 0.3 0.3])
                    plot([ReMa(counter7,2),InMa(counter7,2)],[ReMa(counter7,1),InMa(counter7,1)],'color',[0.3 0.3 0.3])
                    
                end
            end
            
            %Ponit 1 - Negative charactristic
            y_1 = ReMa(counter7,1);
            x_1 = ReMa(counter7,2);
            u_1 = ReMa(counter7,3);
            v_1 = ReMa(counter7,4);
            
            [Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1);
            
            counter = counter + 1;
            counter7 = counter7 + 1;
            
            ReMa(counter7,1:8) = Array(1,1:8);
            DataPoints(counter,1:8) = Array(1,1:8);
            
            NAxisPoints = NAxisPoints + 1;
            AxisPoints(NAxisPoints,1:8) = Array(1,1:8);
            
            M_axis = Array(1,5);
            
            plot([ReMa(counter7,2),ReMa(counter7-1,2)],[ReMa(counter7,1),ReMa(counter7-1,1)],'color',[0.3 0.3 0.3])
            
            Betta = Betta + 1;
            
            InMa = ReMa;
            counter7 = 1;
            
            %Ponit 1
            y_3 = InMa(1,1);
            x_3 = InMa(1,2);
            u_3 = InMa(1,3);
            v_3 = InMa(1,4);
            
            %Ponit 3
            y_1 = InMa(2,1);
            x_1 = InMa(2,2);
            u_1 = InMa(2,3);
            v_1 = InMa(2,4);
            
            counter_2 = counter_2 + 1;
            
        end
        
        NKernell = counter;
        
        NArc = NWallPoints;
        
        hold off
        
    end
    
end

DataPoints2 = DataPoints;
counter2 = counter;
WallPoints2 = WallPoints;
NWallPoints2 = NWallPoints;
NAxisPoints2 = NAxisPoints;
AxisPoints2 = AxisPoints;

%%
for q=1:2
    %%
    %FINAL CHARACTERISTIC LINE
    
    figure(q);
    
    if q==2
        
        subplot(2,1,1);
        
    end
    
    WallPoints = WallPoints2;
    NWallPoints = NWallPoints2;
    
    counter = NKernell;
    
    DataPoints = DataPoints2;
    
    u_last = DataPoints(counter,3);
    
    [y_last] = LastPoint(u_last, M_axis, gamma, R, p_st, T_st, MassFlow_IVL);
    
    X1 = DataPoints(counter,2);
    Y1 = DataPoints(counter-1,1);
    
    [N, A_exit, Dy, FinalChar1, MassFlow_FC, Thrust_FC, counter, alpha_D, y_last2] = FinalChar(Y1, X1, gamma, R, y_last, p_st, T_st, u_last, counter, M_axis, MassFlow_IVL);
            
    DataPoints = [DataPoints; FinalChar1];
    
    Nlastchar = counter;
    
    for i = NKernell:Nlastchar-1   % 429 - counter at Kernel final point when M = M_axis > M_D
        
        hold on
        plot([DataPoints(i,2),DataPoints(i+1,2)],[DataPoints(i,1),DataPoints(i+1,1)],'color',[0.3 0.3 0.3])
        
    end
    
    %%
    % EXPANSION WAVES REGION
    
    InMa = flipud(DataPoints((NKernell - 2*counter_2 - counter_7 + 1):(NKernell - 1),1:8));
    counter4 = 0;
    ReMa = zeros(1,8);    
    counter3 = 1;
    
    mass_FC2 = [];
    
    for i = NKernell + 1:Nlastchar - 1
        
        %Ponit 1 - Negative charactristic
        y_1 = DataPoints(i,1);
        x_1 = DataPoints(i,2);
        u_1 = DataPoints(i,3);
        v_1 = DataPoints(i,4);
        
        %Ponit 2 - Positive characteristic
        y_2 = InMa(counter4 + 1,1);
        x_2 = InMa(counter4 + 1,2);
        u_2 = InMa(counter4 + 1,3);
        v_2 = InMa(counter4 + 1,4);
        
        [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
        
        counter = counter + 1;
        DataPoints(counter,1:8) = Array(1,1:8);
        
        counter4 = counter4 + 1;
        ReMa(counter4,1:8) = Array(1,1:8);
        
        pl3 = plot([InMa(counter4,2),ReMa(counter4,2)],[InMa(counter4,1),ReMa(counter4,1)],'color',[0.3 0.3 0.3]);
        pl4 = plot([DataPoints(i,2),ReMa(counter4,2)],[DataPoints(i,1),ReMa(counter4,1)],'color',[0.3 0.3 0.3]);
        
        %Mass Flow rate calculation along the Final Characteristic
        
        m_fc2 = (Dy/2)*(2*pi*DataPoints(i,7)*u_last)*(DataPoints(i,1) + DataPoints(i-1,1));
        mass_FC2(counter3) = m_fc2;
        counter3 = counter3 + 1;
        f0 = sum(mass_FC2);
        
        rho_1 = ReMa(1,7);
        y_1 = ReMa(1,1);
        V_1 = sqrt(ReMa(1,3)^2 + ReMa(1,4)^2);
        theta_1 = round(atand(ReMa(1,4)/ReMa(1,3)),3);
        
        rho_0 = DataPoints(i,7);
        y_0 = DataPoints(i,1);
        V_0 = sqrt(DataPoints(i,3)^2 + DataPoints(i,4)^2);
        theta_0 = round(atand(DataPoints(i,4)/DataPoints(i,3)),3);
        
        H = ReMa(1,1) - DataPoints(i,1);
        phi2 = atand((ReMa(1,1) - DataPoints(i,1))/(DataPoints(i,2) - ReMa(1,2)));
        phi = 180 - phi2;
        m_T = (H/2)*((2*pi*rho_0*V_0*y_0*((sind(phi - theta_0))/(sind(phi))))+(2*pi*rho_1*V_1*y_1*((sind(phi - theta_1))/(sind(phi)))));
        
        MassFlow_LC = f0 + m_T;
        
        while MassFlow_LC < MassFlow_IVL
            
            %Ponit 1 - Negative charactristic
            y_1 = ReMa(counter4,1);
            x_1 = ReMa(counter4,2);
            u_1 = ReMa(counter4,3);
            v_1 = ReMa(counter4,4);
            
            %Ponit 2 - Positive characteristic
            y_2 = InMa(counter4 + 1,1);
            x_2 = InMa(counter4 + 1,2);
            u_2 = InMa(counter4 + 1,3);
            v_2 = InMa(counter4 + 1,4);
            
            [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
            
            counter = counter +1;
            DataPoints(counter,1:8) = Array(1,1:8);
            
            counter4 = counter4 + 1;
            ReMa(counter4,1:8) = Array(1,1:8);
            
            pl1 = plot([InMa(counter4,2),ReMa(counter4,2)],[InMa(counter4,1),ReMa(counter4,1)],'color',[0.3 0.3 0.3]);
            pl2 = plot([ReMa(counter4 - 1,2),ReMa(counter4,2)],[ReMa(counter4-1,1),ReMa(counter4,1)],'color',[0.3 0.3 0.3]);
            
            counter6 = 0;
            
            for e = 1:counter4
                
                counter6 = counter6 + 1;
                
                if e == 1
                    
                    rho_1 = ReMa(e,7);
                    y_1 = ReMa(e,1);
                    V_1 = sqrt(ReMa(e,3)^2 + ReMa(e,4)^2);
                    theta_1 = round(atand(ReMa(e,4)/ReMa(e,3)),3);
                    
                    rho_0 = DataPoints(i,7);
                    y_0 = DataPoints(i,1);
                    V_0 = sqrt(DataPoints(i,3)^2 + DataPoints(i,4)^2);
                    theta_0 = round(atand(DataPoints(i,4)/DataPoints(i,3)),3);
                    
                    H = ReMa(e,1) - DataPoints(i,1);
                    phi2 = atand((ReMa(e,1) - DataPoints(i,1))/(DataPoints(i,2) - ReMa(e,2)));
                    phi = 180 - phi2;
                    m_T = (H/2)*((2*pi*rho_0*V_0*y_0*((sind(phi - theta_0))/(sind(phi))))+(2*pi*rho_1*V_1*y_1*((sind(phi - theta_1))/(sind(phi)))));
                    
                    mass_LC(counter6) = m_T;
                    
                else
                    
                    rho_1 = ReMa(e,7);
                    y_1 = ReMa(e,1);
                    V_1 = sqrt(ReMa(e,3)^2 + ReMa(e,4)^2);
                    theta_1 = round(atand(ReMa(e,4)/ReMa(e,3)),3);
                    
                    rho_0 = ReMa(e-1,7);
                    y_0 = ReMa(e-1,1);
                    V_0 = sqrt(ReMa(e-1,3)^2 + ReMa(e-1,4)^2);
                    theta_0 = round(atand(ReMa(e-1,4)/ReMa(e-1,3)),3);
                    
                    H = ReMa(e,1) - ReMa(e-1,1);
                    phi2 = atand((ReMa(e,1) - ReMa(e-1,1))/(ReMa(e-1,2) - ReMa(e,2)));
                    phi = 180 - phi2;
                    m_T = (H/2)*((2*pi*rho_0*V_0*y_0*((sind(phi - theta_0))/(sind(phi))))+(2*pi*rho_1*V_1*y_1*((sind(phi - theta_1))/(sind(phi)))));
                    
                    mass_LC(counter6) = m_T;
                    
                end
                
            end
            
            MassFlow_LC = f0 + sum(mass_LC);
            
        end
        
        if MassFlow_LC ~= MassFlow_IVL
            
            NWallPoints = NWallPoints + 1;
            num_guesses2 = 0;
            
            if i == Nlastchar - 1
                
                delete(pl3)
                delete(pl4)
                
                y_1low = DataPoints(i,1);
                x_1low = DataPoints(i,2);
                u_1low = DataPoints(i,3);
                v_1low = DataPoints(i,4);
                
                y_1high = ReMa(counter4,1);
                x_1high = ReMa(counter4,2);
                u_1high = ReMa(counter4,3);
                v_1high = ReMa(counter4,4);
                
            else
                
                delete(pl1)
                delete(pl2)
                
                y_1low = ReMa(counter4-1,1);
                x_1low = ReMa(counter4-1,2);
                u_1low = ReMa(counter4-1,3);
                v_1low = ReMa(counter4-1,4);
                
                y_1high = ReMa(counter4,1);
                x_1high = ReMa(counter4,2);
                u_1high = ReMa(counter4,3);
                v_1high = ReMa(counter4,4);
                
            end
            
            Tol3 = 0.12;
            
            Nmf = numel(InMa(:,1));
            
            while abs(MassFlow_LC - MassFlow_IVL) > Tol3
                
                num_guesses2 = num_guesses2 + 1;
                
                y_1 = (y_1low + y_1high)/2;
                x_1 = (x_1low + x_1high)/2;
                u_1 = u_1high - (y_1high - y_1)/(y_1high - y_1low)*(u_1high - u_1low);
                v_1 = v_1high - (y_1high - y_1)/(y_1high - y_1low)*(v_1high - v_1low);
                V_1 = sqrt(u_1^2 + v_1^2);
                theta_1 = atand(v_1/u_1);
                
                Vg_4 = V_1;
                
                [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                
                Array = [y_1, x_1, u_1, v_1, Ma t, rho, p];
                
                Array_x = Array;
                
                rho_1 = rho;
                
                DataPoints(counter,1:8) = Array(1,1:8);
                
                ReMa(counter4,1:8) = Array(1,1:8);
                
                WallPoints(NWallPoints,1:8) = Array(1,1:8);
                
                if i == Nlastchar - 1
                    
                    rho_0 = DataPoints(i,7);
                    y_0 = DataPoints(i,1);
                    V_0 = sqrt(DataPoints(i,3)^2 + DataPoints(i,4)^2);
                    theta_0 = atand(DataPoints(i,4)/DataPoints(i,3));
                    
                    H = ReMa(counter4,1) - DataPoints(i,1);
                    phi2 = atand((ReMa(counter4,1) - DataPoints(i,1))/(DataPoints(i,2) - ReMa(counter4,2)));
                    phi = 180 - phi2;
                    m_T = (H/2)*((2*pi*rho_0*V_0*y_0*((sind(phi - theta_0))/(sind(phi))))+(2*pi*rho_1*V_1*y_1*((sind(phi - theta_1))/(sind(phi)))));
                    
                else
                    
                    rho_0 = ReMa(counter4-1,7);
                    y_0 = ReMa(counter4-1,1);
                    V_0 = sqrt(ReMa(counter4-1,3)^2 + ReMa(counter4-1,4)^2);
                    theta_0 = atand(ReMa(counter4-1,4)/ReMa(counter4-1,3));
                    
                    H = ReMa(counter4,1) - ReMa(counter4-1,1);
                    phi2 = atand((ReMa(counter4,1) - ReMa(counter4-1,1))/(ReMa(counter4-1,2) - ReMa(counter4,2)));
                    phi = 180 - phi2;
                    m_T = (H/2)*((2*pi*rho_0*V_0*y_0*((sind(phi - theta_0))/(sind(phi))))+(2*pi*rho_1*V_1*y_1*((sind(phi - theta_1))/(sind(phi)))));
                    
                end
                
                mass_LC(counter6) = m_T;
                
                MassFlow_LC = f0 + sum(mass_LC);
                
                if MassFlow_LC - MassFlow_IVL > Tol3
                    
                    y_1high = y_1;
                    x_1high = x_1;
                    u_1high = u_1;
                    v_1high = v_1;
                    
                elseif MassFlow_LC - MassFlow_IVL < -Tol3
                    
                    y_1low = y_1;
                    x_1low = x_1;
                    u_1low = u_1;
                    v_1low = v_1;
                    
                end
            end
            
            if i==NKernell + 1
                
                x_2 = InMa(counter4,2);
                y_2 = InMa(counter4,1);
                u_2 = InMa(counter4,3);
                v_2 = InMa(counter4,4);
                                
                [Array] = WallPointXX(gamma, R, delta, p_st, T_st, a_0, x_2, y_2, u_2, v_2, DataPoints, InMa, Nmf, counter);
                
%                 DataPoints(counter,1:8) = Array(1,1:8);                
%                 WallPoints(NWallPoints,1:8) = Array(1,1:8);
%                 
%                 NWallPoints = NWallPoints + 1;
%                 counter = counter + 1;
%                 DataPoints(counter,1:8) = Array_x(1,1:8);
%                 WallPoints(NWallPoints,1:8) = Array_x(1,1:8);
%                 
%                 plot([InMa(Nmf - 1,2),DataPoints(counter - 1,2)],[InMa(Nmf - 1,1),DataPoints(counter - 1,1)],'color',[0.3 0.3 0.3]);
%                 
            end
            
        end
        
        
        if i == Nlastchar - 1
            
            plot([DataPoints(i,2),ReMa(counter4,2)],[DataPoints(i,1),ReMa(counter4,1)],'color',[0.3 0.3 0.3]);
            plot([InMa(Nmf,2),ReMa(counter4,2)],[InMa(Nmf,1),ReMa(counter4,1)],'k','LineWidth',2);
            
        else
            
            plot([ReMa(counter4 - 1,2),ReMa(counter4,2)],[ReMa(counter4-1,1),ReMa(counter4,1)],'color',[0.3 0.3 0.3]);
            plot([InMa(Nmf,2),ReMa(counter4,2)],[InMa(Nmf,1),ReMa(counter4,1)],'k','LineWidth',2);
            
        end
        
        InMa = [];
        InMa = ReMa;
        ReMa = zeros(1,8);
        counter4 = 0;
        mass_LC = [];
              
    end
    
    plot([DataPoints(counter,2),DataPoints(Nlastchar,2)],[DataPoints(counter,1),DataPoints(Nlastchar,1)],'k','LineWidth',2);
    
    NWallPoints = NWallPoints + 1;
    WallPoints(NWallPoints,1:8) = DataPoints(Nlastchar,1:8);
    
    hold off    
end
%%

NAxisPoints = NAxisPoints2;
AxisPoints = AxisPoints2;

DataPoints3 = DataPoints;

figure(f2)
subplot(2,1,1);

% COMPELTE EXPANSION TO OBTAIN CURVES 
%Ponit 1 - Negative charactristic

ReMa = zeros(1,8);

y_1 = FinalChar1(1,1);
x_1 = FinalChar1(1,2);
u_1 = FinalChar1(1,3);
v_1 = FinalChar1(1,4);

[Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1);

counter = counter +1;
counter0 = 1;

DataPoints(counter,1:8) = Array(1,1:8);
ReMa(1,1:8) = Array(1,1:8);

NAxisPoints = NAxisPoints + 1;
AxisPoints(NAxisPoints,1:8) = Array(1,1:8);

Noo = numel(ReMa(:,1));

hold on
plot([FinalChar1(1,2),DataPoints(counter,2)],[FinalChar1(1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3])

for o=2:numel(FinalChar1(:,1))-1
    
    if o==22
        ddd = 0;
    end
    
    L_1 = o;
    L_2 = 2*L_1 - 1;
    
    for l=L_1:L_2

        if l==L_1
            
            %Ponit 1 - Negative charactristic
            y_1 = FinalChar1(o,1);
            x_1 = FinalChar1(o,2);
            u_1 = FinalChar1(o,3);
            v_1 = FinalChar1(o,4);
            
            %Ponit 2 - Positive characteristic
            y_2 = DataPoints(counter - (Noo - 1),1);
            x_2 = DataPoints(counter - (Noo - 1),2);
            u_2 = DataPoints(counter - (Noo - 1),3);
            v_2 = DataPoints(counter - (Noo - 1),4);
            
            [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
            
            counter = counter +1;
            
            DataPoints(counter,1:8) = Array(1,1:8);
            ReMa(counter0,1:8) = Array(1,1:8);
            
            g5 = plot([FinalChar1(o,2),DataPoints(counter,2)],[FinalChar1(o,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
            g6 = plot([DataPoints(counter - (Noo - 1) - 1,2),DataPoints(counter,2)],[DataPoints(counter - (Noo - 1) - 1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
            
            Tol = 0.001;
            
            Array_1 = Array;
            
            if Array(1,2) > FinalChar1(numel(FinalChar1(:,1)),2)
                
                delete(g5)
                delete(g6)
                
                %Positive characteristic
                
                y_low = DataPoints(counter - (Noo - 1)-1,1);
                x_low = DataPoints(counter - (Noo - 1)-1,2);
                u_low = DataPoints(counter - (Noo - 1)-1,3);
                v_low = DataPoints(counter - (Noo - 1)-1,4);
                
                y_high = DataPoints(counter,1);
                x_high = DataPoints(counter,2);
                u_high = DataPoints(counter,3);
                v_high = DataPoints(counter,4);
                
                while abs(Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                    
                    y_4 = (y_low + y_high)/2;
                    x_4 = (x_low + x_high)/2;
                    u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                    v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                    Vg_4 = sqrt(u_4^2 + v_4^2);
                    
                    [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                    
                    Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                    
                    if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                        
                        y_high = y_4;
                        x_high = x_4;
                        u_high = u_4;
                        v_high = v_4;
                        
                    elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                        
                        y_low = y_4;
                        x_low = x_4;
                        u_low = u_4;
                        v_low = v_4;
                        
                    end
                    
                end
                
                DataPoints(counter,1:8) = Array(1,1:8);
                ReMa(counter0,1:8) = Array(1,1:8);
                
                plot([DataPoints(counter - (Noo - 1)-1,2),DataPoints(counter,2)],[DataPoints(counter - (Noo - 1)-1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
                
                %Negative characteristic
                
                y_low = Array_1(1,1);
                x_low = Array_1(1,2);
                u_low = Array_1(1,3);
                v_low = Array_1(1,4);
                
                y_high = FinalChar1(o,1);
                x_high = FinalChar1(o,2);
                u_high = FinalChar1(o,3);
                v_high = FinalChar1(o,4);
                
                while abs(Array_1(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                    
                    y_4 = (y_low + y_high)/2;
                    x_4 = (x_low + x_high)/2;
                    u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                    v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                    Vg_4 = sqrt(u_4^2 + v_4^2);
                    
                    [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                    
                    Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                    Array_1 = Array;
                    
                    if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                        
                        y_low = y_4;
                        x_low = x_4;
                        u_low = u_4;
                        v_low = v_4;
                        
                    elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                        
                        y_high = y_4;
                        x_high = x_4;
                        u_high = u_4;
                        v_high = v_4;
                        
                    end
                    
                end
                
                counter = counter + 1;
                counter0 = counter0 + 1;
                DataPoints(counter,1:8) = Array(1,1:8);
                ReMa(counter0,1:8) = Array(1,1:8);
       
                plot([FinalChar1(o,2),DataPoints(counter,2)],[FinalChar1(o,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
                
                if DataPoints(counter,1)==DataPoints(counter-1,1)&&DataPoints(counter,2)==DataPoints(counter-1,2)
                    
                    DataPoints(counter,:) = [];
                    ReMa(counter0,:) = [];
                    counter = counter - 1;
                    counter0 = counter0 - 1;
                    
                end
                
            end
            
            if abs(DataPoints(counter,2) - FinalChar1(numel(FinalChar1(:,1)),2)) <= Tol
                
                break
                
            end
            
        elseif l==L_2
            
            %Ponit 1 - Negative charactristic
            y_1 = DataPoints(counter,1);
            x_1 = DataPoints(counter,2);
            u_1 = DataPoints(counter,3);
            v_1 = DataPoints(counter,4);
            
            [Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1);
            
            counter = counter + 1;
            counter0 = counter0 + 1;
            DataPoints(counter,1:8) = Array(1,1:8);
            ReMa(counter0,1:8) = Array(1,1:8);
            
            NAxisPoints = NAxisPoints + 1;
            AxisPoints(NAxisPoints,1:8) = Array(1,1:8);
            
            g7 = plot([DataPoints(counter - 1,2),DataPoints(counter,2)],[DataPoints(counter - 1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
            
            if Array(1,2) > FinalChar1(numel(FinalChar1(:,1)),2)
                
                delete(g7)
                
                %Negative characteristic
                
                y_low = DataPoints(counter,1);
                x_low = DataPoints(counter,2);
                u_low = DataPoints(counter,3);
                v_low = DataPoints(counter,4);
                
                y_high = DataPoints(counter - 1,1);
                x_high = DataPoints(counter - 1,2);
                u_high = DataPoints(counter - 1,3);
                v_high = DataPoints(counter - 1,4);
                
                while abs(Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                    
                    y_4 = (y_low + y_high)/2;
                    x_4 = (x_low + x_high)/2;
                    u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                    v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                    Vg_4 = sqrt(u_4^2 + v_4^2);
                    
                    [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                    
                    Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                    
                    if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                        
                        y_low = y_4;
                        x_low = x_4;
                        u_low = u_4;
                        v_low = v_4;
                        
                    elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                        
                        y_high = y_4;
                        x_high = x_4;
                        u_high = u_4;
                        v_high = v_4;
                        
                    end
                    
                end
                
                DataPoints(counter,1:8) = Array(1,1:8);
                ReMa(counter0,1:8) = Array(1,1:8);
                
                plot([DataPoints(counter-1,2),DataPoints(counter,2)],[DataPoints(counter-1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);

            end
            
            if abs(DataPoints(counter,2) - FinalChar1(numel(FinalChar1(:,1)),2)) <= Tol
                
                break
                
            end
            
        else
            
            %Ponit 1 - Negative charactristic
            y_1 = DataPoints(counter,1);
            x_1 = DataPoints(counter,2);
            u_1 = DataPoints(counter,3);
            v_1 = DataPoints(counter,4);
            
            %Ponit 2 - Positive characteristic
            y_2 = DataPoints(counter - (Noo - 1),1);
            x_2 = DataPoints(counter - (Noo - 1),2);
            u_2 = DataPoints(counter - (Noo - 1),3);
            v_2 = DataPoints(counter - (Noo - 1),4);
            
            [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
            
            counter = counter +1;
            counter0 = counter0 + 1;
            
            ReMa(counter0,1:8) = Array(1,1:8);
            DataPoints(counter,1:8) = Array(1,1:8);
            
            g8 = plot([DataPoints(counter - 1,2),DataPoints(counter,2)],[DataPoints(counter - 1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
            g9 = plot([DataPoints(counter - (Noo - 1) - 1,2),DataPoints(counter,2)],[DataPoints(counter - (Noo - 1) - 1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
            
            Array_1 = Array;
            
            if Array(1,2) > FinalChar1(numel(FinalChar1(:,1)),2)
                
                delete(g8)
                delete(g9)
                
                %Positive characteristic
                
                y_low = DataPoints(counter - (Noo - 1)-1,1);
                x_low = DataPoints(counter - (Noo - 1)-1,2);
                u_low = DataPoints(counter - (Noo - 1)-1,3);
                v_low = DataPoints(counter - (Noo - 1)-1,4);
                
                y_high = DataPoints(counter,1);
                x_high = DataPoints(counter,2);
                u_high = DataPoints(counter,3);
                v_high = DataPoints(counter,4);
                
                while abs(Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                    
                    y_4 = (y_low + y_high)/2;
                    x_4 = (x_low + x_high)/2;
                    u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                    v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                    Vg_4 = sqrt(u_4^2 + v_4^2);
                    
                    [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                    
                    Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                    
                    if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                        
                        y_high = y_4;
                        x_high = x_4;
                        u_high = u_4;
                        v_high = v_4;
                        
                    elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                        
                        y_low = y_4;
                        x_low = x_4;
                        u_low = u_4;
                        v_low = v_4;
                        
                    end
                    
                end
                
                DataPoints(counter,1:8) = Array(1,1:8);
                ReMa(counter0,1:8) = Array(1,1:8);
                
                plot([DataPoints(counter - (Noo - 1)-1,2),DataPoints(counter,2)],[DataPoints(counter - (Noo - 1)-1,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
                
                %Negative characteristic
                
                y_low = Array_1(1,1);
                x_low = Array_1(1,2);
                u_low = Array_1(1,3);
                v_low = Array_1(1,4);
                
                y_high = DataPoints(counter-1,1);
                x_high = DataPoints(counter-1,2);
                u_high = DataPoints(counter-1,3);
                v_high = DataPoints(counter-1,4);
                
                while abs(Array_1(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                    
                    y_4 = (y_low + y_high)/2;
                    x_4 = (x_low + x_high)/2;
                    u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                    v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                    Vg_4 = sqrt(u_4^2 + v_4^2);
                    
                    [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                    
                    Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                    Array_1 = Array;
                    
                    if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                        
                        y_low = y_4;
                        x_low = x_4;
                        u_low = u_4;
                        v_low = v_4;
                        
                    elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                        
                        y_high = y_4;
                        x_high = x_4;
                        u_high = u_4;
                        v_high = v_4;
                        
                    end
                    
                end
                
                counter = counter + 1;
                counter0 = counter0 + 1;
                DataPoints(counter,1:8) = Array(1,1:8);
                ReMa(counter0,1:8) = Array(1,1:8);
       
                plot([DataPoints(counter-2,2),DataPoints(counter,2)],[DataPoints(counter-2,1),DataPoints(counter,1)],'color',[0.3 0.3 0.3]);
                
                if DataPoints(counter,1)==DataPoints(counter-1,1)&&DataPoints(counter,2)==DataPoints(counter-1,2)
                    
                    DataPoints(counter,:) = [];
                    ReMa(counter0,:) = [];
                    counter = counter - 1;
                    counter0 = counter0 - 1;
                    
                end
                
            end
            
            if abs(DataPoints(counter,2) - FinalChar1(numel(FinalChar1(:,1)),2)) <= Tol
                
                break
                
            end
            
        end
                
    end
    
    Noo = numel(ReMa(:,1));
    ReMa = zeros(1,8);
    counter0 = 1;
    
end

hold off
%%

subplot(2,1,2);

% Polynomial fit to the nozzle shape

X = WallPoints(NArc:NWallPoints,2);
Y = WallPoints(NArc:NWallPoints,1);
YY = WallPoints(:,1);
XX = WallPoints(:,2);

syms a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15
eqn1 = WallPoints(NArc,1) == a0 + a1*WallPoints(NArc,2) + a2*WallPoints(NArc,2)^2 + a3*WallPoints(NArc,2)^3 + a4*WallPoints(NArc,2)^4 + a5*WallPoints(NArc,2)^5 + a6*WallPoints(NArc,2)^6 + a7*WallPoints(NArc,2)^7 + a8*WallPoints(NArc,2)^8 + a9*WallPoints(NArc,2)^9 + a10*WallPoints(NArc,2)^10 + a11*WallPoints(NArc,2)^11 + a12*WallPoints(NArc,2)^12 + a13*WallPoints(NArc,2)^13 + a14*WallPoints(NArc,2)^14 + a15*WallPoints(NArc,2)^15;
eqn2 = WallPoints(NArc,4)/WallPoints(NArc,3) == a1 + 2*a2*WallPoints(NArc,2) + 3*a3*WallPoints(NArc,2)^2 + 4*a4*WallPoints(NArc,2)^3 + 5*a5*WallPoints(NArc,2)^4 + 6*a6*WallPoints(NArc,2)^5 + 7*a7*WallPoints(NArc,2)^6 + 8*a8*WallPoints(NArc,2)^7 + 9*a9*WallPoints(NArc,2)^8 + 10*a10*WallPoints(NArc,2)^9 + 11*a11*WallPoints(NArc,2)^10 + 12*a12*WallPoints(NArc,2)^11 + 13*a13*WallPoints(NArc,2)^12 + 14*a14*WallPoints(NArc,2)^13 + 15*a15*WallPoints(NArc,2)^14;
eqn3 = WallPoints(NArc + 6,1) == a0 + a1*WallPoints(NArc + 6,2) + a2*WallPoints(NArc + 6,2)^2 + a3*WallPoints(NArc + 6,2)^3 + a4*WallPoints(NArc + 6,2)^4 + a5*WallPoints(NArc + 6,2)^5 + a6*WallPoints(NArc + 6,2)^6 + a7*WallPoints(NArc + 6,2)^7 + a8*WallPoints(NArc + 6,2)^8 + a9*WallPoints(NArc + 6,2)^9 + a10*WallPoints(NArc + 6,2)^10 + a11*WallPoints(NArc + 6,2)^11 + a12*WallPoints(NArc + 6,2)^12 + a13*WallPoints(NArc + 6,2)^13 + a14*WallPoints(NArc + 6,2)^14 + a15*WallPoints(NArc + 6,2)^15;
eqn4 = WallPoints(NArc + 6,4)/WallPoints(NArc + 6,3) == a1 + 2*a2*WallPoints(NArc + 6,2) + 3*a3*WallPoints(NArc + 6,2)^2 + 4*a4*WallPoints(NArc + 6,2)^3 + 5*a5*WallPoints(NArc + 6,2)^4 + 6*a6*WallPoints(NArc + 6,2)^5 + 7*a7*WallPoints(NArc + 6,2)^6 + 8*a8*WallPoints(NArc + 6,2)^7 + 9*a9*WallPoints(NArc + 6,2)^8 + 10*a10*WallPoints(NArc + 6,2)^9 + 11*a11*WallPoints(NArc + 6,2)^10 + 12*a12*WallPoints(NArc + 6,2)^11 + 13*a13*WallPoints(NArc + 6,2)^12 + 14*a14*WallPoints(NArc + 6,2)^13 + 15*a15*WallPoints(NArc + 6,2)^14;
eqn5 = WallPoints(NArc + 12,1) == a0 + a1*WallPoints(NArc + 12,2) + a2*WallPoints(NArc + 12,2)^2 + a3*WallPoints(NArc + 12,2)^3 + a4*WallPoints(NArc + 12,2)^4 + a5*WallPoints(NArc + 12,2)^5 + a6*WallPoints(NArc + 12,2)^6 + a7*WallPoints(NArc + 12,2)^7 + a8*WallPoints(NArc + 12,2)^8 + a9*WallPoints(NArc + 12,2)^9 + a10*WallPoints(NArc + 12,2)^10 + a11*WallPoints(NArc + 12,2)^11 + a12*WallPoints(NArc + 12,2)^12 + a13*WallPoints(NArc + 12,2)^13 + a14*WallPoints(NArc + 12,2)^14 + a15*WallPoints(NArc + 12,2)^15;
eqn6 = WallPoints(NArc + 12,4)/WallPoints(NArc + 12,3) == a1 + 2*a2*WallPoints(NArc + 12,2) + 3*a3*WallPoints(NArc + 12,2)^2 + 4*a4*WallPoints(NArc + 12,2)^3 + 5*a5*WallPoints(NArc + 12,2)^4 + 6*a6*WallPoints(NArc + 12,2)^5 + 7*a7*WallPoints(NArc + 12,2)^6 + 8*a8*WallPoints(NArc + 12,2)^7 + 9*a9*WallPoints(NArc + 12,2)^8 + 10*a10*WallPoints(NArc + 12,2)^9 + 11*a11*WallPoints(NArc + 12,2)^10 + 12*a12*WallPoints(NArc + 12,2)^11 + 13*a13*WallPoints(NArc + 12,2)^12 + 14*a14*WallPoints(NArc + 12,2)^13 + 15*a15*WallPoints(NArc + 12,2)^14;
eqn7 = WallPoints(NArc + 18,1) == a0 + a1*WallPoints(NArc + 18,2) + a2*WallPoints(NArc + 18,2)^2 + a3*WallPoints(NArc + 18,2)^3 + a4*WallPoints(NArc + 18,2)^4 + a5*WallPoints(NArc + 18,2)^5 + a6*WallPoints(NArc + 18,2)^6 + a7*WallPoints(NArc + 18,2)^7 + a8*WallPoints(NArc + 18,2)^8 + a9*WallPoints(NArc + 18,2)^9 + a10*WallPoints(NArc + 18,2)^10 + a11*WallPoints(NArc + 18,2)^11 + a12*WallPoints(NArc + 18,2)^12 + a13*WallPoints(NArc + 18,2)^13 + a14*WallPoints(NArc + 18,2)^14 + a15*WallPoints(NArc + 18,2)^15;
eqn8 = WallPoints(NArc + 18,4)/WallPoints(NArc + 18,3) == a1 + 2*a2*WallPoints(NArc + 18,2) + 3*a3*WallPoints(NArc + 18,2)^2 + 4*a4*WallPoints(NArc + 18,2)^3 + 5*a5*WallPoints(NArc + 18,2)^4 + 6*a6*WallPoints(NArc + 18,2)^5 + 7*a7*WallPoints(NArc + 18,2)^6 + 8*a8*WallPoints(NArc + 18,2)^7 + 9*a9*WallPoints(NArc + 18,2)^8 + 10*a10*WallPoints(NArc + 18,2)^9 + 11*a11*WallPoints(NArc + 18,2)^10 + 12*a12*WallPoints(NArc + 18,2)^11 + 13*a13*WallPoints(NArc + 18,2)^12 + 14*a14*WallPoints(NArc + 18,2)^13 + 15*a15*WallPoints(NArc + 18,2)^14;
eqn9 = WallPoints(NArc + 24,1) == a0 + a1*WallPoints(NArc + 24,2) + a2*WallPoints(NArc + 24,2)^2 + a3*WallPoints(NArc + 24,2)^3 + a4*WallPoints(NArc + 24,2)^4 + a5*WallPoints(NArc + 24,2)^5 + a6*WallPoints(NArc + 24,2)^6 + a7*WallPoints(NArc + 24,2)^7 + a8*WallPoints(NArc + 24,2)^8 + a9*WallPoints(NArc + 24,2)^9 + a10*WallPoints(NArc + 24,2)^10 + a11*WallPoints(NArc + 24,2)^11 + a12*WallPoints(NArc + 24,2)^12 + a13*WallPoints(NArc + 24,2)^13 + a14*WallPoints(NArc + 24,2)^14 + a15*WallPoints(NArc + 24,2)^15;
eqn10 = WallPoints(NArc + 24,4)/WallPoints(NArc + 24,3) == a1 + 2*a2*WallPoints(NArc + 24,2) + 3*a3*WallPoints(NArc + 24,2)^2 + 4*a4*WallPoints(NArc + 24,2)^3 + 5*a5*WallPoints(NArc + 24,2)^4 + 6*a6*WallPoints(NArc + 24,2)^5 + 7*a7*WallPoints(NArc + 24,2)^6 + 8*a8*WallPoints(NArc + 24,2)^7 + 9*a9*WallPoints(NArc + 24,2)^8 + 10*a10*WallPoints(NArc + 24,2)^9 + 11*a11*WallPoints(NArc + 24,2)^10 + 12*a12*WallPoints(NArc + 24,2)^11 + 13*a13*WallPoints(NArc + 24,2)^12 + 14*a14*WallPoints(NArc + 24,2)^13 + 15*a15*WallPoints(NArc + 24,2)^14;
eqn11 = WallPoints(NArc + 30,1) == a0 + a1*WallPoints(NArc + 30,2) + a2*WallPoints(NArc + 30,2)^2 + a3*WallPoints(NArc + 30,2)^3 + a4*WallPoints(NArc + 30,2)^4 + a5*WallPoints(NArc + 30,2)^5 + a6*WallPoints(NArc + 30,2)^6 + a7*WallPoints(NArc + 30,2)^7 + a8*WallPoints(NArc + 30,2)^8 + a9*WallPoints(NArc + 30,2)^9 + a10*WallPoints(NArc + 30,2)^10 + a11*WallPoints(NArc + 30,2)^11 + a12*WallPoints(NArc + 30,2)^12 + a13*WallPoints(NArc + 30,2)^13 + a14*WallPoints(NArc + 30,2)^14 + a15*WallPoints(NArc + 30,2)^15;
eqn12 = WallPoints(NArc + 30,4)/WallPoints(NArc + 30,3) == a1 + 2*a2*WallPoints(NArc + 30,2) + 3*a3*WallPoints(NArc + 30,2)^2 + 4*a4*WallPoints(NArc + 30,2)^3 + 5*a5*WallPoints(NArc + 30,2)^4 + 6*a6*WallPoints(NArc + 30,2)^5 + 7*a7*WallPoints(NArc + 30,2)^6 + 8*a8*WallPoints(NArc + 30,2)^7 + 9*a9*WallPoints(NArc + 30,2)^8 + 10*a10*WallPoints(NArc + 30,2)^9 + 11*a11*WallPoints(NArc + 30,2)^10 + 12*a12*WallPoints(NArc + 30,2)^11 + 13*a13*WallPoints(NArc + 30,2)^12 + 14*a14*WallPoints(NArc + 30,2)^13 + 15*a15*WallPoints(NArc + 30,2)^14;
eqn13 = WallPoints(NArc + 33,1) == a0 + a1*WallPoints(NArc + 33,2) + a2*WallPoints(NArc + 33,2)^2 + a3*WallPoints(NArc + 33,2)^3 + a4*WallPoints(NArc + 33,2)^4 + a5*WallPoints(NArc + 33,2)^5 + a6*WallPoints(NArc + 33,2)^6 + a7*WallPoints(NArc + 33,2)^7 + a8*WallPoints(NArc + 33,2)^8 + a9*WallPoints(NArc + 33,2)^9 + a10*WallPoints(NArc + 33,2)^10 + a11*WallPoints(NArc + 33,2)^11 + a12*WallPoints(NArc + 33,2)^12 + a13*WallPoints(NArc + 33,2)^13 + a14*WallPoints(NArc + 33,2)^14 + a15*WallPoints(NArc + 33,2)^15;
eqn14 = WallPoints(NArc + 33,4)/WallPoints(NArc + 33,3) == a1 + 2*a2*WallPoints(NArc + 33,2) + 3*a3*WallPoints(NArc + 33,2)^2 + 4*a4*WallPoints(NArc + 33,2)^3 + 5*a5*WallPoints(NArc + 33,2)^4 + 6*a6*WallPoints(NArc + 33,2)^5 + 7*a7*WallPoints(NArc + 33,2)^6 + 8*a8*WallPoints(NArc + 33,2)^7 + 9*a9*WallPoints(NArc + 33,2)^8 + 10*a10*WallPoints(NArc + 33,2)^9 + 11*a11*WallPoints(NArc + 33,2)^10 + 12*a12*WallPoints(NArc + 33,2)^11 + 13*a13*WallPoints(NArc + 33,2)^12 + 14*a14*WallPoints(NArc + 33,2)^13 + 15*a15*WallPoints(NArc + 33,2)^14;
eqn15 = WallPoints(NWallPoints,1) == a0 + a1*WallPoints(NWallPoints,2) + a2*WallPoints(NWallPoints,2)^2 + a3*WallPoints(NWallPoints,2)^3 + a4*WallPoints(NWallPoints,2)^4 + a5*WallPoints(NWallPoints,2)^5 + a6*WallPoints(NWallPoints,2)^6 + a7*WallPoints(NWallPoints,2)^7 + a8*WallPoints(NWallPoints,2)^8 + a9*WallPoints(NWallPoints,2)^9 + a10*WallPoints(NWallPoints,2)^10 + a11*WallPoints(NWallPoints,2)^11 + a12*WallPoints(NWallPoints,2)^12 + a13*WallPoints(NWallPoints,2)^13 + a14*WallPoints(NWallPoints,2)^14 + a15*WallPoints(NWallPoints,2)^15;
eqn16 = WallPoints(NWallPoints,4)/WallPoints(NWallPoints,3) == a1 + 2*a2*WallPoints(NWallPoints,2) + 3*a3*WallPoints(NWallPoints,2)^2 + 4*a4*WallPoints(NWallPoints,2)^3 + 5*a5*WallPoints(NWallPoints,2)^4 + 6*a6*WallPoints(NWallPoints,2)^5 + 7*a7*WallPoints(NWallPoints,2)^6 + 8*a8*WallPoints(NWallPoints,2)^7 + 9*a9*WallPoints(NWallPoints,2)^8 + 10*a10*WallPoints(NWallPoints,2)^9 + 11*a11*WallPoints(NWallPoints,2)^10 + 12*a12*WallPoints(NWallPoints,2)^11 + 13*a13*WallPoints(NWallPoints,2)^12 + 14*a14*WallPoints(NWallPoints,2)^13 + 15*a15*WallPoints(NWallPoints,2)^14;

sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10, eqn11, eqn12, eqn13, eqn14, eqn15, eqn16], [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15]);
a0 = double(sol.a0);
a1 = double(sol.a1);
a2 = double(sol.a2);
a3 = double(sol.a3);
a4 = double(sol.a4);
a5 = double(sol.a5);
a6 = double(sol.a6);
a7 = double(sol.a7);
a8 = double(sol.a8);
a9 = double(sol.a9);
a10 = double(sol.a10);
a11 = double(sol.a11);
a12 = double(sol.a12);
a13 = double(sol.a13);
a14 = double(sol.a14);
a15 = double(sol.a15);

p = [a15,a14,a13,a12,a11,a10,a9,a8,a7,a6,a5,a4,a3,a2,a1,a0];
y_fit = polyval(p,X);

p_matrix = p;

EA = Y - y_fit;
EP = (EA./Y).*100;

Tabla = [Y,y_fit,EA,EP];

hold on
plot(X,y_fit,'k','LineWidth',2);

ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
xlabel('x/y_t','fontweight','bold')  
ylabel('y/y_t','fontweight','bold') 
grid on
grid minor

%% 
%POLYNOMIAL EXPANSION WAVES REGION
    
InMa = DataPoints((NKernell - 2*counter_2 - counter_7 + 1):(NKernell),1:8);
ReMa = zeros(1,8);
counter3 = 1;
counter4 = 0;

Xwall = 0;

while Xwall <= FinalChar1(numel(FinalChar1(:,1)),2)
    
    p = p_matrix;

    E = 2;
    y_2 = InMa(E,1);
    x_2 = InMa(E,2);
    u_2 = InMa(E,3);
    v_2 = InMa(E,4);
    
    [Array] = WallPoint(gamma, R, delta, p_st, T_st, a_0, x_2, y_2, u_2, v_2, p);
    
    counter2 = counter2 + 1;
    DataPoints2(counter2,1:8) = Array(1,1:8);
    
    NWallPoints2 = NWallPoints2 + 1;
    WallPoints2(NWallPoints2,1:8) = Array(1,1:8);
    
    ReMa(1,1:8) = Array(1,1:8);

    hold on
    g1 = plot([InMa(E,2),DataPoints2(counter2,2)],[InMa(E,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
        
    Xwall = DataPoints2(counter2,2);
    
    if Xwall > FinalChar1(numel(FinalChar1(:,1)),2)
        
        delete(g1)
        
        x_4 = FinalChar1(numel(FinalChar1(:,1)),2);
        y_4 = y_fit(numel(y_fit));
        theta_4 = atand(15*p(1,1)*x_4^14 + 14*p(1,2)*x_4^13 + 13*p(1,3)*x_4^12 + 12*p(1,4)*x_4^11 + 11*p(1,5)*x_4^10 + 10*p(1,6)*x_4^9 + 9*p(1,7)*x_4^8 + 8*p(1,8)*x_4^7 + 7*p(1,9)*x_4^6 + 6*p(1,10)*x_4^5 + 5*p(1,11)*x_4^4 + 4*p(1,12)*x_4^3 + 3*p(1,13)*x_4^2 + 2*p(1,14)*x_4 + p(1,15));
        
        x_3 = DataPoints2(counter2 - No,2);
        y_3 = DataPoints2(counter2 - No,1);
        u_3 = DataPoints2(counter2 - No,3);
        v_3 = DataPoints2(counter2 - No,4);
        
        x_1 = DataPoints2(counter2 - No + 1,2);
        y_1 = DataPoints2(counter2 - No + 1,1);
        u_1 = DataPoints2(counter2 - No + 1,3);
        v_1 = DataPoints2(counter2 - No + 1,4);
        
        [Array_2, Array_4] = InverseWallPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_3, y_3, u_3, v_3, x_4, y_4, theta_4);
        
        DataPoints2(counter2,1:8) = Array_2(1,1:8);     
        
        counter2 = counter2 + 1;
        DataPoints2(counter2,1:8) = Array_4(1,1:8);
        
        ReMa(1,1:8) = Array_4(1,1:8);
        
        E = 1;
        
        plot([DataPoints2(counter2-1,2),DataPoints2(counter2,2)],[DataPoints2(counter2-1,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
 
    end
    
    Nmf = numel(InMa(:,1));
    
    counter4 = E;
    
    for i = E + 1: Nmf
        
        %Ponit 1 - Negative charactristic
        y_1 = DataPoints2(counter2,1);
        x_1 = DataPoints2(counter2,2);
        u_1 = DataPoints2(counter2,3);
        v_1 = DataPoints2(counter2,4);
        
        %Ponit 2 - Positive characteristic
        counter4 = counter4 +1;
        y_2 = InMa(counter4,1);
        x_2 = InMa(counter4,2);
        u_2 = InMa(counter4,3);
        v_2 = InMa(counter4,4);
        
        [Array] = InnerPoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1, x_2, y_2, u_2, v_2);
        
        counter2 = counter2 + 1;
        DataPoints2(counter2,1:8) = Array(1,1:8);
        
        counter3 = counter3 + 1;
        ReMa(counter3,1:8) = Array(1,1:8);
        
        g2 = plot([InMa(counter4,2),DataPoints2(counter2,2)],[InMa(counter4,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
        g3 = plot([DataPoints2(counter2 - 1,2),DataPoints2(counter2,2)],[DataPoints2(counter2 - 1,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
        
        Array_1 = Array;
        
        Tol = 0.001;
        
        if Array(1,2) > FinalChar1(numel(FinalChar1(:,1)),2)
            
            delete(g2)
            delete(g3)
            
            %Positive characteristic
            
            y_low = InMa(counter4,1);
            x_low = InMa(counter4,2);
            u_low = InMa(counter4,3);
            v_low = InMa(counter4,4);
            
            y_high = DataPoints2(counter2,1);
            x_high = DataPoints2(counter2,2);
            u_high = DataPoints2(counter2,3);
            v_high = DataPoints2(counter2,4);
          
            while abs(Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                
                y_4 = (y_low + y_high)/2;
                x_4 = (x_low + x_high)/2;
                u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                Vg_4 = sqrt(u_4^2 + v_4^2);               

                [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                
                Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                
                if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                    
                    y_high = y_4;
                    x_high = x_4;
                    u_high = u_4;
                    v_high = v_4;
                    
                elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                    
                    y_low = y_4;
                    x_low = x_4;
                    u_low = u_4;
                    v_low = v_4;
                    
                end
                
            end
            
            DataPoints2(counter2,1:8) = Array(1,1:8);
            
            plot([InMa(counter4,2),DataPoints2(counter2,2)],[InMa(counter4,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
            
            %Negative characteristic
            
            y_low = Array_1(1,1);
            x_low = Array_1(1,2);
            u_low = Array_1(1,3);
            v_low = Array_1(1,4);
            
            y_high = DataPoints2(counter2-1,1);
            x_high = DataPoints2(counter2-1,2);
            u_high = DataPoints2(counter2-1,3);
            v_high = DataPoints2(counter2-1,4);
            
            while abs(Array_1(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                
                y_4 = (y_low + y_high)/2;
                x_4 = (x_low + x_high)/2;
                u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                Vg_4 = sqrt(u_4^2 + v_4^2);
                
                [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                
                Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                Array_1 = Array;
                
                if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                    
                    y_low = y_4;
                    x_low = x_4;
                    u_low = u_4;
                    v_low = v_4;
                    
                elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                    
                    y_high = y_4;
                    x_high = x_4;
                    u_high = u_4;
                    v_high = v_4;
                    
                end
                
            end
            
            counter2 = counter2 + 1;
            DataPoints2(counter2,1:8) = Array(1,1:8);
            
            ReMa(counter3,1:8) = Array(1,1:8);
            
            plot([DataPoints2(counter2-2,2),DataPoints2(counter2,2)],[DataPoints2(counter2-2,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
            
        end
        
    end
    
    if abs(DataPoints2(counter2,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
        
        y_1 = DataPoints2(counter2,1);
        x_1 = DataPoints2(counter2,2);
        u_1 = DataPoints2(counter2,3);
        v_1 = DataPoints2(counter2,4);
        
        [Array] = Axipoint(gamma, R, delta, p_st, T_st, a_0, x_1, y_1, u_1, v_1);
        
        counter2 = counter2 + 1;
        DataPoints2(counter2,1:8) = Array(1,1:8);
        
        NAxisPoints2 = NAxisPoints2 + 1;
        AxisPoints2(NAxisPoints2,1:8) = Array(1,1:8);
        
        counter3 = counter3 + 1;
        ReMa(counter3,1:8) = Array(1,1:8);
        
        g4 = plot([DataPoints2(counter2 - 1,2),DataPoints2(counter2,2)],[DataPoints2(counter2 - 1,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
        
        if Array(1,2) > FinalChar1(numel(FinalChar1(:,1)),2)
            
            delete(g4)
            
            y_low = Array(1,1);
            x_low = Array(1,2);
            u_low = Array(1,3);
            v_low = Array(1,4);
            
            y_high = DataPoints2(counter2-1,1);
            x_high = DataPoints2(counter2-1,2);
            u_high = DataPoints2(counter2-1,3);
            v_high = DataPoints2(counter2-1,4);
            
            while abs(Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2)) > Tol
                
                y_4 = (y_low + y_high)/2;
                x_4 = (x_low + x_high)/2;
                u_4 = u_high - (y_high - y_4)/(y_high - y_low)*(u_high - u_low);
                v_4 = v_high - (y_high - y_4)/(y_high - y_low)*(v_high - v_low);
                Vg_4 = sqrt(u_4^2 + v_4^2);
                
                [Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st);
                
                Array = [y_4, x_4, u_4, v_4, Ma t, rho, p];
                
                if Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) > Tol
                    
                    y_low = y_4;
                    x_low = x_4;
                    u_low = u_4;
                    v_low = v_4;
                    
                elseif Array(1,2) - FinalChar1(numel(FinalChar1(:,1)),2) < -Tol
                    
                    y_high = y_4;
                    x_high = x_4;
                    u_high = u_4;
                    v_high = v_4;
                    
                end
                
            end
            
            DataPoints2(counter2,1:8) = Array(1,1:8);
            
            ReMa(counter3,1:8) = Array(1,1:8);
            
            plot([DataPoints2(counter2-1,2),DataPoints2(counter2,2)],[DataPoints2(counter2-1,1),DataPoints2(counter2,1)],'color',[0.3 0.3 0.3]);
            
            
        end

    end
    
    No = numel(InMa(:,1));
    InMa = [];
    InMa = ReMa;
    ReMa = zeros(1,8);
    counter4 = 0;
    counter3 = 1;
    
end

title('16th Degree Polynomial Contour')
hold off

subplot(2,1,1);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
xlabel('x/y_t','fontweight','bold')  
ylabel('y/y_t','fontweight','bold') 
title('Ideal Contour')

grid on
grid minor

axis([SP21 SP22],[-1.637 DataPoints2(counter2,2) 0 2])

figure(f1)

ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([-1.637 FinalChar1(numel(FinalChar1(:,1)),2)+0.2 0 2])
xlabel('x/y_t','fontweight','bold')  
ylabel('y/y_t','fontweight','bold') 
title('CD Ideal Nozzle Design')

grid on
grid minor

%% 

%Graphic analysis

YY2 = WallPoints2(:,1);
XX2 = WallPoints2(:,2);

x_SS = [];
y_SS = [];
M_SS = [];
t_SS = [];
rho_SS = [];
p_SS = [];

x_SS2 = [];
y_SS2 = [];
M_SS2 = [];
t_SS2 = [];
rho_SS2 = [];
p_SS2 = [];

x_SS3 = [];
y_SS3 = [];
M_SS3 = [];
t_SS3 = [];
rho_SS3 = [];
p_SS3 = [];

for i=1:numel(YY)
    
    x = XX(i);
    y = YY(i);
    A = pi*y^2;
    
    syms z
    eqn = A/A_t == 1/z*((2/(gamma + 1))*(1 + ((gamma - 1)/2)*z^2))^((gamma + 1)/(2*(gamma - 1)));
    sol = vpasolve(eqn, z, [1 M_axis]);
    M = double(sol);
    
    t = T_st/(1 + (gamma - 1)/2*M^2);                 % Static temperature
    p = p_st/(1 + (gamma - 1)/2*M^2)^(gamma/(gamma - 1));
    rho = p/(R*t);
    
    x_SS = [x_SS;x];
    y_SS = [y_SS;y];
    M_SS = [M_SS;M];
    t_SS = [t_SS;t];
    rho_SS = [rho_SS;rho];
    p_SS = [p_SS;p];
    
end

SuperSonic = [y_SS,x_SS,M_SS,t_SS,rho_SS,p_SS];

for i=1:numel(YY)
    
    x = XX(i);
    y = YY(i);
    A = pi*y^2;
    
    syms z
    eqn = A/A_t == 1/z*((2/(gamma + 1))*(1 + ((gamma - 1)/2)*z^2))^((gamma + 1)/(2*(gamma - 1)));
    sol = vpasolve(eqn, z, [0 1]);
    M = double(sol);
    
    t = T_st/(1 + (gamma - 1)/2*M^2);                 % Static temperature
    p = p_st/(1 + (gamma - 1)/2*M^2)^(gamma/(gamma - 1));
    rho = p/(R*t);
    
    x_SS2 = [x_SS2;x];
    y_SS2 = [y_SS2;y];
    M_SS2 = [M_SS2;M];
    t_SS2 = [t_SS2;t];
    rho_SS2 = [rho_SS2;rho];
    p_SS2 = [p_SS2;p];
    
end

SuperSonic2 = [y_SS2,x_SS2,M_SS2,t_SS2,rho_SS2,p_SS2];

for i=1:numel(YY2)
    
    x = XX2(i);
    y = YY2(i);
    A = pi*y^2;
    
    syms z
    eqn = A/A_t == 1/z*((2/(gamma + 1))*(1 + ((gamma - 1)/2)*z^2))^((gamma + 1)/(2*(gamma - 1)));
    sol = vpasolve(eqn, z, [1 M_axis]);
    M = double(sol);
    
    t = T_st/(1 + (gamma - 1)/2*M^2);                 % Static temperature
    p = p_st/(1 + (gamma - 1)/2*M^2)^(gamma/(gamma - 1));
    rho = p/(R*t);
    
    x_SS3 = [x_SS3;x];
    y_SS3 = [y_SS3;y];
    M_SS3 = [M_SS3;M];
    t_SS3 = [t_SS3;t];
    rho_SS3 = [rho_SS3;rho];
    p_SS3 = [p_SS3;p];
    
end

SuperSonic5 = [y_SS3,x_SS3,M_SS3,t_SS3,rho_SS3,p_SS3];

figure(f3)

subsonic3 = [subsonic(1:numel(subsonic(:,1)),1), subsonic(1:numel(subsonic(:,1)),2), subsonic(1:numel(subsonic(:,1)),3)/M_axis, subsonic(1:numel(subsonic(:,1)),4)/T_st, subsonic(1:numel(subsonic(:,1)),5)/rho_st, subsonic(1:numel(subsonic(:,1)),6)/p_st];
WallPoints3  = [WallPoints(1:numel(WallPoints(:,1)),1), WallPoints(1:numel(WallPoints(:,1)),2), WallPoints(1:numel(WallPoints(:,1)),3), WallPoints(1:numel(WallPoints(:,1)),4), WallPoints(1:numel(WallPoints(:,1)),5)/M_axis, WallPoints(1:numel(WallPoints(:,1)),6)/T_st, WallPoints(1:numel(WallPoints(:,1)),7)/rho_st, WallPoints(1:numel(WallPoints(:,1)),8)/p_st];
SuperSonic3 = [SuperSonic(1:numel(SuperSonic(:,1)),1), SuperSonic(1:numel(SuperSonic(:,1)),2), SuperSonic(1:numel(SuperSonic(:,1)),3)/M_axis, SuperSonic(1:numel(SuperSonic(:,1)),4)/T_st, SuperSonic(1:numel(SuperSonic(:,1)),5)/rho_st, SuperSonic(1:numel(SuperSonic(:,1)),6)/p_st];
SuperSonic4 = [SuperSonic2(1:numel(SuperSonic(:,1)),1), SuperSonic2(1:numel(SuperSonic(:,1)),2), SuperSonic2(1:numel(SuperSonic(:,1)),3)/M_axis, SuperSonic2(1:numel(SuperSonic(:,1)),4)/T_st, SuperSonic2(1:numel(SuperSonic(:,1)),5)/rho_st, SuperSonic2(1:numel(SuperSonic(:,1)),6)/p_st];
AxisPoints3 = [AxisPoints(1:numel(AxisPoints(:,1)),1), AxisPoints(1:numel(AxisPoints(:,1)),2), AxisPoints(1:numel(AxisPoints(:,1)),3), AxisPoints(1:numel(AxisPoints(:,1)),4), AxisPoints(1:numel(AxisPoints(:,1)),5)/M_axis, AxisPoints(1:numel(AxisPoints(:,1)),6)/T_st, AxisPoints(1:numel(AxisPoints(:,1)),7)/rho_st, AxisPoints(1:numel(AxisPoints(:,1)),8)/p_st];

WallPoints4  = [WallPoints2(1:numel(WallPoints2(:,1)),1), WallPoints2(1:numel(WallPoints2(:,1)),2), WallPoints2(1:numel(WallPoints2(:,1)),3), WallPoints2(1:numel(WallPoints2(:,1)),4), WallPoints2(1:numel(WallPoints2(:,1)),5)/M_axis, WallPoints2(1:numel(WallPoints2(:,1)),6)/T_st, WallPoints2(1:numel(WallPoints2(:,1)),7)/rho_st, WallPoints2(1:numel(WallPoints2(:,1)),8)/p_st];
AxisPoints4 = [AxisPoints2(1:numel(AxisPoints2(:,1)),1), AxisPoints2(1:numel(AxisPoints2(:,1)),2), AxisPoints2(1:numel(AxisPoints2(:,1)),3), AxisPoints2(1:numel(AxisPoints2(:,1)),4), AxisPoints2(1:numel(AxisPoints2(:,1)),5)/M_axis, AxisPoints2(1:numel(AxisPoints2(:,1)),6)/T_st, AxisPoints2(1:numel(AxisPoints2(:,1)),7)/rho_st, AxisPoints2(1:numel(AxisPoints2(:,1)),8)/p_st];
SuperSonic6 = [SuperSonic5(1:numel(SuperSonic5(:,1)),1), SuperSonic5(1:numel(SuperSonic5(:,1)),2), SuperSonic5(1:numel(SuperSonic5(:,1)),3)/M_axis, SuperSonic5(1:numel(SuperSonic5(:,1)),4)/T_st, SuperSonic5(1:numel(SuperSonic5(:,1)),5)/rho_st, SuperSonic5(1:numel(SuperSonic5(:,1)),6)/p_st];

plot(subsonic3(1:numel(subsonic(:,1)),2),subsonic3(1:numel(subsonic(:,1)),6),'k','LineWidth',1.5)
hold on
plot(subsonic3(1:numel(subsonic(:,1)),2),subsonic3(1:numel(subsonic(:,1)),3),':k','LineWidth',1.5)
plot(subsonic3(1:numel(subsonic(:,1)),2),subsonic3(1:numel(subsonic(:,1)),4),'--k','LineWidth',1.5)
plot(subsonic3(1:numel(subsonic(:,1)),2),subsonic3(1:numel(subsonic(:,1)),5),'-.k','LineWidth',1.5)
plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),6),'k','LineWidth',1.5)
plot(SuperSonic4(1:numel(SuperSonic(:,1)),2),SuperSonic4(1:numel(SuperSonic(:,1)),6),'k','LineWidth',1.5)
plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),3),':k','LineWidth',1.5)
plot(SuperSonic4(1:numel(SuperSonic(:,1)),2),SuperSonic4(1:numel(SuperSonic(:,1)),3),':k','LineWidth',1.5)
plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),4),'--k','LineWidth',1.5)
plot(SuperSonic4(1:numel(SuperSonic(:,1)),2),SuperSonic4(1:numel(SuperSonic(:,1)),4),'--k','LineWidth',1.5)
plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),5),'-.k','LineWidth',1.5)
plot(SuperSonic4(1:numel(SuperSonic(:,1)),2),SuperSonic4(1:numel(SuperSonic(:,1)),5),'-.k','LineWidth',1.5)

lgd = legend({'Pressure','Mach','Temperature','Density'},'Location','best');
title(lgd,'Properties')
xlabel('x/y_t','fontweight','bold')  
ylabel('p/p_0, T/T_0, \rho/\rho_0, M/M_{Exit}','fontweight','bold') 
title('Analysis of the flow field - One-dimensional')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([-1.637 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

%%

% Figure F4

figure(f4)

SP41 = subplot(2,2,1);

plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),6),'k','LineWidth',1.5)
hold on
plot(WallPoints3(1:NWallPoints,2),WallPoints3(1:NWallPoints,8),'--k','LineWidth',1.5)
plot(AxisPoints3(1:NAxisPoints,2),AxisPoints3(1:NAxisPoints,8),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('p/p_0','fontweight','bold') 
title('Analysis of Pressure')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

SP42 = subplot(2,2,2);

plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),4),'k','LineWidth',1.5)
hold on
plot(WallPoints3(1:NWallPoints,2),WallPoints3(1:NWallPoints,6),'--k','LineWidth',1.5)
plot(AxisPoints3(1:NAxisPoints,2),AxisPoints3(1:NAxisPoints,6),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('T/T_0','fontweight','bold') 
title('Analysis of Temperature')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

SP43 = subplot(2,2,3);

plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),5),'k','LineWidth',1.5)
hold on
plot(WallPoints3(1:NWallPoints,2),WallPoints3(1:NWallPoints,7),'--k','LineWidth',1.5)
plot(AxisPoints3(1:NAxisPoints,2),AxisPoints3(1:NAxisPoints,7),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('\rho / \rho_0','fontweight','bold') 
title('Analysis of Density')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

SP44 = subplot(2,2,4);

plot(SuperSonic3(1:numel(SuperSonic(:,1)),2),SuperSonic3(1:numel(SuperSonic(:,1)),3),'k','LineWidth',1.5)
hold on
plot(WallPoints3(1:NWallPoints,2),WallPoints3(1:NWallPoints,5),'--k','LineWidth',1.5)
plot(AxisPoints3(1:NAxisPoints,2),AxisPoints3(1:NAxisPoints,5),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('M/M_Exit','fontweight','bold') 
title('Analysis of Mach Number')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

%%

% Figure F5

figure(f5)

SP51 = subplot(2,2,1);

plot(SuperSonic6(1:numel(SuperSonic6(:,1)),2),SuperSonic6(1:numel(SuperSonic6(:,1)),6),'k','LineWidth',1.5)
hold on
plot(WallPoints4(1:NWallPoints2,2),WallPoints4(1:NWallPoints2,8),'--k','LineWidth',1.5)
plot(AxisPoints4(1:NAxisPoints2,2),AxisPoints4(1:NAxisPoints2,8),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('p/p_0','fontweight','bold') 
title('Analysis of Pressure')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

SP52 = subplot(2,2,2);

plot(SuperSonic6(1:numel(SuperSonic6(:,1)),2),SuperSonic6(1:numel(SuperSonic6(:,1)),4),'k','LineWidth',1.5)
hold on
plot(WallPoints4(1:NWallPoints2,2),WallPoints4(1:NWallPoints2,6),'--k','LineWidth',1.5)
plot(AxisPoints4(1:NAxisPoints2,2),AxisPoints4(1:NAxisPoints2,6),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('T/T_0','fontweight','bold') 
title('Analysis of Temperature')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

SP53 = subplot(2,2,3);

plot(SuperSonic6(1:numel(SuperSonic6(:,1)),2),SuperSonic6(1:numel(SuperSonic6(:,1)),5),'k','LineWidth',1.5)
hold on
plot(WallPoints4(1:NWallPoints2,2),WallPoints4(1:NWallPoints2,7),'--k','LineWidth',1.5)
plot(AxisPoints4(1:NAxisPoints2,2),AxisPoints4(1:NAxisPoints2,7),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('\rho/\rho_0','fontweight','bold') 
title('Analysis of Density')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

SP54 = subplot(2,2,4);

plot(SuperSonic6(1:numel(SuperSonic6(:,1)),2),SuperSonic6(1:numel(SuperSonic6(:,1)),3),'k','LineWidth',1.5)
hold on
plot(WallPoints4(1:NWallPoints2,2),WallPoints4(1:NWallPoints2,5),'--k','LineWidth',1.5)
plot(AxisPoints4(1:NAxisPoints2,2),AxisPoints4(1:NAxisPoints2,5),'-.k','LineWidth',1.5)
lgd = legend({'One-dimensional','2D - Wall','2D - Axis'},'Location','best');
title(lgd,'Type of Analysis')
xlabel('x/y_t','fontweight','bold')  
ylabel('M/M_{Exit}','fontweight','bold') 
title('Analysis of Mach Number')
ax = gca; 
outerpos = ax.OuterPosition; 
ti = ax.TightInset;  
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2); 
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4); 
ax.Position = [left bottom ax_width ax_height];
axis([0 FinalChar1(numel(FinalChar1(:,1)),2) 0 1])
grid on
grid minor
hold off

%%

% NOZZLE SHAPE EXPORT

NozzleShape = subsonic(1:numel(subsonic(:,1))-1,1:2);
NozzleShape = [NozzleShape;WallPoints(1:numel(WallPoints(:,1)),1:2)];

dlmwrite('NozzleShape.txt',NozzleShape)
xlswrite('WallPoints.xlsx',WallPoints);
xlswrite('subsonic.xlsx',subsonic);
xlswrite('SuperSonic.xlsx',SuperSonic);