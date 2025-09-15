function [A_t, Delta_y, CriPro, epsilon, alpha, TableVariables, MassFlow_IVL, Thrust_IVL, MassFlow_1D, C_D, Thrust_1D, lambda] = SonicLine(gamma, R, delta, p_st, T_st, y_t, rho_tup, Npoints, a_0)

%Sauer's analysis - Gas Dynamics, Multidimensional flow
%This function calculate:
%a) Location of the sonic line
%b) Location of the v~ = 0
%c) The fluid properties u~ = V, M, p, rho and t
%d) The mass flow rate and the thrust along the line v~ = 0
%e) The one-dimensional sonic-values of the mass flow and the thrust
%Assume the flow to be irrotational, axisymmetric, and the ambient pressure to be zero

%Geometrical features
Delta_y = y_t/(Npoints-1);               %Delta of vertical distance along v~ = 0
A_t = pi*y_t^2;                          %Throat area 

%LOCATION OF THE ORIGIN OF THE COORDINATE SYSTEM IN SAUER'S ANALYSIS

%epsilon - Distance from the origin point of the global coordinate system to the
%origin of the coordinate system of the Sauer's analysis.

epsilon = -y_t/(2*(3 + delta))*(((gamma + 1)*(1 + delta))/(rho_tup/y_t))^(1/2);

%critical properties
a_cri = ((2/(gamma + 1))*a_0^2)^(1/2);
T_cri = (2*T_st)/(gamma + 1);
p_cri = p_st/((gamma + 1)/2)^(gamma/(gamma - 1));
rho_cri = p_cri/(R*T_cri);

CriPro = [a_cri,T_cri,p_cri,rho_cri];

%Coefficient of the linear axial perturbation velocity
alpha = ((1 + delta)/((gamma + 1)*rho_tup*y_t))^(1/2);

% Vectors definition
x_sonic = [];           %Sonic Line in the Sauer's coordinate system     
x_sonic2 = [];          %Sonic Line in the Cartesian coordinate System
x_v = [];               %(v~=0) Line - Initial-Value Line in MoC - Sauer
x_v2 = [];              %(v~=0) Line - Initial-Value Line in MoC - Cartesian
u_p = [];               %Nondimensional perturbation velocity
u_v = [];               %Dimensional velocity
v_v = [];
a_l = [];               %Local acustic speed
Mach = [];              %Mach number
Point = 0;              
Points = [];
y_div = [];
te = [];                %Temperature
pre = [];               %Pressure
rho1 = [];              %Local density
mass = [];              
Thrust = [];

for y = 0:Delta_y:y_t
    
    Point = Point + 1;
    
    %Sonic and Initial-value lines 
    x_1 = -(((gamma + 1)*alpha)/(2*(1 + delta)))*y^2;  %Sonic line
    x_11 = x_1 - epsilon;                              %Origin
    x_2 =-((gamma + 1)*alpha*y^2)/(2*(3 + delta));     %v~ = 0 line
    x_22 = x_2 - epsilon;                              %Origin 2
    
    %Flowfield properties calculation
    u_prima = alpha*x_2 + ((gamma + 1)*alpha^2*y^2)/(2*(delta + 1));
    u_vi = a_cri*(1 + u_prima);
    v_vi = 0;
    a = (a_0^2 - (gamma - 1)/2*u_vi^2)^(1/2);          
    Ma = u_vi/a;                                       % Mach Number
    
    %Thermodynamics properties 
    t = T_st/(1 + (gamma - 1)/2*Ma^2);                 % Static temperature
    p = p_st/(1 + (gamma - 1)/2*Ma^2)^(gamma/(gamma - 1));
    rho = p/(R*t);
       
    %Mass Flow rate calculation along the initial-value line
    if y == 0
        m = 0;  
    else
        m = (Delta_y/2)*(2*pi)*(rho*u_vi*y + rho1(Point - 1)*u_v(Point - 1)*y_div(Point - 1)); 
    end
    
    %Thrust rate calculation along the initial-value line
    if y == 0
        T_f = 0;  
    else
        T_f = (Delta_y/2)*(2*pi)*((p + rho*u_vi^2)*y + (pre(Point - 1) + rho1(Point - 1)*u_v(Point - 1)^2)*y_div(Point - 1)); 
    end
    
    %Results saved in the array and rounded to 3 decimals
    x_sonic = round([x_sonic;x_1],3);
    x_sonic2 = round([x_sonic2;x_11],3);
    x_v = round([x_v;x_2],3);
    x_v2 = round([x_v2;x_22],3);
    u_p = round([u_p;u_prima],3);
    u_v = round([u_v;u_vi],3);
    v_v = [v_v;v_vi];
    Points = [Points;Point];
    y_div = [y_div;y];
    a_l = [a_l;a];
    Mach = [Mach;Ma];
    te = [te;t];
    pre = [pre;p];
    rho1 = [rho1;rho];
    mass = [mass;m];
    Thrust = [Thrust;T_f];
end

TableVariables = table(Points,y_div,x_sonic,x_sonic2,x_v,x_v2,u_v,v_v,a_l,Mach,te,rho1,pre)
MassFlow_IVL = sum(mass)
Thrust_IVL = sum(Thrust)

%One-dimensional Mass Flow in the Throat
MassFlow_1D = rho1(1)*a_l(1)*A_t

%Discharge coefficient
C_D = MassFlow_IVL/MassFlow_1D

%One-dimensional Sonic Thrust in the Throat
Thrust_1D = pre(1)*A_t + MassFlow_1D*a_l(1)

%Thrust efficient or Thrust ratio
lambda = Thrust_IVL/Thrust_1D
end