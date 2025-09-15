function[Ma, t, p, rho] = TERMO(Vg_4, gamma, R, p_st, T_st)

%This function calculate the termodynamics properties based on the Velocity
%for a perfect gas. It was done by a independent routine in case of the
%need to be changed the function.

Cp = R*gamma/(gamma - 1);
t = T_st - Vg_4^2/(2*Cp);                            % Static temperature
C = (gamma*R*t)^(1/2);                               % Velocity of sound for t
Ma = Vg_4/C;                                         % Mach Number
p = p_st/(T_st/t)^(gamma/(gamma - 1));               % Static pressure
rho = p/(R*t);                                       % Density

end