function [N, A_exit, Dy, FinalChar1, MassFlow_FC, Thrust_FC, counter, alpha_D, y_last2] = FinalChar(Y1, X1, gamma, R, y_last, p_st, T_st, u_last, counter, M_axis, MassFlow_IVL)
%% Predictor - In this part of the function it's got N that represent the number of points in which the final char is divided and Dy that represent the first Guess
%Geometrical features
N = ceil(y_last/Y1);
if rem(N,2) == 0
    N = floor(y_last/Y1);
end
Dy = y_last/(N - 1);                           %Delta of vertical distance along the final characteristic
alpha_D = asind(1/M_axis);

y_FC = [];
x_FC = [];
u_FC = [];
v_FC = [];
M_FC = [];
t_FC = [];
p_FC = [];
rho_FC = [];
mass_FC = [];
Thrust_FC = [];

counter3 = 0;

for y = Dy:Dy:y_last
    
    counter3 = counter3 +1;
    
    t_fc = T_st/(1 + (gamma - 1)/2*M_axis^2);                       % Static temperature
    p_fc = p_st/(1 + (gamma - 1)/2*M_axis^2)^(gamma/(gamma - 1));   % Static pressure
    rho_fc = p_fc/(R*t_fc);                                         % Density
    
    %Mass Flow rate calculation along the initial-value line
    m_fc = (Dy/2)*(2*pi*rho_fc*u_last)*(y + (y - Dy));
    
    %Mass Flow rate calculation along the initial-value line
    T_fc = (Dy/2)*(2*pi)*(p_fc + rho_fc*u_last^2)*(y + (y - Dy));
    
    mass_FC(counter3) = m_fc;
    Thrust_FC(counter3) = T_fc;
    
end

MassFlow_FC = sum(mass_FC);
Thrust_FC = sum(Thrust_FC);

%% Corrector to make the continuity princliple works

if MassFlow_FC ~= MassFlow_IVL
        
    num_guesses = 0;
    low = 0;
    high = 1;
    Tol3 = 0.03;
    
    counter3 = 0;
    
    while abs(MassFlow_FC - MassFlow_IVL) > Tol3
        
        Guess = (high + low)/2;
        num_guesses = num_guesses + 1;
        
        for y = Guess:Guess:N*Guess
            
            counter3 = counter3 +1;
            
            t_fc = T_st/(1 + (gamma - 1)/2*M_axis^2);                       % Static temperature
            p_fc = p_st/(1 + (gamma - 1)/2*M_axis^2)^(gamma/(gamma - 1));   % Static pressure
            rho_fc = p_fc/(R*t_fc);                                         % Density
            
            %Mass Flow rate calculation along the initial-value line
            
            m_fc = (Guess/2)*(2*pi*rho_fc*u_last)*(y + (y - Guess));
            
            %Mass Flow rate calculation along the initial-value line
            T_fc = (Guess/2)*(2*pi)*(p_fc + rho_fc*u_last^2)*(y + (y - Guess)); 
            
            mass_FC(counter3) = m_fc;
            Thrust_FC(counter3) = T_fc;
            
        end
        
        MassFlow_FC = sum(mass_FC);
        Thrust_FC = sum(Thrust_FC);
        
        if MassFlow_FC - MassFlow_IVL > Tol3
            high = Guess;
            counter3 = 0;
        elseif MassFlow_FC - MassFlow_IVL < -Tol3
            low = Guess;
            counter3 = 0;
        end
    end
end

Dy = Guess;
y_last2 = N*Dy;

for y = Dy:Dy:y_last2
    
    counter = counter + 1;
    
    x_fc = X1 + y/tand(alpha_D);
    u_fc = u_last;
    v_fc = 0;
    t_fc = T_st/(1 + (gamma - 1)/2*M_axis^2);                       % Static temperature
    p_fc = p_st/(1 + (gamma - 1)/2*M_axis^2)^(gamma/(gamma - 1));   % Static pressure
    rho_fc = p_fc/(R*t_fc);                                         % Density
    
    y_FC = [y_FC; y];
    x_FC = [x_FC; x_fc];
    u_FC = [u_FC; u_fc];
    v_FC = [v_FC; v_fc];
    M_FC = [M_FC; M_axis];
    t_FC = [t_FC; t_fc];
    p_FC = [p_FC; p_fc];
    rho_FC = [rho_FC; rho_fc];
    
end

A_exit = pi*y_last2^2;                          %Exitt area
FinalChar1 = [y_FC, x_FC, u_FC, v_FC, M_FC, t_FC, rho_FC, p_FC];

end