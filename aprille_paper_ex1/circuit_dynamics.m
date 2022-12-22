function dxdt = circuit_dynamics(t, x)
    C1 = 1e-6;
    C2 = 1e-3;
    C3 = 1e-3;
    R1 = 5;
    L1 = 0.1;
    RL = 1e3;
    Vamp = 10;
    freq = 60;
    Isat = 1e-6;
    kD = 40;
    id = Isat*(exp(kD*x(1))-1);
    Vin = Vamp*sin(2*pi*freq*t);
    dxdt = zeros(4,1);

%     dxdt(1) = 1/C1*(1/R1*(-x(1)-x(2)+Vin) - id);
%     dxdt(2) = 1/C2*(1/R1*(-x(1)-x(2)+Vin) - x(3));
%     dxdt(3) = 1/L1*(x(2)-x(4));
%     dxdt(4) = 1/C3*(x(3)-x(4)/RL);
    dxdt(1) = 1/1e-6*(1/5*(-x(1)-x(2)+10*sin(120*pi*t)) - 1e-6*(exp(40*x(1))-1));
    dxdt(2) = 1/1e-3*(1/5*(-x(1)-x(2)+10*sin(120*pi*t)) - x(3));
    dxdt(3) = 1/0.1*(x(2)-x(4));
    dxdt(4) = 1/1e-3*(x(3)-x(4)/1000);
end