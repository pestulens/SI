tspan = [0 50];
x0 = [0; 0; 0];
u = [10; 10];
[t,x] = ode45(@(t,x) odefun(t,x,u),tspan,x0);
x = x';
My = [1 -1 -1; -1/2 0 0];
Ky = [0 0; 1/2 0];
y = My*x + Ky*u;
plot (t,y,'green');
grid on;
function dxdt = odefun(~,x,u)
        M = [-1/6 0 -1/3; 0 0 1; 1/2 -1/2 -1/2];
        K = [1/6 1/3; 0 0; 0 0];
        
        dxdt = M*x + K*u;
end