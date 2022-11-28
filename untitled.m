
num_of_cycles = 20;
amp = 1; %source amplitude
freq = 2; %source frequency
w_init = 0.5; %the initial state condition [0:1] 

P_coeff = 2;
J = 1;

beta = 9;
a = 4;
c = 0.01;
n = 14;
q = 13;
g = 4;
alpha = 7;

points = 40000;
tspan = [0 num_of_cycles/freq];
t = linspace(tspan(1),tspan(2),points);
delta_t = t(2) - t(1);
V = amp*sin(freq*2*pi*t);
x = (1:points)/10000;
x_c = 0.01;
W = zeros(size((t)));
Wp = zeros(size((t)));
W_dot = zeros(size((t)));
curr = zeros(size((t)));
U = zeros(size((t)));
Up = ones(size((t)));
W(1) = w_init;

for i=2:points
    %Biolek window
    W_dot(i) = a*V(i)^q;
    W(i) = W(i-1)+W_dot(i)*delta_t*(1-(W(i-1)-heaviside(-V(i-1)))^(2*P_coeff));

    W_dot(i) = a*V(i)^q;
    Wp(i) = Wp(i-1)+W_dot(i)*delta_t*(1-(Wp(i-1)-heaviside(-V(i-1)))^(2*P_coeff));
    
    % correct the w vector according to bounds [0 D]
    if W(i) < 0
        W(i) = 0;
        W_dot(i) = 0;
    elseif W(i) > 1
        W(i) = 1;
        W_dot(i) = 0;
    end
    
  curr(i) = W(i)^n*beta*sinh(alpha*V(i))+c*(exp(g*V(i))-1);
  U(i) = -W(i)*(exp(-(x(i)+x_c)^2)/V(i)^2*curr(i)^2)+exp((-(x(i)-x_c)^2)/V(i)^2*curr(i)^2);
  Up(i) = W(i)/2*sin(2*pi*x(i)/(V(i)/curr(i)));
  
  
    
    
end


figure(1);
plot(V(30e3:end),curr(30e3:end));
title('I-V curve');
xlabel('V[volt]');
ylabel('I[amp]');

figure(2);
plot(x(30e3:end),U(30e3:end)+Up(30e3:end),'r');
title('U-x');
xlabel('Position[x]');
legend('U');

