%
%  by Alex Kaiser
%

dx = 1e-4; 

x = -3:dx:3; 
y = zeros(size(x)); 

N = length(x); 

for i = 1:N
    y(i) = delta_5_smooth(x(i)); 
end 

fig = figure; 
plot(x,y)
title('Unscaled delta function')

% printfig(fig, 'delta_5_smooth'); 

% abuse zero pads at end 

minus_1  = [1,1:N-1]; 
minus_2  = [1,1,1:N-2]; 
plus_1 = [2:N,N]; 
plus_2 = [3:N,N,N]; 

deriv = (y(plus_1) - y(minus_1)) / (2 * dx); 


figure(2); 
plot(x, deriv,'r'); hold on;
second_deriv = (y(plus_1) - 2*y + y(minus_1))/(dx^2); 
plot(x, second_deriv,'g--'); 
third_deriv = (y(plus_2) - 2*y(plus_1) + 2*y(minus_1) - y(minus_2))/(2*dx^3); 
plot(x, third_deriv,'b:')
l2 = legend('\phi''(r)', '\phi"(r)', '\phi^{(3)}(r)');
set(l2, 'FontSize',FS)
xlabel('$r$', 'Interpreter','LaTeX')
box on;
tightfig
print('-dpng', 'new5pt_derivatives.png','-r300')
print('-depsc2','new5pt_derivatives.eps')



x = -.5:dx:.5; 
KK = (38 - sqrt(69))/60; 
f = @(r) 3123 - 6840*KK + 3600*KK.^2 - 12440*r.^2 + 25680*KK*r.^2 - 12600*KK.^2*r.^2 + 8080*r.^4 - 8400*KK*r.^4 - 1400*r.^6; 

fig = figure; 
plot(x,f(x)); 
title('Term in root, K = (38 - sqrt(69))/60')
% printfig(fig, 'term_in_root_works'); 

x = -.5:dx:.5; 
KK = (38 + sqrt(69))/60; 
f = @(r) 3123 - 6840*KK + 3600*KK.^2 - 12440*r.^2 + 25680*KK*r.^2 - 12600*KK.^2*r.^2 + 8080*r.^4 - 8400*KK*r.^4 - 1400*r.^6; 

fig = figure; 
plot(x,f(x)); 
title('Term in root, K = (38 + sqrt(69))/60')
% printfig(fig, 'term_in_root_bad'); 




