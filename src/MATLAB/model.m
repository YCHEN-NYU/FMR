%% ========================================================================

% single peak Lorentz profile (Derivative of Lorentzian)

% p(1):a1, p(2):w, Full width at Half Maximum (Tesla), p(3):theta, p(4):x0, Hres(Tesla), p(5):a2

% y_real = a1*(0.5*w*cos(theta)+(x-x0)*sin(theta))./((w/2)^2+(x-x0).^2) + c1 + c2*x + c3*x.^2
% y_imag = a2*(0.5*w*cos(theta+pi/2)+(x-x0)*sin(theta+pi/2))./((w/2)^2+(x-x0).^2)+ c4 + c5*x + c6*x.^2

%% ========================================================================

function y = model(p,x)

y = zeros(length(x),2); % allocate yout
% p(1):a1, p(2):w, Full width at Half Maximum (Tesla), p(3):theta, p(4):x0, Hres(Tesla), p(5):a2

y(:,1) = p(1)*(0.5*p(2)*cos(p(3))+(x-p(4))*sin(p(3)))./(0.25*p(2)^2+(x-p(4)).^2)+p(6)+p(7)*x+p(8)*x.^2;

y(:,2) = p(5)*(0.5*p(2)*cos(p(3)+pi/2)+(x-p(4))*sin(p(3)+pi/2))./(0.25*p(2)^2+(x-p(4)).^2)+p(9)+p(10)*x+p(11)*x.^2;

end