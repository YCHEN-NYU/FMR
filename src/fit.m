%%=========================================================================
% Fit x-y(complex) with Single Lorentzian Model Defined by
% p(1):a1, p(2):w, Full width at Half Maximum (Tesla), p(3):theta, p(4):x0, Hres(Tesla), p(5):a2

% y_real = a1*(0.5*w*cos(theta)+(x-x0)*sin(theta))./((w/2)^2+(x-x0).^2) + c1 + c2*x + c3*x.^2
% y_imag = a2*(0.5*w*cos(theta+pi/2)+(x-x0)*sin(theta+pi/2))./((w/2)^2+(x-x0).^2)+ c4 + c5*x + c6*x.^2

% Fitting Method: least-square curve fitting

% Input Format: x, y = y_real+i*y_imag, frequency(GHz)
% Output Format: FWHM = fitpara(2), Hres(Oe) = fitpara(4)
% 95% confident interval: [FWHM_LowerBound,FWHM_UpperBound] = [min(fitconfint(2,:)),max(fitconfint(2,:))]
%%=========================================================================
function [fitpara, fitconfint, fig] = fit(x,y,frequency)

% set output figure size, position & background color
fig = figure();
set(fig, 'Position', [200, 100, 1000, 600])
set(gcf,'color','w');

% S12 Real
subplot(2,1,1);
plot(x,real(y),'ro','markersize',10);
ylabel('S12 Real','FontSize',36,'FontWeight','bold') 
set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.2g'));

deltax = 0.05*(max(x)- min(x));
xlim([min(x)-deltax,max(x)+deltax]);
deltay_real = 0.05*(max(real(y))- min(real(y)));
ylim([min(real(y))-deltay_real,max(real(y))+deltay_real]);

% S12 Imaginary
subplot(2,1,2);
plot(x,imag(y),'bs','markersize',10);
xlabel('H(T)','FontSize',36,'FontWeight','bold')
ylabel('S12 Img','FontSize',36,'FontWeight','bold') 
set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.2g'));

deltax = 0.05*(max(x)- min(x));
xlim([min(x)-deltax,max(x)+deltax]);
deltay_imag = 0.05*(max(imag(y))- min(imag(y)));
ylim([min(imag(y))-deltay_imag,max(imag(y))+deltay_imag]);


% set a common title on the top of the figure
set(gcf,'NextPlot','add');
axes; 
set(gca,'Visible','off'); 
set(gca,'Fontsize',30,'Linewidth',3);
h = title(['f = ' num2str(frequency) 'GHz'],'fontsize',40,'fontweight','b');
set(h,'Visible','on');

%%=========================================================================
% To select the region for fitting
fprintf('Click on [LEFT, RIGHT] Region!\n')
[xcenter, ~] = ginput(1);
xleft = xcenter - 0.02;
xright = xcenter + 0.02;

fprintf('xleft = %f', xleft);
fprintf('xright = %f', xright);

% swap xleft and xright if xleft > xright
if(xleft > xright)
    xtemp = xright;
    xright = xleft;
    xleft = xtemp;
end
            
% Exclude unselected data 
ind=zeros(length(x),1);
for i =1:1:length(x)
    if (x(i)<xleft || x(i)>xright)
        	ind(i)=i;
    end              
end
            
x=x(setdiff(1:length(x),ind));
y=y(setdiff(1:length(y),ind));
            
% reset the figure and plot with the region just selected. 
clf(fig,'reset');
hold on;

% S12 Real
subplot(2,1,1);
h1 = plot(x,real(y),'ro','markersize',10);

deltax = 0.05*(max(x)- min(x));
xlim([min(x)-deltax,max(x)+deltax]);
deltay_real = 0.05*(max(real(y))- min(real(y)));
ylim([min(real(y))-deltay_real,max(real(y))+deltay_real]);

ylabel('S12 Real','FontSize',36,'FontWeight','bold') 
set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
%set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.1g'));

% S12 Imaginary
subplot(2,1,2);
h2 = plot(x,imag(y),'bs','markersize',10);
xlabel('H(T)','FontSize',36,'FontWeight','bold')
ylabel('S12 Img','FontSize',36,'FontWeight','bold') 
set(gca,'Fontsize',30,'Linewidth',3,'fontweight','bold');

deltax = 0.05*(max(x)- min(x));
xlim([min(x)-deltax,max(x)+deltax]);
deltay_imag = 0.05*(max(imag(y))- min(imag(y)));
ylim([min(imag(y))-deltay_imag,max(imag(y))+deltay_imag]);

set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
%set(gca, 'XTickLabel', num2str(get(gca,'XTick')','%0.1g'));
set(gcf,'NextPlot','add');
axes; 
set(gca,'Visible','off'); 
h = title(['f = ' num2str(frequency) 'GHz'],'fontsize',40,'fontweight','b');
set(h,'Visible','on');

%%
% just checking testx and testy are finite (wich are the plot we want to fit)
% exclude data exclude points from testx,testy through the logic ok_ variable.
% if ok_ is zero we exclude that point
ok_ = isfinite(x) & isfinite(real(y)) & isfinite(imag(y)); %&ok1&ok2&ok3; 
x=x(ok_); 
y=y(ok_);
% ==================================
% store data to be fitted inside xdata (n,1) and ydata (n,2)
xdata = x;
ydata = [real(y),imag(y)];
            %%
%Initialization points construction
%theta = n*pi, symmetric; theta = (n+1/2)*pi, asymmetric
%==================================
% Left Point with FWHM 
fprintf('Click on LEFT!\n');
[xleft,~] = ginput(1);
% find the closest point by xleft position
dleft = abs(xdata - xleft*ones(size(xdata)));
[~,ileft] = min(dleft);
yleft =ydata(ileft);
            
%==================================     
% Left Point with FWHM 
fprintf('Click on RIGHT!\n');
[xright,~] = ginput(1);
% find the closest point by xright position
dright = abs(xdata-xright*ones(size(xdata)));
[~,iright] = min(dright);
yright =ydata(iright);
            
%==================================
% Peak Point with FWHM
fprintf('Click on PEAK!\n');
[xpeak,~] = ginput(1);
% find the closest point by xpeak position
dpeak = abs(xdata-xpeak*ones(size(xdata)));
[~,ipeak] = min(dpeak);
ypeak =ydata(ipeak);
%             
%==================================
% %%% Construct starting parameters (11 of them)
st_a1 = abs(real(ypeak-yleft))+abs(real(yright-ypeak));
st_w = xright-xleft;
st_theta = 0;% start with real symmetric and imaginary asymmetric
st_x0 = xpeak;
st_a2 = abs(imag(ypeak-yleft))+abs(imag(yright-ypeak));
              
st_c10 = mean(ydata(:,1));
st_c11= (ydata(length(x),1)-ydata(1,1))/(xdata(length(x))-xdata(1));
st_c12 = 0;
              
st_c20 = mean(ydata(:,2));
st_c21= (ydata(length(x),2)-ydata(1,2))/(xdata(length(x))-xdata(1));
st_c22 = 0;
st_ = [st_a1, st_w, st_theta,st_x0, st_a2, st_c10,st_c11,st_c12,st_c20,st_c21,st_c22];
%==================================


% set fitting options
options = optimset('Tolx',1e-10,'TolFun',1e-10,'DiffMinChange', 1e-10,'MaxIter',15000,'MaxFunEvals',20000);
% define lower bounds and upper bounds of fitting parameters
x0_lb = [-1,0,-2*pi,0,-1,-1,-1,-1,-1,-1,-1];
x0_ub = [1,0.5,2*pi,5,1,1,1,1,1,1,1];

%==================================
% fit curve by least square non-linear model defined in cplx_fun.m
% p: fitting parameters
[p,~,residuals,~,~,~,jacobian] = lsqcurvefit(@model,st_,xdata,ydata,x0_lb,x0_ub,options);

% ==================================
% Plot out fitting curves for both real and imaginary parts
xmesh = linspace(min(x),max(x),1000);
yfit1 = p(1)*(.5*p(2)*cos(p(3))+(xmesh-p(4))*sin(p(3)))./(.25*p(2)^2+(xmesh-p(4)).^2)+p(6)+p(7)*xmesh+p(8)*xmesh.^2;
yfit2 = p(5)*(.5*p(2)*cos(p(3)+pi/2)+(xmesh-p(4))*sin(p(3)+pi/2))./(.25*p(2)^2+(xmesh-p(4)).^2)++p(9)+p(10)*xmesh+p(11)*xmesh.^2;

Hres_temp=abs(p(4)); % resonant field (T)
w_temp=abs(p(2)); % w_temp is the Full width half maximum(T)
N = 2.5;
subplot(2,1,1);

xlim([Hres_temp-N*w_temp, Hres_temp+N*w_temp]);

h3 = line(xmesh,yfit1,'linewidth',2,'color','r');
lgd1 = legend([h1,h3],'real','real-fit','location','northeast');
set(lgd1,'fontsize',14);
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
line([Hres_temp,Hres_temp], [min(yfit1), max(yfit1)],'Linewidth', 2)
line([Hres_temp - w_temp/2,Hres_temp - w_temp/2], [min(yfit1), max(yfit1)],'linewidth', 2)
line([Hres_temp + w_temp/2,Hres_temp + w_temp/2], [min(yfit1), max(yfit1)],'linewidth', 2)

subplot(2,1,2);
xlim([Hres_temp-N*w_temp, Hres_temp+N*w_temp]);

line([Hres_temp,Hres_temp], [min(yfit2), max(yfit2)],'Linewidth', 2)
line([Hres_temp - w_temp/2,Hres_temp - w_temp/2], [min(yfit2), max(yfit2)],'linewidth', 2)
line([Hres_temp + w_temp/2,Hres_temp + w_temp/2], [min(yfit2), max(yfit2)],'linewidth', 2)


h4 = line(xmesh,yfit2,'linewidth',2,'color','b');
set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%1.3e'));
lgd2 = legend([h2,h4],'imag','imag-fit','location','northeast');
set(lgd2,'FontSize',14);

% set a common title for the figure;
set(gcf,'NextPlot','add');
axes; 
set(gca,'Visible','off'); 
h5 = title(['f = ' num2str(frequency) 'GHz'],'fontsize',40,'fontweight','b');
set(h5,'Visible','on');

% ==================================
% return fitting parameters and 95% confindent intervals
fitpara = p;
% Nonlinear regression parameter confidence intervals
fitconfint = nlparci(p,residuals,'jacobian',jacobian);

end

