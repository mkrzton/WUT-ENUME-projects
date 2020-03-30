%----------------TASK 1----------------%
%generating figure of function presented in task 1


x=linspace(-1,1,1000);
y = sqrt(1-x.^2).*exp(x-1/3);

hold on

figure(1);
plot(x,y,'k');
    movegui(1,'northwest');
    title('Dependance of f(x)');
    xlabel('x');
    ylabel('y');
    
%----------------TASK 2----------------%
%creating 3 slopes representing use of approximating function and applying Fi method   

fi = Fi(250,10);
f2 = figure;

figure(f2);
plot(fi);
    movegui(f2,'north');
    title (' Dependence of Fi');
    xlabel('N');
    ylabel('K');

%N=10 and K=9

N=10;
K=9;

x=linspace(-1,1,1000);
y = sqrt(1-x.^2).*exp(x-1/3);
x2 = linspace(-1,1,K);
y2 = sqrt(1-x2.^2).*exp(x2-1/3);
x3 = linspace(-1,1,1000);
y3 = LSAproximation(N,K);

hold on

figure(3);
plot(x,y,'k',x2,y2,'*b',x3,y3,'r');
    movegui(3,'southwest');
    title (' Approximation function for N=10 , K=9')
    legend({'original', 'estimating points', 'approximation'}, 'Location', 'Northwest')
    xlabel('x');
    ylabel('y');
  
%N=20 and K=18

N=20;
K=18;

y = sqrt(1-x.^2).*exp(x-1/3);
x2 = linspace(-1,1,K);
y2 = sqrt(1-x2.^2).*exp(x2-1/3);
y3 = LSAproximation(N,K);

hold on

figure(4);
plot(x,y,'k',x2,y2,'*b',x3,y3,'r');
    movegui(4,'south');
    title (' Approximation function for N=20 , K=18');
    legend({'original', 'estimating points', 'approximation'}, 'Location', 'Northwest');
    xlabel('x');
    ylabel('y');
 
%N=30 and K=27
    
N=30;
K=27;

y = sqrt(1-x.^2).*exp(x-1/3);
x2 = linspace(-1,1,K);
y2 = sqrt(1-x2.^2).*exp(x2-1/3);
y3 = LSAproximation(N,K);

hold on

figure(5);
plot(x,y,'k',x2,y2,'*b',x3,y3,'r');
    movegui(5,'southeast');
    title (' Approximation function for N=30 , K=27')
    legend({'original', 'estimating points', 'approximation'}, 'Location', 'Northwest')
    xlabel('x');
    ylabel('y');
    
pause();
close all;

%----------------TASK 3----------------%
% generating 3D sketches of maximal error and root mean square error

x=linspace(-1,1,1000)';
y1=sqrt(1-x.^2).*exp(x-1/3);

for N=5:50
   for K=4:N-2
       
       y2 = LSAproximation(N,K);
       RMSE(N,K) = (norm((y2-y1))/norm(y1));
       
   end
end

figure(6);
surf(log10(RMSE));

    movegui(6,'northwest');
    title ('3D sketch of root mean square error');
    xlabel('K');
    ylabel('N');
    zlabel('error');

x=linspace(-1,1,1000)';
y1=sqrt(1-x.^2).*exp(x-1/3);

for N=5:50
   for K=4:N-2
       
       y2 = LSAproximation(N,K);
       Maxerror(N,K)=(norm((y2-y1),inf))/norm(y1, inf);
       
   end
end

figure(7);
surf(log10(Maxerror));

    movegui(7,'northeast');
    title ('3D sketch of maximal error');
    xlabel('K');
    ylabel('N');
    zlabel('error');

pause();
close all;

%----------------TASK 4----------------%
%drawing sketch of a function but with \sigma

x = linspace(-1,1,1000);
y = sqrt(1-x.^2).*exp(x-1/3);
x2 = linspace(-1,1,18);
y2 = sqrt(1-x2.^2).*exp(x2-1/3);
y2 = y2+(randn(size(y2))*0.1);
deviation = LSApproximation_sigma(20,18,0.1);
hold on

figure(8);
plot(x,y,'k',x2,y2,'*b',x,deviation,'r');
    movegui(8,'northwest');
    title('Function approximation with standard deviation \sigma');
    legend({'original', 'deviated estimating points', 'approximation'}, 'Location', 'Northwest')
    xlabel('x');
    ylabel('y');

x = linspace(-1,1,1000)';
y1 = sqrt(1-x.^2).*exp(x-1/3);
sigma = logspace(-5,-1,25)';

for a = 1:25
    for N = 5:50
        for K = 4:N-2
       
             y2 = LSApproximation_sigma(N,K,sigma(a));
             error(N,K) = norm(y2-y1)/norm(y1);
             smallest(a) = min(min(error(error>0)));
             
            %sometimes it appears that program does not compile here if so
            %delete (a) from smallest then build after building add (a)
            %that was removed
        end
    end
end

%drawing dependancies of minimal error values of different sigma's

z = log10(sigma);
y = (smallest)';
z1 = linspace(-5,-1,1000);
p = polyfit(z,y,7);
pol_fit = polyval(p,z1);
 
hold on
 
figure(9);
loglog(z,y,'*b',z1,pol_fit,'r');
    movegui(9,'north');
    title('Dependance of minimum RMS error of approximation with \sigma');
    legend({'minima', 'polyfit'}, 'Location', 'Northwest');
    xlabel('\sigma');
    ylabel('RMS error value');
    
 
x = linspace(-1,1,1000)';
y1 = sqrt(1-x.^2).*exp(x-1/3);
sigma = logspace(-5,-1,25)';

for a = 1:25
    for N = 5:50
        for K = 4:N-2
            
             y2 = LSApproximation_sigma(N,K,sigma(a));
             error(N,K) = norm((y2-y1),inf)/norm(y1,inf);
             smallest(a) = min(min(error(error>0)));
       
        end
    end
end

z = log10(sigma);
y = (smallest)';
z1 = linspace(-5,-1,1000);
p = polyfit(z,y,7);
pol_fit = polyval(p,z1);
 
hold on
 
figure(10);
loglog(z,y,'*b',z1,pol_fit,'r');
    movegui(10,'northeast');
    title('Dependance of minimum maximal error of approximation with \sigma');
    legend({'minima', 'polyfit'}, 'Location', 'Southwest');
    xlabel('\sigma');
    ylabel('Maximal Error value');
    

pause();
close all;

% function responsible for generating rectangular matrix with given values
% of rows and columns

function[matrix]=Fi(N,K)

matrix=zeros(N,K);

x=linspace(-1,1,N)';
xk=linspace(-1,1,K)';

    for k=1:K
        
        position = find((abs(x-xk(k)))<1/6);
        matrix(position,k) = cos(3*pi*(x(position)-xk(k))).^2;
    end
end

% function returns vector of 1000 values of approximating function

function[app]=LSAproximation(N,K)

    x = linspace(-1,1,N)';
    y = sqrt(1-x.^2).*exp(x-1/3);
    FI = Fi(N,K);
    FIapp = Fi(1000,K);
    P =(FI'*FI)\(FI'*y);
    app = FIapp*P;
    
end

% function returns vector of 1000 values of approximating function with
% inclusion of standard deviation - sigma

function[app]=LSApproximation_sigma(N,K,sigma)

    x = linspace(-1,1,N)';
    y = sqrt(1-x.^2).*exp(x-1/3);
    sigmay = y+randn(size(y))*sigma;
    FI = Fi(N,K);
    Fiapp = Fi(1000,K);
    P =(FI'*FI)\(FI'*sigmay);
    app = Fiapp*P;
    
end