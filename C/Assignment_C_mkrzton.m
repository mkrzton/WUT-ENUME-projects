%----------------TASK 1----------------%
 h=0.01;
 t=0:h:10;
 y0=[0;2];
 A=[0,1;(-10/9),(-2/3)];
 
 y=ODE(h,y0);
 ode=y(:,1);
 y=Gauss_Legendre(A,y0,h);
 gauss=y(1,:);
 y=Euler(A,y0,h);
 euler=y(1,:);
 
 figure(14)
         plot(t,ode,'r',t,gauss,'k');
         movegui(14,'northwest');
         legend({'ode45', 'Gauss-Legendre'}, 'Location', 'Southeast');
         title ('Dependence of ode45 and Gauss-Legendre method');
         xlabel('t');
         ylabel('y');
 
 figure(15)
         plot(t,ode,'r',t,euler,'g');
         movegui(15,'north');
         legend({'ode45', 'Euler'}, 'Location', 'Southeast');
         title ('Dependence of ode45 and Euler method');
         xlabel('t');
         ylabel('y');

 figure(16)
         plot(t,gauss,'k',t,euler,'g');
         movegui(16,'northeast');
         legend({'Gauss-Legendre', 'Euler'}, 'Location', 'Southeast');
         title ('Dependence of Gauss-Legendre and Euler method');
         xlabel('t');
         ylabel('y');
 pause();
 
%----------------TASK 2 & 3----------------%
 
 h=logspace(-3,0,10)';
 h=[h;2;3]; 
 
 %Generating slopes of all methods for different integration step values
 
 for i=1:numel(h)
     
     t=0:h(i):10;
     
     y=ODE(h(i),y0);
     ode=y(:,1);
     y=Gauss_Legendre(A,y0,h(i));
     gauss=y(1,:);
     y=Euler(A,y0,h(i));
     euler=y(1,:);
 
    figure(i);
            plot(t,ode,'r',t,gauss,'k',t,euler,'g');
            legend('ode45','Gauss-Legendre','Euler', 'Location', 'Southeast');
            title('Dependence of Gauss-Legendre, ode45, Euler methods on h');
            xlabel('t');
            ylabel('y');
 end
 
 pause();

 h=logspace(-4,0,30)';
 h=[h;2;3;4]; 
  
 %calculating values of RMS and Max error for different integration step
 %values in order to create dependence of those errors
 
  for i=1:numel(h)
      
      y=ODE(h(i),y0);
      ode=y(:,1);
      y=Gauss_Legendre(A,y0,h(i));
      gauss=y(1,:);
      y=Euler(A,y0,h(i));
      euler=y(1,:);
      
      gauss_RMS(i)=(norm((gauss'-ode)))/(norm(ode));
      euler_RMS(i)=(norm((euler'-ode)))/(norm(ode));
      
      gauss_Max(i)=(norm((gauss'-ode),inf))/(norm(ode,inf));
      euler_Max(i)=(norm((euler'-ode),inf))/(norm(ode,inf));
 end
 
 figure(4)
         loglog(h,euler_RMS,'g',h,gauss_RMS,'k');
         movegui(4,'northwest');
         legend({'Euler', 'Gauss-Legendre'}, 'Location', 'Southeast');
         title ('Dependence of RMS error on h');
         xlabel('h');
         ylabel('RMS Error');
 
 figure(5)
         loglog(h,euler_Max,'g',h,gauss_Max,'k');
         movegui(5,'north');
         legend({'Euler', 'Gauss-Legendre'}, 'Location', 'Southeast');
         title ('Dependence of maximal error on h');
         xlabel('h');
         ylabel('Maximal Error');

pause();
close all


function diff = differential_equation(t,y)

    diff = [y(2);(-6*y(2)-10*y(1))/9];  
end

% function which solves my Differential equation using ode45

function y = ODE(h,y0)

    t = 0:h:10;
    step_size = odeset('RelTol',2.22045e-14,'AbsTol',1e-15); 
    [t,y] = ode45(@differential_equation,t,y0,step_size); 
end

% function which solves differential equation using implicit Gauss-Legendre order 5 method 
function y = Gauss_Legendre(A,y0,h)

t=0:h:10;
y=zeros(2,numel(t));
y(:,1)=y0;
I=eye(2);

     Left=[I-(5/36)*h*A, (-2/9+sqrt(15)/15)*h*A, (-5/36+sqrt(15)/30)*h*A;
          (-5/36-sqrt(15)/24)*h*A,I-2/9*h*A,(-5/36+sqrt(15)/24)*h*A;
          (-5/36+sqrt(15)/30)*h*A, (-2/9+sqrt(15)/15)*h*A, I-5/36*h*A];
     
    for n=2:numel(t)
        
        Right=zeros(6,1);
        Right(1:2)=A*y(:,n-1);
        Right(3:4)=A*y(:,n-1);
        Right(5:6)=A*y(:,n-1);
        
        F=Left\Right;
        f1=F(1:2);
        f2=F(3:4);
        f3=F(5:6);
        y(:,n)=(y(:,n-1)+h*((-5/6)*f1+(8/3)*f2+(-5/6)*f3));  
    end
end

% function which solves differential equation using Euler method
 function y = Euler(A,y0,h)

t=0:h:10;   
y=zeros(2,numel(t));
y(:,1)=y0;

    for n=2:numel(t)
        
        y(:,n)=y(:,n-1)+h*A*y(:,n-1);
    end
end







