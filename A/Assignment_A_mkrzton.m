

%================ TASK 2 - Slopes of cond() and det() ============================

aN=linspace(0.99,1.01,1000);
for n=2:1000
   
     slope(1,n)= det(generate_matrix(20,aN(n)));
     slope(2,n)= det(generate_matrix(10,aN(n)));
     slope(3,n)= det(generate_matrix(3,aN(n)));    
end

    f1=figure(1);

    semilogy(aN,slope);
    semilogy(aN,slope);
    semilogy(aN,slope);
    
    movegui(f1,'north');
    title (' Dependence of det(A) on \alpha')
    legend({'20x20', '10x10', '3x3'}, 'Location', 'Southeast')
    xlabel('\alpha');
    ylabel('det(A)');

for n=2:1000
   
     slope(1,n)= cond(generate_matrix(20,aN(n)));
     slope(2,n)= cond(generate_matrix(10,aN(n)));
     slope(3,n)= cond(generate_matrix(3,aN(n))); 
end

    f2=figure(2);
    
    semilogy(aN,slope);
    semilogy(aN,slope);
    semilogy(aN,slope);
    
    movegui(f2,'northeast');
    title (' Dependence of cond(A) on \alpha')
    legend({'20x20', '10x10', '3x3'}, 'Location', 'Northeast')
    xlabel('\alpha');
    ylabel('cond(A)');

pause();

%========= Task 3-Generating inverse Matrices for x=k^2/300=============

%3x3 using LU method

for k=0:21
    
    x=2^k/300;
    information = sprintf('3x3 matrix, x=%.2f, k=%d',x,k);
    disp(information)
   
    matrice=generation_for_task_4(3,k);
    
    LU_function=LU_for_task_4(3,k)
    inv_function=inv(matrice)
    
    disp("=======================================================================================================")
end

pause();
%10x10 using LU method

for k=0:21
  
    x=2^k/300;
    information = sprintf('10x10 matrix, x=%.2f, k=%d',x,k);
    disp(information)
    
    matrice=generation_for_task_4(10,k);
    
    LU_function=LU_for_task_4(10,k)
    inv_function=inv(matrice)
    
    disp("=======================================================================================================")
end

pause();

%20x20 using LU method

for k=0:21
    
    x=2^k/300;
    information = sprintf('20x20 matrix, x=%.2f, k=%d',x,k);
    disp(information)
  
    matrice=generation_for_task_4(20,k);
    
    LU_function=LU_for_task_4(20,k)
    inv_function=inv(matrice)
    
    disp("=======================================================================================================")
end

pause();

%3x3 using LLT method

for k=0:21
    
    x=2^k/300;
    information = sprintf('3x3 matrix, x=%.2f, k=%d',x,k);
    disp(information)
    
    matrice=generation_for_task_4(3,k);
    
    LLT_function=LLT_for_task_4(3,k)
    inv_function=inv(matrice)
    
    disp("=======================================================================================================")
end

pause();

%10x10 using LLT method\

for k=0:21
    
    x=2^k/300;
    information = sprintf('10x10 matrix, x=%.2f,  k=%d',x,k);
    disp(information)
   
    matrice=generation_for_task_4(10,k);
    
    LLT_function=LLT_for_task_4(10,k)
    inv_function=inv(matrice)
    
    disp("=======================================================================================================")
end

pause();

%20x20 using LLT method

for k=0:21
    
    x=2^k/300;
    information = sprintf('20x20 matrix, x=%.2f,  k=%d',x,k);
    disp(information)
   
    matrice=generation_for_task_4(20,k);
    
    LLT_function=LLT_for_task_4(20,k)
    inv_function=inv(matrice)
    
    disp("=======================================================================================================")
end

pause();

%3x3 matrix dependence of maximal error

f3 = figure;
figure(f3);

k=linspace(0,21,22);
for x = 1:22
    a(x) = LU_max_error(3,k(x))+eps;
    b(x) = LLT_max_error(3,k(x));
    c(x) = inv_max_error(3,k(x));
     hold on
end
hold off

    loglog(k,a,'r',k,b,'g',k,c,'b');
    
    movegui(f3,'southwest');
    title (' Dependence of maximum error on 3x3 matrice')
    xlabel('k');
    ylabel('Maximum error value');
    legend('LU','LLT','INV', 'Location', 'Northwest');

%10x10 matrix dependence of maximal error

f4 = figure;
figure(f4);

k=linspace(0,21,22);
for x = 1:22
    
    a(x) = LU_max_error(10,k(x));
    b(x) = LLT_max_error(10,k(x));
    c(x) = inv_max_error(10,k(x));
   
end


    loglog(k,a,'r',k,b,'g',k,c,'b');
    
    movegui(f4,'south');
    title (' Dependence of maximum error on 10x10 matrice')
    xlabel('k');
    ylabel('Maximum error value');
    legend('LU','LLT','INV', 'Location', 'Northwest');

%20x20 matrix dependence of maximal error

f5 = figure;
figure(f5);

k=linspace(0,21,22);
for x = 1:22
    
    a(x) = LU_max_error(20,k(x));
    b(x) = LLT_max_error(20,k(x));
    c(x) = inv_max_error(20,k(x));
    
end

    loglog(k,a,'r',k,b,'g',k,c,'b');
    
    movegui(f5,'southeast');
    title (' Dependence of maximum error on 20x20 matrice')
    xlabel('k');
    ylabel('Maximum error value');
    legend('LU','LLT','INV', 'Location', 'Northwest');

%3x3 matrix dependence of RMS value

f6=figure;
figure(f6);

k=linspace(0,21,22);
for x=1:22
    
    a(x)=LU_RMS(3,k(x))+eps;
    b(x)=LLT_RMS(3,k(x));
    c(x)=inv_RMS(3,k(x));
   
end

    loglog(k,a,'r',k,b,'g',k,c,'b');
    
    movegui(f6,'northwest');
    title (' Dependence of RMS error on 3x3 matrice')
    xlabel('k');
    ylabel('RMS error value');
    legend('LU','LLT','INV', 'Location', 'Northwest');

%10x10 matrix dependence of RMS value

f7=figure;
figure(f7);

k=linspace(0,21,22);
for x=1:22
    
    a(x)=LU_RMS(10,k(x));
    b(x)=LLT_RMS(10,k(x));
    c(x)=inv_RMS(10,k(x));
   
end

    loglog(k,a,'r',k,b,'g',k,c,'b');
    
    movegui(f7,'north');
    title (' Dependence of RMS error on 10x10 matrice')
    xlabel('k');
    ylabel('RMS error value');
    legend('LU','LLT','INV', 'Location', 'Northwest');

%20x20 matrix dependence of RMS value 
    
f8=figure;
figure(f8);

k=linspace(0,21,22);
for x=1:22
   
    a(x)=LU_RMS(20,k(x));
    b(x)=LLT_RMS(20,k(x));
    c(x)=inv_RMS(20,k(x));
  
end
   
    loglog(k,a,'r',k,b,'g',k,c,'b');
    
    movegui(f8,'northeast');
    title (' Dependence of RMS error on 20x20 matrice')
    xlabel('k');
    ylabel('RMS error value');
    legend('LU','LLT','INV', 'Location', 'Northwest');
    
pause();

%========Task 5 - comapring values of norm() and LLT and LU functions===== 

A=generation_for_task_4(3,5);
a=norm(A,2);
b=LU_RMS(3,5);
c=LU_max_error(3,5);
d=LLT_RMS(3,5);
e=LLT_max_error(3,5);
information = sprintf('matrix 3x3, k=5 || norm=%d || LU_RMS=%d || LU_ME=%d || LLT_RMS=%d || LLT_M=%d',a,b,c,d,e);
    disp(information)
    disp("=======================================================================================================")

A=generation_for_task_4(10,15);
a=norm(A,2);
b=LU_RMS(10,15);
c=LU_max_error(10,15);
d=LLT_RMS(3,5);
e=LLT_max_error(10,15);
information = sprintf('matrix 10x10, k=15 || norm=%d || LU_RMS=%d || LU_ME=%d || LLT_RMS=%d || LLT_M=%d',a,b,c,d,e);
    disp(information)
    disp("=======================================================================================================")

A=generation_for_task_4(20,21);
a=norm(A,inf);
b=LU_RMS(20,21);
c=LU_max_error(20,21);
LLT_rms=LLT_RMS(20,21);
e=LLT_max_error(20,21);
information = sprintf('matrix 20x20, k=21 || norm=%d || LU_RMS=%d || LU_ME=%d || LLT_RMS=%d || LLT_M=%d',a,b,c,d,e);
    disp(information)

%============= TASK 1- fucntion for generating matrix ==============

function [matrice]=generate_matrix(N,alpha)

matrice = ones(N,N);
x=log(alpha);

for m = 1:N
    for n = 1:N
		
      if ((m==n)&&(m==1)&&(n==1))
                
           matrice(m,:)=matrice(m,:)*x;
           matrice(:,n)=matrice(:,n)*x;
            
           matrice(m,n+1:N)=matrice(m,n+1:N)*(2/3);
           matrice(m+1:N,n)=matrice(m+1:N,n)*(2/3);
      end
		
      if ((m>1)&&(n>1)&&(m==n))
                
           matrice(m,n:N) = matrice(m,n:N)*(((m)*4)/9);
           matrice(m+1:N,n) = matrice(m+1:N,n)*(((m)*4)/9);
      end
   end
end
end

%========= TASK 4-function that generates matrix for x=2^k/300 ========

function [matrice]=generation_for_task_4(N,k)

matrice = ones(N,N);
x=2^k/300;

for m = 1:N
    for n = 1:N
		
      if ((m==n)&&(m==1)&&(n==1))
                
           matrice(m,:)=matrice(m,:)*x;
           matrice(:,n)=matrice(:,n)*x;
            
           matrice(m,n+1:N)=matrice(m,n+1:N)*(2/3);
           matrice(m+1:N,n)=matrice(m+1:N,n)*(2/3);
      end
		
      if ((m>1)&&(n>1)&&(m==n))
                
           matrice(m,n:N) = matrice(m,n:N)*(((m)*4)/9);
           matrice(m+1:N,n) = matrice(m+1:N,n)*(((m)*4)/9);
      end
   end
end
end

%LU - Generate inverse matrix A^-1 performing LU, TASK 4

function [B] = LU_for_task_4(N,k)

A = generation_for_task_4(N,k);

[L,U,P] = lu(A);

A = zeros(N);
B = zeros(N);

for i=1:N
     for j = 1:N
         
         A(j,i)=(P(j,i)-L(j,1:j-1)*A(1:j-1,i))/L(j,j);
     end
end

for i = N:-1:1
     for j = N:-1:1
         
         B(j,i)=(A(j,i)-U(j,j+1:N)*B(j+1:N,i))/U(j,j);
     end
end
end

% LLT - Generate inverse matrix using LLT, TASK 4

function[M]=LLT_for_task_4(N,x)

A=generation_for_task_4(N,x);

ID = eye(N);
K = zeros(N);
M = zeros(N);

L=chol(A,"lower");
LT = L';

for ColID = 1:N
    for RowL = 1:N
        for ColL = 1:RowL
            
          if ColL == RowL
              K(RowL,ColID)=ID(RowL,ColID)/L(RowL,ColL);
              
          elseif ColL < RowL
              ID(RowL,ColID)=ID(RowL,ColID)-(L(RowL,ColL)*K(ColL,ColID));
          end
       end
    end
    
    for RowLT = N:-1:1
        for ColT = N:-1:RowLT
            
            if ColT == RowLT
                M(RowLT,ColID) = K(RowLT,ColID)/LT(RowLT,ColT);
                
            elseif ColT > RowLT
                K(RowLT,ColID)=K(RowLT,ColID)-(LT(RowLT,ColT)*M(ColT,ColID));
                
            end
        end
    end
end
end

% Task 5 - calculation of RMS using LLT method
function [llt_rms] = LLT_RMS(N,k)

    mat=generation_for_task_4(N,k);
    inv_mat=LLT_for_task_4(N,k);
    I=eye(N);
    
    Y=(mat*inv_mat)-I;
    
    llt_rms=max((eig(Y*Y')).^(1/2));


end

%Task 5 - calculation of RMS using LU method
function [lu_rms] = LU_RMS(N,k)

    mat=generation_for_task_4(N,k);
    inv_mat=LU_for_task_4(N,k);
    I=eye(N);
    
    Y=(mat*inv_mat)-I;
    
    lu_rms=max((eig(Y*Y')).^(1/2));


end

%Task 5 - calculation of RMS using inverse function
function [inv_rms] = inv_RMS(N,k)

    mat=generation_for_task_4(N,k);
    inv_mat=inv(mat);
    I=eye(N);
    
    Y=(mat*inv_mat)-I;
    
    inv_rms=max((eig(Y*Y')).^(1/2));


end

%Task 5 - calculation of ME using LLT method
function [llt_me] = LLT_max_error(N,k)

    mat=generation_for_task_4(N,k);
    inv_mat=LLT_for_task_4(N,k);
    I=eye(N);
    
    Y=(mat*inv_mat)-I;

    llt_me=max(sum(abs(Y),2));

end

%Task 5 - calculation of ME using LU method
function [lu_me] = LU_max_error(N,k)

    mat=generation_for_task_4(N,k);
    inv_mat=LU_for_task_4(N,k);
    I=eye(N);
    
    Y=(mat*inv_mat)-I;

    lu_me=max(sum(abs(Y),2));

end

%Task 5 - calculation of RMS using inverse function
function [inv_me] = inv_max_error(N,k)

    mat=generation_for_task_4(N,k);
    inv_mat=inv(mat);
    I=eye(N);
    Y=(mat*inv_mat)-I;

    inv_me=max(sum(abs(Y),2));

end