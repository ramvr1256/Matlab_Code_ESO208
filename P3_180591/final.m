clc
strfile=input('Enter the file name','s');
file = fopen(strfile,'r');
i=1;
while ~feof(file)
    x(i) = fscanf(file,'%f',1);
    y(i) = fscanf(file,'\t%f',1);
    i=i+1;
end
n=i-1;
fclose(file); 
str2 =input('What do you want to do?\n1.Fit a least square polynomial.\n2.Fit a lagrange interpolation polynomial \n3.Fit a newton interpolation polynomial\n4.Fit cubic splines\n');
if str2==1
    Average = mean(x);
    Variance=0;
    for i=1:n
        Variance = Variance + (y(i)-Average)^2;
    end
    order=input('Enter the Size of polynomial');
    innerprod = zeros(order+1,order+1);
        funprod = zeros(order+1,1);
        for j=1:order+1
            for k=1:order+1
                kkkk=0;
                for i=1: n
                    kkkk = kkkk + x(i)^(j+k-2);
                end
                innerprod(j,k) = kkkk;
            end
            kkkk=0;
            for i=1: n
                kkkk = kkkk + y(i)* (x(i)^(j-1));
            end
            funprod(j) = kkkk;
        end
        
        coeff = fliplr(transpose(innerprod\funprod));  %A* transpose(C) = F
        p = poly2sym(coeff);
        for i=1:n
            yprime(i) = polyval(coeff,x(i));
        end
        epsilon_sqaure=0;
        for i=1:n
            epsilon_sqaure = epsilon_sqaure + (y(i)-yprime(i))^2;
        end
        rSquared = (Variance - epsilon_sqaure)/Variance;
        file = fopen('output.txt','w');
        fprintf(file,'The equation of the polynomial is:\n');
        fprintf(file,'%s',poly2sym(coeff));
        fprintf(file,'\nCoefficient of regression is: ');
        fprintf(file,'%f',rSquared);
        fclose(file);
        figure
        scatter(x,y);
        hold on
        fplot(p,[x(1) x(n)]);
        hold off
    
    
end
if str2==2
    sum=0;
    for i=1:n
        p=1;
        for j=1:n
            if j~=i
                cof1 = poly(x(j))/(x(i)-x(j)); % coeff of polynomial whose roots are x(j)
                p = conv(p,cof1); % coeff of polynomial which is multiplication of polynomials with coefficients p and c
            end
        end
        term = p*y(i);
        sum= sum + term;
    end
    polynomial= poly2sym(sum);
    file = fopen('output.txt','w');
    fprintf(file,'polynomial is:\n');
    fprintf(file,'%s',polynomial);
    fclose(file);
    figure
    scatter(x,y);
    hold on
    fplot(polynomial,[x(1) x(n)]);
    hold off
    
end
if str2==3
    A=zeros(1,n);
    B=A.';
    D=x(1:n);
    D=[D;y(1:n)];
    New=D.';
    New=[New B];
    n;
    New(4,2);
    
    for j=3:n+2
        for i=1:n+2-j
            New(i,j)=(New(i+1,j-1)-New(i,j-1))/(x(i+j-2)-x(i));
        end
    end
    New;
    coeff=zeros(n,n+1);
    Basis=zeros(n,n+1);
    
    Basis(1,n+1)=x(1)*-1;
    Basis(1,n)=1;
    i=2;
    hg=ones(1,n);
    gh=hg.';
    New(:,1);
    Newf=New(:,1)*-1;
    gh=[gh Newf(:,1)];
    
    for k=1:n-2
        Basis(k+1,n-k:n+1)=conv(gh(k+1,:),Basis(k,n-k+1:n+1));
    end
    Basis;
    for k=1:n-1
        coeff(k,:)=New(1,k+2)*Basis(k,:);
    end
    coeff(n,n+1)=y(1);
    coeff
    for k=1:n
        coff(k)=sum(coeff(:,k));
    end
    coff
    p=coff;
    polynomial= poly2sym(p);
    file = fopen('output.txt','w');
    fprintf(file,'polynomial is:\n');
    fprintf(file,'%s',polynomial);
    fclose(file);
    figure
    scatter(x,y);
    hold on
    fplot(polynomial,[x(1) x(n)]);
    hold off
    
    
end
if str2==4
    
    h = zeros(n-1,1);
    divDiff = zeros(n-1,1);
    for i=1:n-1
        h(i) = x(i+1)-x(i);
        divDiff(i) = (y(i+1) - y(i))/(x(i+1)- x(i));
    end
    H = zeros(n-2,n);
    A = zeros(n-2);
    bigDivDMatrix = zeros(n-2,1);
    soln = [];
    for i=1:n-2
        H(i,i) = h(i);
        H(i,i+1) = 2 * (h(i+1) + h(i));
        H(i,i+2) = h(i+1);
        bigDivDMatrix(i) = 6* divDiff(i+1) - 6* divDiff(i);
    end
    str3=input('1. Linear spline\n2. Quadratic spline\n3.NAtural spline \n4. Not-a knot\n5. Periodic\n6.Clamped spline')
    if str3==1
        %Linear Spline
        file = fopen('output.txt','w') ;
        
        for i = 1:n-1
            
            cof1 = -divDiff(i);
            cof0 = y(i)-divDiff(i)*x(i);
            
            fprintf(file,'The Interval for Interpolation %f to %f \n',x(i),x(i+1)) ;
            
            
            fprintf(file,'a1 = %f\n',cof1);
            fprintf(file,'a0 = %f\n',cof0);
            for i = 1:n-1
                y_val = zeros(1,(x(i+1)-x(i))*100 + 1) ; %to compensate the number of x values
                k = 1 ;
                for x_val = x(i):0.01:x(i+1)
                    
                    f=(divDiff(i)*(x_val-x(i)))+y(i);
                    y_val(k) = f  ;
                    k = k+1 ;
                end
                x_val = x(i):0.01:x(i+1) ; % interval length 0.01
                plot(x_val,y_val,'b') ;
                hold on
            end
            scatter(x,y);  %plot x(i)
            hold off
        end
        hold off
        fclose(file) ;
        
    end
    
    if str3 == 2
        u(1)=1;
        file = fopen('output.txt', 'w');
        for i =2:n
            u(i)=(2*(y(i)-y(i-1))/(x(i)-x(i-1)))-u(i-1);
        end
        for i = 1:n-1
            q(i,:)=(u(i+1)/(2*(x(i+1)-x(i))))*poly([x(i) x(i)]) - (u(i)/(2*(x(i+1)-x(i))))*poly([x(i+1) x(i+1)]) + [0 0 (y(i)+(u(i)*(x(i+1)-x(i)))/2)];
            
            x2 = linspace(x(i),x(i+1));
            plot(x2,polyval(q(i,:),x2),'-')
            xlabel('x')
            ylabel('y')
            hold on;
            plot (x,y,'o')
            fprintf(file,'The Interval for Interpolation %f to %f \n',x(i),x(i+1)) ;
            
            
            fprintf(file,'a2 = %f\n',q(i,1));
            fprintf(file,'a1 = %f\n',q(i,2));
            fprintf(file,'a0 = %f\n',q(i,3));
            fprintf(file,'The Value of the first derivative at first node is %4.4f and at second node is %4.4f \n',2*q(i,1)*x(i) + q(i,2),2*q(i,1)*x(i+1) + q(i,2))  ;
            
            fprintf(file, '\n')
        end
        q;
        disp(q(1,:))
    end
    
    if str3 == 3
        % natural spline
        %v(0) = v(n) = 0
        A(:,1) = H(:,1);
        A(:,n-2) = H(:,n-1);
        A(:,2:n-3) = H(:,3:n-2);
        soln= A\bigDivDMatrix;
        soln= [0;soln;0];
        
        file = fopen('output.txt','w') ;
        
        for i = 1:n-1
            cof3 = (soln(i+1)-soln(i))/(6*h(i)) ;
            cof2 = (soln(i+1)*x(i)-soln(i)*x(i+1)) /(2*h(i)) ;
            cof1 = soln(i)*(h(i)/6 - x(i+1)^2/(2*h(i))) + soln(i+1)*(-h(i)/6 + x(i)^2/(2*h(i))) + divDiff(i) ;
            cof0 =  soln(i+1)*(-x(i)^3/h(i) + h(i)*x(i))/6 + soln(i)*(x(i+1)^3/h(i) - h(i)*x(i+1))/6 + (y(i)*x(i+1)-y(i+1)*x(i))/h(i) ;
            
            fprintf(file,'The Interval for Interpolation %f to %f \n',x(i),x(i+1)) ;
            
            fprintf(file,'a3 = %f\n',cof3);
            fprintf(file,'a2 = %f\n',cof2);
            fprintf(file,'a1 = %f\n',cof1);
            fprintf(file,'a0 = %f\n',cof0);
            fprintf(file,'The Value of the first derivative at first node is %4.4f and at second node is %4.4f \n',3*cof3*x(i)^2 + 2*cof2*x(i) + cof1,3*cof3*x(i+1)^2 + 2*cof2*x(i+1) + cof1)  ;
            fprintf(file,'The Value of the second derivative at first node is %4.4f and at second node is %4.4f \n\n',soln(i),soln(i+1)) ;
            %         curve = poly2sym([coeff3 coeff2 coeff1 coeff0]);
            %         fplot(curve,[x(1) x(n)]);
            %         hold on
        end
        hold off
        fclose(file) ;
        
    elseif str3 == 4
        % Not-a-knot
        % (v(1) - v(0))/h(0) = (v(2) - v(1))/h(1) & (v(n-1) -
        % v(n-2))/h(n-2) = (v(n) - v(n-1))/h(n-1)
        A(:,1) = ((h(2)+h(1))/h(2))*H(:,1) + H(:,2) ;
        A(:,2) = -1*(h(1)/h(2))*H(:,1) + H(:,3) ;
        A(:,n-2) = H(:,n-1) + ((h(n-1)+h(n-2))/h(n-2))*H(:,n) ;
        A(:,n-3) = H(:,n-2) -(h(n-1)/h(n-2))*H(:,n) ;
        A(:,3:n-4) = H(:,4:n-3) ;
        
        soln = A\bigDivDMatrix ;
        v0 = soln(1) - (soln(2)-soln(1))*(h(1)/h(2)) ;
        vn = soln(n-2) + (soln(n-2)-soln(n-3))*(h(n-1)/h(n-2)) ;
        
        soln = [v0;soln;vn] ;
        
        file = fopen('output.txt','w') ;
        
        for i = 1:n-1
            cof3 = (soln(i+1)-soln(i))/(6*h(i)) ;
            cof2 = (soln(i+1)*x(i)-soln(i)*x(i+1)) /(2*h(i)) ;
            cof1 = soln(i)*(h(i)/6 - x(i+1)^2/(2*h(i))) + soln(i+1)*(-h(i)/6 + x(i)^2/(2*h(i))) + divDiff(i) ;
            cof0 =  soln(i+1)*(-x(i)^3/h(i) + h(i)*x(i))/6 + soln(i)*(x(i+1)^3/h(i) - h(i)*x(i+1))/6 + (y(i)*x(i+1)-y(i+1)*x(i))/h(i) ;
            
            fprintf(file,'The Interval for Interpolation %f to %f \n',x(i),x(i+1)) ;
            
            fprintf(file,'a3 = %f\n',cof3);
            fprintf(file,'a2 = %f\n',cof2);
            fprintf(file,'a1 = %f\n',cof1);
            fprintf(file,'a0 = %f\n',cof0);
            fprintf(file,'The Value of the first derivative at first node is %4.4f and at second node is %4.4f \n',3*cof3*x(i)^2 + 2*cof2*x(i) + cof1,3*cof3*x(i+1)^2 + 2*cof2*x(i+1) + cof1)  ;
            fprintf(file,'The Value of the second derivative at first node is %4.4f and at second node is %4.4f \n\n',soln(i),soln(i+1)) ;
            %         curve = poly2sym([coeff3 coeff2 coeff1 coeff0]);
            %         fplot(curve,[x(1) x(n)]);
            %         hold on
        end
        hold off
        fclose(file) ;
    elseif str3 == 5
        % Periodic
        %         v(1) = v(n-1);
        %         v(2) = v(n);
        A(:,1) = H(:,2) + H(:,n) ;
        A(:,n-2) = H(:,n-1) + H(:,1) ;
        A(:,2:n-3) = H(:,3:n-2) ;
        
        soln = A\bigDivDMatrix ;
        v0 = soln(n-2);
        vn = soln(1);
        
        soln = [v0;soln;vn] ;
        
        file = fopen('output.txt','w') ;
        
        for i = 1:n-1
            cof3 = (soln(i+1)-soln(i))/(6*h(i)) ;
            cof2 = (soln(i+1)*x(i)-soln(i)*x(i+1)) /(2*h(i)) ;
            cof1 = soln(i)*(h(i)/6 - x(i+1)^2/(2*h(i))) + soln(i+1)*(-h(i)/6 + x(i)^2/(2*h(i))) + divDiff(i) ;
            cof0 =  soln(i+1)*(-x(i)^3/h(i) + h(i)*x(i))/6 + soln(i)*(x(i+1)^3/h(i) - h(i)*x(i+1))/6 + (y(i)*x(i+1)-y(i+1)*x(i))/h(i) ;
            
            fprintf(file,'The Interval for Interpolation %f to %f \n',x(i),x(i+1)) ;
            
            fprintf(file,'a3 = %f\n',cof3);
            fprintf(file,'a2 = %f\n',cof2);
            fprintf(file,'a1 = %f\n',cof1);
            fprintf(file,'a0 = %f\n',cof0);
            fprintf(file,'The Value of the first derivative at first node is %4.4f and at second node is %4.4f \n',3*cof3*x(i)^2 + 2*cof2*x(i) + cof1,3*cof3*x(i+1)^2 + 2*cof2*x(i+1) + cof1)  ;
            fprintf(file,'The Value of the second derivative at first node is %4.4f and at second node is %4.4f \n\n',soln(i),soln(i+1)) ;
            %         curve = poly2sym([coeff3 coeff2 coeff1 coeff0]);
            %         fplot(curve,[x(1) x(n)]);
            %         hold on
        end
        hold off
        fclose(file) ;
    elseif str3 == 6
        % Clamped spline
        %u(0) = alpha & u(n) = beta
        u0 = input('1st derivative value at the starting point\n');
        un = input('1st derivative value at the end point\n');
        col2 = bigDivDMatrix - ((3*(divDiff(1)-u0))/h(1))*H(:,1) ;
        A(:,1) = H(:,2) - (H(:,1)/2) ;
        A(:,n-2) = H(:,n-1) - (H(:,n)/2) ;
        A(:,2:n-3) = H(:,3:n-2) ;
        col2 = col2 -((3*(un - divDiff(n-1)))/h(n-1))*H(:,n) ;
        
        soln = A\col2 ;
        
        v0 = (((6*(divDiff(1)-u0))/h(1))-soln(1))/2 ;
        vn = (((6*(un - divDiff(n-1)))/h(n-1))-soln(n-2))/2 ;
        
        soln = [v0;soln;vn] ;
        
        file = fopen('output.txt','w') ;
        
        for i = 1:n-1
            cof3 = (soln(i+1)-soln(i))/(6*h(i)) ;
            cof2 = (soln(i+1)*x(i)-soln(i)*x(i+1)) /(2*h(i)) ;
            cof1 = soln(i)*(h(i)/6 - x(i+1)^2/(2*h(i))) + soln(i+1)*(-h(i)/6 + x(i)^2/(2*h(i))) + divDiff(i) ;
            cof0 =  soln(i+1)*(-x(i)^3/h(i) + h(i)*x(i))/6 + soln(i)*(x(i+1)^3/h(i) - h(i)*x(i+1))/6 + (y(i)*x(i+1)-y(i+1)*x(i))/h(i) ;
            
            fprintf(file,'The Interval for Interpolation %f to %f \n',x(i),x(i+1)) ;
            
            fprintf(file,'a3 = %f\n',cof3);
            fprintf(file,'a2 = %f\n',cof2);
            fprintf(file,'a1 = %f\n',cof1);
            fprintf(file,'a0 = %f\n',cof0);
            fprintf(file,'The Value of the first derivative at first node is %4.4f and at second node is %4.4f \n',3*cof3*x(i)^2 + 2*cof2*x(i) + cof1,3*cof3*x(i+1)^2 + 2*cof2*x(i+1) + cof1)  ;
            fprintf(file,'The Value of the second derivative at first node is %4.4f and at second node is %4.4f \n\n',soln(i),soln(i+1)) ;
            %         curve = poly2sym([coeff3 coeff2 coeff1 coeff0]);
            %         fplot(curve,[x(1) x(n)]);
            %         hold on
        end
        hold off
        fclose(file) ;
    end
    if str3>2
    for i = 1:n-1
        y_val = zeros(1,(x(i+1)-x(i))*100 + 1) ; %to compensate the number of x values
        k = 1 ;
        for x_val = x(i):0.01:x(i+1)
            f = (soln(i+1)/6)*(( power(x_val-x  (i),3))/h(i) - h(i)*(x_val-x(i))) + (soln(i)/6)*( (-power(x_val-x(i+1),3))/h(i) + h(i)*(x_val-x(i+1))) + (y(i+1)/h(i))*(x_val-x(i)) - (y(i)/h(i))*(x_val-x(i+1)) ;
            y_val(k) = f  ;
            k = k+1 ;
        end
        x_val = x(i):0.01:x(i+1) ; % interval length 0.01
        plot(x_val,y_val,'b') ;
        hold on
    end
    scatter(x,y);  %plot x(i)
    hold off
    end
    end
    
    
