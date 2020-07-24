
str=input('A. Solve a System of Equation \nB. Perform a LU decomposition \nC. Perform a Matrix Inversion ', 's');
if(str=='A')
    str1=input('If the matrix is tri-diagonal ,Enter Y ...else N\n','s');
    if(str1=='Y')
        strfile=input('Enter the file name','s');
        fileID = fopen(strfile,'r');
        formatSpec = '%f';
        sizeA = [1 Inf];
        A = fscanf(fileID,formatSpec,sizeA);
        n=A(1,1);
        l=A(2:n);
        d=A(n+1:2*n);
        u=A(2*n+1:3*n-1);
        b=A(3*n:4*n-1);
        a1=zeros(1,n);
        b1=zeros(1,n);
        x=zeros(1,n);
        a1(1)=d(1);
        b1(1)=b(1);
        for i= 2:n
             a1(i)=d(i)-(l(i-1)/a1(i-1))*u(i-1);
             b1(i)=b(i)-(l(i-1)/a1(i-1))*b1(i-1);
        end
        x(n)=b1(n)/a1(n);
        for i=1:n-1
            x(n-i)=(b1(n-i)-u(n-i)*x(n-i+1))/a1(n-i);
        end
        x(1:n);
        y=x.';
        
        fileID = fopen('output.txt', 'w');
        fprintf(fileID,'x\n');
        fprintf(fileID,'%f\n', y);
        fprintf(fileID,'\n');
        fclose(fileID);
                
    end
    if(str1=='N')
        %
        fileID = fopen('nums2.txt','r');
        formatSpec = '%d';
        sizeA = [1 Inf];
        A = fscanf(fileID,formatSpec,sizeA);
        n=A(1,1);
        M=A(2:n+1);
        for i = 1:n-1
            temp=A((n*i+2):n*(i+1)+1);
            M = [M; temp];
        end
        B=A(n*n+2:n*n+n+1);
        C = B.';
        %
        Aug=[M C];

        for i=1:n
        [~,x] = max(Aug(i:end,i));
        pivot=x+i-1;
        temp=Aug(pivot,:);
        Aug(pivot,:)=Aug(i,:);
        Aug(i,:)=temp;
        Aug
            for j=i+1:n
                 Aug(j,:)=Aug(j,:)-Aug(i,:)*Aug(j,i)/Aug(i,i);
            end
        end
        Aug
        X=zeros(1,n);
        sum=0;
        for i=n:-1:1
            
            for j =n:-1:i+1
               sum= sum+Aug(i,j)*X(j);
            end
            X(i)= (Aug(i,n+1)-sum)/Aug(i,i);
            sum=0;
        end
        X
        y=X.';
        fileID = fopen('output.txt', 'w');
        fprintf(fileID,'x\n');
        fprintf(fileID,'%f\n', y);
        fprintf(fileID,'\n');
        fclose(fileID);
        
    end
end
if(str=='B')
    disp('Convergence criterion for Function')
    str1=input('If the matrix is symmetric and positive definite ,Enter Y ...else N\n','s');
    if(str1=='Y')
        strfile=input('Enter the file name','s');
        fileID = fopen(strfile,'r');
        formatSpec = '%f';
        sizeA = [1 Inf];
        A = fscanf(fileID,formatSpec,sizeA);
        n=A(1,1);
        M=A(2:n+1);
        for i = 1:n-1
            temp=A((n*i+2):n*(i+1)+1);
            M = [M; temp];
        end
        M;
        
        fileID = fopen('output.txt', 'w');
        for i=1:n
            max=-Inf;
            for j=i:n
                for k=i:n
                    if(abs(M(j,k))>max)
                        i1=j;
                        j1=k;
                        max=M(j,k);
                    end
                end
            end
            fprintf(fileID,'At iteration ');
            fprintf(fileID,'%d\n', i);
            temp=M(i,:);
            M(i,:)=M(i1,:);
            M(i1,:)=temp;
            
            temp=M(:,i);
            M(:,i)=M(:,j1);
            M(:,j1)=temp;
            for i4=1:n
                for j4=1:n 
                    fprintf(fileID,'%f  ', M(i4,j4));
                end
                fprintf(fileID,'\n');
            end
            fprintf(fileID,'\n\n');
            M
        end
        M;
        L=zeros(n);
        for j=1:n
            for i=j:n
                sum=0;
                if(i==j)
                    sum=M(i,i);
                    for k=1:j-1
                        sum=sum-L(j,k)*L(j,k);
                    end
                    L(i,j)=sqrt(sum);
                end
                if(i~=j)
                    sum=M(i,j);
                    for k=1:j-1
                        sum=sum-L(i,k)*L(j,k);
                    end
                    L(i,j)=sum/L(j,j);
                end
            end
        end
         L
         fprintf(fileID,'L matrix is \n');
         for i2=1:n
                for j2=1:n 
                    fprintf(fileID,'%f  ', L(i2,j2));
                end
                fprintf(fileID,'\n');
            end
         %L*L.'
        
    end
    
    
    if(str1=='N')
        fileID = fopen('nums2.txt','r');
        formatSpec = '%f';
        sizeA = [1 Inf];
        A = fscanf(fileID,formatSpec,sizeA);
        n=A(1,1);
        M=A(2:n+1);
        for i = 1:n-1
            temp=A((n*i+2):n*(i+1)+1);
            M = [M; temp];
        end
        M;
        fileID = fopen('output.txt', 'w');
        for i=1:n
            max=-Inf;
            for j=i:n
                for k=i:n
                    if(M(j,k)>max)
                        i1=j;
                        j1=k;
                        max=M(j,k);
                    end
                end
            end
            fprintf(fileID,'At iteration ');
            fprintf(fileID,'%d\n', i);
            temp=M(i,:);
            M(i,:)=M(i1,:);
            M(i1,:)=temp;
            
            temp=M(:,i);
            M(:,i)=M(:,j1);
            M(:,j1)=temp;
            
            for i4=1:n
                for j4=1:n 
                    fprintf(fileID,'%f  ', M(i4,j4));
                end
                fprintf(fileID,'\n');
            end
            fprintf(fileID,'\n\n');
            
            M;
        end
        M;
        str2=input('Enter 1 for the Doolittle method and 2 for crout method');
        if(str2==1)
            L=zeros(n);
            U=zeros(n);
            for i=1:n
                L(i,i)=1;
            end
            for i=1:n
                for j=i:n
                    sum=M(i,j);
                    for k=1:i-1
                        sum=sum-L(i,k)*U(k,j);
                    end
                    U(i,j)=sum;
                    
                    sum=M(j,i);
                    for k=1:j-1
                        sum=sum-L(j,k)*U(k,i);
                    end
                    L(j,i)=sum/U(i,i);
                    
                end
            end
            fprintf(fileID,'L matrix is \n');
            for i2=1:n
                for j2=1:n 
                    fprintf(fileID,'%f  ', L(i2,j2));
                end
                fprintf(fileID,'\n');
            end
            fprintf(fileID,'\n\n');
            fprintf(fileID,'U matrix is \n');
            for i2=1:n
                for j2=1:n 
                    fprintf(fileID,'%f  ', L(i2,j2));
                end
                fprintf(fileID,'\n');
            end
            L;
            U;
            L*U;
        end
        if(str2==2)
            L=zeros(n);
            U=zeros(n);
            for i=1:n
                U(i,i)=1;
            end
            for i=1:n
                for j=i:n
                    sum=M(j,i);
                    for k=1:i-1
                        sum=sum-L(j,k)*U(k,i);
                    end
                    L(j,i)=sum;
                    
                    sum=M(i,j);
                    for k=1:i-1
                        sum=sum-L(i,k)*U(k,j);
                    end
                    U(i,j)=sum/L(i,i);   
                end
            end
            for i2=1:n
                for j2=1:n 
                    fprintf(fileID,'%f  ', L(i2,j2));
                end
                fprintf(fileID,'\n');
            end
            fprintf(fileID,'\n\n');
            fprintf(fileID,'U matrix is \n');
            for i2=1:n
                for j2=1:n 
                    fprintf(fileID,'%f  ', L(i2,j2));
                end
                fprintf(fileID,'\n');
            end
            L;
            U;
            L*U;
        end
    end
end

%FINDING INVERSE FUNCTION
if(str=='C')
    strfile=input('Enter the file name','s');
    fileID = fopen(strfile,'r');% File Matrix Scanning
    formatSpec = '%d';
    sizeA = [1 Inf];
    A = fscanf(fileID,formatSpec,sizeA);
    n=A(1,1);
    M=A(2:n+1);
    for i = 1:n-1
        temp=A((n*i+2):n*(i+1)+1);
        M = [M; temp];
    end
    I=eye(n);  %Identity Matrix
    Aug=[M I]; %augmented matrix
    for i=1:n  %Operations steps
        Aug(i,:)=Aug(i,:)/Aug(i,i);
            for j=i+1:n
                 Aug(j,:)=Aug(j,:)-Aug(i,:)*Aug(j,i)/Aug(i,i);
            end
    end
    for i=n:-1:1
            for j=i-1:-1:1
                 Aug(j,:)=Aug(j,:)-Aug(i,:)*Aug(j,i);
            end
    end
    
        Aug;
        In = Aug(1:n,n+1:2*n);
        In
        fileID = fopen('output.txt', 'w');
        fprintf(fileID,'Inverse Matrix\n');
        for i=1:n
            for j=1:n 
                fprintf(fileID,'%f ', In(i,j));
            end
            fprintf(fileID,'\n');
        end
        fprintf(fileID,'\n');
        fclose(fileID);
end


