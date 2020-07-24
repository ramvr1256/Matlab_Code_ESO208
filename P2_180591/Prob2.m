% Name:-Ramveer
% Roll no:-180591

clc
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
        Giverr=A(n*n+2);
str=input('Enter L for the largest eigenvalue or\n Enter A for all eigenvalues', 's');
    if(str=='L')
        b=zeros(1,n);
        b(1)=1;
        y=b.';
        err=Inf;
        for i=1:100
            y=M*y;
            sum=0;
            for j=1:n
                sum=sum+y(j)*y(j);
            end
            Lpnorm=sqrt(sum);
            y;
            y=y/Lpnorm;
            y;
            Lpnorm;
            if(i>1)
                err=abs(100*(Lpnorm-eigenvalue)/Lpnorm);
            end
            eigenvalue=Lpnorm;
            if(err<Giverr)
                iteration=i;
                break
            end
            eigenvalue
            y;
            
            Giverr
        end
        err
        iteration
        eigenvalue;
        
        fileID = fopen('output.txt', 'w');
        fprintf(fileID,'Eigenvalues=\n');
        fprintf(fileID,'%f\n', eigenvalue);
        fprintf(fileID,'\n\n');
        fprintf(fileID,'Eigenvector=\n');
        fprintf(fileID,'%f\n', y);
        fprintf(fileID,'\n\n');
        fprintf(fileID,'Iterations\n');
        fprintf(fileID,'%d\n', iteration);
        fclose(fileID);
        
    end
    
    
    
    
    if(str=='A')
        for i=1:100
            Q=zeros(n);
            R=zeros(n);
            for j=1:n
                x=M(:,j);
                z=x;
            for d=1:j-1
                z=z-((x.')*Q(:,d))*Q(:,d);
            end
            sum=0;
                for k=1:n
                   sum=sum+z(k)*z(k);
                end
            Lpnorm=sqrt(sum);
            Q(:,j)=z/Lpnorm;
            end
            Q;
            for u=1:n
                for v=1:n
                    R(u,v)=(Q(:,u).')*M(:,v);
                end
            end
            R;
            %Giverr=0.000001;
            err=Inf;
            max=0;
            if(i>1)
                for g=1:n
                    temp=abs((R(g,g)-Rold(g,g)));
                    if(temp>max)
                        max=temp;
                    end
                    err=max*100;
                end
            end
            err;
            if(err<Giverr)
                iteration=i;
                break
            end
            Rold=R;
            M=R*Q;
        end
        M;
        iteration
        eigenvalues=zeros(1,n);
        for i=1:n
            eigenvalues(i)=M(i,i);
        end
        Eigenvalues=eigenvalues.';
        iteration;
        
        fileID = fopen('output.txt', 'w');
        fprintf(fileID,'Eigenvalues\n');
        fprintf(fileID,'%f\n', Eigenvalues);
        fprintf(fileID,'\n\n');
        fprintf(fileID,'Iterations\n');
        fprintf(fileID,' =%d\n', iteration);
        fclose(fileID);
    end