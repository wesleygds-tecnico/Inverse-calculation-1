function [] = CinvEGA(A)
D = A;
n = length(A);
X = zeros(n);
if (det(A)~=0)
    for k1 = 1:n
        A = D;
        B = zeros(n,1);
        B(k1) = 1;
        for k = 1:(n-1)
            for i = (n-1):-1:k
                if(A(i+1,k)>A(i,k))
                    C = A(i,:);
                    A(i,:) = A(i+1,:);
                    A(i+1,:) = C;
                    C = B(i);
                    B(i) = B(i+1);
                    B(i+1) = C;
                end
            end
        end
        for j = 1:(n-1)
            for i = (j+1):n
                m = A(i,j)/A(j,j);
                A(i,j) = A(i,j) - A(j,j)*m;
                for k = (j+1):n
                    A(i,k) = A(i,k) - A(j,k)*m; 
                end
                B(i) = B(i) - B(j)*m;
            end
        end
        fprintf('O sistema equivalente 1 é definido por: \n');
        disp(A)
        disp(B)
        x = zeros(n,1);
        x(end) = B(end)/(A(end,end));
        for i = (n-1):-1:1 
            soma = 0;
            for j = (i+1):n
                soma = soma + x(j)*A(i,j);
            end
            x(i) = (B(i) - soma)/A(i,i);
        end
        X(:,k1) = x;
    end
end
fprintf('O resultado do sistema é: \n')
disp(X)
fprintf('Verificação \n')
C = D*X;
disp(C)
end