function [] = L6Ex1(A)
n = length(A);
D = A;
b = zeros(n,1);
for k = 1:n
    e = zeros(k);       %Cria matrizes a serem substituidas pelas matrizes principais
    for i = 1:k         %Fixa linha
        for j = 1:k     %itera coluna
            e(i,j) = A(i,j);    %Seleciona matrizes principais
        end
    end
    a = det(e);         %Calcula determinante das M Principais
    b(k,1) = a;         %Seleciona elementos da matriz anterior
end
for i1 = 1:n     %Seleciona elementos da matriz de determinantes
    if(abs(b(i1,1))<10^-5)  %Comparação
       a = 1;       %Desautoriza algoritmo LU
       fprintf('Não é possível decompor LU a matriz \n')
       break 
    end
end
a = 0;          %Autoriza algoritmo LU
if (a~=1)
    X = zeros(n);   %Matriz a se tornar inversa
        L = eye(n);     %Matriz a se tornar L
        for j = 1:n-1       %Seleciona colunas
            for i = j+1:n   %Seleciona linha
                m = A(i,j)/A(j,j);  %Coeficientes da matriz
                L(i,j) = m;         %Guarda elementos em L
                A(i,j) = A(i,j) - A(j,j)*m;   %Zera elementos abaixo da diagonal principal
                for k = (j+1):n
                    A(i,k) = A(i,k) - A(j,k)*m;     %Realiza a soma para as outras colunas da mesma linha
                end
            end
        end
    U = A;  %Seleciona a matriz U, obtida anteriormente
    for k1 = 1:n    %Realiza os n sistemas
        B = zeros(n,1); %Gera matriz de 1s
        B(k1,1) = 1;    %Aloca o 1 
        y = zeros(n,1);     %Guarda respostas do primeiro sistema
        y(1,1) = B(1,1);    %Primeira solução
        %Soluciona sistema triangular inferior
        for i = 2:n     %Linha
            soma = 0;
            for j = 1:i-1   %Coluna
                soma = soma + (y(j)*L(i,j));
            end
            y(i) = (B(i)-soma);
        end
        %Soluciona sistema triangular superior
        x = zeros(n,1); %Gera matriz a ser solução do sistema 2
        x(n,1) = y(n,1)/U(n,n); %Primera solução
        for i = n-1:-1:1                %Linhas
            soma = 0;
            for j = (i+1):n             %Colunas
                soma = soma + (x(j)*U(i,j));
            end
            x(i) = (y(i)-soma)/(U(i,i));
        end
    X(:,k1) = x;        %Guarda Matriz resultado
    end
end
fprintf('Verificação \n')
C = L*U;
disp(C)
fprintf('Matrizes resultado \n')
fprintf('L \n')
disp(L)
fprintf('U \n')
disp(U)
fprintf('A matriz inversa é: \n')
E = D*X;        %Verificação
disp(X)
fprintf('Verificação da inversa: \n')
disp(E)
end