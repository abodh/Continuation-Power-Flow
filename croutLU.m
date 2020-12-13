function error_VD = croutLU(A, b)
L = zeros(length(A));
U = zeros(length(A));

% the first column of lower traingular is the same as matrix A
L(1:length(A),1) = A(1:length(A),1);

% setting the diagonals of upper traingular to 1
for i = 1:length(A)
    U(i,i) = 1;
end

% calculating the first row of upper traingular
U(1,2:length(A)) = A(1,2:length(A))/L(1,1);

for j = 2 : length(A)
    for k = j : length(A)
        % updating the columns of L
        L(k,j) = A(k,j) - L(k,1:j-1)*U(1:j-1,j);
    end
    for k = j+1 : length(A)
        % updating the rows of U
        U(j,k) = (A(j,k) - (L(j,1:j-1)*U(1:j-1,k)))/L(j,j);
    end
end

% forward substitution
y = zeros(length(A),1);
for i = 1:length(A)
        if i == 1
            y(1) = b(1)/L(1,1);
        else
            y(i) = (b(i)-L(i,1:i-1)*y(1:i-1))/L(i,i);
        end
end

% backward substitution
x = zeros(length(A),1);
for i = length(A):-1:1   
    if i == length(A)
        x(length(A)) = y(length(A))/U(length(A),length(A));
    else
        x(i) = (y(i)-U(i,i+1:length(A))*x(i+1:length(A)))/U(i,i);
    end
end
error_VD = x;
end