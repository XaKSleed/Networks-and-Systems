close all 
clear, clc
%%g_x = [1 0 1 1]; %%x^3 + x + 1
%%g_x = [1 1 0 1]; %%x^3 + x^2 + 1
%%g_x = [1 0 0 1 0]; %%x^4 + x + 1
g_x = [1 0 0 1 0 1]; %% x^5 + x^2 + 1
%%g_x = [1 0 0 0 1 1 1 0 1]; %%x^8 + x^4 + x^3 + X^2 + 1
%%g_x = [1 0 0 0 0 1 0 0 1 1]; %% x^9 + x^4 + x + 1
d = 0; %%code distance
l = 4; %%length message
r = length(g_x) - 1; %% deg g(x)
messages = zeros(2^l,l); %%messages
code_words = zeros(2^l, l+r); %%code words
A_i = zeros(1,r + l); %% array of weights
x_r = zeros(1, length(g_x));
x_r(1) = 1;
%%disp("x^r");
%%disp(x_r);
p = 0:0.01:1; %%chanel error prob
Pe = zeros(1, length(p));
Pe_h = zeros(1, length(p));
n = l + r;
%%creating messages
for i = 1:2^l
    number = i-1;
    for j = 1:l
        messages(i, l - j + 1) = rem(number, 2);
        number = fix(number/2);
    end
end
%%disp(messages);

     %%Encoding algorytm
 for v = 1:2^l
    m = messages(v,:);
    m_x = product(m, x_r);
    c_x = newDivision(m_x, g_x);
    for i = 1:length(c_x)
        if c_x(i) >= 2
            while c_x(i) >= 2
                c_x(i) = c_x(i) - 2;
            end
        end
        if c_x(i) < 0
            while c_x(i) < 0
                c_x(i) = c_x(i) + 2;
            end
        end
    end
    
    a_x = zeros(1, length(m_x));
    k = 1;
    for i = 1:length(m_x)
        a_x(k) = m_x(i);
        k = k + 1;
    end
    for i = 1:length(c_x)
        a_x(i) = a_x(i) + c_x(i);
        k = k + 1;
    end
    code_words(v, :) = a_x;
 end   
    %%counting weights
 for i = 1:2^l
     s_x = sum(code_words(i, :));
     if s_x ~= 0
        A_i(s_x ) = A_i(s_x)+1;
     end
 end
 %%searching minimum code distance
 for i = 2:length(A_i)
    if A_i(i)~= 0
        d = i;
%         disp("code distance is ");
%         disp(i - 1);
        break;
    end
end
 %%right P
 for i = 1:length(p)
    sum = 0;
    for j = d:length(A_i)
        sum = sum + (A_i(j)*p(i)^(j))*((1 - p(i))^(n-(j)));
    end
    Pe(i) = sum;
end
%%hight P
for i = 1:length(p)
    sum = 0;
    for j = 0:2
        C = factorial(n)/(factorial(j)*factorial(n - j));
        sum = sum + (C *(p(i)^j)*((1 - p(i))^(n-(j)))); 
    end
    Pe_h(i) = 1 - sum;
end
figure(1);
hold on;
grid on;
plot(p, Pe, p, Pe_h);
title("x^5 + x^2 + 1");

function polynom = product(polynom1, polynom2) %%polinoms product function
    matrix = zeros(length(polynom1), ((length(polynom1) - 1) + (length(polynom2) - 1)) + 1);
    for i = 1:length(polynom1)
        for j = 1:length(polynom2)
            matrix(i,(((i - 1) + (j - 1)) + 1)) = polynom1(i)*polynom2(j);
        end
    end
    for i = 1:((length(polynom1) - 1) + (length(polynom2) - 1)) + 1
        for j = 2:length(polynom1)
            matrix(1,i) = matrix(1,i) + matrix(j,i);
        end
    end
    polynom = matrix(1,:);
end
function polynom = newDivision(polynom1, polynom2) %%polinoms division function

qpoly = zeros(1, length(polynom1));
res_poly = zeros(1, length(polynom1));
if length(polynom1) < length(polynom2)
    polynom = qpoly;
    return;
else
    for i=1:length(polynom2)
        if polynom2(i) ~= 0
            k = i;
            break;
        end
    end
   
   sub_poly2 = zeros(1, length(polynom2) - k - 1);
    for i=1:length(polynom2)- k + 1
        sub_poly2(i) = polynom2(k);
        k = k+1;
    end
    k = 0;
    
   while true
    
    for i = 1:length(polynom1)
        if polynom1(i) ~=0
            k = i;
            break;
        end
    end
     if k == 0
        polynom = res_poly;
        return;
     end
     if length(polynom1) - k + 1 < length(sub_poly2)
        polynom = polynom1;
        break;
    end
        for j = 1:length(sub_poly2)
            polynom1(k) = polynom1(k) - sub_poly2(j);
            if polynom1(k) < 0
                polynom1(k) = polynom1(k) + 2;
            end
            k = k + 1;
        end
   
   end
end

end