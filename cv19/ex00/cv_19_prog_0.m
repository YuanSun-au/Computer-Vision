
%% 0.1

A = [1 0; 0 4]
A = [1 -2; -2 4]
A = [1 -5; -5 1]

[Q,l] = eig(A)


%% 0.2

a = 1+2i
b = 1-i

a*b
a/b

%% 0.3

norm(a)
norm(b)

angle(a) == atan( imag(a) / real(a) )
angle(b) == atan( imag(b) / real(b) )
