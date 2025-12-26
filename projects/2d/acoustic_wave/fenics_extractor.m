


%%%%% Program to manipulate arrays from fenics %%%%%%%%

delta = 0.01;

u_exact = load('/Users/sudhipv/documents/fenics/study/Chaitanya/acoustic/u_exact_tseries.dat');
u_fenics = load('/Users/sudhipv/documents/fenics/study/Chaitanya/acoustic/u_tseries.dat');
tt = delta:delta:1;
plot(tt,u_fenics(568,:),tt,u_exact(568,:))