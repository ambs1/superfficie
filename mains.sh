f77 -c VCs3_3body.f 
f77 -c anal_3at.f
f77 anal_3at.o VCs3_3body.o
time ./a.out