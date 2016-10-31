library(Rmpfr)
(one <- mpfr(1, 120))
x <- exp(709)
y <- mpfr(x,25)
y
test <- gamma(1000)
test2 <- mpfr(test,1000);test2


#This is the code that works
mpfr(exp(750),128)
gamma(as(3000,"mpfr"))