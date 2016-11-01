library(pwr)
#Calculate power for a one sample, 2-sided t-test
pwr.t.test(d=0.05,n=3150,sig.level=0.05,type="one.sample",alternative="two.sided")
