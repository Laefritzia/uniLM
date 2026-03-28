
# fitted vs. res
#plot(ehat,yhat)
# fitted vs. std res
#plot(r_std,yhat)
# fitted vs. studentized res
#plot(r_t,yhat)
#plot(sqrt(abs(r_t)),yhat)
#lines(lowess(sqrt(abs(r_t)),yhat))

# studentized res vs leverage
#plot(diag(H), r_t)

# cook plot (red line is treshhold for outlier)
#barplot(t(cook),
#        ylim=c(0,max(0.21,max(t(cook)))*1.2),
#        names.arg=1:n)
#abline(h=4/n,col="red")

# QQ: empirical vs theoretical distribution
#plot(qnorm((1:n)/n - 0.01), sort(ehat))
#abline(mean(ehat), sd(ehat))





