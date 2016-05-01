Muts = c(10, 10, 10, 50, 200, 200, 100, 10)
y = Muts[1:length(Muts)-1]
x = Muts[2:length(Muts)]
z = seq(0,length(Muts)-2)
j = (z + 1)^2*x^2/y^2 *(1/x + 1/y)

MutVals=seq(0, length(Muts)-1)

maxColorValue=max(j)
palette <- colorRampPalette(c("gray15","gray95"))(maxColorValue)
cols=palette[cut(log(j), maxColorValue)]

cols=c(cols, "#FCFCFC")

cols

xx<-barplot(Muts, 
        col=cols,
        xlab=expression(italic(x[i])),
        ylab=expression(italic(mP(x[i]))),
        names.arg=MutVals, 
        ylim=c(0,250))

#        cex.lab=1.25, cex.axis=1.2, cex.main=1.25, cex.sub=1.25

#We report standard deviations in red
text(x = xx,y=Muts,label = c(round(sqrt(j), digits=2),NA), pos = 3, cex = .9, col = "red")


a=gray.colors(n=length(unique(j)),start = min(j)/max(j), end=1)
b=sort(unique(j))
