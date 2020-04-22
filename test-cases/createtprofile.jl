using Plots
f1(x,L)=exp(-x/L)
L=2
x=collect(range(1, 10, length=1000))
y=collect(range(1, 101, length=1000))
phi=( (f1.(x,L*0.5)*20).+8 ) .- ((f1.(x,L*0.5)*5).^2.5 );
plot( phi ,-y*40,label=false )
