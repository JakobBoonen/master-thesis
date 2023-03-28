using Thesis
 
function go()
    nx = 1020
    # nx = 100
    # op=OscillatorParam(K=1,kval=2,bval=0.1,ϵ=0.01)
    op=OscillatorParam(K=1,kval=2,bval=0.1,ϵ=100) #reciprocal ϵ
    #lx is length in micrometer
    sp=SimParam(lx=400,nx=nx,time=100,dt=0.001,save_step=100)
    ip=InitParam(u0=12.57083907022486,v0=43.65729831013529,Pi=9.5,Po=10.,pacemaker=nx/100)#nx/100
    # mp=ModelParam(period=19.924,Du=.1,Dv=0.01)
    mp=ModelParam(period=19.924,Du=1,Dv=0.1,α=2)
    plot(op, sp, ip, mp);
end

go()