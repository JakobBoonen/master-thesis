#V is the type of the returned arrays
function initialize_arrays_2D(ip,sp,V)
    u0,v0,Pi,Po,pacemaker = getvalues(ip)
    lx,nx,ntime,dt,save_step = getvalues(sp)
    #Matrix initializations
    u = Float64[u0 for i ∈ 1:nx, j ∈ 1:nx] 
    v = Float64[v0 for i ∈ 1:nx, j  ∈ 1:nx]
    # u = readdlm("last_u.csv", ',', Float64)   
    # v = readdlm("last_v.csv", ',', Float64)   
    # u = reshape(u,nx,nx)
    # v = reshape(v,nx,nx)
    scale = Float64[ifelse((n > nx/4 - pacemaker/2 && n < nx/4 + pacemaker/2 && m > nx/4 - pacemaker/2 && m < nx/4 + pacemaker/2) ||
                           (n > 3*nx/4 - pacemaker/2 && n < 3*nx/4 + pacemaker/2 && m > 3*nx/4 - pacemaker/2 && m < 3*nx/4 + pacemaker/2)
                            , Pi, Po) for n=1:nx, m=1:nx]
    # u,v,scale=V.((cu,cv,cscale))
    unew=copy(u)
    vnew=copy(v)
    uu = Array{Float64}(undef, nx, nx, Int(ntime/dt/save_step)) #returned array
    vv = Array{Float64}(undef, nx, nx, Int(ntime/dt/save_step)) #returned array
    u,v,unew,vnew,scale,uu,vv
end

function plot_wave_2D(op, sp, ip, mp, V)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx  #gridsize

    u, v = numsolPDE_2D(op, sp, ip, mp, V)
    # writedlm("u2D.csv",  u, ',')
    # u = readdlm("u2D.csv", ',', Float64)   
    # u = reshape(u,nx,nx,Int(ntime/dt/save_step))

    xs = -lx/2+0.5*dx:dx:lx/2-0.5*dx
    color_range = (minimum(u),maximum(u))


    fig = Figure()
    ax = Axis(fig[1, 1],aspect = AxisAspect(1))
    sl = Slider(fig[2, 1], range = 1:size(u,3), horizontal = true, startvalue = 1)
    du = lift(sl.value) do t
        u[:,:,t]
    end
    hm = heatmap!(ax, xs, xs, du, colormap = :blues,colorrange = color_range)
    Colorbar(fig[1, 2], hm, height=Relative(1))
    display(fig)    

    # time = Observable(1)
    # us = @lift(u[:,:,$time])
    # fig = Figure()
    # ax = Axis(fig[1, 1],aspect = AxisAspect(1))
    # hm = heatmap!(ax, xs, xs, us, colormap = :blues,colorrange = color_range)
    # Colorbar(fig[1, 2], hm, height=Relative(1))
    # timestamps = 1:size(u,3)
    # GLMakie.record(fig, "figures/animation.mp4", timestamps; framerate = 60) do t
    #     time[] = t
    # end
end

function numsolPDE_2D(op, sp, ip, mp,V,)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    _dx = (nx/lx)^2 #gridsize
    u,v,unew,vnew,scale,uu,vv = initialize_arrays_2D(ip,sp,V)
    period,DuParam,DvParam=getvalues(mp)
    Du = Float64[DuParam for x ∈ 1:nx, x ∈ 1:nx]
    Dv = Float64[DvParam for x ∈ 1:nx, x ∈ 1:nx]
    #Time loop
    t = 0.0; it = 0
    while it<ntime/dt
        reaction_diffusion_step_2D!(u,unew,v,vnew,_dx,dt,scale,op,period,Du,Dv)
        u, unew = unew, u
        v, vnew = vnew, v
        if it % save_step == 0
            uu[:,:,Int(it/save_step)+1] = unew
            vv[:,:,Int(it/save_step)+1] = vnew
        end
        t += dt; it+=1
    end
    # writedlm("last_u.csv",  u, ',')
    # writedlm("last_v.csv",  v, ',')
    return uu, vv
end

function reaction_diffusion_step_2D!(u::Array{Float64,2},unew::Array{Float64,2},v::Array{Float64,2},vnew::Array{Float64,2},_dx,dt,scale,op,period,Du,Dv)
    ns = size(u,1)
    #Four corners
    du, dv = model(u[1,1],v[1,1],op)
    unew[1,1] = u[1,1] + (_dx*Du[1,1]*(u[2,1] - 4*u[1,1] + u[ns,1] + u[1,2] + u[1,ns]) + period/scale[1,1]*du)*dt
    vnew[1,1] = v[1,1] + (_dx*Dv[1,1]*(v[2,1] - 4*v[1,1] + v[ns,1] + v[1,2] + v[1,ns]) + period/scale[1,1]*dv)*dt
    du, dv = model(u[1,ns],v[1,ns],op)
    unew[1,ns] = u[1,ns] + (_dx*Du[1,ns]*(u[1,1] - 4*u[1,ns] + u[1,ns-1] + u[2,ns] + u[ns,ns]) + period/scale[1,ns]*du)*dt
    vnew[1,ns] = v[1,ns] + (_dx*Dv[1,ns]*(v[1,1] - 4*v[1,ns] + v[1,ns-1] + v[2,ns] + v[ns,ns]) + period/scale[1,ns]*dv)*dt
    du, dv = model(u[ns,1],v[ns,1],op)
    unew[ns,1] = u[ns,1] + (_dx*Du[ns,1]*(u[1,1] - 4*u[ns,1] + u[ns-1,1] + u[ns,2] + u[ns,ns]) + period/scale[ns,1]*du)*dt
    vnew[ns,1] = v[ns,1] + (_dx*Dv[ns,1]*(v[1,1] - 4*v[ns,1] + v[ns-1,1] + v[ns,2] + v[ns,ns]) + period/scale[ns,1]*dv)*dt
    du, dv = model(u[ns,ns],v[ns,ns],op)
    unew[ns,ns] = u[ns,ns] + (_dx*Du[ns,ns]*(u[1,ns] - 4*u[ns,ns] + u[ns-1,ns] + u[ns,1] + u[ns,ns-1]) + period/scale[ns,ns]*du)*dt
    vnew[ns,ns] = v[ns,ns] + (_dx*Dv[ns,ns]*(v[1,ns] - 4*v[ns,ns] + v[ns-1,ns] + v[ns,1] + v[ns,ns-1]) + period/scale[ns,ns]*dv)*dt

    @inbounds for j in 2:ns-1
        if u[1,j] != u[1,j-1] || v[1,j] != v[1,j-1]
            du, dv = model(u[1,j],v[1,j],op)
        end    
        unew[1,j] = u[1,j] + (_dx*Du[1,j]*(u[2,j] - 4*u[1,j] + u[ns,j] + u[1,j+1] + u[1,j-1]) + period/scale[1,j]*du)*dt
        vnew[1,j] = v[1,j] + (_dx*Dv[1,j]*(v[2,j] - 4*v[1,j] + v[ns,j] + v[1,j+1] + v[1,j-1]) + period/scale[1,j]*dv)*dt
    end

    @inbounds for j in 2:ns-1
        if u[ns,j] != u[ns,j-1] || v[ns,j] != v[ns,j-1]
            du, dv = model(u[ns,j],v[ns,j],op)
        end    
        unew[ns,j] = u[ns,j] + (_dx*Du[ns,j]*(u[1,j] - 4*u[ns,j] + u[ns-1,j] + u[ns,j+1] + u[ns,j-1]) + period/scale[ns,j]*du)*dt
        vnew[ns,j] = v[ns,j] + (_dx*Dv[ns,j]*(v[1,j] - 4*v[ns,j] + v[ns-1,j] + v[ns,j+1] + v[ns,j-1]) + period/scale[ns,j]*dv)*dt
    end

    @inbounds for i in 2:ns-1   
        if u[i,1] != u[i-1,1] || v[i,1] != v[i-1,1]
            du, dv = model(u[i,1],v[i,1],op)
        end    
        unew[i,1] = u[i,1] + (_dx*Du[i,1]*(u[i+1,1] - 4*u[i,1] + u[i-1,1] + u[i,2] + u[i,ns]) + period/scale[i,1]*du)*dt
        vnew[i,1] = v[i,1] + (_dx*Dv[i,1]*(v[i+1,1] - 4*v[i,1] + v[i-1,1] + v[i,2] + v[i,ns]) + period/scale[i,1]*dv)*dt
    end
    
    @inbounds for i in 2:ns-1   
        if u[i,ns] != u[i-1,ns] || v[i,ns] != v[i-1,ns]
            du, dv = model(u[i,ns],v[i,ns],op)
        end    
        unew[i,ns] = u[i,ns] + (_dx*Du[i,ns]*(u[i+1,ns] - 4*u[i,ns] + u[i-1,ns] + u[i,1] + u[i,ns-1]) + period/scale[i,ns]*du)*dt
        vnew[i,ns] = v[i,ns] + (_dx*Dv[i,ns]*(v[i+1,ns] - 4*v[i,ns] + v[i-1,ns] + v[i,1] + v[i,ns-1]) + period/scale[i,ns]*dv)*dt
    end
    
    du, dv = model(u[2,2],v[2,2],op)
    @inbounds for i in 2:ns-1   
        @inbounds for j in 2:ns-1
            if u[i,j] != u[i-1,j-1] || v[i,j] != v[i-1,j-1]
                du, dv = model(u[i,j],v[i,j],op)
            end    
            unew[i,j] = u[i,j] + (_dx*Du[i,j]*(u[i+1,j] - 4*u[i,j] + u[i-1,j] + u[i,j+1] + u[i,j-1]) + period/scale[i,j]*du)*dt
            vnew[i,j] = v[i,j] + (_dx*Dv[i,j]*(v[i+1,j] - 4*v[i,j] + v[i-1,j] + v[i,j+1] + v[i,j-1]) + period/scale[i,j]*dv)*dt
        end
    end
end