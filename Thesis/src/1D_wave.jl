@inline function APC(x::Float64)
    a1=0.1
    b1=0.9
    # K1=30.
    K1 = 14348907000000000000000
    n1= 15
    @fastmath xpow = x^n1
    return a1 + b1 * xpow/ (K1 + xpow)
end

@inline function Cdc25(x::Float64)
    a2=0.2
    b2=0.8
    # K2=30.
    K2=590490000000000
    n2=10
    @fastmath xpow = x^n2
    return a2 + b2 * xpow/(K2 + xpow)
end

@inline function Wee1(x::Float64)
    a3=0.1
    b3=0.4
    # K3=30.
    K3=24300000
    n3=5
    @fastmath xpow = x^n3
    return a3 + b3 * K3/(K3 + xpow)
end

@inline function model(u::Float64,v::Float64, op)
    K, kval, bval, _e = getvalues(op)
    return @fastmath (Cdc25(u) * (v - u) - Wee1(u) * u) * K *_e, K *(kval - bval * APC(u) * v)
end

@inline function model(u::Vector{Float64},v::Vector{Float64}, op)
    K, kval, bval, _e = getvalues(op)
    return @fastmath @. (Cdc25(u) * (v - u) - Wee1(u) * u) * K *_e, K *(kval - bval * APC(u) * v)
end

function numsolODE(dt,ntime,v0,u0,op)
    u = zeros(0)
    v = zeros(0)
    push!(u, u0)
    push!(v, v0)
    for i = 1:Int(ntime/dt)
        du, dv = model(u[i], v[i], op)
        unew = u[i] .+ du*dt
        vnew = v[i] .+  dv*dt
        push!(u, unew[1])
        push!(v, vnew[1])
    end
    return u, v
end

function numsolPDE(op, sp, ip, mp, Du, Dv, period)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    _dx = (nx/lx)^2 #reciprocal gridsize
    u,v,unew,vnew,scale,uu,vv = initialize_arrays(ip,sp)

    #get non spatial coupled solution
    # u0, v0 = numsolODE(dt,ntime,43.65729831013529,12.57083907022486,op) 

    #Time loop
    it = 1
    while it<ntime/dt+dt
        reaction_diffusion_step!(u,unew,v,vnew,_dx,dt,scale,op,period,Du,Dv)
        u, unew = unew, u
        v, vnew = vnew, v    
        if it % save_step == 0
            uu[:,div(it,save_step)] = unew
            vv[:,div(it,save_step)] = vnew
        end
        it += 1
    end
    return uu, vv;
end

function numsolPDE2(op, sp, ip, mp, Du, Dv, period, usolved )
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = (nx/lx)^2 #reciprocal gridsize
    u,v,unew,vnew,scale,uu,vv,Duu = initialize_arrays(ip,sp)


    umaxofmin = maximum([ifelse(x < 25, x,0) for x in usolved[1,:]])

    Du0 = copy(Du)
    Dv0 = copy(Dv)

    #Time loop
    it = 0
    while it<ntime/dt
        Du .= Du0.*([ifelse(us<umaxofmin, 2*us/umaxofmin,1) for us in u])
        Dv .= Dv0.*([ifelse(us<umaxofmin, 2*us/umaxofmin,1) for us in u])
        reaction_diffusion_step!(u,unew,v,vnew,dx,dt,scale,op,period,Du,Dv)
        u, unew = unew, u
        v, vnew = vnew, v
        if it % save_step == 0
            uu[:,div(it,save_step)+1] = unew
            vv[:,div(it,save_step)+1] = vnew
            Duu[:,div(it,save_step)+1] = Du
        end
        it += 1
    end
    return uu, vv, Duu;
end

function reaction_diffusion_step!(u,unew,v,vnew,_dx,dt,scale,op,period,Du,Dv)
    du, dv = model(u[1],v[1],op)
    unew[1] = u[1] + (_dx*Du[1]*(u[2] - 2*u[1] + u[end]) + period/scale[1]*du)*dt
    vnew[1] = v[1] + (_dx*Dv[1]*(v[2] - 2*v[1] + v[end]) + period/scale[1]*dv)*dt
    @inbounds for i in 2:length(u)-1
        nscale = scale[i]
        if u[i] != u[i-1] || v[i] != v[i-1]
            du, dv = model(u[i],v[i],op)
        end
        @fastmath unew[i] = u[i] + (_dx*Du[i]*(u[i+1] - 2*u[i] + u[i-1]) + period/nscale*du)*dt
        @fastmath vnew[i] = v[i] + (_dx*Dv[i]*(v[i+1] - 2*v[i] + v[i-1]) + period/nscale*dv)*dt
    end
    if u[end] != u[end-1] || v[end] != v[end-1]
        du, dv = model(u[end],v[end],op)
    end
    unew[end] = u[end] + (_dx*Du[end]*(u[1] - 2*u[end] + u[end-1]) + period/scale[end]*du)*dt
    vnew[end] = v[end] + (_dx*Dv[end]*(v[1] - 2*v[end] + v[end-1]) + period/scale[end]*dv)*dt
end

function initialize_arrays(ip,sp)
    u0,v0,Pi,Po,pacemaker = getvalues(ip)
    lx,nx,ntime,dt,save_step = getvalues(sp)
    ΔL = nx
    #Array initializations
    cu = fill(u0, Int(nx))
    cu = cu
    cv = fill(v0, Int(nx))
    cv = cv
    cscale = Float64[ifelse(n%(ΔL) > ΔL/2 - pacemaker/2 && n%(ΔL) < ΔL/2 + pacemaker/2, Pi, Po) for n=1:nx]

    unew=copy(u)
    vnew=copy(v)

    uu = Array{Float64}(undef, nx, Int((ntime/save_step)/dt)) #returned array
    vv = Array{Float64}(undef, nx, Int((ntime/save_step)/dt)) #returned array

    Duu = Array{Float64}(undef, nx, Int((ntime/save_step)/dt)) #returned array
    return u,v,unew,vnew,scale,uu,vv,Duu
end

function reaction_diffusion_step!(u::AbstractArray{Float64,1},unew::AbstractArray{Float64,1},v::AbstractArray{Float64,1},vnew::AbstractArray{Float64,1},dx,dt,Du,Dv,period,scale::AbstractArray{Float64,1},du::AbstractArray{Float64,1},dv::AbstractArray{Float64,1})
    ns = length(u)
    r0 = 1:ns
    unew[r0] .= u[r0] .+ (dx^(-2)*Du*(u[circshift(r0,-1)] .- 2*u[r0] .+ u[circshift(r0,1)]) .+ period./scale[r0].*du[r0])*dt
    vnew[r0] .= v[r0] .+ (dx^(-2)*Dv*(v[circshift(r0,-1)] .- 2*v[r0] .+ v[circshift(r0,1)]) .+ period./scale[r0].*dv[r0])*dt
end

function getperiod(mp, sp)
    lx, nx, time, pacemaker, Du, Dv = getvalues(sp)

    dt = 0.001
    steps = Int(time/dt)

    #Init
    v0 = 0.1
    u0 = 0.1
    usol, vsol = numsolODE(dt,steps,u0,v0,mp)
    u0 = usol[length(usol)]
    v0 = vsol[length(vsol)]
    usol, vsol = numsolODE(dt,steps,u0,v0,mp)

    #GET MINIMUM
    umin = u0
    for i in 1:steps-1
        if usol[i] < umin
            umin = usol[i]
            v0 = vsol[i]
            u0 = umin
        end
    end
    usol, vsol = numsolODE(dt,steps,v0,u0,mp)

    #calculate period
    period = 0.
    for i in 2:length(usol)
        if abs(vsol[i] - v0) < 0.001 && abs(usol[i] - u0) < 0.001
            period = i
            break
        end
    end
    return period*dt, usol[1], vsol[1]
end

function dist(x1, t1, x2, t2)
    return sqrt((x1 - x2)^2 + (t1 - t2)^2)
end

@views function getmaxprofile_top(u,sp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx
    xv = -lx/2+0.5*dx:dx:lx/2-0.5*dx
    tv = 0:dt*save_step:ntime-dt

    uv = u[Int(nx/2),:]
    duv = uv[2:end]-uv[1:end-1]
    lmax = []
    j = 1
    while j < length(duv)
        if duv[j] > 0 && duv[j + 1] < 0
            push!(lmax,j)  # Save the index, which is a time, of the maximum
        end
        j += 1
    end
    # List containing for each maximum an array with the (t, x) profile
    lmaxes = [[(tv[lm], xv[Int(nx/2)])] for lm in lmax]

    @inbounds for i in Int(nx/2)+1:nx-1
        uv = u[i,:]

        # Detect the local maxima
        lmax = []  # start a new lmax
        duv = uv[2:end]-uv[1:end-1]

        j = 1
        while j < length(duv)
            if duv[j] > 0 && duv[j + 1] < 0
                push!(lmax,j)
            end
            j += 1
        end

        for lm in lmax
            t, x = tv[lm], xv[i]

            lastts = [l[end][1] for l in lmaxes]  # get the times
            lastxs = [xv[i - 1] for l in lmaxes]

            distances = [dist(x, t, lastxs[i], lastts[i]) for i in 1:length(lastxs)]
            # Find the smallest distance
            minind = argmin(distances)  # returns the corresponding index
            push!(lmaxes[minind],(t, x))
        end
    end

    return lmaxes
end

@views function getmaxprofile_bottom(u,sp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx
    xv = 0:dx:lx
    tv = dt:dt*save_step:ntime

    #select the bottom halve of the propagation wave
    uv = u[Int(ceil(nx/2)),:]
    duv = uv[2:end]-uv[1:end-1]
    lmax = []
    j = 1
    while j < length(duv)
        if duv[j] > 0 && duv[j + 1] < 0 && duv[j+2] < 0
            push!(lmax,j)  # Save the index, which is a time, of the maximum
        end
        j += 1
    end
    # List containing for each maximum an array with the (t, x) profile
    lmaxes = [[(tv[lm], xv[Int(nx/2)])] for lm in lmax]
    @inbounds for i in Int(nx/2)-1:-1:1
        uv = u[i,:]

        # Detect the local maxima
        lmax = []  # start a new lmax
        duv = uv[2:end] - uv[1:end-1]
        j = 1
        while j < length(duv)
            if duv[j] > 0 && duv[j + 1] < 0 && duv[j + 2] < 0
                push!(lmax,j)
            end
            j += 1
        end
        for lm in lmax
            t, x = tv[lm], xv[i]
            lastts = [l[end][1] for l in lmaxes]  # get the times
            lastxs = [xv[i + 1] for l in lmaxes]

            distances = [dist(x, t, lastxs[i], lastts[i]) for i in 1:length(lastxs)]
            # Find the smallest distance
            minind = argmin(distances)  # returns the corresponding index
            push!(lmaxes[minind],(t, x))
        end
    end
    
    # remove max lines that don't reach x=0
    i = 1
    while i < size(lmaxes,1)+1
        if lmaxes[i][end][2] != 0
            deleteat!(lmaxes,i)
            i -= 1
        end
        i += 1
    end
    return lmaxes
end

function get_wavespeed(u,sp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx  #gridsize
    max = getmaxprofile_bottom(u,sp)
    trf = .97
    wavespeed = zeros(length(max),5)
    # Detect point of the envelope for each profile
    for i=1:length(max)
        max_t_coord = Float64[t for (t,x) in max[i]]
        max_x_coord = Float64[x for (t,x) in max[i]]

        # fit t = a + b(tip_x_location - x)^beta through middle part
        # first select an appropriate part
        amp = maximum(max_t_coord) - minimum(max_t_coord)

        j = 1
        while max_t_coord[j] < minimum(max_t_coord) + (1-trf) * amp && j < length(max_t_coord)
            j += 1
        end
        m1 = j

        j = 1
        while max_t_coord[j] < minimum(max_t_coord) + trf * amp && j < length(max_t_coord)
            j += 1
        end
        m2 = j

        middle_part_t = max_t_coord[m1:m2]
        middle_part_x = max_x_coord[m1:m2]
        ts = (m2-m1)

        cost(t,p) = p[1] .+ p[2] * (t .- max_t_coord[m1]).^p[3]
        p0 = [200., -20, 1]
        fit = curve_fit(cost, middle_part_t, middle_part_x, p0)
        C, D, γ = fit.param
        wavespeed[i,:] .= C, D, γ, max_t_coord[m1], max_t_coord[m2]
    end
    return wavespeed
end

function get_envelope_points(max)
    x = Float64[]
    t = Float64[]
    trf = .95
    # Detect point of the envelope for each profile
    push!(x,200.)
    for i=1:length(max)
        max_t_coord = Float64[t for (t,x) in max[i]]
        max_x_coord = Float64[x for (t,x) in max[i]]

        # fit t = a + b(tip_x_location - x)^beta through middle part
        # first select an appropriate part
        amp = maximum(max_t_coord) - minimum(max_t_coord)

        j = 1
        while max_t_coord[j] < minimum(max_t_coord) + (1-trf) * amp && j < length(max_t_coord)
            j += 1
        end
        m1 = j

        j = 1
        while max_t_coord[j] < minimum(max_t_coord) + trf * amp && j < length(max_t_coord)
            j += 1
        end
        m2 = j

        middle_part_t = max_t_coord[m1:m2]
        middle_part_x = max_x_coord[m1:m2]

        cost(t,p) = p[1] .+ p[2] * (t).^p[3]

        p0 = [minimum(max_t_coord), -0.5, 1.5]
        fit = curve_fit(cost, middle_part_t, middle_part_x, p0)
        a, b, beta = fit.param

        time_coordinate = max_t_coord[end]
        x_intersection =  a + b*time_coordinate^beta
        if x_intersection < 1.
            break
        end
        if x_intersection < x[end]
            push!(x,x_intersection)
            push!(t,time_coordinate)
        end
    end
    popfirst!(x)
    return t, x
end

function get_envelope_speed(u,sp)
    max = getmaxprofile_bottom(u,sp)
    envelope_t, envelope_x = get_envelope_points(max)

    cost(t,p) = p[1] .+ p[2] * (t).^p[3]
    p0 = [200, -0.5,  1]
    fit = curve_fit(cost, envelope_t[2:end-1], envelope_x[2:end-1], p0)
    C, D, γ = fit.param

    return C, D, γ
end

# Solve the system with a fractional laplacian
function numsolfrac(op, sp, ip, mp, Du, Dv, period,α)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx
    u,v,unew,vnew,scale,uu,vv = initialize_arrays(ip,sp)
    nscale = period./scale

    # α = 1.1
    # for nx in 10:200
    #     dx = lx/nx
    #     A = get_A_matrix(nx)
    #     B = get_B_matrix(dx,nx)
    #     T = inv(A)*B
    #     Λ, P = eigen(T)
    #     println(nx," doesn't work")
    #     if (isreal(P) && isreal(Λ))
    #         println(nx," would work")
    #     end
    # end
    # println("end")

    A = get_A_matrix(nx)
    B = get_B_matrix(dx,nx)
    T = inv(A)*B
    Λ, P = eigen(T)
    Λ = Diagonal(Λ)
    Λ = chop.(Λ)
    if (isreal(P) == false && isreal(Λ) == false)
        uu = 0
        vv = 0
        return uu,vv
    end

    Pinv = inv(P)
    Lu = dt*Du[1].*(P)*(Λ.^(α/2)*Pinv)
    Lv = dt*Dv[1].*(P)*(Λ.^(α/2)*Pinv)
    # Lu = dt*Du[1]*Λ^(α/2)
    # Lv = dt*Dv[1]*Λ^(α/2)

    # Ru = exp(-Lu)
    # Qu = exp(-Lu/2)
    Ru = ((24*I - 6*Lu)*inv(24*I + 18*Lu + 6*Lu^2 + Lu^3))
    # Qu = diag(24*(8*I - Lu)*inv(192*I + 72*Lu + 12*Lu^2 + Lu^3))
    # ϕu = diag((96*I + 12*Lu + Lu^2)*inv(192*I + 72*Lu + 12*Lu^2 + Lu^3)*dt)
    # ϕ1u = diag((4*I - Lu)*inv(24*I + 18*Lu + 6*Lu^2 + Lu^3)*dt)
    # ϕ2u = diag(2*(4*I + Lu)*inv(24*I + 18*Lu + 6*Lu^2 + Lu^3)*dt)
    # ϕ3u = diag((4*I + 3*Lu + Lu^2)*inv(24*I + 18*Lu + 6*Lu^2 + Lu^3)*dt)

    # Rv = exp(-Lv)
    # Qv = exp(-Lv/2)
    Rv = ((24*I - 6*Lv)*inv(24*I + 18*Lv + 6*Lv^2 + Lv^3))
    # Qv = diag(24*(8*I - Lv)*inv(192*I + 72*Lv + 12*Lv^2 + Lv^3))
    # ϕv = diag((96*I + 12*Lv + Lv^2)*inv(192*I + 72*Lv + 12*Lv^2 + Lv^3)*dt)
    # ϕ1v = diag((4*I - Lv)*inv(24*I + 18*Lv + 6*Lv^2 + Lv^3)*dt)
    # ϕ2v = diag(2*(4*I + Lv)*inv(24*I + 18*Lv + 6*Lv^2 + Lv^3)*dt)
    # ϕ3v = diag((4*I + 3*Lv + Lv^2)*inv(24*I + 18*Lv + 6*Lv^2 + Lv^3)*dt)

    # fftu = fft(u)
    # fftv = fft(v)
    function ETDRK4!(u,unew,v,vnew,op)
        fu, fv = model(real(ifft(u)),real(ifft(v)),op)
        Nfu = fft(fu.*nscale)
        Nfv = fft(fv.*nscale)
        au = Qu.*fftu + ϕu.*Nfu
        av = Qv.*fftv + ϕv.*Nfv
        fau, fav = model(real(ifft(au)),real(ifft(av)),op)
        Nfau = fft(fau.*nscale)
        Nfav = fft(fav.*nscale)
        bu = Qu.*fftu + ϕu.*Nfau
        bv = Qv.*fftv + ϕv.*Nfav
        fbu, fbv = model(real(ifft(bu)),real(ifft(bv)),op)
        Nfbu = fft(fbu.*nscale)
        Nfbv = fft(fbv.*nscale)
        cu = Qu.*au+ϕu.*(2*Nfbu-Nfu)
        cv = Qv.*av+ϕv.*(2*Nfbv-Nfv)
        fcu, fcv = model(real(ifft(cu)),real(ifft(cv)),op)
        Nfcu = fft(fcu.*nscale)
        Nfcv = fft(fcv.*nscale)
        fftu .= Ru.*fftu + (ϕ1u.*Nfu + ϕ2u.*(Nfau+Nfbu) + ϕ3u.*Nfcu)
        fftv .= Rv.*fftv + (ϕ1v.*Nfv + ϕ2v.*(Nfav+Nfbv) + ϕ3v.*Nfcv)
    end

    function ETD!(u,unew,v,vnew,op)
        fu, fv = model(u,v,op)
        mul!(unew, Ru, u)
        mul!(vnew, Rv, v)
        unew .+= nscale .*fu*dt;
        vnew .+= nscale .*fv*dt;
    end

    it = 1
    while it<ntime/dt+dt
        ETD!(u,unew,v,vnew,op)
        u, unew = unew, u
        v, vnew = vnew, v
        if it % save_step == 0
            # println(it*dt)
            uu[:,div(it,save_step)] = unew
            vv[:,div(it,save_step)] = vnew
            # uu[:,div(it,save_step)] = real(ifft(fftu))
            # vv[:,div(it,save_step)] = real(ifft(fftv))
        end
        it += 1
    end
    return uu, vv;
end

function get_P_matrix(nx)
    P = [exp(-im*(i-1)*(j-1)*2*π/nx) for j in 1:nx, i in 1:nx]
    # P = [cos(i*j*pi/nx) for i in 1:nx, j in 1:nx]
    return P
end

function get_Λ_matrix(dx,nx)
    Λ = [(4*sin((n-1)*π/(2*nx))^2)/(dx^2 * (1 - 1/3*sin((n-1)*π/(2*nx))^2)) for n in 1:nx]
    return Λ
end

function get_A_matrix(nx)
    A = zeros(nx,nx)
    for i in 1:nx
        A[i,i] = 5/6
        A[i,i%nx+1]=1/12
        A[i,(i+nx-2)%nx+1]=1/12
    end
    return A
end

function get_B_matrix(dx,nx)
    A = spzeros(nx, nx)
    for i in 1:nx
        A[i,i] = 2
        A[i,i%nx+1]=-1
        A[i,(i+nx-2)%nx+1]=-1
    end
    return A.*dx^(-2)
end

function Base.chop(x::Float64)
    tol = 1000eps()
    if abs(x)< tol
        return 0.
    else
        return x
    end
end

function Base.chop(x::ComplexF64)
    tol = 10000000eps()
    if imag(x)< tol
        return chop(real(x))
    else
        return x
    end
end