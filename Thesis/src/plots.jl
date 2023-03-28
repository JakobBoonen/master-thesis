function plot(op, sp, ip, mp)
    # envelope_exponent(op, sp, ip, mp)
    # envelope_for_different_α(op, sp, ip, mp)
    wavespeed_for_different_α(op, sp, ip, mp)
    # wave(op, sp, ip, mp)
    # wavespeed(op,sp,ip,mp)
end

function envelope_exponent(op, sp, ip, mp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx
    period,DuParam,DvParam,α =getvalues(mp)

    Du = Float32[DuParam for x ∈ 1:nx]
    Dv = Float32[DvParam for x ∈ 1:nx]

    α = 1
    u, v = numsolfrac(op, sp, ip, mp, Du, Dv, period,α)
    # u, v = numsolPDE(op, sp, ip, mp, Du, Dv, period)

    D, C, gamma = get_envelope_speed(u,sp)

    fig = Figure()
    xs = dx/2:dx:lx-dx/2
    t = dt:dt*save_step:ntime

    ax1 = Axis(fig[1, 1], xlabel = "t", ylabel="x")
    limits!(ax1,0,ntime,0,lx)
    hm = heatmap!(ax1,t,xs,transpose(u),colormap = :blues)
    Colorbar(fig[1, 2], hm, height=Relative(1))

    lines!(ax1,t,D.+C.*t.^gamma)

    display(fig) 
end

function wave(op, sp, ip, mp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx
    period,DuParam,DvParam,α =getvalues(mp)
    # ΔL = nx/3
    # ΔD = 1
    # fraction_of_interior = 0.5
    # Du = Float64[ifelse((x%(ΔL) > (0.5*(1-fraction_of_interior)*ΔL) && x%(ΔL) < ((1-0.5*(1-fraction_of_interior))*ΔL)), ΔD * DuParam, 0.001*DuParam) for x ∈ 1:nx]
    # Dv = Float64[ifelse((x%(ΔL) > (0.5*(1-fraction_of_interior)*ΔL) && x%(ΔL) < ((1-0.5*(1-fraction_of_interior))*ΔL)), ΔD * DvParam, 0.001*DvParam) for x ∈ 1:nx]

    Du = Float32[DuParam for x ∈ 1:nx]
    Dv = Float32[DvParam for x ∈ 1:nx]

    # u, v = numsolPDE(op, sp, ip, mp, Du, Dv, period)
    u, v = numsolfrac(op, sp, ip, mp, Du, Dv, period,α)

    if (u == 0 && v == 0)
        println("No real solution")
    else

    # PLOT
    fig = Figure()
    xs = dx/2:dx:lx-dx/2
    t = dt:dt*save_step:ntime

    ax1 = Axis(fig[1, 1], xlabel = "t", ylabel="x")
    limits!(ax1,0,ntime,0,lx)
    hm = heatmap!(ax1,t,xs,transpose(u),colormap = :blues)
    Colorbar(fig[1, 2], hm, height=Relative(1))

    #Plot diffusion parameter
    # ax0 = Axis(fig[1, 0], xlabel = "Du", ylabel = "x",xlabelsize=20,ylabelsize=20)
    # limits!(ax0,0,2,0,lx)
    # linkyaxes!(ax0, ax1)
    # colsize!(fig.layout, 1, Aspect(1, 0.2))
    # lines!(ax0,Du,xs)

    #Subdiffusion
    # lines!(ax1,t.+6.42,200 .+ 15 .*t.^1,color=:orange,linewidth=3)

    # Plot the envelope points
    # max = getmaxprofile_bottom(u,sp)
    # envelope_t, envelope_x = get_envelope_points(max)
    # scatter!(ax1,envelope_t,envelope_x)

    # Fit wavespeed
    # fit_waves = get_wavespeed(u,sp)
    # for (D,C,gamma,t0,t_end) in eachrow(fit_waves)
    #     ts = 0:0.001:t_end-t0
    #     lines!(ax1,t0.+ts,D .+ C * ts.^gamma,color=:red,linewidth=3,label="min")
    # end

    # Plot the maxima of the profile
    # for max in getmaxprofile_bottom(u,sp)
    #     max_t_coord = Float64[t for (t,x) in max]
    #     max_x_coord = Float64[x for (t,x) in max]
    #     scatter!(ax1,max_t_coord,max_x_coord)
    # end
    
    # axislegend(labelsize=32)
    # save("D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/figures/subdiffusionM5.png",fig)
    display(fig) 
    end
end

function wavespeed(op, sp, ip, mp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    period,DuParam,DvParam,α =getvalues(mp)
    dx = lx/nx

    u, v = numsolfrac(op, sp, ip, mp, DuParam, DvParam, period,α)

    fig = Figure()
    xs = dx/2:dx:lx-dx/2
    t = dt:dt*save_step:ntime

    ax1 = Axis(fig[1, 1], xlabel = "t", ylabel="x")
    limits!(ax1,0,ntime,0,lx)
    hm = heatmap!(ax1,t,xs,transpose(u),colormap = :blues)
    Colorbar(fig[1, 2], hm, height=Relative(1))
  
    # Fit envelope
    for wave in eachrow(get_wavespeed(u,sp))
        ts = 0:0.01:(wave[end]-wave[end-1])
        lines!(ax1,ts.+wave[end-1],wave[1] .+ wave[2] * ts.^wave[3],color=:orange,linewidth=3,label="min")
    end
    # axislegend(labelsize=32)
    display(fig)    
end

function envelope_for_different_α(op,sp,ip,mp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx  #gridsize
    period,DuParam,DvParam=getvalues(mp)
    xs = 0:dx:lx-dx
    t = 0:dt*save_step:ntime-dt

    Du = Float32[DuParam for x ∈ 1:nx]
    Dv = Float32[DvParam for x ∈ 1:nx]

    envelope = Vector{Float64}[]
    for α in 1:0.01:1.02
        u, v = numsolfrac(op, sp, ip, mp, Du, Dv, period,α)
        C, D, γ = get_envelope_speed(u,sp)
        push!(envelope,[abs(D), abs(γ), 2^(α/2)])
        println(length(envelope),"%")
    end

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Diffusion mode α",ylabel="Envelope exponent γ",xticks = 1.4:0.1:2)
    limits!(ax,1.4,2,0.4,1.2)
    scatter!(ax, [x[3] for x in envelope], [x[2] for x in envelope])
    # save("D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/figures/envelope for different alpha.svg",fig)
    display(fig)

    # writedlm( "D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/data/envelope for different alpha.csv",  envelope, ',')
end

function plot_wavespeed_for_different_cell_sizes(op,sp,ip,mp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx  #gridsize
    period,DuParam,DvParam=getvalues(mp)
    xs = 0:dx:lx-dx
    t = 0:dt*save_step:ntime-dt
    ΔL = nx/5
    wavespeed = []
    u=0
    Du = 0

    # 1<i < 10, smaller i is bigger cell with higher diffusion, larger i means bigger boundary with low diffusion
    for ΔD = 1:0.2:3
        fig = Figure()
        for i = 1:9
            interior_size = (0.05*i*ΔL)
            Du = Float64[ifelse((x%(ΔL) > (interior_size) && x%(ΔL) < ((1-0.05*i)*ΔL)), ΔD*DuParam, DuParam) for x ∈ 1:nx]
            Dv = Float64[ifelse((x%(ΔL) > (interior_size) && x%(ΔL) < ((1-0.05*i)*ΔL)), ΔD*DvParam, DvParam) for x ∈ 1:nx]
            u, v = numsolPDE(op, sp, ip, mp, Du, Dv, period);
            wavespeed = get_wavespeed(u,sp)
            ax = Axis(fig[div(i-1,3)+1, mod(i-1,3)+1], xlabel = "D", ylabel = "γ", title = "fraction of interior: $(round((1-2*interior_size/ΔL),digits=2))")
            scatter!(ax,abs.(wavespeed[1:end-1,2]),wavespeed[1:end-1,3],color = 1:length(wavespeed[1:end-1,3]))
            limits!(ax,10,30,0.5,1.5)
        end
        Label(fig[0, :], text = "Diffusion difference: $(round(((ΔD-1)*DuParam),digits=2))", fontsize = 12)
        save("D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/figures/wavespeed/D=$ΔD.svg",fig)
        # display(fig)    
    end
    # envolope_speed = zeros(1,3)
    # envolope_speed[:] .= get_envelope_speed(u,sp)

    # ax1 = Axis(fig[2, :][1, 1][div(i-1,3)+1, mod(i-1,3)+1][1, 2])
    # ax2 = Axis(fig[2, :][1, 1][div(i-1,3)+1, mod(i-1,3)+1][1, 1])
    # lines!(ax2,Du,xs)

    # lines!(ax1, t, A .+ B * t.^γ, linewidth=3, label="$i")
    # heatmap!(ax1, t, xs, transpose(u), colormap = :blues)
    # limits!(ax1,0,ntime,0,lx)
    # limits!(ax2,0,1.5,0,lx)
    

    # axleft =  Axis(fig[1, 1], xlabel = "D", ylabel = "x")
    # axright = Axis(fig[1, 2], xlabel = "t", yticklabelsvisible=false)
    # hm = heatmap!(axright,t,xs,transpose(u),colormap = :blues)
    # Colorbar(fig[1, 3], hm, height=Relative(1))
    # colsize!(fig.layout, 2, Auto(5))

    # limits!(axleft, 0, 1.5, 0, lx)
    # lines!(axleft,Du,xs)
    # hidedecorations!(axright,grid=false)

    # # #plot wavespeed fit
    # ts = 0:0.01:10
    # for wave in wavespeed
    #     for i=8:(size(wave,1))-2
    #         ts = 0:0.01:0.5*i
    #         lines!(axright, ts.+wave[i,1], wave[i,3] .+ wave[i,4].*ts.^wave[i,5],color=:red,linewidth=5)
    #     end
    # end

    # #plot envelope fit
    # lines!(axright, t, envolope_speed[1] .+ envolope_speed[2].*t.^envolope_speed[3])

    # # axislegend(ax,labelsize=32)
    # limits!(axright,0,ntime,0,lx)

    # save("D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/figures/wavepattern_with_sin_diffusion.png",fig)
    # axbottomleft = Axis(fig[2, 1], xlabel = "C", ylabel = "γ")
    # for i=1:length(wavespeed)
    #     # scatter!(axbottomleft,wavespeed[i][16:end-1,1],abs.(wavespeed[i][16:end-1,4]))

    #     scatter!(axbottomleft,wavespeed[i][16:end-2,4],wavespeed[i][16:end-2,5],color=16:length(wavespeed[i])-1)
    # end
    # axtopleft = Axis(fig[2, 1], xlabel = "C", ylabel = "γ")

    # # for i=1:length(wavespeed)
    # #     scatter!(axtopleft,wavespeed[i][16:end-1,1],wavespeed[i][16:end-1,5])
    # # end

    # axbottomright = Axis(fig[2, 2], xlabel = "t", ylabel = "speed (μm/s)")
    # for i=1:length(wavespeed)
    #     lines!(axbottomright,wavespeed[i][16:end-2,1],(abs.(wavespeed[i][16:end-2,4]).*10 .^wavespeed[i][16:end-2,5])/10)
    #     # lines!(axbottomright,wavespeed[i][8:end-1,1],abs.(wavespeed[i][8:end-1,4]).*1 .^wavespeed[i][8:end-1,5])
    # end
    # # linkxaxes!(axright, axbottomright)
    # # linkyaxes!(axright, axleft)
end

function plot_envelope_speed_for_different_cell_sizes(op,sp,ip,mp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx  #gridsize
    period,DuParam,DvParam=getvalues(mp)
    xs = 0:dx:lx-dx
    t = 0:dt*save_step:ntime-dt
    wavespeed = zeros(10,3)
    Du = Float32[DuParam for x ∈ 1:nx]
    Dv = Float32[DvParam for x ∈ 1:nx]
    fig = Figure()
    xs = dx/2:dx:lx-dx/2
    t = dt:dt*save_step:ntime
    for i in 1:9
        α = 2*i/18+1
        u, v = numsolfrac(op, sp, ip, mp, Du, Dv, period, α)
        ax = Axis(fig[div(i-1,3)+1, mod(i-1,3)+1], xlabel = "t", ylabel = "x", title = "α: $(round(α, digits = 2))")
        hm = heatmap!(ax,t,xs,transpose(u),colormap = :blues)

        D, C, gamma = get_envelope_speed(u,sp)
        lines!(ax,t,D .+ C * t.^gamma,color=:orange,linewidth=3,label="gamma = $(round(gamma, digits = 3))")
    
        axislegend(labelsize=10)

        limits!(ax,0,ntime,0,lx)
    end

    # fig = Figure()
    # ax = Axis(fig[1, 1], xlabel = "D",ylabel=L"\gamma")
    # hm = scatter!(ax, wavespeed[:,1], wavespeed[:,2], color=wavespeed[:,3],markersize=5.)
    # Colorbar(fig[1, 2], hm)
    # limits!(ax,0,4,0.5,1.2)
    save("D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/figures/different α.png",fig)
    display(fig)    
end

function wavespeed_for_different_α(op,sp,ip,mp)
    lx, nx, ntime, dt, save_step = getvalues(sp)
    dx = lx/nx  #gridsize
    period,DuParam,DvParam,alfa=getvalues(mp)
    xs = 0:dx:lx-dx
    t = 0:dt*save_step:ntime-dt

    Du = Float32[DuParam for x ∈ 1:nx]
    Dv = Float32[DvParam for x ∈ 1:nx]

    wavespeed = Vector{Float64}[]
    for α in 1:0.01:1.2
        u, v = numsolfrac(op, sp, ip, mp, Du, Dv, period,α)
        wave = get_wavespeed(u,sp)
        push!(wavespeed,[abs(wave[end,:][2]), abs(wave[end,:][3]), 2^(α/2)])
        println(round(length(wavespeed)/21*100,digits=0),"%")
    end

    # wavespeed = readdlm("D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/data/wavespeed for different alpha.csv",',')

    writedlm( "D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/data/wavespeed for different small alpha.csv",  wavespeed, ',')

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Diffusion mode α",ylabel="Wave speed",xticks = 1.4:0.1:2)
    limits!(ax,1.4,2,0,15)
    scatter!(ax, [x[3] for x in eachrow(wavespeed)], [x[1] for x in eachrow(wavespeed)])
    save("D:/Jakob/Documents/Fysica/MA/Thesis/Thesis/figures/wavespeed for different alpha.svg",fig)
    display(fig)

end