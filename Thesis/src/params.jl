Base.@kwdef struct OscillatorParam
    K::Int64
    kval::Int64
    bval::Float64
    ϵ::Int64
end
getvalues(p::OscillatorParam) = p.K,p.kval,p.bval,p.ϵ

Base.@kwdef struct SimParam
    lx::Int64           #Domain size
    nx::Int64           #number of grid points
    time::Float64       #total simulation time
    dt::Float64
    save_step::Int64
end
getvalues(p::SimParam) = p.lx,p.nx,p.time,p.dt,p.save_step

Base.@kwdef struct InitParam
    u0::Float64
    v0::Float64
    Pi::Float64
    Po::Float64
    pacemaker::Float64  #pacemaker size
end
getvalues(p::InitParam) = p.u0,p.v0,p.Pi,p.Po,p.pacemaker

Base.@kwdef struct ModelParam
    period::Float64
    Du::Float64         #Diffusion const
    Dv::Float64         #Diffusion const
    α::Float64
end
getvalues(p::ModelParam) = p.period,p.Du,p.Dv,p.α