### CUSTOM PLOT FUNCTIONS

"""
    # plotSG(f,x,u,solver,dim::DIM,fig_name::String = "Shooting Grid")

    shooting grid plot of the state and input trajectory x::TimeSeries and u::TimeSeries 
    over the time horizon t ∈ (0,dim.t).
    the solution inbetween shooting intervals for the system funciton f is plotted by
    solving the individual initial value problems with the specified ODE solver.
"""
function plotSG(f,x,u,solver,dim::DIM,fig_name::String = "Shooting Grid")

    pl_handle =  plotGRID(x,u,dim)

    # plot points inbetween shooting intervals
    custcol = RGBA(0/255,102/255,128/255,1)   

    for k in 1:dim.N
        tt,xx = IVPt(f,x[k],u[k],solver,dim.dt,dim.dτ)
        NN = length(tt)
        tt .+= (k-1)*dim.dt
        for i in 1:dim.nx
            plot!(pl_handle[i],tt[2:end],getTS(xx,i,NN)[2:end], color = custcol, label = "")
        end
        for i in 1:dim.nu
            plot!(pl_handle[i+dim.nx],tt[[2,end]],ones(Float64,2).*u[k][i], color = custcol, label = "")
        end
    end

    plotAR(pl_handle,fig_name)
    return pl_handle
end
"""
    # plotDMS(x,u,dim::DIM,prob::DMSprob,fig_name::String = "DMS")

    shooting grid plot of the augmented state and input trajectory
    x::TimeSeries and u::TimeSeries over the time horizon t ∈ (0,dim.t).
    the solution inbetween shooting intervals for the system funciton f is plotted by
    solving the individual initial value problems.
"""
function plotDMS(x,u,dim::DIM,prob::DMSprob,fig_name::String = "DMS")

    pl_handle =  plotGRID(x,u,dim)

    # plot points inbetween shooting intervals
    custcol = RGBA(0/255,102/255,128/255,1)   

    for k in 1:dim.N
        tt,xx = IVPt(AugmentedSystem!,x[k],[u[k],prob],prob.solver,dim.dt,dim.dτ)
        NN = length(tt)
        tt .+= (k-1)*dim.dt
        for i in 1:dim.nx+1
            plot!(pl_handle[i],tt[2:end],getTS(xx,i,NN)[2:end], color = custcol, label = "")
        end
        for i in 1:dim.nu
            plot!(pl_handle[i+dim.nx+1],tt[[2,end]],ones(Float64,2).*u[k][i], color = custcol, label = "")
        end
    end

    plotAR(pl_handle,fig_name)
    return pl_handle
end
"""
    # plotGRID(x,u,dim)
    plot TimeSeries x and u on grid dim.
"""
function plotGRID(x,u,dim)
    mrksz = 3 # marker size for shooting grid points
    nx = length(x[1]) #state dimension 
    nu = length(u[1]) #input dimension

    # plot points on shooting grid
    time = collect(0:dim.dt:dim.t)
    pl_handle = Array{Any,1}(undef,nx+nu)
    for i in 1:nx
        pl_handle[i] = scatter(time,getTS(x,i,dim.N+1); label = "x$i", markersize = mrksz, markercolor = :black)
    end
    for i in 1:nu
        pl_handle[i+nx] = scatter(time,getTS(u,i,dim.N+1); label = "u$i", markersize = mrksz, markercolor = :black)
    end

    return pl_handle
end
"""
    # plotAR(pl,name::String) 
    create and show subplots.
"""
function plotAR(pl,name::String) 
    npl = length(pl)
    if npl == 6
        plot(pl[1],pl[2],pl[3],pl[4],pl[5],pl[6], layout = (6,1), link = :x, window_title = name)
    elseif npl == 5
        plot(pl[1],pl[2],pl[3],pl[4],pl[5], layout = (5,1), link = :x, window_title = name)
    end
    gui()
end
"""
    getTS(ts,ind,N)

    get variable x_ind ∈ R  of time series ts of length N.
"""
function getTS(ts,ind,N)
    return collect(ts[k][ind] for k in 1:N)
end
