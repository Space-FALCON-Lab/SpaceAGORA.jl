"""
Modules to find the minimum number of helper satellites required to achieve a specified Δv
on a target satellite within a given time frame using open-cavity laser propulsion.
"""

"""
Helper function that can be used to build a constellation like Main6 with N helpers evenly spaced + 1 target

    inputs:
        N_helpers: Number of helper satellites (Int)
        R_EARTH: Earth radius [m]
        helper_alt_m: Altitude of helper satellites [m]
        target_alt_m: Altitude of target satellite [m]
        cavity_B: Cavity gain parameter [unitless]
        cavity_Pin: Input laser power to cavity [W]

    outputs:
        oe: orbital elements tuple array
        Pm: power matrix (single-pass) zeros
        cav: cavity dictionary for open-cavity links
        helper_num: number of helper satellites
"""
function build_constellation(N_helpers::Int;
    R_EARTH::Float64 = Main.R_EARTH,
    helper_alt_m::Float64 = 300e3,
    target_alt_m::Float64 = 2000e3,
    cavity_B::Float64 = 100.0, #2000
    cavity_Pin::Float64 = 1e4) #1e6

    # helpers
    helper_oe = [
        (a_m=R_EARTH+helper_alt_m, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=(360.0 / N_helpers) * (i - 1))
        for i in 1:N_helpers
    ]
    # target
    target_orbit = (a_m=R_EARTH+target_alt_m, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=1.0)
    oe = vcat(helper_oe, [target_orbit])

    # power matrix (single-pass) default 0
    Pm = zeros(length(oe), length(oe))
    # open-cavity dictionary: helpers -> target
    cav = Dict{Tuple{Int, Int}, Dict{Symbol, Any}}()
    for i in 1:N_helpers
        cav[(i, N_helpers + 1)] = Dict(:B => cavity_B, :Pin => cavity_Pin)
    end

    return (; oe, Pm, cav, helper_num = N_helpers)
end

"""
Run a single scenario (N helpers) and check if Δv target is achieved within T_req

    inputs:
        N_helpers: Number of helper satellites (Int)
        dv_req: Required Δv on target satellite [m/s]
        T_req: Time limit to achieve Δv [s]
        dv_mode: mode for Δv calculation (:component or :total) (default: :component)
        dv_axis: axis for Δv calculation (:R, :T, :N) (default: :R)

    outputs:
        achieved: Bool indicating if Δv was achieved within T_req
        sol: solution object from the simulation
        p: parameter dictionary used in the simulation
        logs: diagnostic logs from the simulation
"""
function can_achieve_dv(N_helpers::Int, dv_req::Float64, T_req::Float64;
    build = build_constellation,
    mass_kg::Float64 = 227.0,
    use_los::Bool = true,
    helper_alt_m::Float64 = 300e3,
    target_alt_m::Float64 = 2000e3,
    cavity_B::Float64 = 100.0,
    cavity_Pin::Float64 = 1e4,
    min_range::Float64 = 0.0,
    max_range::Float64 = 3_000_000.0,
    dv_mode::Symbol = :component,
    dv_axis::Symbol = :R,
    img_dir::String = "Kuang's Prototype Code/output/images/",
    verbose::Bool = false)

    constel = build(N_helpers, helper_alt_m=helper_alt_m, target_alt_m=target_alt_m, cavity_B=cavity_B, cavity_Pin=cavity_Pin)
    oe, Pm, cav, helper_num = constel.oe, constel.Pm, constel.cav, constel.helper_num

    # Solve up to T_req or stop early at dv_req
    sol, p, logs = run_open_cavity_multi(oe;
        mass_kg = mass_kg,
        Pm = Pm,
        cavity = cav,
        use_los = use_los,
        min_range = min_range, max_range = max_range,
        stop_on_dv = true,
        dv_target_sat = helper_num + 1,
        dv_target_value = dv_req,
        dv_target_mode = dv_mode,
        dv_target_axis = dv_axis,
        T_seconds = T_req,
        verbose = verbose, result_plots = false,
        IMG_DIR = img_dir,
        helper_num = helper_num
    )

    # Success heuristic: callback terminated before hard end-time
    # If you set save_positions=(true,true) in the callback, sol.t[end] is the event time, 
    # which is the time from start to when dv_req is met
    achieved = sol.t[end] < T_req - 1e-9
    return achieved, sol, p, logs
end

# --------- Feasibility: RTN-R requirement ---------
function can_achieve_r(N_helpers::Int, r_req_m::Float64, T_req::Float64;
    build = build_constellation,
    mass_kg::Float64 = 227.0,
    use_los::Bool = true,
    helper_alt_m::Float64 = 300e3,
    target_alt_m::Float64 = 2000e3,
    cavity_B::Float64 = 100.0,
    cavity_Pin::Float64 = 1e4,
    min_range::Float64 = 0.0,
    max_range::Float64 = 3_000_000.0,
    img_dir::String = "Kuang's Prototype Code/output/images/",
    verbose::Bool = false)

    constel = build(N_helpers, helper_alt_m=helper_alt_m, target_alt_m=target_alt_m, cavity_B=cavity_B, cavity_Pin=cavity_Pin)
    oe, Pm, cav, helper_num = constel.oe, constel.Pm, constel.cav, constel.helper_num

    # Solve up to T_req or stop early when R-component reaches r_req_m
    sol, p, logs = run_open_cavity_multi(oe;
        mass_kg = mass_kg,
        Pm = Pm,
        cavity = cav,
        use_los = use_los,
        min_range = min_range, max_range = max_range,
        stop_on_dv = false,
        stop_on_r = true,
        r_target_sat = helper_num + 1,
        r_target_value = r_req_m,
        r_target_mode = :component,
        r_target_axis = :R,
        T_seconds = T_req,
        verbose = verbose, result_plots = false,
        IMG_DIR = img_dir,
        helper_num = helper_num
    )

    # Success heuristic: callback terminated before hard end-time
    achieved = sol.t[end] < T_req - 1e-9
    return achieved, sol, p, logs
end

# function can_achieve_dv(N_helpers::Int, dv_req::Float64, T_req::Float64;
    #                         dv_mode::Symbol = :component, dv_axis::Symbol = :R)

    #     constel = build_constellation(N_helpers)
    #     oe, Pm, cav, helper_num = constel.oe, constel.Pm, constel.cav, constel.helper_num

    #     # Solve up to T_req or stop early at dv_req
    #     sol, p, logs = run_open_cavity_multi(oe;
    #         mass_kg=227,
    #         Pm=Pm,
    #         cavity=cav,
    #         use_los=true,
    #         min_range=0.0, max_range=3000e3,
    #         stop_on_dv=true,
    #         dv_target_sat=helper_num + 1,
    #         dv_target_value=dv_req,
    #         dv_target_mode=dv_mode,
    #         dv_target_axis=dv_axis,
    #         T_seconds=T_req,
    #         verbose=false, result_plots=false,
    #         IMG_DIR="Kuang's Prototype Code/output/images/",
    #         helper_num=helper_num
    #     )

    #     # Success heuristic: callback terminated before hard end-time
    #     # If you set save_positions=(true,true) in the callback, sol.t[end] is the event time, 
    #     # which is the time from start to when dv_req is met
    #     achieved = sol.t[end] < T_req - 1e-9 #small tolerance
    #     return achieved, sol, p, logs
    # end


"""
Find the minimum helpers via exponential + binary search:
    First use exponential search to find an upper bound where dv_req is achievable,
    then use binary search between the last insufficient and first sufficient N to find the minimum.

    inputs:
        dv_req: Required Δv on target satellite [m/s]
        T_req: Time limit to achieve Δv [s]
        N_max: Maximum number of helpers to consider (default: 128)
        N_min: Minimum number of helpers to consider (default: 1)
        Other keyword arguments are passed to can_achieve_dv
    outputs:
        If achievable within N_max:
            A named tuple with fields:
                N: minimum number of helpers (Int)
                sol: solution object from the simulation
                p: parameter dictionary used in the simulation
                logs: diagnostic logs from the simulation
        If not achievable within N_max:
            nothing
"""
function min_helpers_for_dv(dv_req::Float64, T_req::Float64;
    N_max::Int = 128,
    N_min::Int = 1,
    build = build_constellation,
    mass_kg::Float64 = 227.0,
    use_los::Bool = true,
    helper_alt_m::Float64 = 300e3,
    target_alt_m::Float64 = 2000e3,
    cavity_B::Float64 = 100.0,
    cavity_Pin::Float64 = 1e4,
    min_range::Float64 = 0.0,
    max_range::Float64 = 3_000_000.0,
    dv_mode::Symbol = :component,
    dv_axis::Symbol = :R,
    img_dir::String = "Kuang's Prototype Code/output/images/",
    verbose::Bool = false)

    # Step 1: Quick check at N_min
    ok, sol_best, p_best, logs_best = can_achieve_dv(N_min, dv_req, T_req;
        build=build, mass_kg=mass_kg, use_los=use_los,
        helper_alt_m=helper_alt_m, target_alt_m=target_alt_m,
        cavity_B=cavity_B, cavity_Pin=cavity_Pin,
        min_range=min_range, max_range=max_range,
        dv_mode=dv_mode, dv_axis=dv_axis,
        img_dir=img_dir, verbose=verbose)
    if ok
        return (N=N_min, sol=sol_best, p=p_best, logs=logs_best) # early return if N_min suffices
    end

    # Step 2: Exponential search for upper bound
    lo = N_min
    hi = min(2*lo, N_max)
    while true # Step 3: loop until a return happens
        ok, sol_hi, p_hi, logs_hi = can_achieve_dv(hi, dv_req, T_req;
            build=build, mass_kg=mass_kg, use_los=use_los,
            helper_alt_m=helper_alt_m, target_alt_m=target_alt_m,
            cavity_B=cavity_B, cavity_Pin=cavity_Pin,
            min_range=min_range, max_range=max_range,
            dv_mode=dv_mode, dv_axis=dv_axis,
            img_dir=img_dir, verbose=verbose) # step 4， 7: check at hi
        if ok # Step 8: if hi is sufficient
            # Binary search in (lo, hi] for minimum N
            best = (N=hi, sol=sol_hi, p=p_hi, logs=logs_hi)
            while lo + 1 < hi
                mid = (lo + hi) >>> 1 #mid = Int(floor((lo + hi) / 2)) # Step 6: take midpoint
                okm, sol_m, p_m, logs_m = can_achieve_dv(mid, dv_req, T_req;
                    build=build, mass_kg=mass_kg, use_los=use_los,
                    helper_alt_m=helper_alt_m, target_alt_m=target_alt_m,
                    cavity_B=cavity_B, cavity_Pin=cavity_Pin,
                    min_range=min_range, max_range=max_range,
                    dv_mode=dv_mode, dv_axis=dv_axis,
                    img_dir=img_dir, verbose=verbose) # Step 9: check at mid
                if okm ##### this if/else is the binary search logic (Part 2) #####
                    hi = mid # move upper bound down
                    best = (N=mid, sol=sol_m, p=p_m, logs=logs_m) # update best
                else
                    lo = mid # move lower bound up
                end
            end
            return best # Step 10: return best found
        else
            lo = hi # Step 5: move lower bound up and keep while looping ##### this is the exponential search logic (Part 1) #####
            hi = min(hi * 2, N_max)   # Step 6: double upper bound, capped at N_max and keep while looping
            if hi == lo || hi > N_max # if cannot increase further
                return nothing  # cannot achieve within N_max
            end
        end
    end
end



# Find minimum number of helpers to achieve R requirement within T_req (binary search)
 function min_helpers_for_r(r_req_m::Float64, T_req::Float64;
     build = build_constellation,
     mass_kg::Float64 = 227.0,
     use_los::Bool = true,
     helper_alt_m::Float64 = 300e3,
     target_alt_m::Float64 = 2000e3,
     cavity_B::Float64 = 100.0,
     cavity_Pin::Float64 = 1e4,
     min_range::Float64 = 0.0,
     max_range::Float64 = 3_000_000.0,
     img_dir::String = "Kuang's Prototype Code/output/images/",
     verbose::Bool = false,
    N_lo::Int = 1,
    N_hi::Int = 200)
 
    # Start from provided lower bound
    lo = N_lo
    # Quick check at lower bound
    achieved_lo, sol_lo, p_lo, logs_lo = can_achieve_r(lo, r_req_m, T_req; build=build, mass_kg=mass_kg, use_los=use_los,
         helper_alt_m=helper_alt_m, target_alt_m=target_alt_m, cavity_B=cavity_B, cavity_Pin=cavity_Pin,
         min_range=min_range, max_range=max_range, img_dir=img_dir, verbose=verbose)
    println("Checking lower bound N = $lo: achieved = $achieved_lo")
    if achieved_lo
        return (N=lo, sol=sol_lo, p=p_lo, logs=logs_lo)
     end

    # Prepare placeholders for best found solution
    bestN = N_hi
    best_sol = nothing
    best_p = nothing
    best_logs = nothing

    # Exponential search to find the first achievable upper bound quickly
    hi = min(max(2 * lo, 2), N_hi)
    while true
        achieved_hi, sol_hi, p_hi, logs_hi = can_achieve_r(hi, r_req_m, T_req; build=build, mass_kg=mass_kg, use_los=use_los,
            helper_alt_m=helper_alt_m, target_alt_m=target_alt_m, cavity_B=cavity_B, cavity_Pin=cavity_Pin,
            min_range=min_range, max_range=max_range, img_dir=img_dir, verbose=verbose)
        println("Checking upper candidate N = $hi: achieved = $achieved_hi")
        if achieved_hi
            # Found achievable upper bound; carry best details as initial best
            bestN = hi
            best_sol = sol_hi
            best_p = p_hi
            best_logs = logs_hi
            break
        else
            lo = hi
            hi = min(hi * 2, N_hi)
            if hi == lo
                return nothing
            end
            if hi >= N_hi && !achieved_hi
                # Final check at N_hi already done; not achievable within bounds
                achieved_cap, _, _, _ = can_achieve_r(N_hi, r_req_m, T_req; build=build, mass_kg=mass_kg, use_los=use_los,
                    helper_alt_m=helper_alt_m, target_alt_m=target_alt_m, cavity_B=cavity_B, cavity_Pin=cavity_Pin,
                    min_range=min_range, max_range=max_range, img_dir=img_dir, verbose=verbose)
                if !achieved_cap
                    return nothing
                end
            end
        end
    end
 
     # Binary search for minimal helpers
    # Continue with binary search between (lo, hi]
     while lo <= hi
        mid = (lo + hi) >>> 1
        ok, sol_m, p_m, logs_m = can_achieve_r(mid, r_req_m, T_req; build=build, mass_kg=mass_kg, use_los=use_los,
             helper_alt_m=helper_alt_m, target_alt_m=target_alt_m, cavity_B=cavity_B, cavity_Pin=cavity_Pin,
             min_range=min_range, max_range=max_range, img_dir=img_dir, verbose=verbose)
        println("Checking N = $mid: achieved = $ok")
        if ok
            bestN = mid
            best_sol = sol_m
            best_p = p_m
            best_logs = logs_m
             hi = mid - 1
        else
             lo = mid + 1
        end
     end
    return (N=bestN, sol=best_sol, p=best_p, logs=best_logs)
 end


#the code below also works:
# # Build a constellation like your Main6: N helpers evenly spaced + 1 target
# function build_constellation(N_helpers::Int;
#     R_EARTH::Float64 = Main.R_EARTH,
#     helper_alt_m::Float64 = 300e3,
#     target_alt_m::Float64 = 2000e3,
#     cavity_B::Float64 = 2000.0,
#     cavity_Pin::Float64 = 1e6)

#     # helpers
#     helper_oe = [
#         (a_m=R_EARTH+helper_alt_m, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=(360.0 / N_helpers) * (i - 1))
#         for i in 1:N_helpers
#     ]
#     # target
#     target_orbit = (a_m=R_EARTH+target_alt_m, e=0.0, i_deg=0.0, Ω_deg=0.0, ω_deg=0.0, ν_deg=0.0)
#     oe = vcat(helper_oe, [target_orbit])

#     # power matrix (single-pass) default 0
#     Pm = zeros(length(oe), length(oe))
#     # open-cavity dictionary: helpers -> target
#     cav = Dict{Tuple{Int, Int}, Dict{Symbol, Any}}()
#     for i in 1:N_helpers
#         cav[(i, N_helpers + 1)] = Dict(:B => cavity_B, :Pin => cavity_Pin)
#     end

#     return (; oe, Pm, cav, helper_num = N_helpers)
# end

# # Run a single scenario (N helpers) and check if Δv target is achieved within T_req
# function can_achieve_dv(N_helpers::Int, dv_req::Float64, T_req::Float64;
#     build = build_constellation,
#     mass_kg::Float64 = 227.0,
#     use_los::Bool = true,
#     min_range::Float64 = 0.0,
#     max_range::Float64 = 3000e3,
#     dv_mode::Symbol = :component,
#     dv_axis::Symbol = :R,
#     img_dir::String = "Kuang's Prototype Code/output/images/",
#     verbose::Bool = false)

#     constel = build(N_helpers)
#     oe, Pm, cav, helper_num = constel.oe, constel.Pm, constel.cav, constel.helper_num

#     # Solve up to T_req or stop early at dv_req
#     sol, p, logs = run_open_cavity_multi(oe;
#         mass_kg=mass_kg,
#         Pm=Pm,
#         cavity=cav,
#         use_los=use_los,
#         min_range=min_range, max_range=max_range,
#         stop_on_dv=true,
#         dv_target_sat=helper_num + 1,
#         dv_target_value=dv_req,
#         dv_target_mode=dv_mode,
#         dv_target_axis=dv_axis,
#         T_seconds=T_req,
#         verbose=verbose, result_plots=false,
#         IMG_DIR=img_dir,
#         helper_num=helper_num
#     )

#     # Success heuristic: callback terminated before hard end-time
#     # If you set save_positions=(true,true) in the callback, sol.t[end] is the event time.
#     achieved = sol.t[end] < T_req - 1e-9
#     return achieved, sol, p, logs
# end

# # Find the minimum helpers via exponential + binary search
# function min_helpers_for_dv(dv_req::Float64, T_req::Float64;
#     N_max::Int = 128,
#     N_min::Int = 1,
#     build = build_constellation,
#     kwargs...)

#     # Quick check at N_min
#     ok, sol_best, p_best, logs_best = can_achieve_dv(N_min, dv_req, T_req; build=build, kwargs...)
#     ok && return (N=N_min, sol=sol_best, p=p_best, logs=logs_best)

#     # Exponential search for upper bound
#     lo = N_min
#     hi = min(2*lo, N_max)
#     while true
#         ok, sol_hi, p_hi, logs_hi = can_achieve_dv(hi, dv_req, T_req; build=build, kwargs...)
#         if ok
#             # Binary search in (lo, hi]
#             best = (N=hi, sol=sol_hi, p=p_hi, logs=logs_hi)
#             while lo + 1 < hi
#                 mid = (lo + hi) >>> 1
#                 okm, sol_m, p_m, logs_m = can_achieve_dv(mid, dv_req, T_req; build=build, kwargs...)
#                 if okm
#                     hi = mid
#                     best = (N=mid, sol=sol_m, p=p_m, logs=logs_m)
#                 else
#                     lo = mid
#                 end
#             end
#             return best
#         else
#             lo = hi
#             hi = min(hi * 2, N_max)
#             if hi == lo || hi > N_max
#                 return nothing  # cannot achieve within N_max
#             end
#         end
#     end
# end
