module NavigationEffectors
    export calcNavigationEffect!

    """
        calcNavigationEffect!(model, u, p, t, sat_idx)

    Hook point for navigation/sensor-estimator models executed by typed periodic callbacks.
    Extend this method for your navigation effector type.
    """
    function calcNavigationEffect!(model, u, p, t::Float64, sat_idx::Int)
        throw(MethodError(calcNavigationEffect!, (model, u, p, t, sat_idx)))
    end
end
