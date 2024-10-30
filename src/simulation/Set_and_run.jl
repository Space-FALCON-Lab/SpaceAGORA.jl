include("../physical_models/Planet_data.jl")
include("../physical_models/Mission.jl")
include("Aerobraking.jl")
include("../config.jl")



function aerobraking_campaign(args, state)
    save_rs = args[:results]

    # Descent towards Mars
    purpose = "Aerobraking around Mars"

    mission = Dict(:Purpose => purpose, :Gravity_Model => args.gravity_model, :)
end