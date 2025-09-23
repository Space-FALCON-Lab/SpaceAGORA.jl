using WGLMakie
using Dash, DashHtmlComponents, DashCoreComponents
using PlotlyJS
using CSV
using DataFrames


function start_viz_dashboard(args,space_object_dict)
    app = dash(external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"])
    # print(space_object_dict[:sc_state_history][1])
    # print(space_object_dict[1].sc_state_history[2])


    app.layout = html_div() do
        html_h1("SpaceAGORA Dashboard"),
        html_div("Dash.jl: Julia interface for Dash"),
        dcc_graph(
            id = "trajectory-3d",
            figure = PlotlyJS.plot(
                [
                    scatter3d(
                        x = [state[1] for state in sc.sc_state_history],
                        y = [state[2] for state in sc.sc_state_history],
                        z = [state[3] for state in sc.sc_state_history],
                        mode = "lines+markers",
                        marker = attr(size=3),
                        line = attr(width=2),
                        name = sc.uid
                    ) for sc in values(space_object_dict)
                ],
                Layout(
                    title = "Space Object Trajectories",
                    scene = attr(
                        xaxis_title = "X",
                        yaxis_title = "Y",
                        zaxis_title = "Z"
                    )
                )
            )
        )
    end

    run_server(app)
end

# function start_viz_dashboard_from_csv(csv_path)
#     df = CSV.read(csv_path, DataFrame)
#     app = dash(external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"])

#     app.layout = html_div() do
#         html_h1("Hello Dash (CSV Data)"),
#         html_div("Dash.jl: Julia interface for Dash, CSV input"),
#         dash_graph(
#             id = "trajectory-3d-csv",
#             figure = PlotlyJS.plot(
#                 [
#                     scatter3d(
#                         x = df[:, 2],
#                         y = df[:, 3],
#                         z = df[:, 4],
#                         mode = "lines+markers",
#                         marker = attr(size=3, color="red"),
#                         line = attr(width=2, color="red"),
#                         name = "Trajectory (CSV)"
#                     )
#                 ],
#                 Layout(
#                     title = "Space Object Trajectory (CSV)",
#                     scene = attr(
#                         xaxis_title = "X",
#                         yaxis_title = "Y",
#                         zaxis_title = "Z"
#                     )
#                 )
#             )
#         )
#     end

#     run_server(app)
# end
# start_viz_dashboard_from_csv("-1_state_data.csv")
# Example main function to call either dashboard