"""This script replicates the simulations shown in figure 13.34 of Katz and Plotkin."""

import pterasoftware as ps
import matplotlib.pyplot as plt
import numpy as np

# These are constants defined by Katz and Plotkin that are specific to the validation simulation.
alpha = 5
num_chordwise_panels = 4
num_spanwise_panels = 13
time_const = 1 / 16
max_num_time_consts = 9

# The original plot shows results for multiple aspect ratios. We will just plot the result for one.
AR = 4

# The particular speed and chord length don't matter because the results are nondimensionalized.
u_inf = 1
c = 1

# Calculate the proper wingspan, time step length, and number of time steps based on Katz and Plotkin's constants.
b = AR * c
delta_time = time_const * c / u_inf
num_steps = int(c * max_num_time_consts / (u_inf * delta_time))

# Define an airplane object that matches the figure's simulations.
airplane = ps.geometry.Airplane(
    name="Plate Airplane",
    wings=[
        ps.geometry.Wing(
            name="Plate Wing",
            symmetric=False,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    chord=c,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=b,
                    z_le=0.0,
                    chord=c,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)

root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=airplane.wings[0].wing_cross_sections[0],
)
tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=airplane.wings[0].wing_cross_sections[1],
)

wing_movement = ps.movement.WingMovement(
    base_wing=airplane.wings[0],
    wing_cross_sections_movements=[
        root_wing_cross_section_movement,
        tip_wing_cross_section_movement,
    ],
)

del root_wing_cross_section_movement
del tip_wing_cross_section_movement

airplane_movement = ps.movement.AirplaneMovement(
    base_airplane=airplane,
    wing_movements=[wing_movement],
)

del wing_movement

operating_point = ps.operating_point.OperatingPoint(
    density=1.225,
    velocity=u_inf,
    alpha=alpha,
    nu=15.06e-6,
)

operating_point_movement = ps.movement.OperatingPointMovement(
    base_operating_point=operating_point,
)

movement = ps.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    num_steps=num_steps,
    delta_time=delta_time,
)

del airplane_movement
del operating_point_movement

problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    unsteady_problem=problem,
)

del problem

solver.run(
    logging_level="Warning",
    prescribed_wake=True,
    calculate_streamlines=False,
)

figure_background_color = "None"
text_color = "#818181"
lift_color = "#EDAD08"
marker_size = 8

x_min = 0.00
x_max = 9.00
y_min = -0.1
y_max = 0.55

final_time_step_time = delta_time * (num_steps - 1)

times = np.linspace(
    0,
    final_time_step_time,
    num_steps,
    endpoint=True,
)

lift_coefficients = np.zeros(num_steps)

for step in range(0, num_steps):
    airplane = solver.steady_problems[step].airplanes[0]

    lift_coefficients[step] = airplane.total_near_field_force_coefficients_wind_axes[2]

figure, axes = plt.subplots()

axes.spines.right.set_visible(False)
axes.spines.top.set_visible(False)
axes.spines.bottom.set_color(text_color)
axes.spines.left.set_color(text_color)
axes.xaxis.label.set_color(text_color)
axes.yaxis.label.set_color(text_color)
axes.tick_params(axis="x", colors=text_color)
axes.tick_params(axis="y", colors=text_color)
figure.patch.set_facecolor(figure_background_color)
axes.set_facecolor(figure_background_color)

axes.plot(
    (times * u_inf) / c,
    lift_coefficients,
    label="AR = " + str(AR),
    color=lift_color,
    marker=".",
    markersize=marker_size,
)

axes.set_xlim(x_min, x_max)
axes.set_ylim(y_min, y_max)
axes.set_xlabel("Time Constants", color=text_color)
axes.set_ylabel("Lift Coefficient", color=text_color)
axes.legend(
    facecolor=figure_background_color,
    edgecolor=figure_background_color,
    labelcolor=text_color,
)
axes.axhline(linewidth=1, linestyle="--", color=text_color)

figure.savefig(
    "Transient Lift Coefficient.png",
    dpi=300,
)
figure.show()

# Uncomment these lines if you'd like to see the classic simulation outputs (plotting all forces and coefficients,
# drawing the airplane on the last time step, and animating the airplane).

# ps.output.plot_results_versus_time(
#     solver=solver,
#     show=True,
#     save=False,
# )
#
# ps.output.draw(
#     solver=solver,
#     scalar_type="lift",
#     show_streamlines=False,
#     show_wake_vortices=False,
#     save=False,
# )
#
# ps.output.animate(
#     solver=solver,
#     scalar_type="lift",
#     show_wake_vortices=True,
#     save=False,
# )
