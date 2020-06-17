
# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp
import tests as tests


# ToDo: Properly document this method.
def make_steady_validation_problem():
    """

    :return:
    """

    steady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_steady_validation_airplane()
    )
    steady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_validation_operating_point()
    )

    steady_validation_problem = asmvp.problems.SteadyProblem(
        airplane=steady_validation_airplane, operating_point=steady_validation_operating_point
    )

    del steady_validation_airplane
    del steady_validation_operating_point

    return steady_validation_problem


# ToDo: Properly document this method.
def make_unsteady_validation_problem():
    """

    :return:
    """

    unsteady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_validation_operating_point()
    )
    unsteady_validation_movement = (
        tests.integration.fixtures.movement_fixtures.make_validation_movement()
    )

    unsteady_validation_problem = asmvp.problems.UnsteadyProblem(
        airplane=unsteady_validation_airplane,
        operating_point=unsteady_validation_operating_point,
        movement=unsteady_validation_movement
    )

    del unsteady_validation_airplane
    del unsteady_validation_operating_point
    del unsteady_validation_movement

    return unsteady_validation_problem
