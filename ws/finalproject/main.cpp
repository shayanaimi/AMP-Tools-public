// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"
#include "hw/HW2.h"

// Include any custom headers you created in your workspace
#include "snowboard.h"

using namespace amp;

int main(int argc, char** argv) {

    // Visualize your potential function
    amp::Problem2D prob = HW2::getWorkspace1();
    MySlope mpf(prob);
    mpf.evenMoguls(8);
    // void amp::Visualizer::makeFigure(const PotentialFunction2D& potential_function, const Problem2D& prob, std::size_t n_grid, bool vector, double u_min, double u_max) {
    // amp::Visualizer::makeFigure(mpf, prob.x_min, prob.x_max, prob.y_min, prob.y_max, 0, 100);
    amp::Visualizer::makeFigure(mpf, prob, 100, false, 0, 100);

    // static void makeFigure(const PotentialFunction2D& potential_function, double x0_min, double x0_max, double x1_min, double x1_max, double u_min = 0.0, double u_max = 100.0)
    Visualizer::showFigures();

  
    return 0;
}