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

    amp::Problem2D prob;
    prob.q_goal = Eigen::Vector2d(0, 0);
    prob.q_init = Eigen::Vector2d(9, 9);
    MySlope mpf(prob);
    mpf.evenMoguls(8);
    std::vector<amp::Obstacle2D> obs;
    for (auto& mog : mpf.moguls) {
        amp::Obstacle2D obs_i;
        double angle = 0;
        double delta_angle = 2 * 3.14159 / 10;
        for (int i = 0; i < 10; i++) {
            double x = mog.getX() + mog.width/2 * cos(angle);
            double y = mog.getY() + mog.width/2 * sin(angle);
            obs_i.verticesCCW().push_back(Eigen::Vector2d(x, y));
            angle += delta_angle;
        }
        obs.push_back(obs_i);
    }
    prob.obstacles = obs;
    MySnowboard board(3, 0.1, .5, 0.5);
    amp::Path2D path = board.plan(prob);
    // void amp::Visualizer::makeFigure(const PotentialFunction2D& potential_function, const Problem2D& prob, std::size_t n_grid, bool vector, double u_min, double u_max) {
    // amp::Visualizer::makeFigure(mpf, prob.x_min, prob.x_max, prob.y_min, prob.y_max, 0, 100);
    amp::Visualizer::makeFigure(mpf, prob, 100, false, 0, 100);
    amp::Visualizer::makeFigure(prob, path);

    // static void makeFigure(const PotentialFunction2D& potential_function, double x0_min, double x0_max, double x1_min, double x1_max, double u_min = 0.0, double u_max = 100.0)
    Visualizer::showFigures();

  
    return 0;
}