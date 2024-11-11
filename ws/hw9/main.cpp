// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"
#include "hw/HW9.h"
#include "hw/HW2.h"
#include "MyKinoRRT.h"

using namespace amp;

// Load problems and map agent for quick testing
std::vector<KinodynamicProblem2D> problems = {HW9::getStateIntProblemWS1(), HW9::getStateIntProblemWS2(), HW9::getFOUniProblemWS1(), HW9::getFOUniProblemWS2(), HW9::getSOUniProblemWS1(), HW9::getSOUniProblemWS2(), HW9::getCarProblemWS1(), HW9::getParkingProblem()};
std::unordered_map<AgentType, std::function<std::shared_ptr<amp::DynamicAgent>()>> agentFactory = {
    {AgentType::SingleIntegrator, []() { return std::make_shared<MySingleIntegrator>(); }},
    {AgentType::FirstOrderUnicycle, []() { return std::make_shared<MyFirstOrderUnicycle>(); }},
    {AgentType::SecondOrderUnicycle, []() { return std::make_shared<MySecondOrderUnicycle>(); }},
    {AgentType::SimpleCar, []() { return std::make_shared<MySimpleCar>(); }}
};

//get the path length for a kino path
double getPathLength(KinoPath path) {
    double length = 0;
    for (int i = 1; i < path.waypoints.size(); i++) {
        length += (path.waypoints[i].head<2>() - path.waypoints[i-1].head<2>()).norm();
    }
    return length;
}

int main(int argc, char** argv) {
    // Select problem, plan, check, and visualize
    int select = 2;
    KinodynamicProblem2D prob = problems[select];
    MyKinoRRT kino_planner;
    // KinoPath path = kino_planner.plan(prob, *agentFactory[prob.agent_type]());
    // HW9::check(path, prob);
    // if (path.valid)
    //     Visualizer::makeFigure(prob, path, false); // Set to 'true' to render animation
    // Visualizer::showFigures();
    // HW9::grade<MyKinoRRT, MySingleIntegrator, MyFirstOrderUnicycle, MySecondOrderUnicycle, MySimpleCar>("shaya.naimi@colorado.edu", argc, argv, std::make_tuple(), std::make_tuple(), std::make_tuple(), std::make_tuple(), std::make_tuple());


    std::list <std::vector<double>> all_computation_times;
    std::list <std::vector<double>> all_path_lengths;
    std::vector<double> all_valid_sols;
    bool bench = true;
    
    if (bench) {
        for (int N : {50000}){
            for (int M: {1, 5, 10, 15}){
            std::cout << "N: " << N << " M: " << M << std::endl;
            
            std::vector<double> path_lengths;
            int valid_sols = 0;
            std::vector<double> computation_times;
            kino_planner.myN = N;
            kino_planner.myM = M;
            for (int i = 0; i < 50; ++i) {
                    // std::cout<< "iteration: " << i << "\n";
                    auto start = std::chrono::high_resolution_clock::now();
                    
                    KinoPath path = kino_planner.plan(prob, *agentFactory[prob.agent_type]());
                    
                    auto stop = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
                    computation_times.push_back(duration.count());
                    double plength = getPathLength(path);
                    
                    if (HW9::check(path, prob, false)) {
                        valid_sols++;
                        path_lengths.push_back(plength);
                    }


                }
                all_computation_times.push_back(computation_times);
                all_path_lengths.push_back(path_lengths);
                all_valid_sols.push_back(valid_sols);
            
        }
    }
    //const std::list<std::vector<double>> values = all_path_lengths;    
    std::string xlabel = "Hyperparameters N and M";
    std::string ylabel = "Computation Time (microseconds)";
    std::string title = "Kino RRT Computation Time (HW2WS1)";
    std::vector<std::string> labels = {"n=50000, M=1", "n=50000, M=5", "n=50000, M=10", "n=50000, M=15"};
    Visualizer::makeBoxPlot(all_computation_times, labels, 
                               title, xlabel, ylabel);  
    xlabel = "Hyperparameters N and M";
    ylabel = "Path Length";
    title = "Kino RRT PathLength (HW2WS1)";
    // labels = {"n=200, r=0.5","n=200, r=1.0","n=200, r=1.5", "n=200, r=2.0", "n=500, r=0.5","n=500, r=1.0","n=500, r=1.5", "n=500, r=2.0"};
    labels = {"n=50000, M=1", "n=50000, M=5", "n=50000, M=10", "n=50000, M=15"};
    Visualizer::makeBoxPlot(all_path_lengths, labels, 
                               title, xlabel, ylabel); 
    xlabel = "Hyperparameters N and M";
    ylabel = "Number of Valid Solutions out of 50 runs";
    title = "Kino RRT Number of Valid Solutions (HW2WS1)";
    // labels = {"n=200, r=0.5","n=200, r=1.0","n=200, r=1.5", "n=200, r=2.0", "n=500, r=0.5","n=500, r=1.0","n=500, r=1.5", "n=500, r=2.0"};
    labels = {"n=50000, M=1", "n=50000, M=5", "n=50000, M=10", "n=50000, M=15"};
    Visualizer::makeBarGraph(all_valid_sols, labels, 
                               title, xlabel, ylabel); 

    Visualizer::showFigures();
    }

    


    return 0;
}