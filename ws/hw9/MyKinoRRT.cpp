#include "MyKinoRRT.h"
#include "AMPCore.h"
#include "hw/HW2.h"
#include <stack>

void MySingleIntegrator::propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) {
    state += dt * control;
};

void MyFirstOrderUnicycle::propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) {
    double r = 0.25;
    state[0] = state[0] + dt * control[0] * cos(state[2]) * r;
    state[1] = state[1] + dt * control[0] * sin(state[2]) * r;
    state[2] = state[2] + dt * control[1];
}

void MySecondOrderUnicycle::propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) {
    double r = 0.25;
    state[0] = state[0] + dt * state[3] * cos(state[2]) * r;
    state[1] = state[1] + dt * state[3] * sin(state[2]) * r;
    state[2] = state[2] + dt * state[4];
    state[3] = state[3] + dt * control[0];
    state[4] = state[4] + dt * control[1];

}

bool MyKinoRRT::inCollision(const Eigen::VectorXd& cspace_state){
    for (amp::Obstacle2D obs : myproblem.obstacles){
            // std::cout<<"checking obstacle" << "\n";
            std::vector<Eigen::Vector2d> points = obs.verticesCCW();
            int k, m, nvert = points.size();
            bool c = false;
            for(k = 0, m = nvert - 1; k < nvert; m = k++) {
                if( ( (points[k][1] >= cspace_state[1] ) != (points[m][1] >= cspace_state[1]) ) &&
                    (cspace_state[0] <= (points[m][0] - points[k][0]) * (cspace_state[1] - points[k][1]) / (points[m][1] - points[k][1]) + points[k][0])
                )   
                c = !c;
                }
                if (c){
                    return true;
            }
        }
        return false;
}

bool MyKinoRRT::inCollision(const Eigen::VectorXd& startpoint, const Eigen::VectorXd& endpoint){
    for (double t = .05; t <= 1; t+= 0.1){
        Eigen::VectorXd point = startpoint + (endpoint - startpoint) * t;
        if (inCollision(point)){
            return true;
        }
    }
    return false;
}

Eigen::VectorXd generateRandomNVector(const std::vector<double> &lower_bounds, const std::vector<double> &upper_bounds)
{

    if (lower_bounds.size() != upper_bounds.size())
    {
        throw std::invalid_argument("Erororororororrr");
    }
    // std::vector<double> lower_bounds_copy = {lower_bounds[0], -3};
    // std::vector<double> upper_bounds_copy = {upper_bounds[0], 3};
    // int n = lower_bounds_copy.size();
    int n = lower_bounds.size();
    Eigen::VectorXd random_vector(n);
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < n; ++i)
    {
        std::uniform_real_distribution<float> dist(lower_bounds[i], upper_bounds[i]);
        random_vector[i] = dist(gen);
    }
    // std::cout << "random vector: " << random_vector.transpose() << "\n";
    return random_vector;
}

std::vector<double> eigenToStdVector(const Eigen::VectorXd& vec){
    std::vector<double> std_vec;
    std_vec.resize(vec.size());
    Eigen::VectorXd::Map(&std_vec[0], vec.size()) = vec;
    return std_vec;
}

bool inGoalRegion(Eigen::VectorXd& state, std::vector<std::pair<double, double>>& q_goal){
    double tolerance = 0;
    if (state[0] >= q_goal[0].first - tolerance && state[0] <= q_goal[0].second + tolerance && state[1] >= q_goal[1].first - tolerance && state[1] <= q_goal[1].second + tolerance){
        return true;
    } else {
        return false;
    }
}   

std::pair<Eigen::VectorXd, Eigen::VectorXd> extractBounds(std::vector<std::pair<double, double>> bounds){
    Eigen::VectorXd firsts(bounds.size());
    Eigen::VectorXd seconds(bounds.size());
    for (int i = 0; i < bounds.size(); i++){
        firsts[i] = bounds[i].first;
        seconds[i] = bounds[i].second;
    }
    return std::make_pair(firsts, seconds);
}

std::pair<std::shared_ptr<amp::Graph<Eigen::VectorXd>>, std::map<amp::Node, Eigen::VectorXd>> MyKinoRRT::makeRRTGraph(int N, double step_size, amp::DynamicAgent& agent){
    std::shared_ptr<amp::Graph<Eigen::VectorXd>> random_graph = std::make_shared<amp::Graph<Eigen::VectorXd>>();  
    //this graph might need to have edges of controls instead of doubles  
    std::map<amp::Node, Eigen::VectorXd> nodes;
    std::vector<Eigen::VectorXd> points;

    Eigen::VectorXd q_init = myproblem.q_init;
    //Goal region bounding each dimension of the state (vector<pair<min, max>>)
    std::vector<std::pair<double, double>> q_goal = myproblem.q_goal;
    Eigen::Vector2d goal_centroid;
    goal_centroid[0] = (q_goal[0].first + q_goal[0].second) / 2.0;
    goal_centroid[1] = (q_goal[1].first + q_goal[1].second) / 2.0;
    auto [qlowerbounds, qupperbounds] = extractBounds(myproblem.q_bounds);
    auto [ulowerbounds, uupperbounds] = extractBounds(myproblem.u_bounds);
    std::cout << "Control lower bounds: " << ulowerbounds.transpose() << std::endl;
    std::cout << "Control upper bounds: " << uupperbounds.transpose() << std::endl;
    points.push_back(q_init);
    nodes[0] = q_init;
    int goalBiasCtr = 0;
    for (int i = 1; i < N; i++){
        Eigen::VectorXd sample = generateRandomNVector(eigenToStdVector(qlowerbounds), eigenToStdVector(qupperbounds));
        if (goalBiasCtr % 20 == 0){ // with a 0.5 chance
            sample = goal_centroid;
        }
        goalBiasCtr = goalBiasCtr + 1;

        while (inCollision(sample)) {
            // std::cout << "can't find a valid sample" << "\n";
            sample = generateRandomNVector(eigenToStdVector(qlowerbounds), eigenToStdVector(qupperbounds));
            //std::cout << "can't find a valid sample" << "\n";
        }

        Eigen::VectorXd closest_point = nodes[0];
        int index_counter = 0;
        int closest_index = 0;
        double min_distance = std::numeric_limits<double>::max(); // initialize the minimum distance to positive infinity
        for (const auto& point : points) { // iterate over all points that have been sampled so far
            double distance = (sample - point).norm(); //this will have to change when we add dimensions
            if (distance < min_distance) {
                min_distance = distance;
                closest_point = point;
                closest_index = index_counter;
            }
            index_counter ++;
        }

        //sample a random control
        Eigen::VectorXd control =  generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
        Eigen::VectorXd mid_point = closest_point;
        agent.propagate(mid_point, control, step_size);
        int failcount = 0;
        Eigen::VectorXd last_try = mid_point;
        bool ditch = false;
        while (inCollision(mid_point, closest_point)) {
            failcount++;
            // std::cout << mid_point.transpose() << "\n";
            control = generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
            // std::cout << control.transpose() << "\n";
            mid_point = closest_point;
            
            agent.propagate(mid_point, control, step_size);
            if (failcount > 10){
                if ((last_try-mid_point).norm() < 0.1){
                    //then ditch this whole sample
                    std::cout << "ditching 1" << "\n";
                    ditch = true;
                    break;
                }
            }
            last_try = mid_point;
        }
        if   (ditch){
            std::cout << "ditching 2" << "\n";
            i = i - 1;
            continue;
        }  
        points.push_back(mid_point);
        nodes[i] = mid_point;        
        //this will change
        double distance_to_sample = (mid_point - closest_point).norm();
        random_graph->connect(closest_index, i, control); // Connect the edges in the graph
        // std::cout << "adding edge to graph from " << closest_index << " to " << i << "\n";
        
        if ((inGoalRegion(mid_point, q_goal))){
            std::cout << "Found goal at midpoint: " << mid_point.transpose() << "\n";
            break;
        }
    }


    std::cout << "returning graph \n";
    return std::make_pair(random_graph, nodes);
}

Eigen::VectorXd getCtrl(amp::Node current, amp::Node came_from, std::shared_ptr<amp::Graph<Eigen::VectorXd>> graph){
    const auto& children = graph->children(came_from);
    const auto& edges = graph->outgoingEdges(came_from);

    int child_index = -1;
    for (int i = 0; i < children.size(); ++i) {
        if (children[i] == current) {
            child_index = i;
            break;
        }
    }

    if (child_index != -1 && child_index < edges.size()) {
        const auto& edge = edges[child_index];
        return edge;
        // Further code to use 'edge' as needed
    }
    // std::cout << "no edge from " << came_from << " to " << current << "\n";
    return Eigen::VectorXd::Zero(1);
}

amp::KinoPath MyKinoRRT::plan(const amp::KinodynamicProblem2D& problem, amp::DynamicAgent& agent) {
    double step_size = 0.5;
    int N = 5000;
    myproblem = problem;
    amp::KinoPath path;
    Eigen::VectorXd state = problem.q_init;
    //MAKE GRAPH
    auto [random_graph, nodes] = makeRRTGraph(N, step_size, agent);

    //GRAPH SEARCH
    amp::Node init_node = 0; //maybe try printing these
    amp::Node goal_node = nodes.size() - 1;  //not sure if this will work
    std::stack<amp::Node> stack;
    std::unordered_map<amp::Node, amp::Node> came_from;
    stack.push(init_node);
    came_from[init_node] = init_node;
    double plength = 0;
    // std::cout << "graph is created, starting breadth-first search" << "\n";
    // std::cout << "is reversible: "  << random_graph->isReversible() << "\n";
    double tol = 1;
    while (!stack.empty()) {   
        amp::Node current = stack.top();
        // std::cout << "current: " << nodes[current] << "\n";
        stack.pop();
        double dist_to_goal = (nodes[current]-nodes[goal_node]).norm();
        // std::cout << "dist_to_goal: " << dist_to_goal << "\n";
        if (current == goal_node ) {
            // std::cout << "found goal near " << nodes[current] << "\n";
            path.waypoints.push_back(nodes[goal_node].head<2>());
            while (current != init_node) {
                if (current == goal_node) {
                    plength = plength + (nodes[current] - nodes[came_from[current]]).norm();
                    current = came_from[current];
                    auto ctrl = getCtrl(current, came_from[current], random_graph);
                    path.controls.push_back(ctrl);
                    path.durations.push_back(step_size);
                    continue;
                }
                auto ctrl = getCtrl(current, came_from[current], random_graph);
                path.controls.push_back(ctrl);
                path.waypoints.emplace_back(nodes[current].head<2>());
                path.durations.push_back(step_size);

                // std::cout << "pushing waypt: " << nodes[current] << "\n";
                plength = plength + (nodes[current] - nodes[came_from[current]]).norm();
                current = came_from[current];
            }
            path.waypoints.push_back(nodes[init_node]); //where hadi's is fucking up

            std::reverse(path.waypoints.begin(), path.waypoints.end());
            std::reverse(path.controls.begin(), path.controls.end());
            // std::cout << "Path1: ";
            // for (const auto& waypoint : path.waypoints) {
            //     std::cout << waypoint.transpose() << " ";
            //     }
            // break;
            // return path;
        }
        // std::cout << "number of children: " << random_graph->children(current).size() << "\n";
        for (const auto& neighbor : random_graph->children(current)) {
            // std::cout << "child: " << nodes[neighbor] << "\n";
            if (came_from.find(neighbor) == came_from.end()) {
                // std::cout << "pushing it" << "\n";
                stack.push(neighbor);
                came_from[neighbor] = current;
                // std::cout << "stack size: " << stack.size() << "\n";
                
            }
        }
        // for (const auto& neighbor : random_graph->parents(current)) {
        //     // std::cout << "parent: " << nodes[neighbor] << "\n";
        //     if (came_from.find(neighbor) == came_from.end()) {
        //         // std::cout << "pushing it" << "\n";
        //         stack.push(neighbor);
        //         came_from[neighbor] = current;
        //         // std::cout << "stack size: " << stack.size() << "\n";
                
        //     }
        // }
    }

    
    // path.waypoints.push_back(state);
    // for (int i = 0; i < 10; i++) {
    //     Eigen::VectorXd control = Eigen::VectorXd::Random(problem.q_init.size());
    //     agent.propagate(state, control, 1.0);
    //     path.waypoints.push_back(state);
    //     path.controls.push_back(control);
    //     path.durations.push_back(1.0);
    // }
    // std::cout << "Path controls: ";
    // for (const auto& control : path.controls) {
    //     std::cout << control.transpose() << " ";
    // }
    std::cout << "\n";
    path.valid = true;
    return path;
}
