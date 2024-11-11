#include "MyKinoRRT.h"
#include "AMPCore.h"
#include "hw/HW2.h"
#include <stack>

Eigen::VectorXd rungeKutta4(const Eigen::VectorXd& x0, const Eigen::VectorXd& control, double dt, Eigen::VectorXd (*func)(const Eigen::VectorXd&, const Eigen::VectorXd&)) {
    Eigen::VectorXd w1 = func(x0, control);
    Eigen::VectorXd w2 = func(x0 + dt / 2.0 * w1, control);
    Eigen::VectorXd w3 = func(x0 + dt / 2.0 * w2, control);
    Eigen::VectorXd w4 = func(x0 + dt * w3, control);
    return x0 + (dt / 6.0) * (w1 + 2.0 * w2 + 2.0 * w3 + w4);
}

void MySimpleCar::HadiRK(Eigen::VectorXd& state, const Eigen::VectorXd& control, double dt) {
    Eigen::VectorXd newState(5), w1(5), w2(5), w3(5), w4(5);

    double x = state(0);
    double y = state(1);
    double theta = state(2);
    double v = state(3);
    double phi = state(4);

    double u1 = control(0);
    double u2 = control(1);

    double der1 = v * cos(theta);
    double der2 = v * sin(theta);
    double der3 = v / 5 * tan(phi);
    double der4 = u1;
    double der5 = u2;

    w1 << der1, der2, der3, der4, der5;

    w2 << w1(0) + 0.5 * dt * w1(0),
          w1(1) + 0.5 * dt * w1(1),
          w1(2) + 0.5 * dt * w1(2),
          w1(3) + 0.5 * dt * w1(3),
          w1(4) + 0.5 * dt * w1(4);

    w3 << w1(0) + 0.5 * dt * w2(0),
          w1(1) + 0.5 * dt * w2(1),
          w1(2) + 0.5 * dt * w2(2),
          w1(3) + 0.5 * dt * w2(3),
          w1(4) + 0.5 * dt * w2(4);

    w4 << w1(0) + dt * w3(0),
          w1(1) + dt * w3(1),
          w1(2) + dt * w3(2),
          w1(3) + dt * w3(3),
          w1(4) + dt * w3(4);

    newState = state + (dt / 6.0) * (w1 + 2 * w2 + 2 * w3 + w4);

    state = newState;
}

void MySimpleCar::HadiRK2(Eigen::VectorXd &state, const Eigen::VectorXd &control, double dt)
{
    Eigen::VectorXd w1(5), w2(5), w3(5), w4(5);

    double x = state(0);
    double y = state(1);
    double theta = state(2);
    double v = state(3);
    double phi = state(4);

    double u1 = control(0);
    double u2 = control(1);

    w1 << v * cos(theta),
        v * sin(theta),
        (v / 5) * tan(phi),
        u1,
        u2;

    Eigen::VectorXd state2 = state + 0.5 * dt * w1;
    w2 << state2(3) * cos(state2(2)),
        state2(3) * sin(state2(2)),
        (state2(3) / 5) * tan(state2(4)),
        u1,
        u2;

    Eigen::VectorXd state3 = state + 0.5 * dt * w2;
    w3 << state3(3) * cos(state3(2)),
        state3(3) * sin(state3(2)),
        (state3(3) / 5) * tan(state3(4)),
        u1,
        u2;

    Eigen::VectorXd state4 = state + dt * w3;
    w4 << state4(3) * cos(state4(2)),
        state4(3) * sin(state4(2)),
        (state4(3) / 5) * tan(state4(4)),
        u1,
        u2;

    state = state + (dt / 6.0) * (w1 + 2 * w2 + 2 * w3 + w4);
}
Eigen::VectorXd SingleIntegratorDynamics(const Eigen::VectorXd& state, const Eigen::VectorXd& control) {
    return control;
}

Eigen::VectorXd FirstOrderUniDynamics(const Eigen::VectorXd& state, const Eigen::VectorXd& control) {
    double r = 0.25;
    Eigen::VectorXd result(3);
    result[0] = control[0] * cos(state[2]) * r;
    result[1] = control[0] * sin(state[2]) * r;
    result[2] = control[1];
    return result;
}

Eigen::VectorXd SecondOrderUniDynamics(const Eigen::VectorXd& state, const Eigen::VectorXd& control) {
    double r = 0.25;
    Eigen::VectorXd result(5);
    result[0] = state[3] * cos(state[2]) * r;
    result[1] = state[3] * sin(state[2]) * r;
    result[2] = state[4];
    result[3] = control[0];
    result[4] = control[1];
    return result;
}

Eigen::VectorXd SimpleCarDynamics(const Eigen::VectorXd& state, const Eigen::VectorXd& control) {
    //x y theta v phi
    double L = 5.0;
    double W = 2.0;
    double x = state[0];
    double y = state[1];
    double theta = state[2];
    double v = state[3];
    double phi = state[4];
    Eigen::VectorXd result(5);
    result[0] = v * cos(theta);
    result[1] = v * sin(theta);
    result[2] = (v / L) * tan(phi);
    result[3] = control[0];
    result[4] = control[1];
    return result;

}

void MySimpleCar::propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) {
    // HadiRK(state, control, dt);
    // state = rungeKutta4(state, control, dt, &SimpleCarDynamics);
    HadiRK(state, control, dt);
}

void MySingleIntegrator::propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) {
    state = rungeKutta4(state, control, dt, &SingleIntegratorDynamics);
}

void MyFirstOrderUnicycle::propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) {
    state = rungeKutta4(state, control, dt, &FirstOrderUniDynamics);
}

void MySecondOrderUnicycle::propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) {
    state = rungeKutta4(state, control, dt, &SecondOrderUniDynamics);

}

std::vector<Eigen::Vector2d> blowUp(std::vector<Eigen::Vector2d> verts, double factor){
    std::vector<Eigen::Vector2d> result;
    Eigen::Vector2d centroid(0,0);
    for (auto p : verts){
        centroid += p;
    }
    centroid /= verts.size();
    for (auto p : verts){
        result.push_back(centroid + (p - centroid)*factor);
    }
    return result;
}

bool MyKinoRRT::inCollision(const Eigen::VectorXd& cspace_state){
    for (amp::Obstacle2D obs : myproblem.obstacles){
            // std::cout<<"checking obstacle" << "\n";
            std::vector<Eigen::Vector2d> points = blowUp(obs.verticesCCW(), 1.2);
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
        if (state.size() == 5) {
            // std::cout << "state size 5, checking " << state << "\n";

            if (state[3] >= q_goal[3].first - tolerance && state[3] <= q_goal[3].second + tolerance && state[4] >= q_goal[4].first - tolerance && state[4] <= q_goal[4].second + tolerance){
            // Handle the case where the number of state variables is 5
            // Add your specific logic here
            return true;
            }
            else { return false;}
            
        }
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

std::pair<std::shared_ptr<amp::Graph<std::pair<Eigen::VectorXd, double>>>, std::map<amp::Node, Eigen::VectorXd>> MyKinoRRT::makeRRTGraph(int N, double step_size, amp::DynamicAgent& agent){
    std::shared_ptr<amp::Graph<std::pair<Eigen::VectorXd, double>>> random_graph = std::make_shared<amp::Graph<std::pair<Eigen::VectorXd, double>>>();  
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
    std::pair<double, double> dt_bounds = myproblem.dt_bounds;
    std::vector<double> dtlower = {dt_bounds.first};
    std::vector<double> dtupper = {dt_bounds.second};
    step_size = generateRandomNVector(dtlower, dtupper)[0] / 10.0;
    // std::cout << "State lower bounds: " << qlowerbounds.transpose() << std::endl;
    // std::cout << "State upper bounds: " << qupperbounds.transpose() << std::endl;
    // std::cout << "Control lower bounds: " << ulowerbounds.transpose() << std::endl;
    // std::cout << "Control upper bounds: " << uupperbounds.transpose() << std::endl;
    points.push_back(q_init);
    nodes[0] = q_init;
    int goalBiasCtr = 0;
    int timesDitched = 0;
    int timesDitchedatGoalBias = 0;
    int actualIterations = 0;
    std::cout << "beginning iterations: " << std::endl;
    for (int i = 1; i <= N; i++) {
        if (nodes.size() < 5){
            std::cout << "Number of nodes: " << nodes.size() << std::endl;
        }
        
        actualIterations += 1;
        if (actualIterations >= 150000){
            std::cout << "exit with failure" << std::endl;
            break;
        }
        if (actualIterations % 10000 == 0){
            std::cout << "Iteration: " << i << std::endl;
        }
       
        Eigen::VectorXd sample = generateRandomNVector(eigenToStdVector(qlowerbounds), eigenToStdVector(qupperbounds));
        if (goalBiasCtr % 20 == 0){ // with a 0.5 chance
            sample = goal_centroid;
        }
        goalBiasCtr = goalBiasCtr + 1;

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
        int m = myM;
        Eigen::VectorXd bestcontrol =  generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
        Eigen::VectorXd bestmid_point = closest_point;
        agent.propagate(bestmid_point, bestcontrol, step_size);

        double min_dist = (sample - bestmid_point).norm();
        for (int j = 0; j < m; j++){
            Eigen::VectorXd control =  generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
            Eigen::VectorXd mid_point = closest_point;
            agent.propagate(mid_point, control, step_size);
            if (inCollision(mid_point, closest_point)){
                continue;
            }
            double distance = (sample - mid_point).norm();
            if (distance < min_dist) {
                min_dist = distance;
                bestcontrol = control;
                bestmid_point = mid_point;
            }
        }
        Eigen::VectorXd control =  bestcontrol;
        Eigen::VectorXd mid_point = bestmid_point;
        // agent.propagate(mid_point, control, step_size);
        // std::cout << "state[2]: " << mid_point[2] << ", state[4]: " << mid_point[4] << std::endl;
        int failcount = 0;
        Eigen::VectorXd last_try = bestmid_point;
        bool ditch = false;

        while (inCollision(mid_point, closest_point)) {
            failcount++;
            if (failcount > 100){
                //basically try this all again but sampling the goal
                sample = goal_centroid;
                // while (inCollision(sample)) {
                //     sample = generateRandomNVector(eigenToStdVector(qlowerbounds), eigenToStdVector(qupperbounds));
                // }
                closest_point = nodes[0];
                index_counter = 0;
                closest_index = 0;
                min_distance = std::numeric_limits<double>::max(); // initialize the minimum distance to positive infinity
                for (const auto& point : points) { // iterate over all points that have been sampled so far
                    double dist = (sample - point).norm(); //this will have to change when we add dimensions
                    if (dist < min_distance) {
                        min_distance = dist;
                        closest_point = point;
                        closest_index = index_counter;
                    }
                    index_counter ++;
                }
                control =  generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
                mid_point = closest_point;
                agent.propagate(mid_point, control, step_size);
                // std::cout << "hacky sampled goal" << "\n";
                // std::cout << "sampled goal" << "\n";
            }
            else{
                control = generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
            // std::cout << control.transpose() << "\n";
                mid_point = closest_point;
                agent.propagate(mid_point, control, step_size);
            }            
            
            if (failcount > 90){
                if ((last_try-mid_point).norm() < 0.1){
                    //then ditch this whole sample
                    // std::cout << "ditching 1" << "\n";
                    timesDitched = timesDitched + 1;
                    if ((goalBiasCtr-1) % 20 == 0){
                        timesDitchedatGoalBias = timesDitchedatGoalBias + 1;
                    }
                    ditch = true;
                    break;
                }
            }
            last_try = mid_point;
        }
        if   (ditch){
            // std::cout << "ditching 2" << "\n";
            i = i - 1;
            continue;
        }    
        points.push_back(mid_point);
        nodes[i] = mid_point;        
        //this will change
        // double distance_to_sample = (mid_point - closest_point).norm();
        random_graph->connect(closest_index, i, std::make_pair(control, step_size)); // Connect the edges in the graph
        // if (!inCollision(mid_point, closest_point)){ //kind of redundant because we already collision checked but whatever
            
        //     points.push_back(mid_point);
        //     nodes[i] = mid_point;        
        //     //this will change
        //     // double distance_to_sample = (mid_point - closest_point).norm();
        //     random_graph->connect(closest_index, i, std::make_pair(control, step_size)); // Connect the edges in the graph
        //     // std::cout << "adding edge to graph from " << closest_index << " to " << i << "\n";
        //     i+=1;
           
        // }
        // else{
        //     continue;
        // }


         if ((inGoalRegion(mid_point, q_goal))){
                std::cout << "Found goal at iteration " << i << "\n";
                break;
            }
    }


    std::cout << "returning graph \n";
    std::cout << "times ditched: " << timesDitched << "\n";
    std::cout << "times ditched at goal bias: " << timesDitchedatGoalBias << "\n";
    return std::make_pair(random_graph, nodes);
}

Eigen::VectorXd getCtrl(amp::Node current, amp::Node came_from, std::shared_ptr<amp::Graph<std::pair<Eigen::VectorXd, double>>> graph){
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
        const auto& edge = edges[child_index].first;
        return edge;
        // Further code to use 'edge' as needed
    }
    // std::cout << "no edge from " << came_from << " to " << current << "\n";
    return Eigen::VectorXd::Zero(1);
}

double getDuration(amp::Node current, amp::Node came_from, std::shared_ptr<amp::Graph<std::pair<Eigen::VectorXd, double>>> graph){
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
        const auto& edge = edges[child_index].second;
        return edge;
        // Further code to use 'edge' as needed
    }
    // std::cout << "no edge from " << came_from << " to " << current << "\n";
    return 0.1;
}

amp::KinoPath MyKinoRRT::plan(const amp::KinodynamicProblem2D& problem, amp::DynamicAgent& agent) {
    //check what agent it is
    std::cout <<  "agent dim length " << agent.agent_dim.length << "\n";
    double step_size;
    if (problem.q_bounds.size() == 2){
        step_size = 0.1;
    }
    else{

        step_size = 0.05;
    } //if it is 5d (and maybe also 3d) then play with random step sizes from 0.01-0.1
    std::cout << "step size can vary between " << myproblem.dt_bounds.first << " and " << myproblem.dt_bounds.second << ". Using " << step_size << "\n";
    step_size = 0.05;
    int N = 50000;
    myproblem = problem;
    amp::KinoPath path;
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
    std::cout << "nodes size: " << nodes.size() << "\n";
    
    double tol = 1;
    while (!stack.empty()) {   
        amp::Node current = stack.top();
        // std::cout << "current: " << nodes[current] << "\n";
        stack.pop();
        // double dist_to_goal = (nodes[current]-nodes[goal_node]).norm();
        // std::cout << "dist_to_goal: " << dist_to_goal << "\n";
        if (current == goal_node ) {
            // std::cout << "found goal near " << nodes[current] << "\n";
            std::cout << "goal node: " << nodes[goal_node] << "\n";
            // path.waypoints.push_back(nodes[goal_node].head<2>());
            while (current != init_node) {
                plength = plength + (nodes[current].head<2>() - nodes[came_from[current]].head<2>()).norm();
                path.waypoints.push_back(nodes[current]);
                path.controls.push_back(getCtrl(current, came_from[current], random_graph));
                path.durations.push_back(getDuration(current, came_from[current], random_graph));
                current = came_from[current];
            }
            
            path.waypoints.push_back(nodes[init_node]); //where hadi's is fucking up

            std::reverse(path.waypoints.begin(), path.waypoints.end());
            std::reverse(path.controls.begin(), path.controls.end());
            std::reverse(path.durations.begin(), path.durations.end());
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
        
    }
    std::cout << "path length: " << plength << "\n";
    std::cout << "\n";
    path.valid = true;
    return path;
}


// while (current != init_node) {
//                 if (current == goal_node) {
//                     plength = plength + (nodes[current].head<2>() - nodes[came_from[current]].head<2>()).norm();
//                     current = came_from[current];
//                     auto ctrl = getCtrl(current, came_from[current], random_graph);
//                     auto dur = getDuration(current, came_from[current], random_graph);
//                     path.controls.push_back(ctrl);
//                     path.durations.push_back(dur);
//                     continue;
//                 }
//                 auto ctrl = getCtrl(current, came_from[current], random_graph);
//                 path.controls.push_back(ctrl);
//                 // path.waypoints.emplace_back(nodes[current].head<2>());
//                 path.waypoints.emplace_back(nodes[current]);
//                 path.durations.push_back(getDuration(current, came_from[current], random_graph));

//                 // std::cout << "pushing waypt: " << nodes[current] << "\n";
//                 plength = plength + (nodes[current] - nodes[came_from[current]]).norm();
//                 current = came_from[current];
//             }



//GRAVEYARD

// while (inCollision(mid_point, closest_point)) {
//             failcount++;
//             if (failcount > 100){
//                 //basically try this all again but sampling the goal
//                 sample = goal_centroid;
//                 // while (inCollision(sample)) {
//                 //     sample = generateRandomNVector(eigenToStdVector(qlowerbounds), eigenToStdVector(qupperbounds));
//                 // }
//                 closest_point = nodes[0];
//                 index_counter = 0;
//                 closest_index = 0;
//                 min_distance = std::numeric_limits<double>::max(); // initialize the minimum distance to positive infinity
//                 for (const auto& point : points) { // iterate over all points that have been sampled so far
//                     double dist = (sample - point).norm(); //this will have to change when we add dimensions
//                     if (dist < min_distance) {
//                         min_distance = dist;
//                         closest_point = point;
//                         closest_index = index_counter;
//                     }
//                     index_counter ++;
//                 }
//                 control =  generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
//                 mid_point = closest_point;
//                 agent.propagate(mid_point, control, step_size);
//                 // std::cout << "hacky sampled goal" << "\n";
//                 // std::cout << "sampled goal" << "\n";
//             }
//             else{
//                 control = generateRandomNVector(eigenToStdVector(ulowerbounds), eigenToStdVector(uupperbounds));
//             // std::cout << control.transpose() << "\n";
//                 mid_point = closest_point;
//                 agent.propagate(mid_point, control, step_size);
//             }            
            
//             if (failcount > 90){
//                 if ((last_try-mid_point).norm() < 0.1){
//                     //then ditch this whole sample
//                     // std::cout << "ditching 1" << "\n";
//                     timesDitched = timesDitched + 1;
//                     if ((goalBiasCtr-1) % 20 == 0){
//                         timesDitchedatGoalBias = timesDitchedatGoalBias + 1;
//                     }
//                     ditch = true;
//                     break;
//                 }
//             }
//             last_try = mid_point;
//         }
//         if   (ditch){
//             // std::cout << "ditching 2" << "\n";
//             i = i - 1;
//             continue;
//         }  