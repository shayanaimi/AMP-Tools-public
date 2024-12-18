#include "snowboard.h"


double euclidian(const Eigen::Vector2d& q, const Eigen::Vector2d& p) {
			return sqrt(pow((q[0]-p[0]), 2) + pow(q[1]-p[1],2));
			// return 1;
		}

Eigen::Vector2d centroid(const amp::Obstacle2D obs)
{

	std::vector<Eigen::Vector2d> vertices = obs.verticesCCW();
	int vertexCount = vertices.size();
    Eigen::Vector2d centroid; centroid << 0, 0; 
    double signedArea = 0.0;
    double x0 = 0.0; // Current vertex X
    double y0 = 0.0; // Current vertex Y
    double x1 = 0.0; // Next vertex X
    double y1 = 0.0; // Next vertex Y
    double a = 0.0;  // Partial signed area

    // For all vertices except last
    int i=0;
    for (i=0; i<vertexCount-1; ++i)
    {
        x0 = vertices[i][0];
        y0 = vertices[i][1];
        x1 = vertices[i+1][0];
        y1 = vertices[i+1][1];
        a = x0*y1 - x1*y0;
        signedArea += a;
        centroid[0] += (x0 + x1)*a;
        centroid[1] += (y0 + y1)*a;
    }
    // Do last vertex separately to avoid performing an expensive
    // modulus operation in each iteration.
    x0 = vertices[i][0];
    y0 = vertices[i][1];
    x1 = vertices[0][0];
    y1 = vertices[0][1];
    a = x0*y1 - x1*y0;
    signedArea += a;
    centroid[0] += (x0 + x1)*a;
    centroid[1] += (y0 + y1)*a;

    signedArea *= 0.5;
    centroid[0] /= (6.0*signedArea);
    centroid[1] /= (6.0*signedArea);

    return centroid;
}

double dist_to_obs(const Eigen::Vector2d& q, const amp::Obstacle2D obs) {
			//TODO: might have to edit this so that it takes into accounts all points on the obstacle instead of just the vertices, but since it's convex I think it's mostly fine
			std::vector<Eigen::Vector2d> verts = obs.verticesCCW();
			double dist = euclidian(q, verts[0]);
			// for (Eigen::Vector2d vert : verts){
			// 	if (euclidian(q, vert) < dist){ dist=euclidian(q, vert); } 
			// }
			for (int i = 0; i < verts.size(); i++){
				Eigen::Vector2d current = verts[i];
				Eigen::Vector2d next;
				if (i == verts.size()-1){
					next = verts[0];
				}
				else{
					next = verts[i + 1];
				}
				for (double t = current[0]; t <= next[0]; t+=.1 ){
					Eigen::Vector2d point; point << current[0] + t, current[1] + t*((next[1]-current[1])/(next[0]-current[0]));
					if (euclidian(point, q) < dist){dist=euclidian(q, point);}
				}
			
		}
		return dist;
}

Eigen::Vector2d closest_vert(const Eigen::Vector2d& q, const amp::Obstacle2D obs) {
			//TODO: might have to edit this so that it takes into accounts all points on the obstacle instead of just the vertices, but since it's convex I think it's mostly fine
			std::vector<Eigen::Vector2d> verts = obs.verticesCCW();
			//double dist = euclidian(q, verts[0]);
			Eigen::Vector2d poop = verts[0];
			//Eigen::Vector2D closest = verts[0];
			// for (Eigen::Vector2d vert : verts){
			// 	if (euclidian(q, vert) < euclidian(q, poop)){ poop=vert; } 
			// }
			for (int i = 0; i < verts.size(); i++){
				Eigen::Vector2d current = verts[i];
				Eigen::Vector2d next;
				if (i == verts.size()-1){
					next = verts[0];
				}
				else{
					next = verts[i + 1];
				}
				for (double t = current[0]; t <= next[0]; t+=.1 ){
					Eigen::Vector2d point; point << current[0] + t, current[1] + t*((next[1]-current[1])/(next[0]-current[0]));
					if (euclidian(point, q) < euclidian(poop, q)){poop = point;}
				}
		
		}
		return poop;
}


Eigen::Vector2d projected_closest(const Eigen::Vector2d& q, const amp::Obstacle2D obs) {
			//TODO: might have to edit this so that it takes into accounts all points on the obstacle instead of just the vertices, but since it's convex I think it's mostly fine
			std::vector<Eigen::Vector2d> verts = obs.verticesCCW();
			//double dist = euclidian(q, verts[0]);
			Eigen::Vector2d closest =  closest_vert(q, obs);
			Eigen::Vector2d cent = centroid(obs);
			//double dist = euclid(closest, cent);
			double x_dist = cent[0] - closest[0];
			double y_dist = cent[1] - closest[1];
			Eigen::Vector2d point; point << closest[0] + x_dist*3/4, closest[0] + y_dist*3/4;
			//find the point that is exactly halfway between the centroid and closest vert
		return point;
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

Eigen::VectorXd rungeKutta4(const Eigen::VectorXd& x0, const Eigen::VectorXd& control, double dt, Eigen::VectorXd (*func)(const Eigen::VectorXd&, const Eigen::VectorXd&)) {
	Eigen::VectorXd w1 = func(x0, control);
    Eigen::VectorXd w2 = func(x0 + dt / 2.0 * w1, control);
    Eigen::VectorXd w3 = func(x0 + dt / 2.0 * w2, control);
    Eigen::VectorXd w4 = func(x0 + dt * w3, control);
    return x0 + (dt / 6.0) * (w1 + 2.0 * w2 + 2.0 * w3 + w4);
}

Eigen::VectorXd founi_propagate(Eigen::VectorXd& state, Eigen::Vector2d& control, double dt) {
    state = rungeKutta4(state, control, dt, &FirstOrderUniDynamics);
	return state;
}

Eigen::VectorXd souni_propagate(Eigen::VectorXd& state, Eigen::Vector2d& control, double dt) {
    state = rungeKutta4(state, control, dt, &SecondOrderUniDynamics);
	return state;
}

amp::KinoPath MySnowboard::souni_plan(const amp::Problem2D& problem) {
	//STATE INTEGRATOR PLANNER
    amp::Path2D path;
	amp::KinoPath kinopath;
	MySlope slope(problem);
	slope.evenMoguls(8);
    path.waypoints.push_back(problem.q_init);
    Eigen::Vector2d q = problem.q_init;
    double epsilon = .25;
	int count = 0;
    // double heading = atan2(problem.q_init[1] - problem.q_goal[1], problem.q_init[0] - problem.q_goal[0]);
	double heading = atan2(problem.q_goal[1] - problem.q_init[1], problem.q_goal[0] - problem.q_init[0]);
    while (euclidian(q, problem.q_goal) > epsilon){
        //ATTRACTIVE FORCE
		Eigen::Vector2d att;
		double euclid = euclidian(q, problem.q_goal);
		if (euclid <= d_star){
			att[0] =  zetta * (q[0] - problem.q_goal[0]);
			att[1] =  zetta * (q[1] - problem.q_goal[1]);
		}
		else{
			att[0] = d_star * zetta * (q[0] - problem.q_goal[0]) / euclid ;
			att[1] = d_star * zetta * (q[1] - problem.q_goal[1]) / euclid ;
		}
		//REPULSIVE FORCE
		
		Eigen::Vector2d rep; rep << 0, 0;
		for (mogul mog: slope.moguls){
			//std::cout<<"checking obs" << "\n";
			double di = slope.dist_to_mog(q, mog);
			if (di <= Q_star){
				//std::cout<<"adding to rep" << "\n";
				// rep += 0.5 * eta * pow(((1 / di)-(1/Q_star)), 2);
				
				rep[0] -= .35*(q[0] - mog.getX())*sin(di*M_PI/mog.getWidth())/di;
				rep[1] -= .35*(q[1] - mog.getY())*sin(di*M_PI/mog.getWidth())/di;
				//mog.getHeight()*cos(di*M_PI/mog.getWidth());

			}
		}

		//check if we are stuck in a local minimum
		Eigen::Vector2d prev_q = q;
		q = q - 0.1 * (att + rep);
		if (euclidian(q, prev_q) < 0.03){
			Eigen::Vector2d random_step;
			random_step[0] = ((double)rand() / RAND_MAX - 0.5) * 0.1; // Small random step in x
			random_step[1] = ((double)rand() / RAND_MAX - 0.5) * 0.1; // Small random step in y
			q += random_step;
		}

		count++;
		if (count == 10000
		){
			std::cout<<"too many iterations" << "\n";
			break;
		}
			
			double max_steering_change = .1; // Maximum change in steering angle in radians
			Eigen::Vector2d gradient = -att-rep;
			double gradient_heading = atan2(gradient[1], gradient[0]);
			double heading_diff = heading - gradient_heading;
			while (heading_diff > M_PI){
				heading_diff -= 2 * M_PI;
			}
			while (heading_diff < -M_PI){
				heading_diff += 2 * M_PI;
			}

			if (abs(heading_diff) > max_steering_change){
				// Limit the change in steering angle to max_steering_change
				// This is done to prevent the robot from making too sharp of turns
				// We calculate the new heading by adding the maximum allowed steering change
				// to the current heading, in the direction of the gradient
				double steering_change = max_steering_change * (heading_diff / abs(heading_diff));
				gradient_heading = heading - steering_change;
				
			}
			

			//for first order uni, controls are velocity and steering angle
			// Eigen::Vector2d control; control << log(gradient.norm()+1)*13, gradient_heading;
			Eigen::Vector2d control; control << gradient.norm()*10, gradient_heading; //change this to be steering angle rate, not steering angle
			double stepSize = .2;
			Eigen::VectorXd state(3);
			state << q[0], q[1], heading;
			state = souni_propagate(state, control, stepSize);
			q = state.head<2>();
			heading = gradient_heading;
			// heading = gradient_heading;
			// Add the new position to the waypoints
			path.waypoints.push_back(q);
			kinopath.controls.push_back(control);
			kinopath.durations.push_back(stepSize);
			kinopath.waypoints.push_back(q);
    }
    
	for (int i = 1; i < path.waypoints.size(); i++){
		path.waypoints[i] += Eigen::Vector2d(-.1, -.1);
	}
	
    path.waypoints.push_back(problem.q_goal);
    return kinopath;
}

amp::Path2D MySnowboard::founi_plan(const amp::Problem2D& problem) {
	//STATE INTEGRATOR PLANNER
    amp::Path2D path;
	MySlope slope(problem);
	slope.evenMoguls(8);
    path.waypoints.push_back(problem.q_init);
    Eigen::Vector2d q = problem.q_init;
    double epsilon = .25;
	int count = 0;
    // double heading = atan2(problem.q_init[1] - problem.q_goal[1], problem.q_init[0] - problem.q_goal[0]);
	double heading = atan2(problem.q_goal[1] - problem.q_init[1], problem.q_goal[0] - problem.q_init[0]);
    while (euclidian(q, problem.q_goal) > epsilon){
        //ATTRACTIVE FORCE
		Eigen::Vector2d att;
		double euclid = euclidian(q, problem.q_goal);
		if (euclid <= d_star){
			att[0] =  zetta * (q[0] - problem.q_goal[0]);
			att[1] =  zetta * (q[1] - problem.q_goal[1]);
		}
		else{
			att[0] = d_star * zetta * (q[0] - problem.q_goal[0]) / euclid ;
			att[1] = d_star * zetta * (q[1] - problem.q_goal[1]) / euclid ;
		}
		//REPULSIVE FORCE
		
		Eigen::Vector2d rep; rep << 0, 0;
		for (mogul mog: slope.moguls){
			//std::cout<<"checking obs" << "\n";
			double di = slope.dist_to_mog(q, mog);
			if (di <= Q_star){
				//std::cout<<"adding to rep" << "\n";
				// rep += 0.5 * eta * pow(((1 / di)-(1/Q_star)), 2);
				
				rep[0] -= .35*(q[0] - mog.getX())*sin(di*M_PI/mog.getWidth())/di;
				rep[1] -= .35*(q[1] - mog.getY())*sin(di*M_PI/mog.getWidth())/di;
				//mog.getHeight()*cos(di*M_PI/mog.getWidth());

			}
		}

		//check if we are stuck in a local minimum
		Eigen::Vector2d prev_q = q;
		q = q - 0.1 * (att + rep);
		if (euclidian(q, prev_q) < 0.03){
			Eigen::Vector2d random_step;
			random_step[0] = ((double)rand() / RAND_MAX - 0.5) * 0.1; // Small random step in x
			random_step[1] = ((double)rand() / RAND_MAX - 0.5) * 0.1; // Small random step in y
			q += random_step;
		}

		count++;
		if (count == 10000
		){
			std::cout<<"too many iterations" << "\n";
			break;
		}
			
			double max_steering_change = .1; // Maximum change in steering angle in radians
			Eigen::Vector2d gradient = -att-rep;
			double gradient_heading = atan2(gradient[1], gradient[0]);
			double heading_diff = heading - gradient_heading;
			while (heading_diff > M_PI){
				heading_diff -= 2 * M_PI;
			}
			while (heading_diff < -M_PI){
				heading_diff += 2 * M_PI;
			}

			if (abs(heading_diff) > max_steering_change){
				// Limit the change in steering angle to max_steering_change
				// This is done to prevent the robot from making too sharp of turns
				// We calculate the new heading by adding the maximum allowed steering change
				// to the current heading, in the direction of the gradient
				double steering_change = max_steering_change * (heading_diff / abs(heading_diff));
				gradient_heading = heading - steering_change;
				
			}
			

			//for first order uni, controls are velocity and steering angle
			// Eigen::Vector2d control; control << log(gradient.norm()+1)*13, gradient_heading;
			Eigen::Vector2d control; control << gradient.norm()*13, gradient_heading;
			double stepSize = .1;
			Eigen::VectorXd state(3);
			state << q[0], q[1], heading;
			state = founi_propagate(state, control, stepSize);
			q = state.head<2>();
			heading = gradient_heading;
			// heading = gradient_heading;
			// Add the new position to the waypoints
			path.waypoints.push_back(q);
    }
    
	for (int i = 1; i < path.waypoints.size(); i++){
		// path.waypoints[i] += Eigen::Vector2d(.1, .1);
	}
	
    path.waypoints.push_back(problem.q_goal);
    return path;
}

//**
// * @brief STATE INTEGRATOR PLANNER
amp::Path2D MySnowboard::plan(const amp::Problem2D& problem) {
	//STATE INTEGRATOR PLANNER
    amp::Path2D path;
	MySlope slope(problem);
	slope.evenMoguls(8);
    path.waypoints.push_back(problem.q_init);
    Eigen::Vector2d q = problem.q_init;
    double epsilon = .25;
	int count = 0;
    // double heading = atan2(problem.q_init[1] - problem.q_goal[1], problem.q_init[0] - problem.q_goal[0]);
	double heading = atan2(problem.q_goal[1] - problem.q_init[1], problem.q_goal[0] - problem.q_init[0]);
    while (euclidian(q, problem.q_goal) > epsilon){
        //ATTRACTIVE FORCE
		Eigen::Vector2d att;
		double euclid = euclidian(q, problem.q_goal);
		if (euclid <= d_star){
			att[0] =  zetta * (q[0] - problem.q_goal[0]);
			att[1] =  zetta * (q[1] - problem.q_goal[1]);
		}
		else{
			att[0] = d_star * zetta * (q[0] - problem.q_goal[0]) / euclid ;
			att[1] = d_star * zetta * (q[1] - problem.q_goal[1]) / euclid ;
		}
		//REPULSIVE FORCE
		Eigen::Vector2d rep; rep << 0, 0;
		for (mogul mog: slope.moguls){
			//std::cout<<"checking obs" << "\n";
			double di = slope.dist_to_mog(q, mog);
			if (di <= Q_star){
				//std::cout<<"adding to rep" << "\n";
				// rep += 0.5 * eta * pow(((1 / di)-(1/Q_star)), 2);
				
				rep[0] -= .35*(q[0] - mog.getX())*sin(di*M_PI/mog.getWidth())/di;
				rep[1] -= .35*(q[1] - mog.getY())*sin(di*M_PI/mog.getWidth())/di;
				//mog.getHeight()*cos(di*M_PI/mog.getWidth());

			}
		}

		//check if we are stuck in a local minimum
		Eigen::Vector2d prev_q = q;
		q = q - 0.1 * (att + rep);
		if (euclidian(q, prev_q) < 0.03){
			Eigen::Vector2d random_step;
			random_step[0] = ((double)rand() / RAND_MAX - 0.5) * 0.1; // Small random step in x
			random_step[1] = ((double)rand() / RAND_MAX - 0.5) * 0.1; // Small random step in y
			q += random_step;
		}

		count++;
		if (count == 1000
		){
			std::cout<<"too many iterations" << "\n";
			break;
		}

			double max_steering_change = .3; // Maximum change in steering angle in radians
			Eigen::Vector2d gradient = -att-rep;
			double gradient_heading = atan2(gradient[1], gradient[0]);
			double heading_diff = heading - gradient_heading;
			while (heading_diff > M_PI){
				heading_diff -= 2 * M_PI;
			}
			while (heading_diff < -M_PI){
				heading_diff += 2 * M_PI;
			}

			if (abs(heading_diff) > max_steering_change){
				// Limit the change in steering angle to max_steering_change
				// This is done to prevent the robot from making too sharp of turns
				// We calculate the new heading by adding the maximum allowed steering change
				// to the current heading, in the direction of the gradient
				double steering_change = max_steering_change * (heading_diff / abs(heading_diff));
				gradient_heading = heading - steering_change;
				
			}
			
			//now we need to use the new gradient heading as the control input and propagate that control for a small amount of time
			Eigen::Vector2d control; control << cos(gradient_heading), sin(gradient_heading);
			double stepSize = 0.2;
			q = q + stepSize * control;
			heading = gradient_heading;
			// Add the new position to the waypoints
			path.waypoints.push_back(q);
    }
    
	for (int i = 1; i < path.waypoints.size(); i++){
		path.waypoints[i] += Eigen::Vector2d(.1, .1);
	}
	
    path.waypoints.push_back(problem.q_goal);
    return path;
}


bool mogulCollision(mogul m, std::vector<mogul> moguls){
	for (int i = 0; i < moguls.size(); i++){
		if (euclidian(moguls[i].getPos(), m.getPos()) < m.width/2 + moguls[i].width/2){
			return true;
		}
	}
	return false;
}
void MySlope::randomize_moguls(int n){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 10); //change to problem bounds

	for (int i = 0; i < n; i++){
		double x = dis(gen);
		double y = dis(gen);
		mogul m(x, y);
		m.height = 4;
		m.width = 3;
		if (!mogulCollision(m, moguls)){
			moguls.push_back(m);
		}
	}
	
}

void MySlope::evenMoguls(int n){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 10); //change to problem bounds
	//n rows and n columns
	//calculate mogul widths so that they fill the space
	//mogul width should be width of slope / number of moguls
	double slope_width = x_max - x_min;
	double mwidth = slope_width / 8*1.4;
	
	std::cout << "mogul width: " << mwidth << "\n";
	//space them evenly throughout the space, with each row staggered by 1/2 width

	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			double x = mwidth/2 + j*mwidth;
			double y = mwidth/2 + i*mwidth;
			mogul m(x, y);
			std::uniform_real_distribution<> dis_height(.2, .5);
			m.height = dis_height(gen);
			m.width = mwidth;
			if (!mogulCollision(m, moguls)){
				moguls.push_back(m);
			}
		}
		
	}
	if (!moguls.empty()) {
		moguls.pop_back();
		moguls.erase(moguls.begin());
	}
	std::cout << "moguls size: " << moguls.size() << "\n";
}


double MySlope::operator()(const Eigen::Vector2d& q) const {
			
			double d_star = 3; //CHANGE ALL OF THESE TO PARAMETERS
			double zetta = .1;
			double Q_star = .5;
			double eta = 0.5;

			//ATTRACTIVE FORCE (downhill)
			double att = 0;
			double euclid = euclidian(q, problem.q_goal);
			if (euclid <= .2){
				att = zetta * pow(euclid, 2);
			}
			else{
				att = d_star * zetta * euclid - 0.5 * zetta * pow(d_star, 2);
			}

			//REPULSIVE FORCE (moguls)
			//mogul location: 5, 5
			//mogul height: 4
			//mogul width: 6

			double rep = 0;
			for (mogul mog: moguls){
				//std::cout<<"checking obs" << "\n";
				double di = dist_to_mog(q, mog);
				if (di <= Q_star){
					//std::cout<<"adding to rep" << "\n";
					// rep += 0.5 * eta * pow(((1 / di)-(1/Q_star)), 2);
					rep += mog.getHeight()*cos(di*M_PI/mog.getWidth());

				}
			}
            return att + rep;
        }

