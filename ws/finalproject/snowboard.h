#pragma once

// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"

class mogul {
	public:
		mogul(double x, double y) : x(x), y(y) {}
		double x, y;
		double height;
		double width;
		std::vector<double> xbounds = {0, 10};
		std::vector<double> ybounds = {0, 10};
		double getHeight() { return height; }
		double getWidth() { return width; }
		double getX() { return x; }
		double getY() { return y; }
		Eigen::Vector2d getPos() { return Eigen::Vector2d(x, y); }
		
};

class MySnowboard : public amp::GDAlgorithm {
	public:
		// Consider defining a class constructor to easily tune parameters, for example: 
		MySnowboard(double d_star, double zetta, double Q_star, double eta) :
			d_star(d_star),
			zetta(zetta),
			Q_star(Q_star),
			eta(eta) {}

		// Override this method to solve a given problem.
		virtual amp::Path2D plan(const amp::Problem2D& problem) override;
	private:
		double d_star, zetta, Q_star, eta;
		// Add additional member variables here...
};

class MySlope : public amp::PotentialFunction2D {
    public:
		        MySlope(const amp::Problem2D& prob) : problem(prob) {
            x_max = prob.x_max;
            x_min = prob.x_min;
            y_max = prob.y_max;
            y_min = prob.y_min;

			//problem = prob;
		}
		// Returns the potential function value (height) for a given 2D point. 

        virtual double operator()(const Eigen::Vector2d& q) const override;

		double euclidian(const Eigen::Vector2d& q, const Eigen::Vector2d& p) const{
			return sqrt(pow((q[0]-p[0]), 2) + pow(q[1]-p[1],2));
			// return 1;
		}

		double dist_to_obs(const Eigen::Vector2d& q, const amp::Obstacle2D obs) const{
			//TODO: might have to edit this so that it takes into accounts all points on the obstacle instead of just the vertices, but since it's convex I think it's mostly fine
			std::vector<Eigen::Vector2d> verts = obs.verticesCCW();
			double dist = euclidian(q, verts[0]);
			for (Eigen::Vector2d vert : verts){
				if (euclidian(q, vert) < dist){ dist=euclidian(q, vert); } 
			}
			return dist;
		}

		double dist_to_mog(const Eigen::Vector2d& q,  mogul mog) const{
			double dist = euclidian(q, Eigen::Vector2d(mog.getX(), mog.getY()));
			return dist;
		}
		virtual Eigen::Vector2d getGradient(const Eigen::Vector2d& q) const override {
            return Eigen::Vector2d(q[0] * q[0],  q[1] * q[1]);
        }
		std::vector<mogul> moguls;
		void add_moguls(std::vector<mogul> moguls) { this->moguls = moguls; }
		void randomize_moguls(int n);
		void evenMoguls(int n);
		double x_min, x_max, y_min, y_max;
	private: 
		const amp::Problem2D& problem;
};


