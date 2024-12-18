#pragma once
#include "AMPCore.h"

class MyMACollChecker : public amp::ConfigurationSpace {
    public:
        MyMACollChecker(const Eigen::VectorXd& lower_bounds, const Eigen::VectorXd& upper_bounds)
            : ConfigurationSpace(lower_bounds, upper_bounds) {}

        // Implement the pure virtual method
        virtual bool inCollision(const Eigen::VectorXd& cspace_state) const override;
        bool inCollision(amp::MultiAgentProblem2D problem, const Eigen::VectorXd& cspace_state) const;
        bool inCollision(amp::MultiAgentProblem2D problem, const Eigen::VectorXd& startpoint, const Eigen::VectorXd& endpoint) const;
        void set_problem(amp::MultiAgentProblem2D& problem)  {
            myproblem = problem;
        }
        bool diskCollision(Eigen::Vector2d p1, Eigen::Vector2d p2, double radius) const;
        amp::MultiAgentProblem2D myproblem; //change this to multiagentproblem
        double robot_radius;
        double cautious_radius = 1;
        std::vector<double> agent_radii;
        void set_robot_radius(double radius) {robot_radius = radius;}
        void set_cautious_radius(double radius) {cautious_radius = radius;}
        bool circlePoly(amp::Obstacle2D obs, Eigen::Vector2d point, double radius) const;
        void set_agent_radii(std::vector<double> radii)  {agent_radii = radii;}
        double get_agent_radius(int i) const {return agent_radii[i];}



    private:
        
        
};