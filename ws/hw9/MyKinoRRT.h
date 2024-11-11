#pragma once

// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework headers
#include "hw/HW9.h"

class MyKinoRRT : public amp::KinodynamicRRT {
    public:
        virtual amp::KinoPath plan(const amp::KinodynamicProblem2D& problem, amp::DynamicAgent& agent) override;
        std::pair<std::shared_ptr<amp::Graph<std::pair<Eigen::VectorXd, double>>>, std::map<amp::Node, Eigen::VectorXd>> makeRRTGraph(int N, double step_size, amp::DynamicAgent& agent);
        bool inCollision(const Eigen::VectorXd& cspace_state);
        bool inCollision(const Eigen::VectorXd& startpoint, const Eigen::VectorXd& endpoint);
        int myN = 10000;
        int myM = 5;
    private:
        amp::KinodynamicProblem2D myproblem = amp::HW9::getStateIntProblemWS1();
};  

class MySingleIntegrator : public amp::DynamicAgent {
    public:
        virtual void propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) override;
        double agent_step_size = 0.1;        
};

class MyFirstOrderUnicycle : public amp::DynamicAgent {
    public:
        virtual void propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) override;
};

class MySecondOrderUnicycle : public amp::DynamicAgent {
    public:
        virtual void propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) override;
};

class MySimpleCar : public amp::DynamicAgent {
    public:
        virtual void propagate(Eigen::VectorXd& state, Eigen::VectorXd& control, double dt) override;
        void HadiRK(Eigen::VectorXd& state, const Eigen::VectorXd& control, double dt);
        void HadiRK2(Eigen::VectorXd &state, const Eigen::VectorXd &control, double dt);
};