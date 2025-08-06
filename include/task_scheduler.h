#ifndef TASK_SCHEDULER_H
#define TASK_SCHEDULER_H

#include <cstdlib>
#include <vector>

#include "hypercube.h"
#include "utils.h"
#include "task.h"
#include "fusion.h"
#include "barrier_function.h"


using namespace mosek::fusion;
using namespace monty;


template <int Dim>
class TaskScheduler {

    public:
    // methods
    TaskScheduler(std::vector<Task<Dim>>& input_tasks, int num_agents, double arena_size, double max_input, std::vector<double>& initial_state)
        : tasks(input_tasks), num_agents(num_agents), arena_size(arena_size), max_input(max_input), initial_state(initial_state), global_dim(num_agents * Dim) {

            this->time_optimizer    = new Model(); // this is a custom spart pointer type (like shared pointer)
            this->barrier_optimizer = new Model(); // this is a custom spart pointer type (like shared pointer)
            make_time_schedule();
            create_barrier_functions();
            optimize_barrier_function();
        };

    ~TaskScheduler() {
        this->time_optimizer->dispose();    // cleanup explicitly
        this->barrier_optimizer->dispose(); // cleanup explicitly
    }

    // make time schedule
    void make_time_schedule(){
        // Implementation for creating a time schedule based on the tasks
        // This could involve determining the start and end times for each task
        // and ensuring that they do not overlap.

        std::vector<int> eventually_tasks_index;
        std::vector<int> always_tasks_index;

        // divide the tasks into eventually tasks and always tasks 
        for (size_t i = 0; i < tasks.size(); ++i) {
            if (tasks[i].getType() == TaskType::EVENTUALLY) {
                eventually_tasks_index.push_back(i);
            } else if (tasks[i].getType() == TaskType::ALWAYS) {
                always_tasks_index.push_back(i);
            }
        }

        if (eventually_tasks_index.empty()){
            std::cout << "No tasks of type EVENTUALLY found. Noting to optimize." << std::endl;
            return;
        }

        // create a variable vector for the two tasks types
        auto t_eventually = time_optimizer->variable(eventually_tasks_index.size(), Domain::unbounded()); // eventua;;y task have only one time instant to take case

        // add constraints on the eventually time variable (must be contained in the interval of the eventually operator)
        for (size_t i : eventually_tasks_index) {
            Task<Dim>* task = &tasks[i];
            time_optimizer->constraint(t_eventually->index(i), Domain::greaterThan(task->start_time));
            time_optimizer->constraint(t_eventually->index(i), Domain::lessThan(task->end_time));
        }
        
        // define now the cost for the optimizaiton

        Expression::t cost = Expr::zeros(1);

        for (size_t i : eventually_tasks_index) {
            
            Task<Dim>* task = &tasks[i];
            // find the task eventually with the closest end_time to the start time of the eventually
            
            int closest_lower_time = task->start_time; // initialize to a large negative number
            int closest_upper_time = task->end_time; // initialize to a large positive number
            
            for (size_t j : always_tasks_index) {
                Task<Dim>* always_task = &tasks[j];
                
                // Some checks on the tasks for now 
                if (!always_task->is_same_edge(*task)) {
                    continue; // skip if the two tasks are not on the same edge
                }

                else if (task->is_contained_in(*always_task)) {
                   throw std::runtime_error("The task of type EVENTUALLY is contained in a task of type ALWAYS. This is not allowed for now. Only intersections are allowed.");
                }
                else if (always_task->is_contained_in(*task)) {
                    throw std::runtime_error("The task of type ALWAYS is contained in a task of type EVENTUALLY. This is not allowed for now. Only intersections are allowed.");
                }


                if ( (always_task->end_time - task->start_time) > closest_lower_time && 
                             always_task->end_time <= task->end_time) {
                    closest_lower_time = always_task->end_time ; // cae of always task interescting from the left
                }

                if (always_task->start_time - task->end_time < closest_upper_time && 
                             always_task->start_time >= task->start_time) {
                    closest_upper_time = always_task->start_time; // case of alwats task intersecting on the right
                }

            }

            // print the bounds
            std::cout << "Task " << i << ": Closest lower time: " << closest_lower_time 
                      << ", Closest upper time: " << closest_upper_time << std::endl;
            
            if (closest_lower_time >= closest_upper_time && always_tasks_index.size() > 0) {
                throw std::runtime_error("No valid time schedule found for the task. At least one task of type eventually can not be satisfied over one edge due to overallping with two tasks of type always");
            }
            else {
                auto yl = this->time_optimizer->variable(1);
                auto yu = this->time_optimizer->variable(1);

                // create cone constraint
                auto delta_l = Expr::sub(closest_lower_time, t_eventually->index(i)); // -(t - t_l)
                auto delta_u = Expr::sub(t_eventually->index(i), closest_upper_time); // -(t_u - t)
                auto xl      = Expr::vstack(yl,1, delta_l); // represents y = e^(-(t + t_l))
                auto xu      = Expr::vstack(yu,1, delta_u); // represents y = e^(-(t_u - t))

                // add cost amd exponential cone constraint
            
                cost = Expr::add(cost, yl);
                cost = Expr::add(cost, yu);

                // exponential cone constraint
                time_optimizer->constraint( xl, Domain::inPExpCone());
                time_optimizer->constraint( xu, Domain::inPExpCone());

            }
        }


        this->time_optimizer->objective("obj", ObjectiveSense::Minimize, cost);
        this->time_optimizer->solve();

        // print solution
        ndarray<double, 1> sol = *t_eventually->level();
        std::cout << "Optimal time schedule for EVENTUALLY tasks: " << std::endl;
        std::cout << sol << std::endl;
        

        for (size_t i : eventually_tasks_index) {

            Task<Dim>& task = tasks[i];
            task.set_tau(sol[i]); // assuming agent_i is the index of the task in the eventually_tasks_index vector
            // print the time
            std::cout << "Task " << i << ": Optimal time: " << task.tau << std::endl;
        }

        std::cout << "Tasks size: " << tasks.size() << std::endl;

        // create time grid
        for (size_t i=0; i < tasks.size(); ++i){
            if (tasks[i].getType() == TaskType::EVENTUALLY) {

                //print tau value
                std::cout << "Task " << i << ": tau = " << tasks[i].tau << std::endl;

                time_grid.push_back(tasks[i].tau);
            } else if (tasks[i].getType() == TaskType::ALWAYS) {
                // print start and end time
                std::cout << "Task " << i << ": start_time = " << tasks[i].start_time 
                          << ", end_time = " << tasks[i].end_time << std::endl;
                time_grid.push_back(tasks[i].start_time);
                time_grid.push_back(tasks[i].end_time);
            }


        }

        time_grid.push_back(0.0); // add 0.0 to the time grid
    
        // sort the time grid
        std::sort(time_grid.begin(), time_grid.end());
        auto last = std::unique(time_grid.begin(), time_grid.end());
        time_grid.erase(last, time_grid.end());

        // Print result
        std::cout << "Time grid: " << std::endl;
        for (double v : time_grid) {
            std::cout << v << " " << std::endl;
        }
        
    };

    void create_barrier_functions() {
        // Implementation for creating barrier functions based on the tasks
        // This could involve setting up the barrier functions for each task
        // and ensuring that they are optimized correctly.

        barriers.reserve(tasks.size());

        for (auto& task : tasks) {
            std::vector<double> time_grid_barrier;

            for (double t : time_grid) {
                if (task.getType() == TaskType::EVENTUALLY) {
                    // For EVENTUALLY tasks, we only consider the time instant tau
                    if (t <= task.tau) {
                        time_grid_barrier.push_back(t);
                    }
                } else if (task.getType() == TaskType::ALWAYS) {
                    // For ALWAYS tasks, we consider the start and end times
                    if (t <= task.end_time) {
                        time_grid_barrier.push_back(t);
                    }
                }
            };

            // create the barrier function
            BarrierFunction<Dim> barrier = BarrierFunction<Dim>( task, time_grid_barrier, barrier_optimizer );
            barriers.push_back(barrier);
        }
    }


    std::vector<int> active_tasks(double t) const {
        // Implementation for determining which tasks are active at a given time t
        // This could involve checking the start and end times of each task
        // and returning a list of active tasks.

        std::vector<int> active_tasks;
        for (size_t i = 0; i < barriers.size(); ++i) {
            if (t<  barriers[i].time_grid.back()) {
                active_tasks.push_back(i);
            }
        }
        return active_tasks;
    };

    
    
    
    void optimize_barrier_function() {
        // Implementation for optimizing the barrier function
        // This could involve setting up a barrier function and optimizing it
        // based on the tasks and their time constraints.
        
        // For now, this is just a placeholder
        std::cout << "Optimizing barrier function..." << std::endl;
        
        int num_vertices      = std::pow(2, global_dim); // number of vertices in the hypercube
        int num_intervals     = time_grid.size() - 1; // number of intervals in the time grid
        auto center           = std::vector<double>(global_dim, 0.0); // center of the hypercube

        // auto slack  = this->barrier_optimizer->variable(Domain::greaterThan(0.0)); // slack variables for the barrier function
        // auto UV     = this->barrier_optimizer->variable(new_array_ptr<int,1>({global_dim, 2*num_vertices*num_intervals})); // variables for the vertices of the hypercube
        // auto k_gain = this->barrier_optimizer->parameter(); // Scalar parameter
        auto x_0    = this->barrier_optimizer->parameter(new_array_ptr<int,1>({global_dim,1})); // initial state of the system
        x_0->setValue(new_array_ptr<double>(initial_state));
    
        // create mosek matrix of the vertices
        // int  repetitions = 2;
        // auto VV         = generate_hypercube_vertices_mosek(global_dim, arena_size, repetitions); // vertices of the hypercube, repeated twice and stack in a matrix

        // Define the problem onstraints and objective 

        // 1: input constraints
        
        // auto BB = generate_hypercube_B_matrix_mosek(global_dim, max_input*2,2*num_vertices*num_intervals);
        // auto A  = generate_hypercube_A_matrix_mosek(global_dim, max_input*2);

        // auto c = this->barrier_optimizer->constraint( Expr::sub(Expr::mul(A,UV), BB) , Domain::lessThan(0.) );

        // // 2 : add forward invariance constraints
        // for (int i = 0; i < num_intervals; ++i) {

        //     double s_j  = time_grid[i];
        //     double s_j1 = time_grid[i+1]-0.00001;
        //     auto   U_j  = UV->slice(new_array_ptr<int,1>({0,i*2*num_vertices}), new_array_ptr<int,1>({global_dim,(i+1)*2*num_vertices}));
        //     auto   dyn  = Expr::add(VV, U_j);



        //     for(int j : this->active_tasks(s_j)) {
        //         BarrierFunction<Dim>& barrier = this->barriers[j];
                
                
        //         auto gamma_sj  = barrier.gamma_vector(s_j);
        //         auto gamma_sj1 = barrier.gamma_vector(s_j1);

        //         // repeat the vector for the number of vertices
        //         auto gamma_sj_rep   = Expr::repeat(gamma_sj,  num_vertices ,1);
        //         auto gamma_sj1_rep  = Expr::repeat(gamma_sj1, num_vertices ,1);
        //         auto gamma          = Expr::hstack(new_array_ptr<Expression::t>({gamma_sj_rep, gamma_sj1_rep}));
                
        //         auto e_sj      = barrier.e_vector(s_j);
        //         auto e_sj1     = barrier.e_vector(s_j1);
        //         auto e_sj_rep  = Expr::repeat(e_sj,  num_vertices ,1);
        //         auto e_sj1_rep = Expr::repeat(e_sj1, num_vertices ,1);

        //         auto e     = Expr::hstack(new_array_ptr<Expression::t>({e_sj_rep, e_sj1_rep}));
        //         auto D     = barrier.get_D_matrix(num_agents);
        //         auto C     = barrier.get_c_vector_stacked(2*num_vertices);
        //         auto ones  = Matrix::ones(2*Dim, 2*num_vertices );


        //         auto exp1  = Expr::mul(D,dyn);
        //         auto exp21 = Expr::mul(Expr::constTerm(D), VV);


        //         auto exp2 = Expr::mul(exp21,k_gain);
        //         auto exp3 = Expr::mul(C,k_gain);
        //         auto exp4 = Expr::mul(gamma,k_gain);
        //         auto exp5 = Expr::mul(ones,slack); 
        //         auto exp6 = e;

        //         auto exps = new_array_ptr<Expression::t,1>({exp1, exp2, exp5});
        //         // print diemnsions of each term

        //        // print shape of experessions
        //         std::cout << "Shape of exp1: " << (*exp1->getShape())[0] << ", " << (*exp1->getShape())[1] << std::endl;
        //         std::cout << "Shape of exp2: " << (*exp2->getShape())[0] << ", " << (*exp2->getShape())[1] << std::endl;
        //         std::cout << "Shape of exp3: " << (*exp3->getShape())[0] << ", " << (*exp3->getShape())[1] << std::endl;
        //         std::cout << "Shape of exp4: " << (*exp4->getShape())[0] << ", " << (*exp4->getShape())[1] << std::endl;
        //         std::cout << "Shape of exp5: " << (*exp5->getShape())[0] << ", " << (*exp5->getShape())[1] << std::endl;
        //         std::cout << "Shape of exp6: " << (*exp6->getShape())[0] << ", " << (*exp6->getShape())[1] << std::endl;


        //         auto constraint = Expr::add(exps);

        //         // some stupid print
        //         std::cout << "Print something for debugging purposes" << std::endl;
                
        //         this->barrier_optimizer->constraint(constraint  , Domain::greaterThan(0.));


        //     }

        // add flat constraint
        for (BarrierFunction<Dim> barrier : barriers) {
            barrier.set_gamma_flat_constraint();
            barrier.set_max_input_constraint(max_input);
        }

        // define the beta vector 
        std::vector<double> betas;
        for (BarrierFunction<Dim> barrier : barriers) {
            betas.push_back(barrier.time_grid.back());
        }
        betas.push_back(0.0); // add 0.0 to the betas vector


        // order the betas vector in ascending order
        std::sort(betas.begin(), betas.end());

        // print the betas vector
        std::cout << "Betas vector: " << std::endl;
        for (double beta : betas) {
            std::cout << beta << " ";
        }
        std::cout << std::endl;

        auto zeta_vars = this->barrier_optimizer->variable(new_array_ptr<int,1>({global_dim,betas.size()})); // zeta variables for the beta vector

        for (int l =1; l < betas.size(); ++l) {
            double beta_l   = betas[l];
            auto   zeta_l   = zeta_vars->slice(new_array_ptr<int>({0,l}), new_array_ptr<int>({global_dim,l+1}));
            
            
            for (int l_tilde : active_tasks(betas[l-1])) {

                auto gamma_at_beta_l = barriers[l_tilde].gamma_vector(beta_l);      // gamma vector at beta_l - EPSILON
                auto D               = barriers[l_tilde].get_D_matrix();  // D matrix for the barrier function
                auto c               = barriers[l_tilde].get_c_vector();            // c vector for the barrier

                auto exps    = Expr::add(new_array_ptr<Expression::t,1>({Expr::mul(D,zeta_l), Expr::constTerm(c), gamma_at_beta_l}));
                this->barrier_optimizer->constraint(exps , Domain::greaterThan(0.));
            }
        }

        // initial condition constraint
        auto zeta_0  = zeta_vars->slice(new_array_ptr<int,1>({0,0}), new_array_ptr<int,1>({global_dim,1})); // pick first column

        for (BarrierFunction<Dim> barrier : barriers) {
            
            auto gamma_0 = barrier.gamma_vector(0.0); // gamma
            auto D_0     = barrier.get_D_matrix(); // D matrix for the barrier function
            auto c_0     = barrier.get_c_vector(); // c vector for the

            auto exps    = Expr::add(new_array_ptr<Expression::t,1>({Expr::mul(D_0,zeta_0), Expr::constTerm(c_0), gamma_0}));
            
            this->barrier_optimizer->constraint(exps, Domain::greaterThan(0.));
        }

        this->barrier_optimizer->constraint(Expr::sub(x_0, zeta_0), Domain::equalsTo(0.));



        // create objective of the optimization problem 
        
        std::vector<Expression::t>  sum_vec;
        double slack_penalty = 1E5; // penalty for the slack variable

        for (auto barrier : barriers){
            sum_vec.push_back(Expr::sum(barrier.r_var));
            sum_vec.push_back(Expr::neg(Expr::mul(barrier.slack_var, slack_penalty))); // add the slack variable with a penalty

        }

        this->barrier_optimizer->objective("obj", ObjectiveSense::Maximize, Expr::add(new_array_ptr<Expression::t>(sum_vec)));
        this->barrier_optimizer->solve();

        // print barriers slopes
        for (size_t i = 0; i < barriers.size(); ++i) {
            std::cout << "Barrier " << i << ": Slopes: " << *barriers[i].slopes_var->level() << std::endl;
            std::cout << "Barrier " << i << ": Gamma0: " << *barriers[i].gamma0_var->level() << std::endl;
            std::cout << "Barrier " << i << ": r: "      << *barriers[i].r_var->level()      << std::endl;
            std::cout << "Barrier " << i << ": Slack: "  << *barriers[i].slack_var->level()  << std::endl;
        }
        
        
    };





    //attributes
    std::vector<Task<Dim>>            tasks;              // List of tasks to be
    std::vector<BarrierFunction<Dim>> barriers;           // list of barriers to be optimized
    std::vector<double>               time_grid;          // Time grid for the tasks
    std::vector<Expr>                 time_constraints;   // time constraints for selecting the time of satisfaction of a task
    std::vector<double>               initial_state;      // Initial state of the system, used for the optimization
   
    
    Model::t                          time_optimizer    = nullptr; // Custom smart pointer type for the time_optimizer
    Model::t                          barrier_optimizer = nullptr; // Custom smart pointer type for the barrier_optimizer
    
    
    const int                         num_agents; // Number of agents in the system
    const double                      arena_size; // Size of the arena for the tasks
    const double                      max_input;  // Maximum input value for the tasks
    const int                         global_dim;  // num_agents * Dim, global dimension of the system

};



#endif