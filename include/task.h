#ifndef TASK_H
#define TASK_H

#include <vector>

#include "hypercube.h"
#include "utils.h"
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;



/**
 * @file task.h
 * @brief Defines the Task class template representing a temporal constraint between two agents.
 *
 * A Task encodes a geometric region (a hypercube) and a time interval during which
 * a relation (either ALWAYS or EVENTUALLY) between two agents must be enforced.
 * The task is used to construct constraints in optimization problems.
 */

/**
 * @brief Represents a temporal task involving two agents in a multi-agent system.
 *
 * Each task is defined by:
 * - A spatial constraint represented by a hypercube
 * - A time interval [start_time, end_time]
 * - A task type: ALWAYS (must hold for the entire interval) or EVENTUALLY (must hold at least once)
 *
 * The class provides methods to generate constraint matrices and utility checks for scheduling.
 *
 * @tparam Dim Dimension of the space in which the hypercube lies (e.g., 2D or 3D).
 */

template <int Dim>
class Task {
public:
    
    // Methods
    Task(  double                             side_length, 
           const std::vector<double>&         center, 
           const double                       start, 
           const double                       end, 
           TaskType                           type,
           const int                          agent_i,
           const int                          agent_j,
           const int                          num_agents)
           
           : hypercube(side_length, center), 
             start_time(start), 
             end_time(end), 
             type(type), 
             agent_i(agent_i), 
             agent_j(agent_j),
             num_agents(num_agents)  {


            D = matrix_D_mosek(Dim,num_agents,agent_i, agent_j);
             };

    
    /**
     * @brief Returns the type of the task.
     * @return Task type (ALWAYS or EVENTUALLY).
     */
    const TaskType& getType() const {
        return type;
    };

    

    bool is_same_edge(const Task<Dim>& other) const {
        return (agent_i == other.agent_i && agent_j == other.agent_j) || 
               (agent_i == other.agent_j && agent_j == other.agent_i);
    };

    bool is_contained_in(const Task<Dim>& other) const {
        return (start_time >= other.start_time && end_time <= other.end_time);
    };

    void set_tau(double tau) {
        this->tau = tau;
    };



    // ====== Attributes ======

    Hypercube<Dim> hypercube;   ///< Spatial constraint region
    const double   start_time;  ///< Start time of the task
    const double   end_time;    ///< End time of the task
    const TaskType type;        ///< Task type (ALWAYS or EVENTUALLY)
    const int      agent_i;     ///< Index of the first agent
    const int      agent_j;     ///< Index of the second agent
    const int      num_agents;  ///< Total number of agents in the system
    Matrix::t      D;           ///< Mosek matrix D for the inequality Dx + c >=0 hypercurbe representation 
                                ///<(stored in the task as it requires the knowledge of the number of agents and the edges over which the task is defined)
    double         tau = -1.0;  ///< Time at which the EVENTUALLY task is satisfied (used in optimization)

                        
};



#endif
