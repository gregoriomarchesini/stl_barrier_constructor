#include "task.h"
#include "fusion.h"
#include "hypercube.h"
#include "utils.h"
#include "task_scheduler.h"

using namespace mosek::fusion;
using namespace monty;


int main() {
    constexpr int dim        = 2;
    constexpr int num_agents = 100;
    double       side_length = 4.0;
    double       arena_size  = 10.0;
    double       max_input   = 10;

    std::vector<double> initial_state(num_agents * dim, 10.0); // Initialize with zeros

    // Create center vector using Eigen
    std::vector<double> center = {0.0, 0.0};

    int agent_i = 1;
    int agent_j = 2;

   
    auto task  = Task<dim>(side_length, center, 0.0, 10.0 , TaskType::EVENTUALLY, agent_i, agent_j, num_agents);
    auto task1 = Task<dim>(side_length, center, 20.0, 40.0, TaskType::EVENTUALLY, agent_i, agent_j, num_agents);
    auto task2 = Task<dim>(side_length, center, 15.0, 17.0, TaskType::ALWAYS    , agent_i, agent_j, num_agents);


    std::vector<Task<dim>> tasks;

    tasks.push_back(task);
    tasks.push_back(task1);
    tasks.push_back(task2);

    TaskScheduler<dim> scheduler(tasks, num_agents, arena_size , max_input, initial_state);


    return 0;
}