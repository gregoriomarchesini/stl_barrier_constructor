#include <stdexcept>
#include <vector>

#include "utils.h"
#include "fusion.h" 

using namespace mosek::fusion;
using namespace monty;

std::vector<std::vector<double>> generate_hypercube_vertices(const int dim, const double size) {
    int num_vertices = 1 << dim;  // 2^dim
    double half_size = size / 2.0;

    std::vector<std::vector<double>> vertices(dim, std::vector<double>(num_vertices)); // allocate vector

    for (int v = 0; v < num_vertices; ++v) {
        for (int d = 0; d < dim; ++d) {
            int bit = (v >> d) & 1;
            vertices[d][v] = (bit == 0) ? -half_size : +half_size;
        }
    }

    return vertices;
}


Matrix::t generate_hypercube_vertices_mosek(const int dim, const double size){
    int num_vertices = 1 << dim;  // 2^dim
    double half_size = size / 2.0;

    std::vector<double> vertices(num_vertices * dim);

    // Fill in column-major order: [x0 x1 ... xn y0 y1 ... yn ...]
    for (int d = 0; d < dim; ++d) {
        for (int v = 0; v < num_vertices; ++v) {
            int bit = (v >> d) & 1;
            vertices[d * num_vertices + v] = (bit == 0) ? -half_size : +half_size;
        }
    }

    return Matrix::dense(dim,std::pow(2,dim) ,new_array_ptr<double>(vertices));
}


Matrix::t generate_hypercube_vertices_mosek(const int dim, const double size, const int repetitions){
    int num_vertices = 1 << dim;  // 2^dim
    double half_size = size / 2.0;

    std::vector<double> vertices(repetitions*num_vertices * dim);

    // Fill in column-major order: [x0 x1 ... xn y0 y1 ... yn ...]
    for (int r = 0; r < repetitions; ++r) {
        for (int d = 0; d < dim; ++d) {
            for (int v = 0; v < num_vertices; ++v) {
                int bit = (v >> d) & 1;
                vertices[r * dim * num_vertices + d * num_vertices + v] = (bit == 0) ? -half_size : +half_size;
            }
        }
    }
    
    return  Matrix::dense(dim,std::pow(2,dim)*repetitions ,new_array_ptr<double>(vertices));
}



Matrix::t output_matrix_C_mosek(const int agents_dim,const int num_agents,const int agent_i, const int agent_j){

    if (agents_dim != 2 && agents_dim != 3) {
        throw std::invalid_argument("agents_dim must be either 2 or 3.");
    }


    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> values;


    for (int i = 0; i < agents_dim; ++i) {
        int row = i;
        int col_agent_i = agents_dim * agent_i + i;
        int col_agent_j = agents_dim * agent_j + i;

        rows.push_back(row);
        cols.push_back(col_agent_i);
        values.push_back(1.0);  // Coefficient for agent i in agent j
        
        if (agent_i != agent_j) { 
            // add the row for agent j
            rows.push_back(row);
            cols.push_back(col_agent_j);
            values.push_back(-1.0); 

        }
        
    }

    // return the output matrix as a sparse matrix    
    return Matrix::sparse(agents_dim, agents_dim*num_agents,
                                           new_array_ptr<int>(rows), 
                                           new_array_ptr<int>(cols), 
                                           new_array_ptr<double>(values));


};


Matrix::t matrix_D_mosek(const int agents_dim,const int num_agents,const int agent_i, const int agent_j){
    // A@C

    if (agents_dim != 2 && agents_dim != 3) {
        throw std::invalid_argument("agents_dim must be either 2 or 3.");
    }


    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> values;


    for (int i = 0; i < agents_dim; ++i) {
        int row = 2*i;
        int col_agent_i = agents_dim * agent_i + i;
        int col_agent_j = agents_dim * agent_j + i;

        rows.push_back(row);
        cols.push_back(col_agent_i);
        values.push_back(-1.0);  // Coefficient for agent i in agent j

        rows.push_back(row+1);
        cols.push_back(col_agent_i);
        values.push_back(1.0);  // Coefficient for agent i in agent j
        
        if (agent_i != agent_j) { 
            // add the row for agent j
            rows.push_back(row);
            cols.push_back(col_agent_j);
            values.push_back(1.0); 

            // add the row for agent j
            rows.push_back(row);
            cols.push_back(col_agent_j);
            values.push_back(-1.0); 

        }
        
    }

    // return the output matrix as a sparse matrix    
    return Matrix::sparse(2*agents_dim, agents_dim*num_agents,
                                           new_array_ptr<int>(rows), 
                                           new_array_ptr<int>(cols), 
                                           new_array_ptr<double>(values));


};



