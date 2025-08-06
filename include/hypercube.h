#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include <iostream>
#include <vector>

#include "utils.h"


template <int Dim>
class Hypercube {
    static_assert(Dim == 2 || Dim == 3, "Hypercube only supports Dim = 2 or 3");
public:
    


    /**
     * @brief Defined an hypercube by its center and its side_length. The hypercube is defined as A@(x-c)<= b' 
     * where A is the matrix of the hypercube, x is the state of the system and c is the center of the hypercube.
     * The class stores the vector b = b' + A@C directly. The matrix A is defined based on the dimension of the hypercube and it is
     * not stire in the class. Namely we have
     * - If Dim = 2 then A = [[1,0],[-1,0],[0,1],[0,-1]]
     * - If Dim = 3 then A = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
     * 
     * @param side_length side length of the hypercube.
     * @param center Center of the hypercube.
     */
    

    
    Hypercube(double side_length, const std::vector<double>& center) 
    : size(side_length), center(center)
    {
        // Initialize A as a 4x2 sparse matrix with given entries
        static_assert(Dim == 2 || Dim == 3, "Vectornd only supports Dim == 2 or 3");
        

        if constexpr (Dim == 2){
            // Initialize b as zero vector of size 4 (or you can set it differently)
            b = new_array_ptr<double>({size / 2 + center[0], size / 2 - center[0], size / 2 + center[1], size / 2 - center[1]});
        }
        else if constexpr (Dim == 3) {
            // Initialize b as zero vector of size 4 (or you can set it differently)
            b = new_array_ptr<double>({ size / 2 + center[0] , size / 2 - center[0] , size / 2 + center[1] ,size / 2 - center[1],  size / 2 + center[2], size / 2 - center[2]});
        }

        A = get_A_matrix();        
    };


    mosek::fusion::Matrix::t get_b_in_matrix_form(int repetitions) const {
        // Create a vector of size 2*Dim*repetitions with the elements of b
        std::vector<double> b_vec;
        for (int i = 0; i < repetitions; ++i) {
            b_vec.insert(b_vec.end(), (*b).begin(), (*b).end());
        }
        
        return mosek::fusion::Matrix::dense(2*Dim, repetitions, new_array_ptr<double>(b_vec));
        
    }



    mosek::fusion::Matrix::t get_A_matrix(){


    // Returns the sparse matrix representing the multiplacation A@C where
    // C is the selction matrix such that x_i-x_j = C x where x is the state of the system
    // and A is the matrix of an hypercube described as A@(x-c)<= b. if the systems is 2D then
    // A = [[1 ,0],
    //      [-1,0],
    //      [0, 0],
    //      [0,-1]]
    // if if the system is 3D then
    // A = [[1,0,0],
    //      [-1,0,0],
    //      [0,1,0],
    //      [0,-1,0],
    //      [0,0,1],
    //      [0,0,-1]]


    std::vector<double> elements;

    if constexpr (Dim == 2) {
        elements = {1, 0, -1, 0, 0, 1, 0, -1};
    } else if constexpr (Dim == 3) {
        elements = {1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
    }
    
    return mosek::fusion::Matrix::dense(2*Dim,Dim, new_array_ptr<double>(elements));
    };

    mosek::fusion::Matrix::t get_output_matrix(const int num_agents, const int agent_i, const int agent_j) const {
        // Returns the output matrix A@C where C is the selection matrix such that x_i-x_j = C x
        // and A is the matrix of an hypercube described as A@(x-c)<= b.
        return output_matrix_C_mosek(Dim, num_agents, agent_i, agent_j);
    };


    
    // Attributes

    std::shared_ptr<ndarray<double>>      b;      // vector b of the hypercube defined as b = b' + A@C directly
    double                               size;   // side length of the hypercube
    std::vector<double>                  center; // center of the hypercube
    Matrix::t                            A;      // matrix A of the hypercube defined as A@(x-c)<= b
};

#endif 
