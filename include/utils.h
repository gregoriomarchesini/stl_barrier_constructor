#ifndef HYPERCUBE_UTILS_H
#define HYPERCUBE_UTILS_H

#include <vector>

#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

/**
 * @file hypercube_utils.h
 * @brief Utility functions for generating hypercube structures and constraint matrices for optimization.
 *
 * This header provides tools to:
 * - Generate vertices of an n-dimensional hypercube
 * - Format those vertices for use with the MOSEK Fusion API
 * - Construct A matrices for task-based constraints
 */

/**
 * @brief Generate all vertices of a hypercube centered at the origin.
 *
 * Each column of the returned matrix represents a vertex of the hypercube.
 *
 * @param dim  The dimension of the hypercube.
 * @param size The length of the hypercube's side.
 * @return A 2D vector where each inner vector is a coordinate dimension,
 *         and each column is a vertex of the hypercube.
 */
std::vector<std::vector<double>> generate_hypercube_vertices(const int dim, const double size);

/**
 * @brief Generate hypercube vertices and return them as a MOSEK-compatible array.
 *
 * The result is ordered such that the matrix can be reshaped as (dim x 2^dim).
 *
 * @param dim  The dimension of the hypercube.
 * @param size The length of the hypercube's side.
 * @return A shared pointer to a 1D MOSEK ndarray containing the flattened hypercube vertices.
 */
Matrix::t generate_hypercube_vertices_mosek(const int dim, const double size);

/**
 * @brief Generate hypercube vertices and return them as a MOSEK-compatible array.
 *
 * The result is ordered such that the matrix can be reshaped as (dim x 2^dim).
 *
 * @param dim  The dimension of the hypercube.
 * @param size The length of the hypercube's side.
 * @param repetitions Number of times to repeat the vertices in the output.
 * @return A shared pointer to a 1D MOSEK ndarray containing the flattened hypercube vertices.
 */
Matrix::t generate_hypercube_vertices_mosek(const int dim, const double size, const int repetitions);

/**
 * @brief Given a constraint of the form A(x_i-x_j) <= b, returns the matrix A@C where x_1-x_j = Cx.
 * so as to apply the constraint A@Cx<= b directly in mosek. The matrix A@C is returned as sparse mosek matrix.
 *
 * @param agents_dim The dimension of each agent's state.
 * @param num_agents The total number of agents.
 * @param index_i    The index of the first agent (corresponding to +1).
 * @param index_j    The index of the second agent (corresponding to -1).
 * @return A sparse Fusion matrix representing the A block row.
 */



Matrix::t output_matrix_C_mosek(const int agents_dim,const int num_agents,const int agent_i, const int agent_j);

Matrix::t matrix_D_mosek(const int agents_dim,const int num_agents,const int agent_i, const int agent_j);



/**
 * @brief Enumeration of temporal task types.
 */
enum class TaskType {
    ALWAYS     = 0, ///< Task must hold at all time steps.
    EVENTUALLY = 1  ///< Task must hold at least once within a given interval.
};

#endif  // HYPERCUBE_UTILS_H