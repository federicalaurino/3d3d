//Useful type definitions (built using the predefined types in Gmm++)
//Special class for small (dim < 16) vectors 
using bgeot::base_small_vector;  
// Geometrical nodes (derived from base_small_vector)
using bgeot::base_node;   	 	 
// Double-precision FP numbers 
using bgeot::scalar_type; 	 	
// Unsigned long integers 
using bgeot::size_type;   	 	 
// Short integers 
using bgeot::short_type;         
// Type for vector of integers 
typedef std::vector<size_type> vector_size_type;
// Type for dense vector
typedef std::vector<scalar_type> vector_type;
// Type for sparse vector (std::vector) 
typedef gmm::rsvector<scalar_type> sparse_vector_type;
// Type for dense matrix
typedef gmm::row_matrix<vector_type> matrix_type;
// Type for sparse matrix
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;    

// Pi
const scalar_type pi = std::atan(1)*4;
