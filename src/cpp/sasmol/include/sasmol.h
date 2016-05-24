#pragma once

#include <vector>
#include "Eigen/Dense"
#include "sasio.h"
#include "sascalc.h"
#include "sasop.h"
#include "sasmath.h"
#include "sassubset.h"

namespace sasmol
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class SasAtom: public sasio::Files, public sascalc::Prop, public sasop::Move, public sasmath::Math, public sassubset::Mask
    {

        public:
            SasAtom() = default ; ///< (Brief description)
            
            int natoms = 0 ; ///< (Brief description)
            int number_of_frames = 0 ; ///< (Brief description)
            float total_mass = 0.0 ; ///< (Brief description)

            std::vector<std::string> atom_record ; ///< (Brief description)
            std::vector<int> atom_index ; ///< (Brief description)
            std::vector<std::string> atom_name ; ///< (Brief description)
            std::vector<std::string> atom_altloc ; ///< (Brief description)
            std::vector<std::string> atom_resname ; ///< (Brief description)
            std::vector<std::string> atom_chain ; ///< (Brief description)
            std::vector<int> atom_resid ; ///< (Brief description)
            std::vector<std::string> atom_icode ; ///< (Brief description)

            Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> x ; ///< (Brief description)
            Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> y ; ///< (Brief description)
            Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> z ; ///< (Brief description)

            std::vector<std::string> atom_occupancy ; ///< (Brief description)
            std::vector<std::string> atom_beta ; ///< (Brief description)
            std::vector<std::string> atom_segname ; ///< (Brief description)
            std::vector<std::string> atom_element ; ///< (Brief description)
            std::vector<std::string> atom_selement ; ///< (Brief description)
            std::vector<std::string> atom_charge ; ///< (Brief description)
            
            std::vector<std::string> unique_record ; ///< (Brief description)
            std::vector<int> unique_index ; ///< (Brief description)
            std::vector<std::string> unique_name ; ///< (Brief description)
            std::vector<std::string> unique_altloc ; ///< (Brief description)
            std::vector<std::string> unique_resname ; ///< (Brief description)
            std::vector<std::string> unique_chain ; ///< (Brief description)
            std::vector<int> unique_resid ; ///< (Brief description)
            std::vector<std::string> unique_icode ; ///< (Brief description)
            std::vector<std::string> unique_occupancy ; ///< (Brief description)
            std::vector<std::string> unique_beta ; ///< (Brief description)
            std::vector<std::string> unique_segname ; ///< (Brief description)
            std::vector<std::string> unique_element ; ///< (Brief description)
            std::vector<std::string> unique_selement ; ///< (Brief description)
            std::vector<std::string> unique_charge ; ///< (Brief description)

            std::vector<float> atom_com ; ///< (Brief description)
            
            Eigen::Matrix<float,3,1> uk ;     // eigenvalues
            Eigen::Matrix3f ak ;        // eigenvectors
        
            Eigen::Array<float,Eigen::Dynamic,1> atom_mass ; ///< (Brief description)
            
            ~SasAtom() = default ; ///< (Brief description) 

        public:
            virtual int & _natoms() { return natoms; } ///< (Brief description)
            virtual int & _number_of_frames() { return number_of_frames; } ///< (Brief description)
            virtual float & _total_mass() { return total_mass; } ///< (Brief description)

            virtual std::vector<std::string> & _atom_record() { return atom_record; } ///< (Brief description)
            virtual std::vector<int> & _atom_index() { return atom_index; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_name() { return atom_name; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_altloc() { return atom_altloc; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_resname() { return atom_resname; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_chain() { return atom_chain; } ///< (Brief description)
            virtual std::vector<int> & _atom_resid() { return atom_resid; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_icode() { return atom_icode; } ///< (Brief description)

            virtual Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> & _x() { return x; } ///< (Brief description)
            virtual Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> & _y() { return y; } ///< (Brief description)
            virtual Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> & _z() { return z; } ///< (Brief description)

            virtual std::vector<std::string> & _atom_occupancy() { return atom_occupancy; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_beta() { return atom_beta; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_segname() { return atom_segname; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_element() { return atom_element; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_selement() { return atom_selement; } ///< (Brief description)
            virtual std::vector<std::string> & _atom_charge() { return atom_charge; } ///< (Brief description)

            virtual std::vector<std::string> & _unique_record() { return unique_record; } ///< (Brief description)
            virtual std::vector<int> & _unique_index() { return unique_index; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_name() { return unique_name; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_altloc() { return unique_altloc; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_resname() { return unique_resname; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_chain() { return unique_chain; } ///< (Brief description)
            virtual std::vector<int> & _unique_resid() { return unique_resid; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_icode() { return unique_icode; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_occupancy() { return unique_occupancy; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_beta() { return unique_beta; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_segname() { return unique_segname; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_element() { return unique_element; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_selement() { return unique_selement; } ///< (Brief description)
            virtual std::vector<std::string> & _unique_charge() { return unique_charge; } ///< (Brief description)

            virtual std::vector<float> & _atom_com() { return atom_com; } ///< (Brief description)
            
            virtual Eigen::Matrix<float,3,1> & _uk() { return uk; } ///< (Brief description)
            virtual Eigen::Matrix3f & _ak() { return ak; } ///< (Brief description)
        
            virtual Eigen::Array<float,Eigen::Dynamic,1> & _atom_mass() { return atom_mass; } ///< (Brief description)
    };

    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class SasMol : public SasAtom{

        public:
            SasMol() = default ;
            SasMol(const int Natoms, const int Nframes,
                   const long long *const  p_index, const std::string * const p_name, const std::string * const p_altloc, const std::string * const p_resname, const std::string * const p_chain, const long long *const p_resid, const std::string * const p_icode,
                   const double *const p_x, const double *const p_y, const double *const p_z,
                   const std::string * const p_occupancy, const std::string * const p_beta, const std::string * const p_segname, const std::string * const p_element, const std::string * const p_charge);
            int none ;
        
    } ;
}


