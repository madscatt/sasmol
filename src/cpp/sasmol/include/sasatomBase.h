#pragma once

#include <vector>
#include "Eigen/Dense"

namespace sasmol
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class SasAtomBase
    {

        protected:
            virtual int & _natoms() = 0;
            virtual int & _number_of_frames() = 0;
            virtual float & _total_mass() = 0;

            virtual std::vector<std::string> & _atom_record() = 0;
            virtual std::vector<int> & _atom_index() = 0;
            virtual std::vector<std::string> & _atom_name() = 0;
            virtual std::vector<std::string> & _atom_altloc() = 0;
            virtual std::vector<std::string> & _atom_resname() = 0;
            virtual std::vector<std::string> & _atom_chain() = 0;
            virtual std::vector<int> & _atom_resid() = 0;
            virtual std::vector<std::string> & _atom_icode() = 0;

            virtual Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> & _x() = 0;
            virtual Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> & _y() = 0;
            virtual Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> & _z() = 0;

            virtual std::vector<std::string> & _atom_occupancy() = 0;
            virtual std::vector<std::string> & _atom_beta() = 0;
            virtual std::vector<std::string> & _atom_segname() = 0;
            virtual std::vector<std::string> & _atom_element() = 0;
            virtual std::vector<std::string> & _atom_selement() = 0;
            virtual std::vector<std::string> & _atom_charge() = 0;

            virtual std::vector<std::string> & _unique_record() = 0;
            virtual std::vector<int> & _unique_index() = 0;
            virtual std::vector<std::string> & _unique_name() = 0;
            virtual std::vector<std::string> & _unique_altloc() = 0;
            virtual std::vector<std::string> & _unique_resname() = 0;
            virtual std::vector<std::string> & _unique_chain() = 0;
            virtual std::vector<int> & _unique_resid() = 0;
            virtual std::vector<std::string> & _unique_icode() = 0;
            virtual std::vector<std::string> & _unique_occupancy() = 0;
            virtual std::vector<std::string> & _unique_beta() = 0;
            virtual std::vector<std::string> & _unique_segname() = 0;
            virtual std::vector<std::string> & _unique_element() = 0;
            virtual std::vector<std::string> & _unique_selement() = 0;
            virtual std::vector<std::string> & _unique_charge() = 0;

            virtual std::vector<float> & _atom_com() = 0;
            
            virtual Eigen::Matrix<float,3,1> & _uk() = 0;
            virtual Eigen::Matrix3f & _ak() = 0;
        
            virtual Eigen::Array<float,Eigen::Dynamic,1> & _atom_mass() = 0;

    };
}

