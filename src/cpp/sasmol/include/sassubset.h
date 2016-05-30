#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>

#include "sasatomBase.h"

/********* declarations   ******************/



namespace sassubset 
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Mask: public sasmol::SasAtomBase
    {
        protected:
            std::vector<boost::dynamic_bitset<>> mask_names; ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> mask_resnames; ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> mask_resids; ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> mask_chains; ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> mask_occupancies; ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> mask_betas; ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> mask_elements; ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> mask_segnames; ///< (Brief description)

        public:
            std::vector<boost::dynamic_bitset<>> get_dihedral_subset_mask(const std::vector<int> & flexible_residues, const int mtype); ///< (Brief description)
            std::vector<int> get_indices_from_mask(const boost::dynamic_bitset<> & index); ///< (Brief description)
            boost::dynamic_bitset<> get_mask_from_indices(const int natoms, const std::vector<int> & index); ///< (Brief description)
	        boost::dynamic_bitset<> get_subset_mask(const std::string & basis_filter); ///< (Brief description)
            void merge_two_molecules(sasmol::SasMol & mol1, sasmol::SasMol & mol2); ///< (Brief description)
            void copy_molecule_using_mask(sasmol::SasMol & other, const boost::dynamic_bitset<> & mask, const size_t frame); ///< (Brief description)
            Eigen::Matrix<float,3,Eigen::Dynamic> get_coor_using_mask(const size_t frame, const boost::dynamic_bitset<> & mask); ///< (Brief description)
            void set_coor_using_mask(sasmol::SasMol & mol, const size_t frame, const boost::dynamic_bitset<> & mask); ///< (Brief description)
            sassubset::Mask & initialize_children(); ///< (Brief description)
	        void duplicate_molecule(sasmol::SasMol & mol, const int number_of_duplicates, const int frame, std::vector<std::vector<float>> & com_coor, const bool flag_rotate=true);
            template <typename T> void set_descriptor_using_mask(const boost::dynamic_bitset<> & mask, std::vector<T> & descriptor, const T value); ///< (Brief description)

        private:
            void get_mask_array(std::vector<boost::dynamic_bitset<>> & farray, int nflexible, int natoms, const std::vector<std::string> & name, const std::vector<int> & resid, const std::vector<int> & flexible_residues, const int nresidues, const int mtype); ///< (Brief description)
        public:
            std::vector<boost::dynamic_bitset<>> & _mask_names() { return mask_names;} ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> & _mask_resnames() { return mask_resnames;} ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> & _mask_resids() { return mask_resids;} ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> & _mask_chains() { return mask_chains;} ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> & _mask_occupancies() { return mask_occupancies;} ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> & _mask_betas() { return mask_betas;} ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> & _mask_elements() { return mask_elements;} ///< (Brief description)
            std::vector<boost::dynamic_bitset<>> & _mask_segnames() { return mask_segnames;} ///< (Brief description)

    };

    /// @note Where is my home?
}

template <typename T>
void
sassubset::Mask::set_descriptor_using_mask(const boost::dynamic_bitset<> & mask, std::vector<T> & descriptor, const T value)
{
    for (int i=0; i<_natoms(); ++i)
    {
        if (mask[i])
        {
            descriptor[i] = value;
        }
    }
}
