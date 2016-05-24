#include <sasmol.h>

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
sasmol::SasMol::
SasMol(const int Natoms, const int Nframes,
       const long long *const  p_index, const std::string * const p_name, const std::string * const p_altloc, const std::string * const p_resname, const std::string * const p_chain, const long long *const p_resid, const std::string * const p_icode,
       const double * p_x, const double * p_y, const double * p_z,
       const std::string * const p_occupancy, const std::string * const p_beta, const std::string * const p_segname, const std::string * const p_element, const std::string * const p_charge)
{
    natoms = Natoms;
    number_of_frames = Nframes;

    _x().setZero(natoms,number_of_frames);
    _y().setZero(natoms,number_of_frames);
    _z().setZero(natoms,number_of_frames);
    for (int i=0; i<natoms; ++i)
    {
        _atom_index().push_back(p_index[i]);
        _atom_name().push_back(p_name[i]);
        _atom_altloc().push_back(p_altloc[i]);
        _atom_resname().push_back(p_resname[i]);
        _atom_chain().push_back(p_chain[i]);
        _atom_resid().push_back(p_resid[i]);
        _atom_icode().push_back(p_icode[i]);
        for (int j=0; j<number_of_frames; ++j)
        {
            _x()(i,j) = p_x[j*natoms + i];
            _y()(i,j) = p_y[j*natoms + i];
            _z()(i,j) = p_z[j*natoms + i];
        }
        _atom_occupancy().push_back(p_occupancy[i]);
        _atom_beta().push_back(p_beta[i]);
        _atom_segname().push_back(p_segname[i]);
        _atom_element().push_back(p_element[i]);
        _atom_charge().push_back(p_charge[i]);
    }
}
