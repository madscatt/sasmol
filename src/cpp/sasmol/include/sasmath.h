#pragma once
#include <iostream>
#include "Eigen/Dense"
#include "sasatomBase.h"

// forward declaration
namespace sasmol
{
    class SasMol;
}

namespace sasmath
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Math: public sasmol::SasAtomBase
    {
        public:
            Eigen::Matrix3f find_u(sasmol::SasMol &mol, int &frame) ; ///< (Brief description)
    };

    /// @note Where is my home?
    float signed_angle(Eigen::Matrix<float,3,1> &coor1, Eigen::Matrix<float,3,1> &coor2, Eigen::Matrix<float,3,1> &coor3) ; ///< (Brief description)
    float dihedral_angle(Eigen::Matrix<float,3,1> &coor1, Eigen::Matrix<float,3,1> &coor2, Eigen::Matrix<float,3,1> &coor3, Eigen::Matrix<float,3,1> &coor4) ; ///< (Brief description)
    float calc_angle(Eigen::Matrix<float,3,1> &coor1, Eigen::Matrix<float,3,1> &coor2, Eigen::Matrix<float,3,1> &coor3) ; ///< (Brief description)
}
