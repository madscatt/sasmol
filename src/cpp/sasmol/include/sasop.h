#pragma once
#include <iostream>
#include <cmath>

#include "sasatomBase.h"

// forward declaration
namespace sasmol
{
    class SasMol;
}

namespace sasop
{
    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Move: public sasmol::SasAtomBase
    {
        public:
            void translate(int &frame, std::vector<float> &value) ; ///< (Brief description)

            void center(int &frame) ; ///< (Brief description)

            void move_to(int &frame, std::vector<float> &value) ; ///< (Brief description)

            void align(sasmol::SasMol &mol2, int &frame) ; ///< (Brief description)

            void rotate(int &frame, std::string &axis, float &theta) ; ///< (Brief description)

            void general_axis_rotate(int &frame, float &theta, float &ux, float &uy, float &uz) ; ///< (Brief description)
    };
}


