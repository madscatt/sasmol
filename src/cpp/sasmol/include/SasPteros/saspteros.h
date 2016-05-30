/*
 * Modified and adpated into SASSIE by Hailiang Zhang in April 2014
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#ifndef SELECTION_H
#define SELECTION_H

#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/signals2.hpp>
#include <boost/shared_ptr.hpp>

#include "saspteros_parser.h"
#include <sasatomBase.h>

// forward declaration
namespace sasmol
{
    class SasMol;
}


// Forward declaration of frient class
class SasPteros_parser;

/** @brief SasPteros class.
*
*   Technically speaking the saspteros
*   is just an array, which contains indexes of selected atoms in particlar sasmol.
*   SasPteros does not hold the copies of the atoms or their coordinates, it
*   just points to them serving like a handy alias for certain subset of atoms.
*   Thus saspteross may overlap arbitrarily.
*   SasPteross are used to perform various operations on the group of selected atoms.
*   The changes become immediately visible to all other saspteross, which point to
*   some of changed atoms.
*   Each saspteros is bound to particular SasMol. There are neither 'parentless' saspteros nor the
*   saspteross, which combine the atoms from different sasmols.
*   SasPteross are created using the syntax, which is very similar to those used in VMD.
*/
class SasPteros {
  /// SasMol and SasPteros are friends because they are closely integrated.
  friend class sasmol::SasMol;
  friend class SasPteros_parser;

  public:
    /// Ensure correct 16-bytes-alignment for Eigen sse2 optimizations
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** Main constructor.
    *   @param mol SasMol pointed by this saspteros
    *   @param str SasPteros string
    *   if with_signal is true SasMol with send signals to saspteros automatically
    */
    SasPteros(sasmol::SasMol& mol, std::string str);

    /** Constructor with delayed parsing.
    *   Associates saspteros with the sasmol @param mol,
        but does not parse saspteros.
    *   SasPteros text should be passed later by overloaded << operator or by
    *   calling modify() function.
    *   if with_signal is true SasMol with send signals to saspteros automatically
    */
    SasPteros(sasmol::SasMol& mol);

    /// Default constructor for absolutely empty saspteros.
    SasPteros();

    /** Constructor, which creates saspteros from the interval of indexes
        instead of saspteros string.
        It is much faster then parsing corresponding string, but is limited
        to contigous interval of indexes.
        if with_signal is true SasMol with send signals to saspteros automatically
     */
    SasPteros(sasmol::SasMol& mol, int ind1, int ind2);

    /// Assignment operator
    SasPteros& operator=(SasPteros);    

    /// Copy constructor
    SasPteros(const SasPteros&);

    /// Equality operators
    bool operator==(const SasPteros &other) const;

    bool operator!=(const SasPteros &other) const {
        return !(*this == other);
    }

    /// Destructor
    ~SasPteros();

    /** Recomputes saspteros without re-parsing saspteros text.
    *   Only makes sense for coordinate-dependent saspteross when the coordinates change.
    *   Called authomatically for coordinate-dependent saspteross by set_frame()
    *   If saspteros is not coordinate-dependent does nothing.
    */
    void apply();

    /// @name Frame functions
    /// @{

    /// Get current frame saspteros is pointing to
    int get_frame() const {return frame;}
    /// Set current frame for saspteros
    void set_frame(int fr);
    /// @}

    /** Clears saspteros and frees memory, but do not delete it.
    *   SasPteros remains registered in parent sasmol, but is cleared from any data.
    *   Subsequent call of modify() may populate it again.
    */
    void clear();

public:
    // Row text of saspteros
    std::string sel_text;
    // Used with << operator
    std::ostringstream ss;
    // Indexes of atoms in saspteros
    std::vector<int> index;
    // Pointer to target sasmol
    sasmol::SasMol* sasmol;

    // Stores current frame
    int frame;

    // Holds an instance of saspteros parser
    boost::shared_ptr<SasPteros_parser> parser;
    void allocate_parser();

    // Private functions for creating saspteros
    void create_internal(sasmol::SasMol& mol, std::string& str);
    void create_internal(sasmol::SasMol& mol, int ind1, int ind2);
    // Private function for deleting saspteros
    void delete_internal();    

    /// Get vector of all resnames in selection
    std::vector<std::string> get_resname() const;
};


#endif /* SELECTION_H */
