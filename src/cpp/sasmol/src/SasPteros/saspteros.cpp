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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <queue>
#include <boost/algorithm/string.hpp> // String algorithms
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <sasmol.h>

#include <SasPteros/saspteros.h>
#include <SasPteros/saspteros_error.h>
#include <SasPteros/saspteros_macro.h>
//#include "grid_search.h"

#include <Eigen/Geometry>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

void SasPteros::allocate_parser(){
    // Parse saspteros here
    // Parser is heavy object, so if saspteros is not persistent
    // we will delete it after parsing is complete
    parser.reset();
    parser = boost::shared_ptr<SasPteros_parser>(new SasPteros_parser);
    parser->create_ast(sel_text);    
    parser->apply(sasmol, frame, index);
    if(!parser->has_coord){
        parser.reset();
    }
}

SasPteros::SasPteros(){
    sasmol = NULL;
    parser.reset();
};

// Aux function, which creates saspteros
void SasPteros::create_internal(sasmol::SasMol& mol, string& str){
    // Set saspteros string
    sel_text = str;
    boost::trim(sel_text);

    // Expand macro-definitions in the string
    for(int i=0;i<Nmacro;++i)
        boost::replace_all(sel_text,macro[2*i],macro[2*i+1]);

    // Add saspteros to mol and make connection
    sasmol = &mol;

    // By default points to frame 0
    frame = 0;

    allocate_parser();
}

// Aux function, which creates saspteros
void SasPteros::create_internal(sasmol::SasMol& mol, int ind1, int ind2){
    // Set saspteros string
    sel_text = "index "+boost::lexical_cast<string>(ind1)+"-"+boost::lexical_cast<string>(ind2);

    // Add saspteros to mol and save self-pointer
    sasmol = &mol;

    // By default points to frame 0
    frame = 0;
    // No parser needed
    parser.reset();

    // Populate saspteros directly
    index.clear();
    for(int i=ind1; i<=ind2; ++i) index.push_back(i);
}


// Main constructor
SasPteros::SasPteros(sasmol::SasMol& mol, string str){
    create_internal(mol, str);
}

// Constructor without immediate parsing
SasPteros::SasPteros(sasmol::SasMol& mol){
    // Add saspteros to mol and save self-pointer
    sasmol = &mol;   

    frame = 0;
    parser.reset();
}

SasPteros::SasPteros(sasmol::SasMol& mol, int ind1, int ind2){
    create_internal(mol, ind1, ind2);
}


// Destructor
SasPteros::~SasPteros(){
    // Delete parser if it is persistent
    parser.reset();
    // All the rest will be destroyed automatically
}

// Free memory used by saspteros.
// SasPteros is still present in parent.
void SasPteros::clear(){
    // Clear index
    index.clear();
    // If parser is present (persistent), delete it
    parser.reset();
}


// Re-apply AST tree for coordinate-dependent saspteross
void SasPteros::apply(){
    if(parser){
        // If parser is persistent, do quick eval using ast tree
        parser->apply(sasmol, frame, index);
    }
}


void SasPteros::set_frame(int fr){
    if(fr>=0 && fr < sasmol->x.rows() ){
        frame = fr;
        // If parser is persistent, do quick update
        // This will only work for coordinate-dependent saspteross
        apply();
    } else {
        throw SasPteros_error("Invalid frame to set!");
    }
}

vector<string> SasPteros::get_resname() const {
    vector<string> res;
    int i,n;
    n = index.size();
    res.resize(n);
    for(i=0; i<n; ++i) res[i] = sasmol->_atom_resname()[index[i]];
    return res;
}
