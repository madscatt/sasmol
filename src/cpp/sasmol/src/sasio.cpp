#include <dcdio.h>
#include <sasio.h>
#include <sasmol.h>
#include <sascalc.h>
#include <util.h>
#include <stdexcept> 

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::string element_from_name(std::string & tmp)
{

    std::size_t found ;
    std::string element_string ;
      
    if(util::has_only_spaces(tmp))
    {
        element_string = "  " ;        
    }
    else
    {
          found = tmp.find_first_not_of(" ");
          element_string = tmp[found] ;
          // std::cout << "first character is " << tmp[found] << std::endl ;
          if(islower(tmp[found+1]))
          {
                 element_string +=  tmp[found+1] ;
                // std::cout << "second character is " << tmp[found+1] << std::endl ;
          }
          std::cout << "ELEMENT STRING = " << element_string << std::endl ;
    }

    return element_string ; 
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::stringstream & operator<<(std::stringstream & ss, float x)
{
   ss.clear();
   ss.str(std::string());
   ss.width(8) ;
   ss.setf(std::ios_base::fixed, std::ios_base::floatfield) ;
   size_t precision;
   float x_abs = abs(x);
   if (x_abs<1.e3) precision=3;
   else if (x_abs>=1.e3 && x_abs<1.e4) precision = (x>=0.0)?3:2;
   else if (x_abs>=1.e4 && x_abs<1.e5) precision = (x>=0.0)?2:1;
   else if (x_abs>=1.e5 && x_abs<1.e6) precision = (x>=0.0)?1:0;
   else if (x>=1.e6 && x<1.e7) precision=0;
   else
   {
       precision = (x>=0.0)?2:1;
       ss.setf(std::ios_base::scientific, std::ios_base::floatfield);
   }
   ss.precision(precision);
   ss << std::right << x;

   return ss;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
resize_array()
{

        int nf = _number_of_frames() ;

        Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> temp_x(_natoms(),nf) ;
        Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> temp_y(_natoms(),nf) ;
        Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic> temp_z(_natoms(),nf) ;

        temp_x = _x() ;
        temp_y = _y() ;
        temp_z = _z() ;

        _x().resize(_natoms(),_number_of_frames()+1) ;
        _y().resize(_natoms(),_number_of_frames()+1) ;
        _z().resize(_natoms(),_number_of_frames()+1) ;

        // see: slice.cpp in retired_ideas/play_with_eigen/
        //    mat1.block(i,j,rows,cols)
        //       a.block(0, 0, 3, 1) = b.block(0,0,3,1) ;

        _x().block(0, 0, _natoms(), nf) = temp_x.block(0,0,_natoms(),nf) ;
        _y().block(0, 0, _natoms(), nf) = temp_y.block(0,0,_natoms(),nf) ;
        _z().block(0, 0, _natoms(), nf) = temp_z.block(0,0,_natoms(),nf) ;

        return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
read_xyz(const std::string &filename)
{

    /* NOT CHECKED FOR MULTIPLE FRAMES */

    std::string line, word ; 

    std::ifstream infile(filename) ;

    std::string local_natoms ;
    infile >> local_natoms ;  

    _natoms() = stoi(local_natoms) ;
    _number_of_frames() = 1 ;    // HARDWIRED

    // #### need to cleanly read the end of the first line
    // ####     and the next line

    getline(infile,line) ;    
    std::istringstream record(line) ;
    while(record >> word) ;
    getline(infile,line) ;    
    while(record >> word) ;

    std::vector<float> vector_x;
    std::vector<float> vector_y;
    std::vector<float> vector_z;

    std::vector<float> coor;
    
    while(getline(infile,line))
    {
        std::istringstream record(line) ;
        record >> word ;
        _atom_name().push_back(word) ;
        while(record >> word)
        {
            coor.push_back(stod(word));
        }
    }

    auto begin = coor.begin() ;
    auto end = coor.end() ;
    std::vector<float>::iterator ii ;

    for(ii = begin ; ii != end ; ii+=3)
    {
          vector_x.push_back(*ii) ; 
          vector_y.push_back(*(ii+1)) ; 
          vector_z.push_back(*(ii+2)); 
    }

    infile.close() ;

    if(_number_of_frames() == 0)
    {
          _x().setZero(_natoms(), 1);
          _y().setZero(_natoms(), 1);
          _z().setZero(_natoms(), 1);
    }
    else
    {
        resize_array() ;
    }

    for(int i = 0 ; i < _natoms() ; ++i)
    {
          _x()(i,_number_of_frames()) = vector_x[i] ;
          _y()(i,_number_of_frames()) = vector_y[i] ;
          _z()(i,_number_of_frames()) = vector_z[i] ;
    }

    _number_of_frames()++ ;
    
    return ;

/*
5
  test 
  N         -3.259000       -5.011000      -14.360000
*/    

} // end of read_xyz

/****** dcdio ******/

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
FILE *
//sasio::Files::
sasio::
open_write_dcd_file(const std::string &filename, int &natoms, int &nset)
{

    FILE *outfile = 0 ;

    char *c_filename = new char[filename.size()+1] ;
    c_filename[filename.size()] = 0 ;
    memcpy(c_filename,filename.c_str(),filename.size()) ;

    outfile = open_dcd_write(c_filename) ;

    int header_result = 0 ;

    header_result = sasio::write_dcd_header(outfile,filename,natoms,nset) ;

    return outfile ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
int
//sasio::Files::
sasio::
write_dcd_header(FILE *outfile, const std::string &filename, int &natoms, int &nset)
{

    char *c_filename = new char[filename.size()+1] ;
    c_filename[filename.size()] = 0 ;
    memcpy(c_filename,filename.c_str(),filename.size()) ;
         
    int istart = 0 , nsavc = 1 ;
    double delta = 1.0 ;
    int header_result = 0 ;

    header_result = write_dcdheader(outfile, c_filename, natoms, nset, \
               istart, nsavc, delta);

    return header_result ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
write_dcd_step(FILE *outfile, int &frame, int &step)
{
    int step_result = 0 ;

    step_result = write_dcdstep(outfile,_natoms(),&_x()(0,frame),&_y()(0,frame),&_z()(0,frame),step+1) ;
    
    return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
write_dcd(const std::string &filename)
{
    /*

    UNTESTED!!!

    */
    std::cout << "writing molecule to dcd file" << std::endl ;

    FILE *outfile = sasio::open_write_dcd_file(filename,_natoms(),_number_of_frames()) ;

    int step_result = 0 ;
    int result = 0 ;
    
    std::cout << "writing dcd frames . . . " << std::endl ;
    
    for(int i = 0 ; i < _number_of_frames() ; ++i)
    {
        write_dcd_step(outfile, i, i);
    }

    result = close_dcd_write(outfile);

} // end of write_dcd 

void
sasio::Files::
write_pdb(const std::string &filename, int &frame)
{
    std::ofstream outfile(filename) ;

    std::string time_and_user_string ;
    time_and_user_string = util::time_and_user();
    std::string temp_remark = "REMARK PDB FILE WRITTEN "; 
    temp_remark += time_and_user_string ;
    std::string remark(temp_remark,0,80) ;        
    std::cout << remark << std::endl ;

    outfile << remark << std::endl;
    //std::string dum1 ="         1         2         3         4         5         6         7         8";
    //std::string dum2 ="12345678901234567890123456789012345678901234567890123456789012345678901234567890";
    //outfile << dum1 << std::endl;
    //outfile << dum2 << std::endl;
    
    std::stringstream line ;
    std::stringstream my_stringstream;

    for(int i = 0 ; i < _natoms() ; ++i)
    { 
        line.str( std::string() );
        line.clear() ;
        std::string tmp(_atom_record()[i],0,6) ;
        line << std::setfill(' ') << std::setw(6) << tmp ;    

        if(_atom_index()[i] > 99999)
        {
            std::string tmp2(std::to_string(99999),0,5) ;
            line << std::setfill(' ') << std::setw(5) << std::right << tmp2 ;    
        }
        else if(_atom_index()[i] < -9999)
        {
            std::string tmp2(std::to_string(-9999),0,5) ;
            line << std::setfill(' ') << std::setw(5) << std::right << tmp2 ;    
        }
        else
        {    
            std::string tmp2(std::to_string(_atom_index()[i]),0,5) ;
            line << std::setfill(' ') << std::setw(5) << std::right << tmp2 ;    
        }    

        std::string tmp0 = " ";
        line << tmp0 ;
    
        std::string tmp3(_atom_name()[i],0,4) ;
        line << std::setfill(' ') << std::setw(4) << std::left << tmp3 ;    

        std::string tmp4(_atom_altloc()[i],0,1) ;
        line << std::setw(1) << tmp4 ;    

        std::string tmp5(_atom_resname()[i],0,3) ;
        line << std::setfill(' ') << std::setw(4) << std::left << tmp5 ;    

        std::string tmp6(_atom_chain()[i],0,1) ;
        line << std::setfill(' ') << std::setw(1) << std::left << tmp6 ;    

        if(_atom_resid()[i] > 9999)
        {
            std::string tmp7(std::to_string(9999),0,4) ;
            line << std::setfill(' ') << std::setw(4) << std::right << tmp7 ;    
        }
        else if(_atom_resid()[i] < -999)
        {
            std::string tmp7(std::to_string(-999),0,4) ;
            line << std::setfill(' ') << std::setw(4) << std::right << tmp7 ;    
        }
        else
        {    
            std::string tmp7(std::to_string(_atom_resid()[i]),0,4) ;
            line << std::setfill(' ') << std::setw(4) << std::right << tmp7 ;    
        }    

        std::string tmp8(_atom_icode()[i],0,1) ;
        line << std::setw(4) << std::left << tmp8 ;    

        //std::string f_str = std::to_string(f);

        //here
        //
        my_stringstream << _x()(i,frame) ;
        //my_stringstream << std::to_string(_x()(i,frame)) ;
        //line << std::setfill(' ') << std::setw(8) << std::right << my_stringstream.str() ;
        line << std::setfill(' ') << std::setw(8) << std::right << my_stringstream.str() ;
        
        my_stringstream << _y()(i,frame) ;
        line << my_stringstream.str() ;

        my_stringstream << _z()(i,frame) ;
        line << my_stringstream.str() ;

        std::string tmp12(_atom_occupancy()[i],0,6) ;
        line << std::setfill(' ') << std::setw(6) << std::right << tmp12 ;    

        std::string tmp13(_atom_beta()[i],0,6) ;
        line << std::setfill(' ') << std::setw(6) << std::right << tmp13 ;    

        std::string tmp14(_atom_segname()[i],0,4) ;
        line << std::setfill(' ') << std::setw(10) << std::right << tmp14 ;    

        std::string tmp15(_atom_element()[i],0,2) ;
        line << std::setfill(' ') << std::setw(2) << std::right << tmp15 ;    

        std::string tmp16(_atom_charge()[i],0,2) ;
        line << std::setfill(' ') << std::setw(2) << std::right << tmp16 ;    
        
        outfile << line.str() << std::endl ;
    }

    outfile << "END" << std::endl;
    outfile.close() ;
}


///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
_set_unique_attributes()
{
    _unique_record() = util::unique_vector_string(_atom_record());
    _unique_index() = util::unique_vector_int(_atom_index());
    _unique_name() = util::unique_vector_string(_atom_name());
    _unique_altloc() = util::unique_vector_string(_atom_altloc());
    _unique_resname() = util::unique_vector_string(_atom_resname());
    _unique_chain() = util::unique_vector_string(_atom_chain());
    _unique_resid() = util::unique_vector_int(_atom_resid());
    _unique_icode() = util::unique_vector_string(_atom_icode());
    _unique_occupancy() = util::unique_vector_string(_atom_occupancy());
    _unique_beta() = util::unique_vector_string(_atom_beta());
    _unique_segname() = util::unique_vector_string(_atom_segname());
    _unique_element() = util::unique_vector_string(_atom_element());
    _unique_selement() = util::unique_vector_string(_atom_selement());
    _unique_charge() = util::unique_vector_string(_atom_charge());
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
read_pdb(const std::string &filename)
{
    std::string line, word ;

    std::string element_string ;

    std::ifstream infile(filename) ;

    int count = 0 ;
    int frame = 0 ;

    //string s(line, index, number_characters) ;

    std::vector<float> vector_x;
    std::vector<float> vector_y;
    std::vector<float> vector_z;
  
    while(getline(infile,line))
    {
         std::istringstream record(line) ;
         record >> word ;

         std::string temp1(word,0,5) ;
             
         if(temp1 != "ATOM" && temp1 != "HETAT")
         {
             std::cout << "excluding: " << word << std::endl ;
         }
         else
         {
             std::string tmp(line,0,6) ;
             _atom_record().push_back(tmp) ;

             std::string tmp2(line,6,5) ;
             _atom_index().push_back(stoi(tmp2)) ;

             std::string tmp3(line,12,4) ;
             _atom_name().push_back(tmp3) ;

             std::string tmp4(line,16,1) ;
             _atom_altloc().push_back(tmp4) ;

             std::string tmp5(line,17,3) ;
             _atom_resname().push_back(tmp5) ;

             std::string tmp6(line,21,1) ;
             _atom_chain().push_back(tmp6) ;

             std::string tmp7(line,22,4) ;
             _atom_resid().push_back(stoi(tmp7)) ;

             std::string tmp8(line,26,1) ;
             _atom_icode().push_back(tmp8) ;

             std::string tmp9(line,30,8) ;
             vector_x.push_back(stof(tmp9)) ;

             std::string tmp10(line,38,8) ;
             vector_y.push_back(stof(tmp10)) ;

             std::string tmp11(line,46,8) ;
             vector_z.push_back(stof(tmp11)) ;

             try
             {
                 std::string tmp12(line,54,6) ;
                 if(util::has_only_spaces(tmp12))
                 {
                     _atom_occupancy().push_back("0.00") ;
                 }
                 else
                 {
                     _atom_occupancy().push_back(tmp12) ;
                 }
             }
             catch(const std::out_of_range& oor)
             {
                 _atom_occupancy().push_back("0.00") ;
                 std::cerr<<"Occupancy: Out of range error: "<< oor.what() <<std::endl ;
             }

             try
             {
                 std::string tmp13(line,60,6) ;
                 if(util::has_only_spaces(tmp13))
                 {
                     _atom_beta().push_back("0.00") ;
                 }
                 else
                 {
                     _atom_beta().push_back(tmp13) ;
                 }
             }
             catch(const std::out_of_range& oor)
             {
                 _atom_beta().push_back("0.00") ;
                 std::cerr<<"Beta: Out of range error: "<< oor.what() <<std::endl ;
             }

             try
             {
                 std::string tmp14(line,72,4) ;
                 if(util::has_only_spaces(tmp14))
                 {
                     _atom_segname().push_back("SEGN") ;
                 }
                 else
                 {        
                     _atom_segname().push_back(tmp14) ;
                 }
             }
             catch(const std::out_of_range& oor)
             {
                 _atom_segname().push_back("SEGN") ;
                 std::cerr<<"Segname: Out of range error: "<< oor.what() <<std::endl ;
             }

             try
             {
                 std::string tmp15(line,76,2) ;
                 if(util::has_only_spaces(tmp15))
                 {
                     std::cout << "Element not found" << std::endl;
                     element_string = element_from_name(tmp3) ;
                     _atom_element().push_back(element_string) ;
                 }
                 else
                 {
                     _atom_element().push_back(tmp15) ;
                 }
             }
             catch(const std::out_of_range& oor)
             {
                 element_string = element_from_name(tmp3) ;
                 _atom_element().push_back(element_string) ;
                 std::cerr<<"Element: Out of range error: "<< oor.what() <<std::endl ;
             }

             try
             {
                 std::string tmp16(line,78,2) ;
                 if(util::has_only_spaces(tmp16))
                 {
                     _atom_charge().push_back(" ") ;    
                 }
                 else
                 {
                     _atom_charge().push_back(tmp16) ;
                 }
             }
             catch(const std::out_of_range& oor)
             {
                 _atom_charge().push_back(" ") ;    
                 std::cerr<<"Charge: Out of range error: "<< oor.what() <<std::endl ;
             }

             count++ ;
         }
    }

    _natoms() = count ;
    infile.close() ;

    int nf = _number_of_frames() ;

    if(_number_of_frames() == 0)
    {
            _x().setZero(_natoms(), 1);
            _y().setZero(_natoms(), 1);
            _z().setZero(_natoms(), 1);
      }
      else
      {
        resize_array() ;
      }

      for(int i = 0 ; i < _natoms() ; ++i)
      {
            _x()(i,_number_of_frames()) = vector_x[i] ;
            _y()(i,_number_of_frames()) = vector_y[i] ;
            _z()(i,_number_of_frames()) = vector_z[i] ;
      }

      std::vector<std::string> s_element ;
      s_element = util::strip_white_space(_atom_element()) ;
      _atom_selement() = s_element ;

      dynamic_cast<sasmol::SasMol*>(this)->calc_mass() ;
      _atom_com() = dynamic_cast<sasmol::SasMol*>(this)->calc_com(frame) ;

      _number_of_frames() += 1 ;

      _set_unique_attributes();

/*
          1         2         3         4         5         6         7         8
012345678901234567890123456789012345678901234567890123456789012345678901234567890
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
0    6    1 2   7
ATOM      1  N   GLY X   1     -21.525 -67.562  86.759  1.00  0.00      GAG  N
ATOM      2  HT1 GLY X   1     -22.003 -68.460  86.892  1.00  0.00      GAG  H
ATOM      3  HT2 GLY X   1     -21.905 -66.929  87.525  1.00  0.00      GAG  H

*/
}


///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
FILE *
//sasio::Files::
sasio::
open_dcd_read(std::string &filename, int &natoms, int &nset, int &reverseEndian) {

    FILE *dcd_infile = 0 ;

    char *c_filename = new char[filename.size()+1] ;
    c_filename[filename.size()] = 0 ;
    memcpy(c_filename,filename.c_str(),filename.size()) ;

    dcd_infile = ::open_dcd_read(c_filename) ; /// @note I have to cite the global declaration

    int num_fixed = 0 ; 
    int charmm = 0 ;
    int result = 1 ;

    int read_header_result, istart, nsavc, namnf ;
    double delta ; 

    read_header_result = read_dcdheader(dcd_infile,&natoms,&nset,&istart,&nsavc,&delta,&namnf,&reverseEndian,&charmm) ;

    return dcd_infile ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
read_dcd_step(FILE *dcd_infile, int &frame, int &natoms, int &reverseEndian)
{

    int charmm = 0 ;
    int num_fixed = 0 ;
    int result ;

    std::vector<float> vector_x(natoms) ;
    std::vector<float> vector_y(natoms) ;
    std::vector<float> vector_z(natoms) ;

    result = read_dcdstep(dcd_infile,natoms,&vector_x[0],&vector_y[0],&vector_z[0],num_fixed,frame,reverseEndian,charmm) ;

    if(_x().cols() < frame)
    {
        std::cout << "resizing array: cols() =  " << _x().cols() << " frame = " << frame  << std::endl ;
        resize_array() ;
    }

    for(int i = 0 ; i < _natoms() ; ++i)
    {
          _x()(i,frame) = vector_x[i] ;
          _y()(i,frame) = vector_y[i] ;
          _z()(i,frame) = vector_z[i] ;
    } 

    _number_of_frames()++ ;

    return ;
}


///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasio::Files::
read_dcd(std::string &dcd_input_file_name)
{

    int input_natoms, input_nset, input_reverseEndian ;
    FILE *dcd_file_pointer ;

    dcd_file_pointer = sasio::open_dcd_read(dcd_input_file_name, input_natoms, input_nset, input_reverseEndian) ;

    int input_frame = 0 ;

    _number_of_frames() = 0 ;

    _x().setZero(_natoms(), input_nset);
    _y().setZero(_natoms(), input_nset);
    _z().setZero(_natoms(), input_nset);

    std::cout << "_x().cols() = " << _x().cols() << std::endl ;
    std::cout << "input_nset = " << input_nset << std::endl ;

    for(int i = 0 ; i < input_nset ; ++i)
    {
        std::cout << "." ;
        read_dcd_step(dcd_file_pointer, i, input_natoms, input_reverseEndian) ;
    }

    int result = close_dcd_read(dcd_file_pointer) ;

    std::cout << std::endl ;

    return ;
}

/*

to do:

    def read_single_dcd_step(self,filename,frame) 
        --> input_filter, nmer_dihedral, dihedral_monte_carlo, open_minimize

    def write_dcd_frames(self, filename, start, end) --> not used

*/

