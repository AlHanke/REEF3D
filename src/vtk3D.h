/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef VTK3D_H_
#define VTK3D_H_

#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>

class lexer;
class fdm;
class fdm_fnpf;
class fdm_nhf;
class ghostcell;

class vtk3D
{
    public:
        virtual void folder(const char*){};
        virtual void offset(lexer*,int*,int&){};
        virtual void structureWrite(lexer*, fdm*, std::vector<char>&, int&){};
        virtual void structureWrite(lexer*, fdm_fnpf*, std::vector<char>&, int&){};
        virtual void structureWrite(lexer*, fdm_nhf*, std::vector<char>&, int&){};
        virtual void extent(lexer*,ghostcell*){};
        virtual void fileName(char*, const char*, const int, const int){};
        virtual void parallelFileName(char*, const char*, const int){};
        virtual void beginning(lexer*, std::stringstream&){};
        virtual void beginningParallel(lexer*, std::ofstream&){};
        virtual void ending(std::stringstream&, const int*, int&){};
        virtual void endingParallel(std::ofstream&, const char*, const int, const int){};
    protected:
        void vtkVersion(std::ofstream &result){result<<"version=\"1.0\" byte_order=\"LittleEndian\">\n";}; // header_type=\"UInt32\"
        void vtkVersion(std::stringstream &result){result<<"version=\"1.0\" byte_order=\"LittleEndian\">\n";}; // header_type=\"UInt32\"
        void xmlVersion(std::ofstream &result){result<<"<?xml version=\"1.0\"?>\n";};
        void xmlVersion(std::stringstream &result){result<<"<?xml version=\"1.0\"?>\n";};
        void timeValue(std::ofstream &result, const double time){result<<"<FieldData>\n<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<std::setprecision(7)<<time<<"\n</DataArray>\n</FieldData>\n";};
        void timeValue(std::stringstream &result, const double time){result<<"<FieldData>\n<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<std::setprecision(7)<<time<<"\n</DataArray>\n</FieldData>\n";};
        void appendData(std::stringstream &result){result<<"<AppendedData encoding=\"raw\">\n_";};
};

#endif