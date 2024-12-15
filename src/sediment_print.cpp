/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"bedshear.h"
#include<sstream>
#include<vector>
#include<cstring>

void sediment_f::name_ParaView_parallel_bedload(lexer *p, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_qbe\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"ST_qb\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"ST_cbe\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"ST_cb\"/>\n";
}

void sediment_f::name_ParaView_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_qbe\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_qb\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_cbe\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_cb\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void sediment_f::name_ParaView_bedload(lexer *p, ghostcell *pgc, stringstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_qbe\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_qb\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_cbe\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_cb\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void sediment_f::offset_ParaView_2D_bedload(lexer *p, ghostcell *pgc, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_ParaView_bedload(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}

void sediment_f::print_3D_bedload(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;

	// qbe
    pgc->gcsl_start4(p,s->qbe,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->qbe));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // qb
    pgc->gcsl_start4(p,s->qb,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->qb));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // cbe
    pgc->gcsl_start4(p,s->cbe,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->cbe));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // cb
    pgc->gcsl_start4(p,s->cbe,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->cb));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
}

void sediment_f::print_2D_bedload(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;

	// qbe
    pgc->gcsl_start4(p,s->qbe,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->qbe));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // qb
    pgc->gcsl_start4(p,s->qb,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->qb));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // cbe
    pgc->gcsl_start4(p,s->cbe,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->cbe));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // cb
    pgc->gcsl_start4(p,s->cbe,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->cb));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
}

void sediment_f::name_ParaView_parallel_bedshear(lexer *p, ghostcell *pgc, ofstream &result)
{
    if(p->P79==1)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_tau_eff\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"ST_tau_crit\"/>\n";
    }
    
    if(p->P79==2)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shearvel_eff\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shearvel_crit\"/>\n";
    }
    
    if(p->P79==3)
    {
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shields_eff\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"ST_shields_crit\"/>\n";
    }
}

void sediment_f::name_ParaView_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    if(p->P79==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_eff\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_crit\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    
    if(p->P79==2)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_eff\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_crit\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    
    if(p->P79==3)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_eff\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_crit\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
}

void sediment_f::name_ParaView_bedshear(lexer *p, ghostcell *pgc, stringstream &result, int *offset, int &n)
{
    if(p->P79==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_tau_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==2)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shearvel_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    
    if(p->P79==3)
    {
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_eff\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_shields_crit\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
}

void sediment_f::offset_ParaView_2D_bedshear(lexer *p, ghostcell *pgc, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_ParaView_bedshear(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}

void sediment_f::print_2D_bedshear(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;
    
    if(p->P79==1)
    {
        
	// tau_eff
    pgc->gcsl_start4(p,s->tau_eff,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->tau_eff));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // tau_crit
    pgc->gcsl_start4(p,s->tau_crit,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->tau_crit));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    
    if(p->P79==2)
    {
	// shearvel_eff
    pgc->gcsl_start4(p,s->shearvel_eff,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_eff));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // shearvel_crit
    pgc->gcsl_start4(p,s->shearvel_crit,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_crit));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    
    if(p->P79==3)
    {
	// shields_eff
    pgc->gcsl_start4(p,s->shields_eff,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shields_eff));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // shields_crit
    pgc->gcsl_start4(p,s->shields_crit,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->shields_crit));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    }
}

void sediment_f::print_3D_bedshear(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;
    
    if(p->P79==1)
    {
        
	// tau_eff
    pgc->gcsl_start4(p,s->tau_eff,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->tau_eff));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // tau_crit
    pgc->gcsl_start4(p,s->tau_crit,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->tau_crit));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    
    if(p->P79==2)
    {
	// shearvel_eff
    pgc->gcsl_start4(p,s->shearvel_eff,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_eff));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // shearvel_crit
    pgc->gcsl_start4(p,s->shearvel_crit,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shearvel_crit));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    
    if(p->P79==3)
    {
	// shields_eff
    pgc->gcsl_start4(p,s->shields_eff,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shields_eff));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // shields_crit
    pgc->gcsl_start4(p,s->shields_crit,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->shields_crit));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    }
}

void sediment_f::print_2D_parameter1(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;
    
    // alpha
    pgc->gcsl_start4(p,s->alpha,1);
	
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->alpha));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    
    // teta
    pgc->gcsl_start4(p,s->teta,1);
	
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->teta));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // gamma
    pgc->gcsl_start4(p,s->gamma,1);
	
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->gamma));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // beta
    pgc->gcsl_start4(p,s->beta,1);
	
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->beta));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // phi
    pgc->gcsl_start4(p,s->phi,1);
	
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->phi));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
}

void sediment_f::print_3D_parameter1(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;
    
    // alpha
    pgc->gcsl_start4(p,s->alpha,1);
	
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->alpha));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    
    // teta
    pgc->gcsl_start4(p,s->teta,1);
	
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->teta));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // gamma
    pgc->gcsl_start4(p,s->gamma,1);
	
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->gamma));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // beta
    pgc->gcsl_start4(p,s->beta,1);
	
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->beta));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // phi
    pgc->gcsl_start4(p,s->phi,1);
	
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->phi));
    ffn*=(180.0/PI);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
}

void sediment_f::name_ParaView_parallel_parameter1(lexer *p, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_alpha\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_teta\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_gamma\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_beta\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_phi\"/>\n";
}

void sediment_f::name_ParaView_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_alpha\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_teta\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_gamma\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_beta\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_phi\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void sediment_f::name_ParaView_parameter1(lexer *p, ghostcell *pgc, stringstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_alpha\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_teta\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_gamma\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_beta\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_phi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void sediment_f::offset_ParaView_2D_parameter1(lexer *p, ghostcell *pgc, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_ParaView_parameter1(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}


void sediment_f::print_2D_parameter2(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;
	
    // dh,bedch,reduce,threshold,slideflag
    
    // dh
    pgc->gcsl_start4(p,s->vz,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->vz));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // bedchange
    pgc->gcsl_start4(p,s->bedch,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->bedch));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // reduce
    pgc->gcsl_start4(p,s->reduce,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->reduce));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // threshold
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=0.0;
    
    if(s->tau_eff(i,j)>s->tau_crit(i,j))
    ffn=1.0;
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // slideflag
    pgc->gcsl_start4(p,s->slideflag,1);
    
	iin=4*(p->pointnum2D);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPSLICELOOP
	{
    ffn=float(p->sl_ipol4(s->slideflag));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
}

void sediment_f::print_3D_parameter2(lexer* p, ghostcell *pgc, std::vector<char> &buffer, int &m)
{	
	float ffn;
	int iin;
	
    // dh,bedch,reduce,threshold,slideflag
    
    // dh
    pgc->gcsl_start4(p,s->vz,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->vz));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // bedchange
    pgc->gcsl_start4(p,s->bedch,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->dtsed*p->sl_ipol4(s->bedch));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // reduce
    pgc->gcsl_start4(p,s->reduce,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->reduce));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // threshold
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=0.0;
    
    if(s->tau_eff(i,j)>s->tau_crit(i,j))
    ffn=1.0;
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    // slideflag
    pgc->gcsl_start4(p,s->slideflag,1);
    
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
    ffn=float(p->sl_ipol4(s->slideflag));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
}

void sediment_f::name_ParaView_parallel_parameter2(lexer *p, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"ST_dh\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_bedch\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_reduce\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_threshold\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"ST_slideflag\"/>\n";
}

void sediment_f::name_ParaView_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_dh\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_bedch\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_reduce\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_threshold\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_slideflag\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void sediment_f::name_ParaView_parameter2(lexer *p, ghostcell *pgc, stringstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"ST_dh\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_bedch\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_reduce\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_threshold\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"ST_slideflag\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void sediment_f::offset_ParaView_2D_parameter2(lexer *p, ghostcell *pgc, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sediment_f::offset_ParaView_parameter2(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}
