/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/
/*
#include"mooring_DGSEM.h"
#include"lexer.h"

void mooring_DGSEM::facenodegen(lexer *p, double xS, double xE)
{
	p->Darray(vx, H+1);

	for (int i = 0; i < (H+1); i++)
	{
		vx[i] = xS + i*(xE-xS)/H;
	}
}

void mooring_DGSEM::getqp(lexer *p)
{
	p->Darray(r, P+1);	

	if (P == 1) 
	{
		r[0] = -1.0;
		r[1] = 1.0;
	}
	else if (P == 2)
	{
		r[0] = -1.0;
		r[1] = 0.0;
		r[2] = 1.0;
	}
	else if (P == 3)
	{
		r[0] = -1.0;
		r[1] = -0.447213595499958;
		r[2] = 0.447213595499958;
		r[3] = 1.0;
	}
	else if (P == 4)
	{
		r[0] = -1.0;
		r[1] = -0.654653670707977;
		r[2] = 0.0;
		r[3] = 0.654653670707977;
		r[4] = 1.0;
	}
	else if (P == 5)
	{
		r[0] = -1.0;
		r[1] = -0.765055323929465;
		r[2] = -0.285231516480645;
		r[3] = 0.285231516480645;
		r[4] = 0.765055323929465;
		r[5] = 1.0;
	}
	else if (P == 6)
	{
		r[0] = -1.0;
		r[1] = -0.830223896278567;
		r[2] = -0.468848793470714;
		r[3] = 0.0;
		r[4] = 0.468848793470714;
		r[5] = 0.830223896278567;
		r[6] = 1.0;
	}
	else 
	{
		cout<<"Add quadrature points"<<endl;
	}
}


void mooring_DGSEM::getDr(lexer *p)
{
	p->Darray(Dr,P+1,P+1);

	if (P == 1) 
	{
		Dr[0][0] = -0.5;
		Dr[0][1] = 0.5;
	}
	else if (P == 2)
	{
		Dr[0][0] = -1.5;
		Dr[0][1] = 2.0;
		Dr[0][2] = -0.5;
		Dr[1][1] = 0.0;
		Dr[1][2] = 0.5;
	}
	else if (P == 3)
	{
		Dr[0][0] = -3.0;
		Dr[0][1] = 4.045084971874735;
		Dr[0][2] = -1.545084971874736;
		Dr[0][3] = 0.5;
		Dr[1][1] = 0.0;
		Dr[1][2] = 1.118033988749895;
		Dr[1][3] = -0.309016994374947;
		Dr[2][2] = 0.0;
		Dr[2][3] = 0.809016994374947;
	}
	else if (P == 4)
	{
		Dr[0][0] = -5.0;
		Dr[0][1] = 6.756502488724239;  
		Dr[0][2] = -2.66666666666666;   
		Dr[0][3] = 1.410164177942428; 
		Dr[0][4] = -0.5;
		Dr[1][1] = 0.0;
		Dr[1][2] = 1.745743121887938;  
		Dr[1][3] = -0.763762615825973;   
		Dr[1][4] = 0.259009746969017;
		Dr[2][2] = 0.0;
		Dr[2][3] = 1.336584577695453;
		Dr[2][4] = -0.375;
		Dr[3][3] = 0.0;
		Dr[3][4] = 1.240990253030984;
	}
	else if (P == 5)
	{
		Dr[0][0] = -7.5;
		Dr[0][1] = 10.141415936319669;  
		Dr[0][2] = -4.036187270305346;  
		Dr[0][3] = 2.244684648176166;    
		Dr[0][4] = -1.349913314190484;   
		Dr[0][5] = 0.5;
		Dr[1][1] = 0.0;
		Dr[1][2] = 2.523426777429454;  
		Dr[1][3] = -1.152828158535929;    
		Dr[1][4] = 0.653547507429800; 
		Dr[1][5] = -0.237781177984231;
		Dr[2][2] = 0.0;
		Dr[2][3] = 1.752961966367866;  
		Dr[2][4] = -0.786356672223241;   
		Dr[2][5] = 0.269700610832039;
		Dr[3][3] = 0.0;
		Dr[3][4] = 1.721256952830234;
		Dr[3][5] = -0.484951047853569;
		Dr[4][4] = 0.0;
		Dr[4][5] = 1.786364948339093;
	}
	else if (P == 6)
	{
		Dr[0][0] = -10.5;
		Dr[0][1] = 14.201576602919813;  
		Dr[0][2] = -5.668985225545485;   
		Dr[0][3] = 3.199999999999988;  
		Dr[0][4] = -2.049964813076737;        
		Dr[0][5] = 1.31737343570243;
		Dr[0][6] = -0.5;
		Dr[1][1] = 0.0;
		Dr[1][2] = 3.455828214294279;  
		Dr[1][3] = -1.598606688098368;   
		Dr[1][4] = 0.961339797288714;     
		Dr[1][5] = -0.602247179635786;
		Dr[1][6] = 0.226611870395445;
		Dr[2][2] = 0.0;
		Dr[2][3] = 2.266698087086001;  
		Dr[2][4] = -1.066441904006376;    
		Dr[2][5] = 0.616390835517580; 
		Dr[2][6] = -0.226099400942575;
		Dr[3][3] = 0.0;
		Dr[3][4] = 2.006969240588754;  
		Dr[3][5] = -0.907544471268821;
		Dr[3][6] = 0.3125;
		Dr[4][4] = 0.0;
		Dr[4][5] = 2.215804283169971;
		Dr[4][6] = -0.625256665515342;
		Dr[5][5] = 0.0;
		Dr[5][6] = 2.44292601424428;
	}
	
	// Use anti-symmetry to fill rest of matrices
	for (int i = 0; i <= P; i++)
	{
		for (int j = 0; j <= P; j++)
		{
			if (j>i)
			{
				Dr[P-i][P-j] = -Dr[i][j];
			}
		}		
	}
	Dr[P][P] = -Dr[0][0];	
}


void mooring_DGSEM::getsurfInt(lexer *p)
{
	p->Darray(sInt, P+1, 2);

	if (P == 1) 
	{
		sInt[0][0] = 2.0;
		sInt[0][1] = -1.0;
		sInt[1][0] = -1.0;
		sInt[1][1] = 2.0;
	}
	else if (P == 2)
	{
		sInt[0][0] = 4.5;
		sInt[0][1] = 1.5;
		sInt[1][0] = -0.75;
		sInt[1][1] = -0.75;
		sInt[2][0] = 1.5;
		sInt[2][1] = 4.5;		
	}
	else if (P == 3)
	{
		sInt[0][0] = 8.0;
		sInt[0][1] = -2.0;
		sInt[1][0] = -0.894427190999917;
		sInt[1][1] = 0.894427190999917;
		sInt[2][0] = 0.894427190999917;
		sInt[2][1] = -0.894427190999917;
		sInt[3][0] = -2.0;
		sInt[3][1] = 8.0;	
	}
	else if (P == 4)
	{
		sInt[0][0] = 12.5;
		sInt[0][1] = 2.5;
		sInt[1][0] = -1.071428571428572;
		sInt[1][1] = -1.071428571428572;
		sInt[2][0] = 0.9375;
		sInt[2][1] = 0.9375;
		sInt[3][0] = -1.071428571428573;
		sInt[3][1] = -1.071428571428573;	
		sInt[4][0] = 2.5;
		sInt[4][1] = 12.5;	
	}
	else if (P == 5)
	{
		sInt[0][0] = 18.0;
		sInt[0][1] = -3.0;
		sInt[1][0] = -1.259090802393857;
		sInt[1][1] = 1.259090802393857;
		sInt[2][0] = 1.039883175166253;
		sInt[2][1] = -1.039883175166253;
		sInt[3][0] = -1.039883175166253;
		sInt[3][1] = 1.039883175166253;	
		sInt[4][0] = 1.259090802393862;
		sInt[4][1] = -1.259090802393862;	
		sInt[5][0] = -3.0;
		sInt[5][1] = 18.0;	
	}
	else if (P == 6)
	{
		sInt[0][0] = 24.5;  
		sInt[0][1] = 3.5;
		sInt[1][0] =  -1.451626611323416;
		sInt[1][1] = -1.451626611323442;
		sInt[2][0] =  1.162370412976334;   
		sInt[2][1] =  1.162370412976333;
		sInt[3][0] = -1.09375;
		sInt[3][1] = -1.09375;
		sInt[4][0] =  1.162370412976333;   
		sInt[4][1] =  1.162370412976333;
		sInt[5][0] = -1.451626611323440;  
		sInt[5][1] = -1.451626611323439;
		sInt[6][0] =  3.5;  
		sInt[6][1] =  24.5;
	}
	else 
	{
		cout<<"Add quadrature points"<<endl;
	}
}


void mooring_DGSEM::getV(lexer *p)
{
	p->Darray(V,P+1,P+1);
	p->Darray(invV,P+1,P+1);

	if (P == 1) 
	{
		V[0][0] = 0.707106781186547;
		V[0][1] = -1.224744871391589;
		V[1][0] = 0.707106781186547;
		V[1][1] = 1.224744871391589;
		
		invV[0][0] = 0.707106781186547;
		invV[0][1] = 0.707106781186547;
		invV[1][0] = -0.408248290463863;
		invV[1][1] = 0.408248290463863;
	}
	else if (P == 2)
	{
		V[0][0] = 0.707106781186547;
		V[0][1] = -1.224744871391589;
		V[0][2] = 1.581138830084190;
		V[1][0] = 0.707106781186547;
		V[1][1] = 0.0;
		V[1][2] = -0.790569415042095;
		V[2][0] = 0.707106781186547;
		V[2][1] = 1.224744871391589;
		V[2][2] = 1.581138830084190;
		
		invV[0][0] = 0.235702260395516;   
		invV[0][1] = 0.942809041582063;   
		invV[0][2] = 0.235702260395516; 
		invV[1][0] = -0.408248290463863;                    
		invV[1][1] = 0.0;   
		invV[1][2] = 0.408248290463863; 
		invV[2][0] = 0.210818510677892;   
		invV[2][1] = -0.421637021355784;    
		invV[2][2] = 0.210818510677892; 
	}
	else if (P == 3)
	{
		V[0][0] = 0.707106781186547;  
		V[0][1] = -1.224744871391589;   
		V[0][2] = 1.581138830084190;  
		V[0][3] = -1.870828693386972;
		V[1][0] = 0.707106781186547;  
		V[1][1] = -0.547722557505166;  
		V[1][2] = -0.316227766016838;  
		V[1][3] = 0.836660026534076;
		V[2][0] = 0.707106781186547;  
		V[2][1] = 0.547722557505166; 
		V[2][2] = -0.316227766016838;  
		V[2][3] = -0.836660026534076;
		V[3][0] = 0.707106781186547;   
		V[3][1] = 1.224744871391589;   
		V[3][2] = 1.581138830084190;   
		V[3][3] = 1.870828693386972;
		
		invV[0][0] = 0.117851130197758;   
		invV[0][1] = 0.589255650988790;    
		invV[0][2] = 0.589255650988789;    
		invV[0][3] = 0.117851130197758; 
		invV[1][0] = -0.204124145231931;   
		invV[1][1] = -0.456435464587639;    
		invV[1][2] = 0.456435464587639;    
		invV[1][3] = 0.204124145231931; 
		invV[2][0] = 0.263523138347365;   
		invV[2][1] = -0.263523138347365;   
		invV[2][2] = -0.263523138347365;    
		invV[2][3] = 0.263523138347365; 
		invV[3][0] = -0.133630620956212;    
		invV[3][1] = 0.298807152333598;   
		invV[3][2] = -0.298807152333598;    
		invV[3][3] = 0.133630620956212; 
	}
	else if (P == 4)
	{
		V[0][0] = 0.707106781186547;  
		V[0][1] = -1.224744871391589;   
		V[0][2] = 1.581138830084190;  
		V[0][3] = -1.870828693386972;   
		V[0][4] = 2.121320343559644;
		V[1][0] = 0.707106781186547;  
		V[1][1] = -0.801783725737273;   
		V[1][2] = 0.225876975726313;   
		V[1][3] = 0.524890659167824;  
		V[1][4] = -0.909137290096990;
		V[2][0] = 0.707106781186547;  
		V[2][1] = 0.0;  
		V[2][2] = -0.790569415042095;   
		V[2][3] = 0.0; 
		V[2][4] = 0.795495128834866;
		V[3][0] = 0.707106781186547;   
		V[3][1] = 0.801783725737273;   
		V[3][2] = 0.225876975726313;  
		V[3][3] = -0.524890659167823;  
		V[3][4] = -0.909137290096990;
		V[4][0] = 0.707106781186547;   
		V[4][1] = 1.224744871391589;   
		V[4][2] = 1.581138830084190;   
		V[4][3] = 1.870828693386972;   
		V[4][4] = 2.121320343559644;
		
		invV[0][0] = 0.070710678118655;    
		invV[0][1] = 0.384980358646009;    
		invV[0][2] = 0.502831488843767;    
		invV[0][3] = 0.384980358646009;    
		invV[0][4] = 0.070710678118655; 
		invV[1][0] = -0.122474487139159;   
		invV[1][1] = -0.436526695123627;    
		invV[1][2] = 0.0;    
		invV[1][3] = 0.436526695123627;    
		invV[1][4] = 0.122474487139159; 
		invV[2][0] = 0.158113883008419;    
		invV[2][1] = 0.122977464562104;   
		invV[2][2] = -0.562182695141045;    
		invV[2][3] = 0.122977464562104;    
		invV[2][4] = 0.158113883008419; 
		invV[3][0] = -0.187082869338697;    
		invV[3][1] = 0.285773803324704;   
		invV[3][2] = 0.0;   
		invV[3][3] = -0.285773803324704;   
		invV[3][4] = 0.187082869338697; 
		invV[4][0] = 0.094280904158206;   
		invV[4][1] = -0.219988776369148;    
		invV[4][2] = 0.251415744421884;   
		invV[4][3] = -0.219988776369148;    
		invV[4][4] = 0.094280904158206; 
	}
	else if (P == 5)
	{
		V[0][0] = 0.707106781186547;  
		V[0][1] = -1.224744871391589;   
		V[0][2] = 1.581138830084190;  
		V[0][3] = -1.870828693386972;   
		V[0][4] = 2.121320343559644;  
		V[0][5] = -2.345207879911717;
		V[1][0] = 0.707106781186547;  
		V[1][1] = -0.936997584313443;   
		V[1][2] = 0.597614304667198;   
		V[1][3] = 0.052565288801477;  
		V[1][4] = -0.681137663582889;   
		V[1][5] = 0.984276557099482;
		V[2][0] = 0.707106781186547;  
		V[2][1] = -0.349335836968916;  
		V[2][2] = -0.597614304667197;   
		V[2][3] = 0.691894769379491;   
		V[2][4] = 0.209733142791857;  
		V[2][5] = -0.812914072195837;
		V[3][0] = 0.707106781186547;   
		V[3][1] = 0.349335836968916;  
		V[3][2] = -0.597614304667197;  
		V[3][3] = -0.691894769379491;   
		V[3][4] = 0.209733142791858;   
		V[3][5] = 0.812914072195837;
		V[4][0] = 0.707106781186547;   
		V[4][1] = 0.936997584313443;   
		V[4][2] = 0.597614304667197;  
		V[4][3] = -0.052565288801479;  
		V[4][4] = -0.681137663582890;  
		V[4][5] = -0.984276557099482;
		V[5][0] = 0.707106781186547;   
		V[5][1] = 1.224744871391589;   
		V[5][2] = 1.581138830084190;   
		V[5][3] = 1.870828693386972;   
		V[5][4] = 2.121320343559644;   
		V[5][5] = 2.345207879911717;
		
		invV[0][0] = 0.047140452079103;    
		invV[0][1] = 0.267622208107490;    
		invV[0][2] = 0.392344120999955;    
		invV[0][3] = 0.392344120999955;    
		invV[0][4] = 0.267622208107490;    
		invV[0][5] = 0.047140452079103; 
		invV[1][0] = -0.081649658092772;   
		invV[1][1] = -0.354630119774219;   
		invV[1][2] = -0.193831915540906;  
		invV[1][3] = 0.193831915540906;    
		invV[1][4] = 0.354630119774219;    
		invV[1][5] = 0.081649658092773; 
		invV[2][0] = 0.105409255338946;    
		invV[2][1] = 0.226182047841886;   
		invV[2][2] = -0.331591303180832;   
		invV[2][3] = -0.331591303180832;    
		invV[2][4] = 0.226182047841886;    
		invV[2][5] = 0.105409255338946; 
		invV[3][0] = -0.124721912892465;    
		invV[3][1] = 0.019894645381923;    
		invV[3][2] = 0.383903608817246;   
		invV[3][3] = -0.383903608817246;   
		invV[3][4] = -0.019894645381923;    
		invV[3][5] = 0.124721912892465; 
		invV[4][0] = 0.141421356237310;   
		invV[4][1] = -0.257793547457352;    
		invV[4][2] = 0.116372191220042;    
		invV[4][3] = 0.116372191220042;   
		invV[4][4] = -0.257793547457352;    
		invV[4][5] = 0.141421356237309; 
		invV[5][0] = -0.071066905451870;    
		invV[5][1] = 0.169329103151465;   
		invV[5][2] = -0.205023719439950;    
		invV[5][3] = 0.205023719439950;   
		invV[5][4] = -0.169329103151464;    
		invV[5][5] = 0.071066905451870; 
	}
	else if (P == 6)
	{
		V[0][0] = 0.707106781186547;  
		V[0][1] = -1.224744871391589;   
		V[0][2] = 1.581138830084190;  
		V[0][3] = -1.870828693386972;   
		V[0][4] = 2.121320343559644;  
		V[0][5] = -2.345207879911717;   
		V[0][6] = 2.549509756796394;
		V[1][0] = 0.707106781186547;  
		V[1][1] = -1.016812459073918;   
		V[1][2] = 0.844182001556940;  
		V[1][3] = -0.346643573228288;  
		V[1][4] = -0.278372647949331;   
		V[1][5] = 0.807538898284274;  
		V[1][6] = -1.057410345369826;
		V[2][0] = 0.707106781186547;  
		V[2][1] = -0.574220155261391;  
		V[2][2] = -0.269222426980870;   
		V[2][3] = 0.833675471702112;  
		V[2][4] = -0.504704283282048;  
		V[2][5] = -0.365166265456773;   
		V[2][6] = 0.846707059684174;
		V[3][0] = 0.707106781186547;   
		V[3][1] = 0.0;  
		V[3][2] = -0.790569415042095;  
		V[3][3] = 0.0;  
		V[3][4] = 0.795495128834866;   
		V[3][5] = 0.0; 
		V[3][6] = -0.796721798998873;
		V[4][0] = 0.707106781186547;   
		V[4][1] = 0.574220155261392;  
		V[4][2] = -0.269222426980870;  
		V[4][3] = -0.833675471702112;  
		V[4][4] = -0.504704283282048;   
		V[4][5] = 0.365166265456773;   
		V[4][6] = 0.846707059684174;
		V[5][0] = 0.707106781186547;   
		V[5][1] = 1.016812459073917;   
		V[5][2] = 0.844182001556939;   
		V[5][3] = 0.346643573228285;  
		V[5][4] = -0.278372647949334;  
		V[5][5] = -0.807538898284277;  
		V[5][6] = -1.057410345369826;
		V[6][0] = 0.707106781186547;   
		V[6][1] = 1.224744871391589;   
		V[6][2] = 1.581138830084190;   
		V[6][3] = 1.870828693386972;   
		V[6][4] = 2.121320343559644;   
		V[6][5] = 2.345207879911717;   
		V[6][6] = 2.549509756796394;
		
		invV[0][0] = 0.033671751485073;  
		invV[0][1] = 0.195745575298432;   
		invV[0][2] = 0.305290086799465;   
		invV[0][3] = 0.344798735207154;   
		invV[0][4] = 0.305290086799465;   
		invV[0][5] = 0.195745575298432;   
		invV[0][6] = 0.033671751485074;
		invV[1][0] = -0.058321184351980;  
		invV[1][1] = -0.281480173953427;  
		invV[1][2] = -0.247916899831716;   
		invV[1][3] = 0.0;   
		invV[1][4] = 0.247916899831716;   
		invV[1][5] = 0.281480173953427;   
		invV[1][6] = 0.058321184351980;
		invV[2][0] = 0.075292325242104;   
		invV[2][1] = 0.233691566744783;  
		invV[2][2] = -0.116235539367100; 
		invV[2][3] = -0.385496705239574;  
		invV[2][4] = -0.116235539367100;   
		invV[2][5] = 0.233691566744783;   
		invV[2][6] = 0.075292325242104;
		invV[3][0] = -0.089087080637475;  
		invV[3][1] = -0.095959970220076;   
		invV[3][2] = 0.359935534335340;  
		invV[3][3] = 0.0;  
		invV[3][4] = -0.359935534335340;   
		invV[3][5] = 0.095959970220076;   
		invV[3][6] = 0.089087080637475;
		invV[4][0] = 0.101015254455221;  
		invV[4][1] = -0.077060799825387;  
		invV[4][2] = -0.217903743183859;   
		invV[4][3] = 0.387898577108049;  
		invV[4][4] = -0.217903743183859;  
		invV[4][5] = -0.077060799825387;   
		invV[4][6] = 0.101015254455221;
		invV[5][0] = -0.111676565710082;   
		invV[5][1] = 0.223547801302750;  
		invV[5][2] = -0.157658848484616;  
		invV[5][3] = 0.0;   
		invV[5][4] = 0.157658848484616;  
		invV[5][5] = -0.223547801302750;   
		invV[5][6] = 0.111676565710082;
		invV[6][0] = 0.056033181468053;  
		invV[6][1] = -0.135100950622134;   
		invV[6][2] = 0.168720859502966;  
		invV[6][3] = -0.179306180697768;   
		invV[6][4] = 0.168720859502966;  
		invV[6][5] = -0.135100950622134;   
		invV[6][6] = 0.056033181468053;
	}
}


void mooring_DGSEM::nodegen(lexer *p)
{
	p->Darray(x,H,P+1);

	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			x[i][j] = vx[i] + 0.5*(r[j] + 1.0)*(vx[i+1] - vx[i]);
		}
	}

	// Scale factor assuming constant delta H
	rx = 0.0;
	for (int j = 0; j < (P+1); j++)
	{	
		rx += Dr[0][j]*x[0][j];
	}
	rx = 1.0/(rx);
}*/
