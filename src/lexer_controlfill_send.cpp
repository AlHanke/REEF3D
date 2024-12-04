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

#include"lexer.h"

void lexer::ctrlsend()
{
    int n;
    int ii,dd;
    
    ii=dd=0;

    for(n=0;n<ctrlsize;++n)
    {
        ictrl[n]=0;
        dctrl[n]=0.0;
    }
    
    ictrl[ii++] = A10;
    ictrl[ii++] = A209;
    ictrl[ii++] = A210;
    ictrl[ii++] = A211;
    ictrl[ii++] = A212;
    ictrl[ii++] = A214;
    ictrl[ii++] = A215;
    ictrl[ii++] = A216;
    ictrl[ii++] = A217;
    ictrl[ii++] = A218;
    ictrl[ii++] = A219;
    ictrl[ii++] = A220;
    ictrl[ii++] = A221;
    dctrl[dd++] = A223;
    ictrl[ii++] = A230;
    ictrl[ii++] = A240;
    ictrl[ii++] = A241;
    ictrl[ii++] = A242;
    ictrl[ii++] = A243;
    ictrl[ii++] = A244;
    dctrl[dd++] = A244_val;
    ictrl[ii++] = A245;
    dctrl[dd++] = A245_val;
    ictrl[ii++] = A246;
    dctrl[dd++] = A247;
    ictrl[ii++] = A248;
    dctrl[dd++] = A249;
    dctrl[dd++] = A250;
    ictrl[ii++] = A251;
    ictrl[ii++] = A260;
    dctrl[dd++] = A261;
    dctrl[dd++] = A262;

    ictrl[ii++] = A310;
    ictrl[ii++] = A311;
    ictrl[ii++] = A312;
    ictrl[ii++] = A313;
    ictrl[ii++] = A320;
    ictrl[ii++] = A321;
    ictrl[ii++] = A322;
    ictrl[ii++] = A323;
    ictrl[ii++] = A329;
    dctrl[dd++] = A340;
    dctrl[dd++] = A341;
    dctrl[dd++] = A342;
    ictrl[ii++] = A343;
    ictrl[ii++] = A344;
    dctrl[dd++] = A344_val;
    ictrl[ii++] = A345;
    dctrl[dd++] = A345_val;
    dctrl[dd++] = A346;
    ictrl[ii++] = A347;
    ictrl[ii++] = A348;
    ictrl[ii++] = A350;
    ictrl[ii++] = A351;
    ictrl[ii++] = A352;
    ictrl[ii++] = A353;
    dctrl[dd++] = A354;
    dctrl[dd++] = A355;
    dctrl[dd++] = A356;
    ictrl[ii++] = A357;
    ictrl[ii++] = A358;
    ictrl[ii++] = A361;
    ictrl[ii++] = A362;
    ictrl[ii++] = A363;
    dctrl[dd++] = A365;
    ictrl[ii++] = A368;

    ictrl[ii++] = A410;
    dctrl[dd++] = A440;
    
    ictrl[ii++] = A501;
    ictrl[ii++] = A509;
    ictrl[ii++] = A510;
    ictrl[ii++] = A511;
    ictrl[ii++] = A512;
    ictrl[ii++] = A514;
    ictrl[ii++] = A515;
    ictrl[ii++] = A516;
    ictrl[ii++] = A517;
    ictrl[ii++] = A518;
    ictrl[ii++] = A519;
    ictrl[ii++] = A520;
    ictrl[ii++] = A521;
    dctrl[dd++] = A522;
    dctrl[dd++] = A523;
    dctrl[dd++] = A531;
    ictrl[ii++] = A540;
    dctrl[dd++] = A541;
    dctrl[dd++] = A542;
    ictrl[ii++] = A543;
    dctrl[dd++] = A544;
    dctrl[dd++] = A545;
    ictrl[ii++] = A550;
    ictrl[ii++] = A551;
    ictrl[ii++] = A552;
    ictrl[ii++] = A553;
    ictrl[ii++] = A560;
    ictrl[ii++] = A570;
    dctrl[dd++] = A571_u;
    dctrl[dd++] = A571_dir;
    ictrl[ii++] = A573;
    ictrl[ii++] = A580;
    dctrl[dd++] = A580_xs;
    dctrl[dd++] = A580_xe;
    dctrl[dd++] = A580_ys;
    dctrl[dd++] = A580_ye;
    ictrl[ii++] = A581;
    ictrl[ii++] = A583;
    ictrl[ii++] = A584;
    ictrl[ii++] = A585;
    ictrl[ii++] = A586;
    ictrl[ii++] = A587;
    ictrl[ii++] = A588;
    ictrl[ii++] = A589;
    ictrl[ii++] = A590;
    ictrl[ii++] = A591;
    dctrl[dd++] = A591_x;
    dctrl[dd++] = A591_y;
    dctrl[dd++] = A591_z;
    ictrl[ii++] = A592;
    dctrl[dd++] = A592_x;
    dctrl[dd++] = A592_y;
    dctrl[dd++] = A592_z;
    ictrl[ii++] = A593;
    dctrl[dd++] = A593_x;
    dctrl[dd++] = A593_y;
    dctrl[dd++] = A593_z;
    dctrl[dd++] = A593_phi;
    dctrl[dd++] = A593_theta;
    dctrl[dd++] = A593_psi;
    ictrl[ii++] = A594;
    
    ictrl[ii++] = B10;
    ictrl[ii++] = B20;
    ictrl[ii++] = B23;
    dctrl[dd++] = B29;
    ictrl[ii++] = B30;
    dctrl[dd++] = B31;
    ictrl[ii++] = B32;
    dctrl[dd++] = B32_x;
    dctrl[dd++] = B32_y;
    dctrl[dd++] = B32_z;
    ictrl[ii++] = B33;
    dctrl[dd++] = B50;
    dctrl[dd++] = B51;
    dctrl[dd++] = B52;
    dctrl[dd++] = B53;
    dctrl[dd++] = B54;
    dctrl[dd++] = B55;
    dctrl[dd++] = B56;
    ictrl[ii++] = B60;
    ictrl[ii++] = B61;
    ictrl[ii++] = B71;
    ictrl[ii++] = B75;
    ictrl[ii++] = B76;
    ictrl[ii++] = B77;
    ictrl[ii++] = B81;
    dctrl[dd++] = B81_1;
    dctrl[dd++] = B81_2;
    dctrl[dd++] = B81_3;
    ictrl[ii++] = B82;
    dctrl[dd++] = B83;
    ictrl[ii++] = B84;
    ictrl[ii++] = B85;
    ictrl[ii++] = B86;
    ictrl[ii++] = B87;
    dctrl[dd++] = B87_1;
    dctrl[dd++] = B87_2;
    dctrl[dd++] = B88;
    ictrl[ii++] = B89;
    ictrl[ii++] = B90;
    ictrl[ii++] = B91;
    dctrl[dd++] = B91_1;
    dctrl[dd++] = B91_2;
    ictrl[ii++] = B92;
    ictrl[ii++] = B93;
    dctrl[dd++] = B93_1;
    dctrl[dd++] = B93_2;
    ictrl[ii++] = B94;
    dctrl[dd++] = B94_wdt;
    dctrl[dd++] = B96_1;
    dctrl[dd++] = B96_2;
    ictrl[ii++] = B98;
    ictrl[ii++] = B99;
    ictrl[ii++] = B101;
    dctrl[dd++] = B102;
    ictrl[ii++] = B105;
    dctrl[dd++] = B105_1;
    dctrl[dd++] = B105_2;
    dctrl[dd++] = B105_3;
    ictrl[ii++] = B106;
    ictrl[ii++] = B107;
    ictrl[ii++] = B110;
    dctrl[dd++] = B110_zs;
    dctrl[dd++] = B110_ze;
    ictrl[ii++] = B108;
    dctrl[dd++] = B111_zs;
    dctrl[dd++] = B111_ze;
    dctrl[dd++] = B112_zs;
    dctrl[dd++] = B112_z2;
    dctrl[dd++] = B112_ze;
    ictrl[ii++] = B115;
    ictrl[ii++] = B116;
    dctrl[dd++] = B117;
    dctrl[dd++] = B120;
    dctrl[dd++] = B122;
    dctrl[dd++] = B123;
    ictrl[ii++] = B125;
    dctrl[dd++] = B125_y;
    ictrl[ii++] = B127;
    ictrl[ii++] = B130;
    dctrl[dd++] = B131;
    dctrl[dd++] = B132_s;
    dctrl[dd++] = B132_e;
    ictrl[ii++] = B133;
    dctrl[dd++] = B134;
    dctrl[dd++] = B135;
    ictrl[ii++] = B136;
    ictrl[ii++] = B138;
    ictrl[ii++] = B138_1;
    ictrl[ii++] = B138_2;
    ictrl[ii++] = B139;
    ictrl[ii++] = B160;
    ictrl[ii++] = B170;
    ictrl[ii++] = B180;
    ictrl[ii++] = B181;
    dctrl[dd++] = B181_1;
    dctrl[dd++] = B181_2;
    dctrl[dd++] = B181_3;
    ictrl[ii++] = B182;
    dctrl[dd++] = B182_1;
    dctrl[dd++] = B182_2;
    dctrl[dd++] = B182_3;
    ictrl[ii++] = B183;
    dctrl[dd++] = B183_1;
    dctrl[dd++] = B183_2;
    dctrl[dd++] = B183_3;
    ictrl[ii++] = B191;
    dctrl[dd++] = B191_1;
    dctrl[dd++] = B191_2;
    dctrl[dd++] = B191_3;
    dctrl[dd++] = B191_4;
    ictrl[ii++] = B192;
    dctrl[dd++] = B192_1;
    dctrl[dd++] = B192_2;
    dctrl[dd++] = B192_3;
    dctrl[dd++] = B192_4;
    dctrl[dd++] = B194_s;
    dctrl[dd++] = B194_e;
    ictrl[ii++] = B240;
    ictrl[ii++] = B241;
    ictrl[ii++] = B242;
    ictrl[ii++] = B243;
    dctrl[dd++] = B260;
    dctrl[dd++] = B264;
    dctrl[dd++] = B267;
    ictrl[ii++] = B269;
    ictrl[ii++] = B270;
    ictrl[ii++] = B274;
    ictrl[ii++] = B281;
    ictrl[ii++] = B282;
    ictrl[ii++] = B291;
    ictrl[ii++] = B295;
    ictrl[ii++] = B307;
    ictrl[ii++] = B308;
    dctrl[dd++] = B309;
    ictrl[ii++] = B310;
    ictrl[ii++] = B321;
    ictrl[ii++] = B322;
    ictrl[ii++] = B411;
    ictrl[ii++] = B412;
    ictrl[ii++] = B413;
    ictrl[ii++] = B414;
    ictrl[ii++] = B415;
    ictrl[ii++] = B416;
    ictrl[ii++] = B417;
    ictrl[ii++] = B418;
    ictrl[ii++] = B421;
    ictrl[ii++] = B422;
    ictrl[ii++] = B440;
    ictrl[ii++] = B441;
    ictrl[ii++] = B442;
    
    
    dctrl[dd++] = C1;
    dctrl[dd++] = C2;
    dctrl[dd++] = C3;
    dctrl[dd++] = C4;
    dctrl[dd++] = C5;
    ictrl[ii++] = C9;
    ictrl[ii++] = C10;
    ictrl[ii++] = C15;
    ictrl[ii++] = C20;
    dctrl[dd++] = C50_1;
    dctrl[dd++] = C50_2;
    dctrl[dd++] = C51;
    dctrl[dd++] = C52;
    dctrl[dd++] = C53;
    dctrl[dd++] = C54;
    dctrl[dd++] = C55;
    dctrl[dd++] = C56;
    dctrl[dd++] = C57_1;
    dctrl[dd++] = C57_2;
    dctrl[dd++] = C57_3;
    dctrl[dd++] = C57_4;
    dctrl[dd++] = C58_1;
    dctrl[dd++] = C58_2;
    dctrl[dd++] = C58_3;
    dctrl[dd++] = C58_4;
    ictrl[ii++] = C75;

    ictrl[ii++] = D10;
    ictrl[ii++] = D11;
    ictrl[ii++] = D20;
    ictrl[ii++] = D21;
    ictrl[ii++] = D30;
    ictrl[ii++] = D31;
    ictrl[ii++] = D33;
    ictrl[ii++] = D37;
    

    ictrl[ii++] = F10;
    ictrl[ii++] = F30;
    ictrl[ii++] = F31;
    ictrl[ii++] = F32;
    dctrl[dd++] = F33;
    ictrl[ii++] = F34;
    ictrl[ii++] = F35;
    ictrl[ii++] = F36;
    dctrl[dd++] = F39;
    ictrl[ii++] = F40;
    dctrl[dd++] = F42;
    dctrl[dd++] = F43;
    ictrl[ii++] = F44;
    dctrl[dd++] = F45;
    ictrl[ii++] = F46;
    ictrl[ii++] = F47;
    ictrl[ii++] = F49;
    ictrl[ii++] = F50;
    ictrl[ii++] = F50_flag;
    dctrl[dd++] = F51;
    dctrl[dd++] = F52;
    dctrl[dd++] = F53;
    dctrl[dd++] = F54;
    dctrl[dd++] = F55;
    dctrl[dd++] = F56;
    dctrl[dd++] = F57_1;
    dctrl[dd++] = F57_2;
    dctrl[dd++] = F57_3;
    dctrl[dd++] = F57_4;
    dctrl[dd++] = F58_1;
    dctrl[dd++] = F58_2;
    dctrl[dd++] = F58_3;
    dctrl[dd++] = F58_4;
    dctrl[dd++] = F59_xm;
    dctrl[dd++] = F59_ym;
    dctrl[dd++] = F59_zs;
    dctrl[dd++] = F59_ze;
    dctrl[dd++] = F59_r;
    dctrl[dd++] = F60;
    dctrl[dd++] = F61;
    dctrl[dd++] = F62;
    dctrl[dd++] = F63;
    ictrl[ii++] = F64;
    dctrl[dd++] = F64_xs;
    dctrl[dd++] = F64_ys;
    dctrl[dd++] = F64_zs;
    dctrl[dd++] = F64_alpha;
    ictrl[ii++] = F70;
    ictrl[ii++] = F71;
    ictrl[ii++] = F72;
    ictrl[ii++] = F80;
    dctrl[dd++] = F84;
    ictrl[ii++] = F85;
    ictrl[ii++] = F150;
    ictrl[ii++] = F151;
    ictrl[ii++] = F300;
    ictrl[ii++] = F305;
    ictrl[ii++] = F310;
    dctrl[dd++] = F321;
    dctrl[dd++] = F322;
    dctrl[dd++] = F323;
    ictrl[ii++] = F350;
    dctrl[dd++] = F360;
    dctrl[dd++] = F361;
    dctrl[dd++] = F362;
    ictrl[ii++] = F369;
    ictrl[ii++] = F370;
    ictrl[ii++] = F371;
    ictrl[ii++] = F374;
    ictrl[ii++] = F375;
    ictrl[ii++] = F378;
    ictrl[ii++] = F379;
    dctrl[dd++] = F380;
    dctrl[dd++] = F381;
    dctrl[dd++] = F382;
    ictrl[ii++] = F390;
    ictrl[ii++] = F391;
    ictrl[ii++] = F394;
    ictrl[ii++] = F395;
    ictrl[ii++] = F398;
    ictrl[ii++] = F399;
    
    
    ictrl[ii++] = G1;
    ictrl[ii++] = G2;
    ictrl[ii++] = G3;
    ictrl[ii++] = G10;
    ictrl[ii++] = G11;
    ictrl[ii++] = G12;
    ictrl[ii++] = G20;
    ictrl[ii++] = G21;
    ictrl[ii++] = G22;
    ictrl[ii++] = G30;
    ictrl[ii++] = G40;
    
    dctrl[dd++] = H1;
    dctrl[dd++] = H2;
    ictrl[ii++] = H3;
    ictrl[ii++] = H4;
    dctrl[dd++] = H4_beta1;
    dctrl[dd++] = H4_beta2;
    ictrl[ii++] = H9;
    ictrl[ii++] = H10;
    ictrl[ii++] = H15;
    dctrl[dd++] = H50_1;
    dctrl[dd++] = H50_2;
    dctrl[dd++] = H51;
    dctrl[dd++] = H52;
    dctrl[dd++] = H53;
    dctrl[dd++] = H54;
    dctrl[dd++] = H55;
    dctrl[dd++] = H56;
    dctrl[dd++] = H57_1;
    dctrl[dd++] = H57_2;
    dctrl[dd++] = H57_3;
    dctrl[dd++] = H57_4;
    dctrl[dd++] = H58_1;
    dctrl[dd++] = H58_2;
    dctrl[dd++] = H58_3;
    dctrl[dd++] = H58_4;
    ictrl[ii++] = H61;
    dctrl[dd++] = H61_T;
    ictrl[ii++] = H62;
    dctrl[dd++] = H62_T;
    ictrl[ii++] = H63;
    dctrl[dd++] = H63_T;
    ictrl[ii++] = H64;
    dctrl[dd++] = H64_T;
    ictrl[ii++] = H65;
    dctrl[dd++] = H65_T;
    ictrl[ii++] = H66;
    dctrl[dd++] = H66_T;

    ictrl[ii++] = I10;
    ictrl[ii++] = I11;
    ictrl[ii++] = I12;
    ictrl[ii++] = I13;
    ictrl[ii++] = I21;
    ictrl[ii++] = I30;
    ictrl[ii++] = I40;
    ictrl[ii++] = I41;
    ictrl[ii++] = I44;
    dctrl[dd++] = I50;
    dctrl[dd++] = I55;
    ictrl[ii++] = I56;
    dctrl[dd++] = I58_1;
    dctrl[dd++] = I58_2;
    ictrl[ii++] = I230;
    dctrl[dd++] = I231;
    dctrl[dd++] = I232;
    dctrl[dd++] = I233;
    ictrl[ii++] = I240;
    dctrl[dd++] = I241;
    
    ictrl[ii++] = M10;

    ictrl[ii++] = N10;
    ictrl[ii++] = N11;
    ictrl[ii++] = N18;
    ictrl[ii++] = N20;
    ictrl[ii++] = N22;
    ictrl[ii++] = N23;
    ictrl[ii++] = N24;
    ictrl[ii++] = N25;
    ictrl[ii++] = N26;
    ictrl[ii++] = N40;
    dctrl[dd++] = N41;
    dctrl[dd++] = N43;
    dctrl[dd++] = N44;
    ictrl[ii++] = N45;
    ictrl[ii++] = N46;
    dctrl[dd++] = N47;
    ictrl[ii++] = N48;
    dctrl[dd++] = N49;
    ictrl[ii++] = N50;
    ictrl[ii++] = N60;
    dctrl[dd++] = N61;
    

    ictrl[ii++] = P10;
    ictrl[ii++] = P11;
    ictrl[ii++] = P12;
    ictrl[ii++] = P15;
    ictrl[ii++] = P16;
    ictrl[ii++] = P20;
    ictrl[ii++] = P21;
    dctrl[dd++] = P22;
    ictrl[ii++] = P23;
    ictrl[ii++] = P24;
    ictrl[ii++] = P25;
    ictrl[ii++] = P26;
    ictrl[ii++] = P27;
    ictrl[ii++] = P28;
    ictrl[ii++] = P29;
    dctrl[dd++] = P30;
    dctrl[dd++] = P34;
    ictrl[ii++] = P35;
    ictrl[ii++] = P40;
    ictrl[ii++] = P41;
    dctrl[dd++] = P42;
    ictrl[ii++] = P43;
    dctrl[dd++] = P43_xs;
    dctrl[dd++] = P43_xe;
    dctrl[dd++] = P43_ys;
    dctrl[dd++] = P43_ye;
    ictrl[ii++] = P44;
    ictrl[ii++] = P45;
    ictrl[ii++] = P46;
    ictrl[ii++] = P46_is;
    ictrl[ii++] = P46_ie;
    ictrl[ii++] = P47;
    dctrl[dd++] = P47_ts;
    dctrl[dd++] = P47_te;
    ictrl[ii++] = P50;
    ictrl[ii++] = P51;
    ictrl[ii++] = P52;
    ictrl[ii++] = P53;
    ictrl[ii++] = P54;
    dctrl[dd++] = P55;
    ictrl[ii++] = P56;
    ictrl[ii++] = P57;
    ictrl[ii++] = P58;
    ictrl[ii++] = P59;
    ictrl[ii++] = P61;
    ictrl[ii++] = P62;
    ictrl[ii++] = P63;
    ictrl[ii++] = P64;
    ictrl[ii++] = P65;
    ictrl[ii++] = P66;
    ictrl[ii++] = P71;
    ictrl[ii++] = P72;
    ictrl[ii++] = P73;
    ictrl[ii++] = P74;
    ictrl[ii++] = P75;
    ictrl[ii++] = P76;
    ictrl[ii++] = P77;
    ictrl[ii++] = P78;
    ictrl[ii++] = P79;
    ictrl[ii++] = P80;
    ictrl[ii++] = P81;
    ictrl[ii++] = P82;
    ictrl[ii++] = P85;
    ictrl[ii++] = P88;
    dctrl[dd++] = P91;
    ictrl[ii++] = P92;
    ictrl[ii++] = P101;
    dctrl[dd++] = P101_xm;
    dctrl[dd++] = P101_ym;
    dctrl[dd++] = P101_zs;
    dctrl[dd++] = P101_ze;
    dctrl[dd++] = P101_r1;
    dctrl[dd++] = P101_r2;
    ictrl[ii++] = P110;
    dctrl[dd++] = P111;
    ictrl[ii++] = P120;
    ictrl[ii++] = P121;
    ictrl[ii++] = P122;
    ictrl[ii++] = P123;
    ictrl[ii++] = P124;
    ictrl[ii++] = P125;
    ictrl[ii++] = P126;
    ictrl[ii++] = P131;
    ictrl[ii++] = P132;
    ictrl[ii++] = P133;
    ictrl[ii++] = P134;
    ictrl[ii++] = P140;
    dctrl[dd++] = P141;
    ictrl[ii++] = P151;
    ictrl[ii++] = P152;
    ictrl[ii++] = P166;
    ictrl[ii++] = P167;
    ictrl[ii++] = P168;
    ictrl[ii++] = P180;
    ictrl[ii++] = P181;
    dctrl[dd++] = P182;
    ictrl[ii++] = P184;
    ictrl[ii++] = P185;
    ictrl[ii++] = P190;
    ictrl[ii++] = P191;
    dctrl[dd++] = P192;
    ictrl[ii++] = P194;
    ictrl[ii++] = P195;
    ictrl[ii++] = P230;
    ictrl[ii++] = P240;
    ictrl[ii++] = P351;
    ictrl[ii++] = P352;


    ictrl[ii++] = Q10;
    ictrl[ii++] = Q11;
    ictrl[ii++] = Q12;
    ictrl[ii++] = Q13;
    dctrl[dd++] = Q14;
    dctrl[dd++] = Q15;
    dctrl[dd++] = Q16;
    dctrl[dd++] = Q17;
    ictrl[ii++] = Q20;
    dctrl[dd++] = Q22;
    dctrl[dd++] = Q23;
    ictrl[ii++] = Q24;
    dctrl[dd++] = Q25;
    ictrl[ii++] = Q29;
    dctrl[dd++] = Q30;
    dctrl[dd++] = Q41;
    ictrl[ii++] = Q43;
    ictrl[ii++] = Q61;
    ictrl[ii++] = Q73;
    ictrl[ii++] = Q101;
    dctrl[dd++] = Q102;
    ictrl[ii++] = Q110;
    ictrl[ii++] = Q111;
    ictrl[ii++] = Q120;
    ictrl[ii++] = Q121;
    ictrl[ii++] = Q122;
    ictrl[ii++] = Q180;
    ictrl[ii++] = Q181;
    dctrl[dd++] = Q182;
    ictrl[ii++] = Q183;
    ictrl[ii++] = Q200;
    ictrl[ii++] = Q201;
    ictrl[ii++] = Q202;
    

    ictrl[ii++] = S10;
    ictrl[ii++] = S11;
    ictrl[ii++] = S12;
    dctrl[dd++] = S13;
    dctrl[dd++] = S14;
    ictrl[ii++] = S15;
    ictrl[ii++] = S16;
    ictrl[ii++] = S17;
    dctrl[dd++] = S19;
    dctrl[dd++] = S20;
    dctrl[dd++] = S21;
    dctrl[dd++] = S22;
    dctrl[dd++] = S23;
    dctrl[dd++] = S24;
    ictrl[ii++] = S25;
    dctrl[dd++] = S26_a;
    dctrl[dd++] = S26_b;
    ictrl[ii++] = S27;
    dctrl[dd++] = S30;
    ictrl[ii++] = S31;
    ictrl[ii++] = S32;
    ictrl[ii++] = S33;
    ictrl[ii++] = S34;
    ictrl[ii++] = S37;
    ictrl[ii++] = S41;
    ictrl[ii++] = S42;
    ictrl[ii++] = S43;
    ictrl[ii++] = S44;
    dctrl[dd++] = S45;
    dctrl[dd++] = S46;
    dctrl[dd++] = S47;
    dctrl[dd++] = S48;
    ictrl[ii++] = S50;
    dctrl[dd++] = S57;
    dctrl[dd++] = S60;
    dctrl[dd++] = S71;
    dctrl[dd++] = S72;
    ictrl[ii++] = S73;
    ictrl[ii++] = S77;
    dctrl[dd++] = S77_xs;
    dctrl[dd++] = S77_xe;
    ictrl[ii++] = S78;
    ictrl[ii++] = S79;
    ictrl[ii++] = S80;
    dctrl[dd++] = S81;
    dctrl[dd++] = S82;
    ictrl[ii++] = S83;
    ictrl[ii++] = S84;
    ictrl[ii++] = S85;
    ictrl[ii++] = S90;
    ictrl[ii++] = S91;
    dctrl[dd++] = S92;
    dctrl[dd++] = S93;
    ictrl[ii++] = S100;
    ictrl[ii++] = S101;

    ictrl[ii++] = T10;
    ictrl[ii++] = T12;
    ictrl[ii++] = T21;
    dctrl[dd++] = T31;
    dctrl[dd++] = T32;
    ictrl[ii++] = T33;
    dctrl[dd++] = T35;
    ictrl[ii++] = T36;
    dctrl[dd++] = T37;
    dctrl[dd++] = T38;
    ictrl[ii++] = T39;
    ictrl[ii++] = T41;
    dctrl[dd++] = T42;
    dctrl[dd++] = T43;
    dctrl[dd++] = T44;
    ictrl[ii++] = T45;
    
    dctrl[dd++] = W1;
    dctrl[dd++] = W2;
    dctrl[dd++] = W3;
    dctrl[dd++] = W4;
    dctrl[dd++] = W5;
    dctrl[dd++] = W6;
    dctrl[dd++] = W7;
    dctrl[dd++] = W10;
    ictrl[ii++] = W11;
    dctrl[dd++] = W11_u;
    dctrl[dd++] = W11_v;
    dctrl[dd++] = W11_w;
    ictrl[ii++] = W12;
    dctrl[dd++] = W12_u;
    dctrl[dd++] = W12_v;
    dctrl[dd++] = W12_w;
    ictrl[ii++] = W13;
    dctrl[dd++] = W13_u;
    dctrl[dd++] = W13_v;
    dctrl[dd++] = W13_w;
    ictrl[ii++] = W14;
    dctrl[dd++] = W14_u;
    dctrl[dd++] = W14_v;
    dctrl[dd++] = W14_w;
    ictrl[ii++] = W15;
    dctrl[dd++] = W15_u;
    dctrl[dd++] = W15_v;
    dctrl[dd++] = W15_w;
    ictrl[ii++] = W16;
    dctrl[dd++] = W16_u;
    dctrl[dd++] = W16_v;
    dctrl[dd++] = W16_w;
    dctrl[dd++] = W20;
    dctrl[dd++] = W21;
    dctrl[dd++] = W22;
    dctrl[dd++] = W29_x;
    dctrl[dd++] = W29_y;
    dctrl[dd++] = W29_z;
    ictrl[ii++] = W30;
    dctrl[dd++] = W31;
    ictrl[ii++] = W41;
    dctrl[dd++] = W50;
    ictrl[ii++] = W50_air;
    ictrl[ii++] = W90;
    dctrl[dd++] = W95;
    dctrl[dd++] = W96;
    dctrl[dd++] = W97;
    dctrl[dd++] = W98;
    ictrl[ii++] = W101;
    dctrl[dd++] = W102_phi;
    dctrl[dd++] = W102_c;
    dctrl[dd++] = W103;
    dctrl[dd++] = W104;
    ictrl[ii++] = W110;
    ictrl[ii++] = W111;
    dctrl[dd++] = W112;
    
    ictrl[ii++] = X10;
    ictrl[ii++] = X11_u;
    ictrl[ii++] = X11_v;
    ictrl[ii++] = X11_w;
    ictrl[ii++] = X11_p;
    ictrl[ii++] = X11_q;
    ictrl[ii++] = X11_r;
    ictrl[ii++] = X12;
    ictrl[ii++] = X14;           
    ictrl[ii++] = X19;
    ictrl[ii++] = X21;
    dctrl[dd++] = X21_d;
    ictrl[ii++] = X22;
    dctrl[dd++] = X22_m;
    ictrl[ii++] = X23;
    dctrl[dd++] = X23_x;
    dctrl[dd++] = X23_y;
    dctrl[dd++] = X23_z;
    ictrl[ii++] = X24;
    dctrl[dd++] = X24_Ix;
    dctrl[dd++] = X24_Iy;
    dctrl[dd++] = X24_Iz;
    dctrl[dd++] = X25_Cp;
    dctrl[dd++] = X25_Cq;
    dctrl[dd++] = X25_Cr;
    dctrl[dd++] = X26_Cu;
    dctrl[dd++] = X26_Cv;
    dctrl[dd++] = X26_Cw;
    ictrl[ii++] = X31;
    ictrl[ii++] = X32;
    ictrl[ii++] = X33;
    ictrl[ii++] = X34;
    ictrl[ii++] = X39;
    ictrl[ii++] = X40;
    dctrl[dd++] = X41;
    dctrl[dd++] = X42;
    dctrl[dd++] = X43;
    dctrl[dd++] = X44;
    ictrl[ii++] = X45;
    ictrl[ii++] = X46;
    ictrl[ii++] = X48;
    ictrl[ii++] = X49;
    ictrl[ii++] = X50;
    ictrl[ii++] = X60;
    ictrl[ii++] = X100;
    dctrl[dd++] = X100_x;
    dctrl[dd++] = X100_y;
    dctrl[dd++] = X100_z;
    ictrl[ii++] = X101;
    dctrl[dd++] = X101_phi;
    dctrl[dd++] = X101_theta;
    dctrl[dd++] = X101_psi;
    ictrl[ii++] = X102;
    dctrl[dd++] = X102_u;
    dctrl[dd++] = X102_v;
    dctrl[dd++] = X102_w;
    ictrl[ii++] = X103;
    dctrl[dd++] = X103_p;
    dctrl[dd++] = X103_q;
    dctrl[dd++] = X103_r;
    ictrl[ii++] = X110;
    ictrl[ii++] = X120;
    dctrl[dd++] = X120_rad;
    dctrl[dd++] = X120_xc;
    dctrl[dd++] = X120_yc;
    dctrl[dd++] = X120_zc;
    ictrl[ii++] = X131;
    dctrl[dd++] = X131_rad;
    dctrl[dd++] = X131_h;
    dctrl[dd++] = X131_xc;
    dctrl[dd++] = X131_yc;
    dctrl[dd++] = X131_zc;
    ictrl[ii++] = X132;
    dctrl[dd++] = X132_rad;
    dctrl[dd++] = X132_h;
    dctrl[dd++] = X132_xc;
    dctrl[dd++] = X132_yc;
    dctrl[dd++] = X132_zc;
    ictrl[ii++] = X133;
    dctrl[dd++] = X133_rad;
    dctrl[dd++] = X133_h;
    dctrl[dd++] = X133_xc;
    dctrl[dd++] = X133_yc;
    dctrl[dd++] = X133_zc;
    ictrl[ii++] = X153;
    dctrl[dd++] = X153_xs;
    dctrl[dd++] = X153_xe;
    dctrl[dd++] = X153_ys;
    dctrl[dd++] = X153_ye;
    dctrl[dd++] = X153_zs;
    dctrl[dd++] = X153_ze;
    ictrl[ii++] = X163;
    ictrl[ii++] = X164;
    ictrl[ii++] = X180;
    ictrl[ii++] = X181;
    dctrl[dd++] = X181_x;
    dctrl[dd++] = X181_y;
    dctrl[dd++] = X181_z;
    ictrl[ii++] = X182;
    dctrl[dd++] = X182_x;
    dctrl[dd++] = X182_y;
    dctrl[dd++] = X182_z;
    ictrl[ii++] = X183;
    dctrl[dd++] = X183_x;
    dctrl[dd++] = X183_y;
    dctrl[dd++] = X183_z;
    dctrl[dd++] = X183_phi;
    dctrl[dd++] = X183_theta;
    dctrl[dd++] = X183_psi;
    ictrl[ii++] = X185;
    dctrl[dd++] = X186;
    ictrl[ii++] = X188;
    ictrl[ii++] = X205;
    ictrl[ii++] = X206;
    dctrl[dd++] = X206_ts;
    dctrl[dd++] = X206_te;
    ictrl[ii++] = X207;
    dctrl[dd++] = X207_ts;
    dctrl[dd++] = X207_te;
    ictrl[ii++] = X210;
    dctrl[dd++] = X210_u;
    dctrl[dd++] = X210_v;
    dctrl[dd++] = X210_w;
    ictrl[ii++] = X211;
    dctrl[dd++] = X211_p;
    dctrl[dd++] = X211_q;
    dctrl[dd++] = X211_r;
    ictrl[ii++] = X240;
    dctrl[dd++] = X241;
    dctrl[dd++] = X242_x;
    dctrl[dd++] = X242_y;
    dctrl[dd++] = X242_z;
    dctrl[dd++] = X243;
    ictrl[ii++] = X310;        
    ictrl[ii++] = X311;
    ictrl[ii++] = X312;
    ictrl[ii++] = X313;
    ictrl[ii++] = X314;
    ictrl[ii++] = X315;
    ictrl[ii++] = X320;
    ictrl[ii++] = X321;
    dctrl[dd++] = X323_m;
    dctrl[dd++] = X323_d;
    dctrl[dd++] = X323_l;
    ictrl[ii++] = X324;
    dctrl[dd++] = X325_dt;
    dctrl[dd++] = X325_relX;
    dctrl[dd++] = X325_relY;
    dctrl[dd++] = X325_relZ;
    ictrl[ii++] = X400;
    dctrl[dd++] = X401_p0;
    dctrl[dd++] = X401_cl;
    dctrl[dd++] = X401_cb;
    dctrl[dd++] = X401_a;
    
    ictrl[ii++] = Y1;
    ictrl[ii++] = Y2;
    ictrl[ii++] = Y3;
    ictrl[ii++] = Y4;
    ictrl[ii++] = Y5;
    ictrl[ii++] = Y40;
    ictrl[ii++] = Y50;
    ictrl[ii++] = Y60;
    ictrl[ii++] = Y71;
    ictrl[ii++] = Y72;
    ictrl[ii++] = Y73;
    ictrl[ii++] = Y74;
    
    ictrl[ii++] = Z10;
    ictrl[ii++] = Z11;
    dctrl[dd++] = Z12_cdx;
    dctrl[dd++] = Z12_cdy;
    dctrl[dd++] = Z12_cdz;
    dctrl[dd++] = Z12_ckx;
    dctrl[dd++] = Z12_cky;
    dctrl[dd++] = Z12_ckz;
    
    
// --------------------------

    for(n=0;n<A581;++n)
    {
    dctrl[dd++] = A581_xs[n];
    dctrl[dd++] = A581_xe[n];
    dctrl[dd++] = A581_ys[n];
    dctrl[dd++] = A581_ye[n];
    dctrl[dd++] = A581_zs[n];
    dctrl[dd++] = A581_ze[n];
    }
    
    for(n=0;n<A583;++n)
    {
    dctrl[dd++] = A583_xc[n];
    dctrl[dd++] = A583_zc[n];
    dctrl[dd++] = A583_ys[n];
    dctrl[dd++] = A583_ye[n];
    dctrl[dd++] = A583_r[n];
    }
    
    for(n=0;n<A584;++n)
    {
    dctrl[dd++] = A584_xc[n];
    dctrl[dd++] = A584_yc[n];
    dctrl[dd++] = A584_zs[n];
    dctrl[dd++] = A584_ze[n];
    dctrl[dd++] = A584_r[n];
    }
    
    for(n=0;n<A585;++n)
    {
    dctrl[dd++] = A585_xm1[n];
    dctrl[dd++] = A585_ym1[n];
    dctrl[dd++] = A585_zm1[n];
    dctrl[dd++] = A585_r1[n];
    dctrl[dd++] = A585_xm2[n];
    dctrl[dd++] = A585_ym2[n];
    dctrl[dd++] = A585_zm2[n];
    dctrl[dd++] = A585_r2[n];
    }
    
    for(n=0;n<A586;++n)
    {
    dctrl[dd++] = A586_xm[n];
    dctrl[dd++] = A586_ym[n];
    dctrl[dd++] = A586_zm[n];
    dctrl[dd++] = A586_r[n];
    }
    
    for(n=0;n<A587;++n)
    {
    dctrl[dd++] = A587_xs[n];
    dctrl[dd++] = A587_xe[n];
    dctrl[dd++] = A587_ys[n];
    dctrl[dd++] = A587_ye[n];
    dctrl[dd++] = A587_zs[n];
    dctrl[dd++] = A587_ze[n];
    }
    
    for(n=0;n<A588;++n)
    {
    dctrl[dd++] = A588_xs[n];
    dctrl[dd++] = A588_xe[n];
    dctrl[dd++] = A588_ys[n];
    dctrl[dd++] = A588_ye[n];
    dctrl[dd++] = A588_zs[n];
    dctrl[dd++] = A588_ze[n];
    }
    
    for(n=0;n<A589;++n)
    {
    dctrl[dd++] = A589_xs[n];
    dctrl[dd++] = A589_xe[n];
    dctrl[dd++] = A589_ys[n];
    dctrl[dd++] = A589_ye[n];
    dctrl[dd++] = A589_zs[n];
    dctrl[dd++] = A589_ze[n];
    }
    
    for(n=0;n<B71;++n)
    {
    dctrl[dd++] = B71_val[n];
    dctrl[dd++] = B71_dist[n];
    dctrl[dd++] = B71_b[n];
    dctrl[dd++] = B71_x[n];
    dctrl[dd++] = B71_y[n];
    }
    
    for(n=0;n<B106;++n)
    {
    dctrl[dd++] = B106_b[n];
    dctrl[dd++] = B106_x[n];
    dctrl[dd++] = B106_y[n];
    }
    
    for(n=0;n<B107;++n)
    {
    dctrl[dd++] = B107_xs[n];
    dctrl[dd++] = B107_xe[n];
    dctrl[dd++] = B107_ys[n];
    dctrl[dd++] = B107_ye[n];
    dctrl[dd++] = B107_d[n];
    }
    
    for(n=0;n<B108;++n)
    {
    dctrl[dd++] = B108_xs[n];
    dctrl[dd++] = B108_xe[n];
    dctrl[dd++] = B108_ys[n];
    dctrl[dd++] = B108_ye[n];
    dctrl[dd++] = B108_d[n];
    }
    
    
    for(n=0;n<B240;++n)
    {
    dctrl[dd++] = B240_C[n];
    dctrl[dd++] = B240_D[n];
    dctrl[dd++] = B240_xs[n];
    dctrl[dd++] = B240_xe[n];
    dctrl[dd++] = B240_ys[n];
    dctrl[dd++] = B240_ye[n];
    dctrl[dd++] = B240_zs[n];
    dctrl[dd++] = B240_ze[n];
    }
    
    for(n=0;n<B270;++n)
    {
    dctrl[dd++] = B270_xs[n];
    dctrl[dd++] = B270_xe[n];
    dctrl[dd++] = B270_ys[n];
    dctrl[dd++] = B270_ye[n];
    dctrl[dd++] = B270_zs[n];
    dctrl[dd++] = B270_ze[n];
    dctrl[dd++] = B270_n[n];
    dctrl[dd++] = B270_d50[n];
    dctrl[dd++] = B270_alpha[n];
    dctrl[dd++] = B270_beta[n];
    }
    
    for(n=0;n<B274;++n)
    {
    dctrl[dd++] = B274_xc[n];
    dctrl[dd++] = B274_yc[n];
    dctrl[dd++] = B274_zs[n];
    dctrl[dd++] = B274_ze[n];
    dctrl[dd++] = B274_r[n];
    dctrl[dd++] = B274_n[n];
    dctrl[dd++] = B274_d50[n];
    dctrl[dd++] = B274_alpha[n];
    dctrl[dd++] = B274_beta[n];
    }
    
    for(n=0;n<B281;++n)
    {
    dctrl[dd++] = B281_xs[n];
    dctrl[dd++] = B281_xe[n];
    dctrl[dd++] = B281_ys[n];
    dctrl[dd++] = B281_ye[n];
    dctrl[dd++] = B281_zs[n];
    dctrl[dd++] = B281_ze[n];
    dctrl[dd++] = B281_n[n];
    dctrl[dd++] = B281_d50[n];
    dctrl[dd++] = B281_alpha[n];
    dctrl[dd++] = B281_beta[n];
    }
    
    for(n=0;n<B282;++n)
    {
    dctrl[dd++] = B282_xs[n];
    dctrl[dd++] = B282_xe[n];
    dctrl[dd++] = B282_ys[n];
    dctrl[dd++] = B282_ye[n];
    dctrl[dd++] = B282_zs[n];
    dctrl[dd++] = B282_ze[n];
    dctrl[dd++] = B282_n[n];
    dctrl[dd++] = B282_d50[n];
    dctrl[dd++] = B282_alpha[n];
    dctrl[dd++] = B282_beta[n];
    }
    
    for(n=0;n<B291;++n)
    {
    dctrl[dd++] = B291_xs[n];
    dctrl[dd++] = B291_xe[n];
    dctrl[dd++] = B291_ys[n];
    dctrl[dd++] = B291_ye[n];
    dctrl[dd++] = B291_zs[n];
    dctrl[dd++] = B291_ze[n];
    dctrl[dd++] = B291_d[n];
    dctrl[dd++] = B291_n[n];
    dctrl[dd++] = B291_d50[n];
    dctrl[dd++] = B291_alpha[n];
    dctrl[dd++] = B291_beta[n];
    }
    
    for(n=0;n<B310;++n)
    {
    dctrl[dd++] = B310_xs[n];
    dctrl[dd++] = B310_xe[n];
    dctrl[dd++] = B310_ys[n];
    dctrl[dd++] = B310_ye[n];
    dctrl[dd++] = B310_zs[n];
    dctrl[dd++] = B310_ze[n];
    dctrl[dd++] = B310_N[n];
    dctrl[dd++] = B310_D[n];
    dctrl[dd++] = B310_Cd[n];
    }
    
    for(n=0;n<B321;++n)
    {
    dctrl[dd++] = B321_xs[n];
    dctrl[dd++] = B321_xe[n];
    dctrl[dd++] = B321_ys[n];
    dctrl[dd++] = B321_ye[n];
    dctrl[dd++] = B321_zs[n];
    dctrl[dd++] = B321_ze[n];
    dctrl[dd++] = B321_N[n];
    dctrl[dd++] = B321_D[n];
    dctrl[dd++] = B321_Cd[n];
    }
    
    for(n=0;n<B322;++n)
    {
    dctrl[dd++] = B322_xs[n];
    dctrl[dd++] = B322_xe[n];
    dctrl[dd++] = B322_ys[n];
    dctrl[dd++] = B322_ye[n];
    dctrl[dd++] = B322_zs[n];
    dctrl[dd++] = B322_ze[n];
    dctrl[dd++] = B322_N[n];
    dctrl[dd++] = B322_D[n];
    dctrl[dd++] = B322_Cd[n];
    }
    
    for(n=0;n<B411;++n)
    {
    ictrl[ii++] = B411_ID[n];
    dctrl[dd++] = B411_Q[n];
    }
    
    for(n=0;n<B412;++n)
    {
    ictrl[ii++] = B412_ID[n];
    dctrl[dd++] = B412_pressBC[n];
    }
    
    for(n=0;n<B413;++n)
    {
    ictrl[ii++] = B413_ID[n];
    dctrl[dd++] = B413_h[n];
    }
    
    for(n=0;n<B414;++n)
    {
    ictrl[ii++] = B414_ID[n];
    dctrl[dd++] = B414_Uio[n];
    }
    
    for(n=0;n<B415;++n)
    {
    ictrl[ii++] = B415_ID[n];
    dctrl[dd++] = B415_U[n];
    dctrl[dd++] = B415_V[n];
    dctrl[dd++] = B415_W[n];
    }
    
    for(n=0;n<B416;++n)
    {
    ictrl[ii++] = B416_ID[n];
    dctrl[dd++] = B416_alpha[n];
    }
    
    for(n=0;n<B417;++n)
    {
    ictrl[ii++] = B417_ID[n];
    dctrl[dd++] = B417_Nx[n];
    dctrl[dd++] = B417_Ny[n];
    dctrl[dd++] = B417_Nz[n];
    }
    
    for(n=0;n<B418;++n)
    {
    ictrl[ii++] = B418_ID[n];
    ictrl[ii++] = B418_pio[n];
    }
    
    for(n=0;n<B421;++n)
    {
    ictrl[ii++] = B421_ID[n];
    ictrl[ii++] = B421_Q[n];
    }
    
    for(n=0;n<B422;++n)
    {
    ictrl[ii++] = B422_ID[n];
    ictrl[ii++] = B422_FSF[n];
    }
    
    for(n=0;n<B440;++n)
    {
    ictrl[ii++] = B440_ID[n];
    ictrl[ii++] = B440_face[n];
    dctrl[dd++] = B440_xs[n];
    dctrl[dd++] = B440_xe[n];
    dctrl[dd++] = B440_ys[n];
    dctrl[dd++] = B440_ye[n];
    }
    
    for(n=0;n<B441;++n)
    {
    ictrl[ii++] = B441_ID[n];
    ictrl[ii++] = B441_face[n];
    dctrl[dd++] = B441_xs[n];
    dctrl[dd++] = B441_xe[n];
    dctrl[dd++] = B441_ys[n];
    dctrl[dd++] = B441_ye[n];
    dctrl[dd++] = B441_zs[n];
    dctrl[dd++] = B441_ze[n];
    }
    
    for(n=0;n<B442;++n)
    {
    ictrl[ii++] = B442_ID[n];
    ictrl[ii++] = B442_face[n];
    dctrl[dd++] = B442_xm[n];
    dctrl[dd++] = B442_ym[n];
    dctrl[dd++] = B442_zm[n];
    dctrl[dd++] = B442_r[n];
    }
    
    for(n=0;n<C75;++n)
    {
    dctrl[dd++] = C75_x[n];
    dctrl[dd++] = C75_z[n];
    dctrl[dd++] = C75_a[n];
    dctrl[dd++] = C75_s[n];
    dctrl[dd++] = C75_l[n];
    dctrl[dd++] = C75_v[n];
    }
    
    for(n=0;n<F70;++n)
    {
    dctrl[dd++] = F70_xs[n];
    dctrl[dd++] = F70_xe[n];
    dctrl[dd++] = F70_ys[n];
    dctrl[dd++] = F70_ye[n];
    dctrl[dd++] = F70_zs[n];
    dctrl[dd++] = F70_ze[n];
    }
    
    for(n=0;n<F71;++n)
    {
    dctrl[dd++] = F71_xs[n];
    dctrl[dd++] = F71_xe[n];
    dctrl[dd++] = F71_ys[n];
    dctrl[dd++] = F71_ye[n];
    dctrl[dd++] = F71_zs[n];
    dctrl[dd++] = F71_ze[n];
    }
    
    for(n=0;n<F72;++n)
    {
    dctrl[dd++] = F72_xs[n];
    dctrl[dd++] = F72_xe[n];
    dctrl[dd++] = F72_ys[n];
    dctrl[dd++] = F72_ye[n];
    dctrl[dd++] = F72_h[n];
    }

    for(n=0;n<F369;++n)
    {
    dctrl[dd++] = F369_x[n];
    dctrl[dd++] = F369_z[n];
    dctrl[dd++] = F369_a[n];
    dctrl[dd++] = F369_s[n];
    dctrl[dd++] = F369_l[n];
    dctrl[dd++] = F369_v[n];
    }
    
    for(n=0;n<F370;++n)
    {
    dctrl[dd++] = F370_xs[n];
    dctrl[dd++] = F370_xe[n];
    dctrl[dd++] = F370_ys[n];
    dctrl[dd++] = F370_ye[n];
    dctrl[dd++] = F370_zs[n];
    dctrl[dd++] = F370_ze[n];
    }
    
    for(n=0;n<F371;++n)
    {
    dctrl[dd++] = F371_xs[n];
    dctrl[dd++] = F371_xe[n];
    dctrl[dd++] = F371_ys[n];
    dctrl[dd++] = F371_ye[n];
    dctrl[dd++] = F371_zs[n];
    dctrl[dd++] = F371_ze[n];
    }
    
    for(n=0;n<F374;++n)
    {
    dctrl[dd++] = F374_xc[n];
    dctrl[dd++] = F374_zc[n];
    dctrl[dd++] = F374_r[n];
    }
    
    for(n=0;n<F375;++n)
    {
    dctrl[dd++] = F375_xc[n];
    dctrl[dd++] = F375_zc[n];
    dctrl[dd++] = F375_r[n];
    }
    
    for(n=0;n<F378;++n)
    {
    dctrl[dd++] = F378_xc[n];
    dctrl[dd++] = F378_yc[n];
    dctrl[dd++] = F378_zc[n];
    dctrl[dd++] = F378_r[n];
    }
    
    for(n=0;n<F379;++n)
    {
    dctrl[dd++] = F379_xc[n];
    dctrl[dd++] = F379_yc[n];
    dctrl[dd++] = F379_zc[n];
    dctrl[dd++] = F379_r[n];
    }
    
    for(n=0;n<F390;++n)
    {
    dctrl[dd++] = F390_xs[n];
    dctrl[dd++] = F390_xe[n];
    dctrl[dd++] = F390_ys[n];
    dctrl[dd++] = F390_ye[n];
    dctrl[dd++] = F390_zs[n];
    dctrl[dd++] = F390_ze[n];
    }
    
    for(n=0;n<F391;++n)
    {
    dctrl[dd++] = F391_xs[n];
    dctrl[dd++] = F391_xe[n];
    dctrl[dd++] = F391_ys[n];
    dctrl[dd++] = F391_ye[n];
    dctrl[dd++] = F391_zs[n];
    dctrl[dd++] = F391_ze[n];
    }
    
    for(n=0;n<F394;++n)
    {
    dctrl[dd++] = F394_xc[n];
    dctrl[dd++] = F394_zc[n];
    dctrl[dd++] = F394_r[n];
    }
    
    for(n=0;n<F395;++n)
    {
    dctrl[dd++]   = F395_xc[n];
    dctrl[dd++] = F395_zc[n];
    dctrl[dd++] = F395_r[n];
    }
    
    for(n=0;n<F398;++n)
    {
    dctrl[dd++] = F398_xc[n];
    dctrl[dd++] = F398_yc[n];
    dctrl[dd++] = F398_zc[n];
    dctrl[dd++] = F398_r[n];
    }
    
    for(n=0;n<F399;++n)
    {
    dctrl[dd++] = F399_xc[n];
    dctrl[dd++] = F399_yc[n];
    dctrl[dd++] = F399_zc[n];
    dctrl[dd++] = F399_r[n];
    }

    for(n=0;n<P35;++n)
    {
    dctrl[dd++]  = P35_ts[n];
    dctrl[dd++]  = P35_te[n];
    dctrl[dd++] = P35_dt[n];
    }

    for(n=0;n<P50;++n)
    {
    dctrl[dd++] = P50_x[n];
    dctrl[dd++] = P50_y[n];
    }
    
    for(n=0;n<P51;++n)
    {
    dctrl[dd++] = P51_x[n];
    dctrl[dd++] = P51_y[n];
    }

    for(n=0;n<P52;++n)
    {
    dctrl[dd++] = P52_y[n];
    }
    
    for(n=0;n<P56;++n)
    {
    dctrl[dd++] = P56_x[n];
    }
    
    for(n=0;n<P58;++n)
    {
    dctrl[dd++] = P58_x[n];
    dctrl[dd++] = P58_y[n];
    dctrl[dd++] = P58_T[n];
    }
    
    for(n=0;n<P61;++n)
    {
    dctrl[dd++] = P61_x[n];
    dctrl[dd++] = P61_y[n];
    dctrl[dd++] = P61_z[n];
    }
    
    for(n=0;n<P62;++n)
    {
    dctrl[dd++] = P62_xs[n];
    dctrl[dd++] = P62_xe[n];
    dctrl[dd++] = P62_ys[n];
    dctrl[dd++] = P62_ye[n];
    dctrl[dd++] = P62_zs[n];
    dctrl[dd++] = P62_ze[n];
    }
    
    for(n=0;n<P63;++n)
    {
    dctrl[dd++] = P63_x[n];
    dctrl[dd++] = P63_y[n];
    }
    
    for(n=0;n<P64;++n)
    {
    dctrl[dd++] = P64_x[n];
    dctrl[dd++] = P64_y[n];
    dctrl[dd++] = P64_z[n];
    }
    
    for(n=0;n<P65;++n)
    {
    dctrl[dd++] = P65_x[n];
    dctrl[dd++] = P65_y[n];
    dctrl[dd++] = P65_z[n];
    }
    
    for(n=0;n<P66;++n)
    {
    dctrl[dd++] = P66_x[n];
    dctrl[dd++] = P66_y[n];
    dctrl[dd++] = P66_z[n];
    }
    
    for(n=0;n<P167;++n)
    {
    dctrl[dd++] = P167_x[n];
    }
    
    for(n=0;n<P168;++n)
    {
    dctrl[dd++] = P168_x[n];
    dctrl[dd++] = P168_zs[n];
    dctrl[dd++] = P168_ze[n];
    }
    
    for(n=0;n<P81;++n)
    {
    dctrl[dd++] = P81_xs[n];
    dctrl[dd++] = P81_xe[n];
    dctrl[dd++] = P81_ys[n];
    dctrl[dd++] = P81_ye[n];
    dctrl[dd++] = P81_zs[n];
    dctrl[dd++] = P81_ze[n];
    }
    
    for(n=0;n<P85;++n)
    {
    dctrl[dd++] = P85_x[n];
    dctrl[dd++] = P85_y[n];
    dctrl[dd++] = P85_r[n];
    dctrl[dd++] = P85_cd[n];
    dctrl[dd++] = P85_cm[n];
    }
    
    for(n=0;n<P88;++n)
    {
    dctrl[dd++] = P88_x[n];
    dctrl[dd++] = P88_y[n];
    }
    
    for(n=0;n<P121;++n)
    {
    dctrl[dd++] = P121_x[n];
    dctrl[dd++] = P121_y[n];
    }
    
    for(n=0;n<P123;++n)
    {
    dctrl[dd++] = P123_y[n];
    }
    
    for(n=0;n<P124;++n)
    {
    dctrl[dd++] = P124_x[n];
    }
    
    for(n=0;n<P125;++n)
    {
    dctrl[dd++] = P125_x[n];
    dctrl[dd++] = P125_y[n];
    }
    
    for(n=0;n<P133;++n)
    {
    dctrl[dd++] = P133_y[n];
    }
    
    for(n=0;n<P134;++n)
    {
    dctrl[dd++] = P134_y[n];
    }
    
    for(n=0;n<P140;++n)
    {
    dctrl[dd++] = P140_x[n];
    dctrl[dd++] = P140_y[n];
    }
    
    for(n=0;n<P184;++n)
    {
    ictrl[ii++] = P184_its[n];
    ictrl[ii++] = P184_ite[n];
    ictrl[ii++] = P184_dit[n];
    }
    
    for(n=0;n<P185;++n)
    {
    dctrl[dd++]  = P185_ts[n];
    dctrl[dd++]  = P185_te[n];
    dctrl[dd++] = P185_dt[n];
    }

    for(n=0;n<P194;++n)
    {
    ictrl[ii++]  = P194_its[n];
    ictrl[ii++]  = P194_ite[n];
    ictrl[ii++] = P194_dit[n];
    }
    
    for(n=0;n<P195;++n)
    {
    dctrl[dd++]  = P195_ts[n];
    dctrl[dd++]  = P195_te[n];
    dctrl[dd++] = P195_dt[n];
    }
    
    for(n=0;n<P230;++n)
    {
    dctrl[dd++] = P230_x[n];
    }
    
    for(n=0;n<P240;++n)
    {
    dctrl[dd++] = P240_x[n];
    }
    
    for(n=0;n<P351;++n)
    {
    dctrl[dd++] = P351_x[n];
    dctrl[dd++] = P351_y[n];
    }
    
    for(n=0;n<P352;++n)
    {
    dctrl[dd++] = P352_x[n];
    dctrl[dd++] = P352_y[n];
    }

    for(n=0;n<Q61;++n)
    {
    dctrl[dd++] = Q61_x[n];
    dctrl[dd++] = Q61_y[n];
    dctrl[dd++] = Q61_z[n];
    ictrl[ii++] = Q61_i[n];
    }

    for(n=0;n<Q73;++n)
    {
    dctrl[dd++] = Q73_val[n];
    dctrl[dd++] = Q73_dist[n];
    dctrl[dd++] = Q73_b[n];
    dctrl[dd++] = Q73_x[n];
    dctrl[dd++] = Q73_y[n];
    }
    
    for(n=0;n<Q110;++n)
    {
    dctrl[dd++] = Q110_xs[n];
    dctrl[dd++] = Q110_xe[n];
    dctrl[dd++] = Q110_ys[n];
    dctrl[dd++] = Q110_ye[n];
    dctrl[dd++] = Q110_zs[n];
    dctrl[dd++] = Q110_ze[n];
    }
    for(n=0;n<Q111;++n)
    {
    dctrl[dd++] = Q111_xs[n];
    dctrl[dd++] = Q111_xe[n];
    dctrl[dd++] = Q111_ys[n];
    dctrl[dd++] = Q111_ye[n];
    dctrl[dd++] = Q111_zs[n];
    dctrl[dd++] = Q111_ze[n];
    }
    
    for(n=0;n<S73;++n)
    {
    dctrl[dd++] = S73_val[n];
    dctrl[dd++] = S73_dist[n];
    dctrl[dd++] = S73_b[n];
    dctrl[dd++] = S73_x[n];
    dctrl[dd++] = S73_y[n];
    }
    
    for(n=0;n<W41;++n)
    {
    dctrl[dd++] = W41_xc[n];
    dctrl[dd++] = W41_yc[n];
    dctrl[dd++] = W41_zs[n];
    dctrl[dd++] = W41_ze[n];
    dctrl[dd++] = W41_vel[n];
    dctrl[dd++] = W41_beta[n];
    }
    
    for(n=0;n<X110;++n)
    {
    dctrl[dd++] = X110_xs[n];
    dctrl[dd++] = X110_xe[n];
    dctrl[dd++] = X110_ys[n];
    dctrl[dd++] = X110_ye[n];
    dctrl[dd++] = X110_zs[n];
    dctrl[dd++] = X110_ze[n];
    }
    
    for(n=0;n<X163;++n)
    {
    dctrl[dd++] = X163_x1[n];
    dctrl[dd++] = X163_y1[n];
    dctrl[dd++] = X163_z1[n];
    dctrl[dd++] = X163_x2[n];
    dctrl[dd++] = X163_y2[n];
    dctrl[dd++] = X163_z2[n];
    dctrl[dd++] = X163_x3[n];
    dctrl[dd++] = X163_y3[n];
    dctrl[dd++] = X163_z3[n];
    dctrl[dd++] = X163_x4[n];
    dctrl[dd++] = X163_y4[n];
    dctrl[dd++] = X163_z4[n];
    dctrl[dd++] = X163_x5[n];
    dctrl[dd++] = X163_y5[n];
    dctrl[dd++] = X163_z5[n];
    dctrl[dd++] = X163_x6[n];
    dctrl[dd++] = X163_y6[n];
    dctrl[dd++] = X163_z6[n];
    }
    
    for(n=0;n<X164;++n)
    {
    dctrl[dd++] = X164_x1[n];
    dctrl[dd++] = X164_y1[n];
    dctrl[dd++] = X164_z1[n];
    dctrl[dd++] = X164_x2[n];
    dctrl[dd++] = X164_y2[n];
    dctrl[dd++] = X164_z2[n];
    dctrl[dd++] = X164_x3[n];
    dctrl[dd++] = X164_y3[n];
    dctrl[dd++] = X164_z3[n];
    dctrl[dd++] = X164_x4[n];
    dctrl[dd++] = X164_y4[n];
    dctrl[dd++] = X164_z4[n];
    dctrl[dd++] = X164_x5[n];
    dctrl[dd++] = X164_y5[n];
    dctrl[dd++] = X164_z5[n];
    dctrl[dd++] = X164_x6[n];
    dctrl[dd++] = X164_y6[n];
    dctrl[dd++] = X164_z6[n];
    dctrl[dd++] = X164_x7[n];
    dctrl[dd++] = X164_y7[n];
    dctrl[dd++] = X164_z7[n];
    dctrl[dd++] = X164_x8[n];
    dctrl[dd++] = X164_y8[n];
    dctrl[dd++] = X164_z8[n];
    }
   
    for(n=0;n<X311;++n)
    {
    dctrl[dd++] = X311_xs[n];
    dctrl[dd++] = X311_xe[n];
    dctrl[dd++] = X311_ys[n];
    dctrl[dd++] = X311_ye[n];
    dctrl[dd++] = X311_zs[n];
    dctrl[dd++] = X311_ze[n];
    dctrl[dd++] = X311_w[n];
    dctrl[dd++] = X311_rho_c[n];
    dctrl[dd++] = X311_EA[n];
    dctrl[dd++] = X311_d[n];
    dctrl[dd++] = X311_l[n];
    dctrl[dd++] = X311_H[n]; 
    dctrl[dd++] = X311_P[n];  
    dctrl[dd++] = X311_facT[n];   
    }

    for(n=0;n<X312;++n)
    {
    dctrl[dd++] = X311_xs[n];
    dctrl[dd++] = X311_xe[n];
    dctrl[dd++] = X311_ys[n];
    dctrl[dd++] = X311_ye[n];
    dctrl[dd++] = X311_zs[n];
    dctrl[dd++] = X311_ze[n];
    dctrl[dd++] = X312_k[n];
    dctrl[dd++] = X312_T0[n]; 
    }

    if (X314 > 0)
    {
        for(n=0;n<X311;++n)
        {
            dctrl[dd++] = X314_T[n]; 
        }
        for(n=0;n<X312;++n)
        {
            dctrl[dd++] = X314_T[n]; 
        }
    }
    if (X315 > 0)
    {
        for(n=0;n<X311;++n)
        {
            dctrl[dd++] = X315_t[n];
        }
        for(n=0;n<X312;++n)
        {
            dctrl[dd++] = X315_t[n];
        }
    }
    
    for(n=0;n<X320;++n)
    {
    ictrl[ii++] = X320_type[n];
    }

    for(n=0;n<X321;++n)
    {
    dctrl[dd++] = X321_Sn[n];
    dctrl[dd++] = X321_d[n];
    dctrl[dd++] = X321_lambda[n];
    dctrl[dd++] = X321_dk[n];
    dctrl[dd++] = X321_rho[n];
    dctrl[dd++] = X321_nd[n];
    dctrl[dd++] = X321_nl[n];
    
    dctrl[dd++] = X322_D[n];
    dctrl[dd++] = X322_L[n];
    dctrl[dd++] = X322_x0[n];
    dctrl[dd++] = X322_y0[n];
    dctrl[dd++] = X322_z0[n];
    dctrl[dd++] = X322_phi[n];
    dctrl[dd++] = X322_theta[n];
    dctrl[dd++] = X322_psi[n];
    }      
    
    for(n=0;n<X324;++n)
    {
    dctrl[dd++] = X324_x[n];
    dctrl[dd++] = X324_y[n];
    dctrl[dd++] = X324_z[n];
    }
    
    for(n=0;n<Z11;++n)
    {
    dctrl[dd++] = Z11_x[n];
    dctrl[dd++] = Z11_y[n];
    dctrl[dd++] = Z11_z[n];
    dctrl[dd++] = Z11_l[n];
    dctrl[dd++] = Z11_w[n];
    dctrl[dd++] = Z11_t[n];
    dctrl[dd++] = Z11_rho[n];
    dctrl[dd++] = Z11_e[n];
    dctrl[dd++] = Z11_ix[n];
    dctrl[dd++] = Z11_iy[n];
    dctrl[dd++] = Z11_iz[n];
    dctrl[dd++] = Z11_nu[n];
    dctrl[dd++] = Z11_n[n];
    }
}


