
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

void lexer::ctrlrecv()
{
    int n;
    
    int ii,dd;
    
    ii=dd=0;




    
    A10 = ictrl[ii++];
    A209 = ictrl[ii++];
    A210 = ictrl[ii++];
    A211 = ictrl[ii++];
    A212 = ictrl[ii++];
    A214 = ictrl[ii++];
    A215 = ictrl[ii++];
    A216 = ictrl[ii++];
    A217 = ictrl[ii++];
    A218 = ictrl[ii++];
    A219 = ictrl[ii++];
    A220 = ictrl[ii++];
    A221 = ictrl[ii++];
    A223 = dctrl[dd++];
    A230 = ictrl[ii++];
    A240 = ictrl[ii++];
    A241 = ictrl[ii++];
    A242 = ictrl[ii++];
    A243 = ictrl[ii++];
    A244 = ictrl[ii++];
    A244_val = dctrl[dd++];
    A245 = ictrl[ii++];
    A245_val = dctrl[dd++];
    A246 = ictrl[ii++];
    A247 = dctrl[dd++];
    A248 = ictrl[ii++];
    A249 = dctrl[dd++];
    A250 = dctrl[dd++];
    A251 = ictrl[ii++];
    A260 = ictrl[ii++];
    A261 = dctrl[dd++];
    A262 = dctrl[dd++];
    
    A310 = ictrl[ii++];
    A311 = ictrl[ii++];
    A312 = ictrl[ii++];
    A313 = ictrl[ii++];
    A320 = ictrl[ii++];
    A321 = ictrl[ii++];
    A322 = ictrl[ii++];
    A323 = ictrl[ii++];
    A329 = ictrl[ii++];
    A340 = dctrl[dd++];
    A341 = dctrl[dd++];
    A342 = dctrl[dd++];
    A343 = ictrl[ii++];
    A344 = ictrl[ii++];
    A344_val = dctrl[dd++];
    A345 = ictrl[ii++];
    A345_val = dctrl[dd++];
    A346 = dctrl[dd++];
    A347 = ictrl[ii++];
    A348 = ictrl[ii++];
    A350 = ictrl[ii++];
    A351 = ictrl[ii++];
    A352 = ictrl[ii++];
    A353 = ictrl[ii++];
    A354 = dctrl[dd++];
    A355 = dctrl[dd++];
    A356 = dctrl[dd++];
    A357 = ictrl[ii++];
    A358 = ictrl[ii++];
    A361 = ictrl[ii++];
    A362 = ictrl[ii++];
    A363 = ictrl[ii++];
    A365 = dctrl[dd++];
    A368 = ictrl[ii++];
    
    A410 = ictrl[ii++];
    A440 = dctrl[dd++];
    
    A501 = ictrl[ii++];
    A509 = ictrl[ii++];
    A510 = ictrl[ii++];
    A511 = ictrl[ii++];
    A512 = ictrl[ii++];
    A514 = ictrl[ii++];
    A515 = ictrl[ii++];
    A516 = ictrl[ii++];
    A517 = ictrl[ii++];
    A518 = ictrl[ii++];
    A519 = ictrl[ii++];
    A520 = ictrl[ii++];
    A521 = ictrl[ii++];
    A522 = dctrl[dd++];
    A523 = dctrl[dd++];
    A531 = dctrl[dd++];
    A540 = ictrl[ii++];
    A541 = dctrl[dd++];
    A542 = dctrl[dd++];
    A543 = ictrl[ii++];
    A544 = dctrl[dd++];
    A545 = dctrl[dd++];
    A550 = ictrl[ii++];
    A551 = ictrl[ii++];
    A552 = ictrl[ii++];
    A553 = ictrl[ii++];
    A560 = ictrl[ii++];
    A570 = ictrl[ii++];
    A571_u = dctrl[dd++];
    A571_dir = dctrl[dd++];
    A573 = ictrl[ii++];
    A580 = ictrl[ii++];
    A580_xs = dctrl[dd++];
    A580_xe = dctrl[dd++];
    A580_ys = dctrl[dd++];
    A580_ye = dctrl[dd++];
    A581 = ictrl[ii++];
    A583 = ictrl[ii++];
    A584 = ictrl[ii++];
    A585 = ictrl[ii++];
    A586 = ictrl[ii++];
    A587 = ictrl[ii++];
    A588 = ictrl[ii++];
    A589 = ictrl[ii++];
    A590 = ictrl[ii++];
    A591 = ictrl[ii++];
    A591_x = dctrl[dd++];
    A591_y = dctrl[dd++];
    A591_z = dctrl[dd++];
    A592 = ictrl[ii++];
    A592_x = dctrl[dd++];
    A592_y = dctrl[dd++];
    A592_z = dctrl[dd++];
    A593 = ictrl[ii++];
    A593_x = dctrl[dd++];
    A593_y = dctrl[dd++];
    A593_z = dctrl[dd++];
    A593_phi = dctrl[dd++];
    A593_theta = dctrl[dd++];
    A593_psi = dctrl[dd++];
    A594 = ictrl[ii++];
    
    
    B10 = ictrl[ii++];
    B20 = ictrl[ii++];
    B23 = ictrl[ii++];
    B29 = dctrl[dd++];
    B30 = ictrl[ii++];
    B31 = dctrl[dd++];
    B32 = ictrl[ii++];
    B32_x = dctrl[dd++];
    B32_y = dctrl[dd++];
    B32_z = dctrl[dd++];
    B33 = ictrl[ii++];
    B50 = dctrl[dd++];
    B51 = dctrl[dd++];
    B52 = dctrl[dd++];
    B53 = dctrl[dd++];
    B54 = dctrl[dd++];
    B55 = dctrl[dd++];
    B56 = dctrl[dd++];
    B60 = ictrl[ii++];
    B61 = ictrl[ii++];
    B71 = ictrl[ii++];
    B75 = ictrl[ii++];
    B76 = ictrl[ii++];
    B77 = ictrl[ii++];
    B81 = ictrl[ii++];
    B81_1 = dctrl[dd++];
    B81_2 = dctrl[dd++];
    B81_3 = dctrl[dd++];
    B82 = ictrl[ii++];
    B83 = dctrl[dd++];
    B84 = ictrl[ii++];
    B85 = ictrl[ii++];
    B86 = ictrl[ii++];
    B87 = ictrl[ii++];
    B87_1 = dctrl[dd++];
    B87_2 = dctrl[dd++];
    B88 = dctrl[dd++];
    B89 = ictrl[ii++];
    B90 = ictrl[ii++];
    B91 = ictrl[ii++];
    B91_1 = dctrl[dd++];
    B91_2 = dctrl[dd++];
    B92 = ictrl[ii++];
    B93 = ictrl[ii++];
    B93_1 = dctrl[dd++];
    B93_2 = dctrl[dd++];
    B94 = ictrl[ii++];
    B94_wdt = dctrl[dd++];
    B96_1 = dctrl[dd++];
    B96_2 = dctrl[dd++];
    B98 = ictrl[ii++];
    B99 = ictrl[ii++];
    B101 = ictrl[ii++];
    B102 = dctrl[dd++];
    B105 = ictrl[ii++];
    B105_1 = dctrl[dd++];
    B105_2 = dctrl[dd++];
    B105_3 = dctrl[dd++];
    B106 = ictrl[ii++];
    B107 = ictrl[ii++];
    B110 = ictrl[ii++];
    B110_zs = dctrl[dd++];
    B110_ze = dctrl[dd++];
    B108 = ictrl[ii++];
    B111_zs = dctrl[dd++];
    B111_ze = dctrl[dd++];
    B112_zs = dctrl[dd++];
    B112_z2 = dctrl[dd++];
    B112_ze = dctrl[dd++];
    B115 = ictrl[ii++];
    B116 = ictrl[ii++];
    B117 = dctrl[dd++];
    B120 = dctrl[dd++];
    B122 = dctrl[dd++];
    B123 = dctrl[dd++];
    B125 = ictrl[ii++];
    B125_y = dctrl[dd++];
    B127 = ictrl[ii++];
    B130 = ictrl[ii++];
    B131 = dctrl[dd++];
    B132_s = dctrl[dd++];
    B132_e = dctrl[dd++];
    B133 = ictrl[ii++];
    B134 = dctrl[dd++];
    B135 = dctrl[dd++];
    B136 = ictrl[ii++];
    B138 = ictrl[ii++];
    B138_1 = ictrl[ii++];
    B138_2 = ictrl[ii++];
    B139 = ictrl[ii++];
    B160 = ictrl[ii++];
    B170 = ictrl[ii++];
    B180 = ictrl[ii++];
    B181 = ictrl[ii++];
    B181_1 = dctrl[dd++];
    B181_2 = dctrl[dd++];
    B181_3 = dctrl[dd++];
    B182 = ictrl[ii++];
    B182_1 = dctrl[dd++];
    B182_2 = dctrl[dd++];
    B182_3 = dctrl[dd++];
    B183 = ictrl[ii++];
    B183_1 = dctrl[dd++];
    B183_2 = dctrl[dd++];
    B183_3 = dctrl[dd++];
    B191 = ictrl[ii++];
    B191_1 = dctrl[dd++];
    B191_2 = dctrl[dd++];
    B191_3 = dctrl[dd++];
    B191_4 = dctrl[dd++];
    B192 = ictrl[ii++];
    B192_1 = dctrl[dd++];
    B192_2 = dctrl[dd++];
    B192_3 = dctrl[dd++];
    B192_4 = dctrl[dd++];
    B194_s = dctrl[dd++];
    B194_e = dctrl[dd++];
    B240 = ictrl[ii++];
    B241 = ictrl[ii++];
    B242 = ictrl[ii++];
    B243 = ictrl[ii++];
    B260 = dctrl[dd++];
    B264 = dctrl[dd++];
    B267 = dctrl[dd++];
    B269 = ictrl[ii++];
    B270 = ictrl[ii++];
    B274 = ictrl[ii++];
    B281 = ictrl[ii++];
    B282 = ictrl[ii++];
    B291 = ictrl[ii++];
    B295 = ictrl[ii++];
    B307 = ictrl[ii++];
    B308 = ictrl[ii++];
    B309 = dctrl[dd++];
    B310 = ictrl[ii++];
    B321 = ictrl[ii++];
    B322 = ictrl[ii++];
    B411 = ictrl[ii++];
    B412 = ictrl[ii++];
    B413 = ictrl[ii++];
    B414 = ictrl[ii++];
    B415 = ictrl[ii++];
    B416 = ictrl[ii++];
    B417 = ictrl[ii++];
    B418 = ictrl[ii++];
    B421 = ictrl[ii++];
    B422 = ictrl[ii++];
    B440 = ictrl[ii++];
    B441 = ictrl[ii++];
    B442 = ictrl[ii++];
    
    C1 = dctrl[dd++];
    C2 = dctrl[dd++];
    C3 = dctrl[dd++];
    C4 = dctrl[dd++];
    C5 = dctrl[dd++];
    C9 = ictrl[ii++];
    C10 = ictrl[ii++];
    C15 = ictrl[ii++];
    C20 = ictrl[ii++];
    C50_1=dctrl[dd++];
    C50_2=dctrl[dd++];
    C51 = dctrl[dd++];
    C52 = dctrl[dd++];
    C53 = dctrl[dd++];
    C54 = dctrl[dd++];
    C55 = dctrl[dd++];
    C56 = dctrl[dd++];
    C57_1 = dctrl[dd++];
    C57_2 = dctrl[dd++];
    C57_3 = dctrl[dd++];
    C57_4 = dctrl[dd++];
    C58_1 = dctrl[dd++];
    C58_2 = dctrl[dd++];
    C58_3 = dctrl[dd++];
    C58_4 = dctrl[dd++];
    C75 = ictrl[ii++];

    D10 = ictrl[ii++];
    D11 = ictrl[ii++];
    D20 = ictrl[ii++];
    D21 = ictrl[ii++];
    D30 = ictrl[ii++];
    D31 = ictrl[ii++];
    D33 = ictrl[ii++];
    D37 = ictrl[ii++];
    
    F10 = ictrl[ii++];
    F30 = ictrl[ii++];
    F31 = ictrl[ii++];
    F32 = ictrl[ii++];
    F33 = dctrl[dd++];
    F34 = ictrl[ii++];
    F35 = ictrl[ii++];
    F36 = ictrl[ii++];
    F39 = dctrl[dd++];
    F40 = ictrl[ii++];
    F42 = dctrl[dd++];
    F43 = dctrl[dd++];
    F44 = ictrl[ii++];
    F45 = dctrl[dd++];
    F46 = ictrl[ii++];
    F47 = ictrl[ii++];
    F49 = ictrl[ii++];
    F50 = ictrl[ii++];
    F50_flag = ictrl[ii++];
    F51 = dctrl[dd++];
    F52 = dctrl[dd++];
    F53 = dctrl[dd++];
    F54 = dctrl[dd++];
    F55 = dctrl[dd++];
    F56 = dctrl[dd++];
    F57_1 = dctrl[dd++];
    F57_2 = dctrl[dd++];
    F57_3 = dctrl[dd++];
    F57_4 = dctrl[dd++];
    F58_1 = dctrl[dd++];
    F58_2 = dctrl[dd++];
    F58_3 = dctrl[dd++];
    F58_4 = dctrl[dd++];
    F59_xm = dctrl[dd++];
    F59_ym = dctrl[dd++];
    F59_zs = dctrl[dd++];
    F59_ze = dctrl[dd++];
    F59_r = dctrl[dd++];
    F60 = dctrl[dd++];
    F61 = dctrl[dd++];
    F62 = dctrl[dd++];
    F63 = dctrl[dd++];
    F64 = ictrl[ii++];
    F64_xs = dctrl[dd++];
    F64_ys = dctrl[dd++];
    F64_zs = dctrl[dd++];
    F64_alpha = dctrl[dd++];
    F70 = ictrl[ii++];
    F71 = ictrl[ii++];
    F72 = ictrl[ii++];
    F80 = ictrl[ii++];
    F84 = dctrl[dd++];
    F85 = ictrl[ii++];
    F150 = ictrl[ii++];
    F151 = ictrl[ii++];
    F300 = ictrl[ii++];
    F305 = ictrl[ii++];
    F310 = ictrl[ii++];
    F321 = dctrl[dd++];
    F322 = dctrl[dd++];
    F323 = dctrl[dd++];
    F350 = ictrl[ii++];
    F360 = dctrl[dd++];
    F361 = dctrl[dd++];
    F362 = dctrl[dd++];
    F369 = ictrl[ii++];
    F370 = ictrl[ii++];
    F371 = ictrl[ii++];
    F374 = ictrl[ii++];
    F375 = ictrl[ii++];
    F378 = ictrl[ii++];
    F379 = ictrl[ii++];
    F380 = dctrl[dd++];
    F381 = dctrl[dd++];
    F382 = dctrl[dd++];
    F390 = ictrl[ii++];
    F391 = ictrl[ii++];
    F394 = ictrl[ii++];
    F395 = ictrl[ii++];
    F398 = ictrl[ii++];
    F399 = ictrl[ii++];
    

    G1  = ictrl[ii++];
    G2  = ictrl[ii++];
    G3  = ictrl[ii++];
    G10 = ictrl[ii++];
    G11 = ictrl[ii++];
    G12 = ictrl[ii++];
    G20 = ictrl[ii++];
    G21 = ictrl[ii++];
    G22 = ictrl[ii++];
    G30 = ictrl[ii++];
    G40 = ictrl[ii++];

    H1 = dctrl[dd++];
    H2 = dctrl[dd++];
    H3 = ictrl[ii++];
    H4 = ictrl[ii++];
    H4_beta1 = dctrl[dd++];
    H4_beta2 = dctrl[dd++];
    H9 = ictrl[ii++];
    H10 = ictrl[ii++];
    H15 = ictrl[ii++];
    H50_1=dctrl[dd++];
    H50_2=dctrl[dd++];
    H51 = dctrl[dd++];
    H52 = dctrl[dd++];
    H53 = dctrl[dd++];
    H54 = dctrl[dd++];
    H55 = dctrl[dd++];
    H56 = dctrl[dd++];
    H57_1 = dctrl[dd++];
    H57_2 = dctrl[dd++];
    H57_3 = dctrl[dd++];
    H57_4 = dctrl[dd++];
    H58_1 = dctrl[dd++];
    H58_2 = dctrl[dd++];
    H58_3 = dctrl[dd++];
    H58_4 = dctrl[dd++];
    H61 = ictrl[ii++];
    H61_T = dctrl[dd++];
    H62 = ictrl[ii++];
    H62_T = dctrl[dd++];
    H63 = ictrl[ii++];
    H63_T = dctrl[dd++];
    H64 = ictrl[ii++];
    H64_T = dctrl[dd++];
    H65 = ictrl[ii++];
    H65_T = dctrl[dd++];
    H66 = ictrl[ii++];
    H66_T = dctrl[dd++];


    I10 = ictrl[ii++];
    I11 = ictrl[ii++];
    I12 = ictrl[ii++];
    I13 = ictrl[ii++];
    I21 = ictrl[ii++];
    I30 = ictrl[ii++];
    I40 = ictrl[ii++];
    I41 = ictrl[ii++];
    I44 = ictrl[ii++];
    I50 = dctrl[dd++];
    I55 = dctrl[dd++];
    I56 = ictrl[ii++];
    I58_1 = dctrl[dd++];
    I58_2 = dctrl[dd++];
    I230 = ictrl[ii++];
    I231 = dctrl[dd++];
    I232 = dctrl[dd++];
    I233 = dctrl[dd++];
    I240 = ictrl[ii++];
    I241 = dctrl[dd++];

    M10 = ictrl[ii++];

    N10 = ictrl[ii++];
    N11 = ictrl[ii++];
    N18 = ictrl[ii++];
    N20 = ictrl[ii++];
    N22 = ictrl[ii++];
    N23 = ictrl[ii++];
    N24 = ictrl[ii++];
    N25 = ictrl[ii++];
    N26 = ictrl[ii++];
    N40 = ictrl[ii++];
    N41 = dctrl[dd++];
    N43 = dctrl[dd++];
    N44 = dctrl[dd++];
    N45 = ictrl[ii++];
    N46 = ictrl[ii++];
    N47 = dctrl[dd++];
    N48 = ictrl[ii++];
    N49 = dctrl[dd++];
    N50 = ictrl[ii++];
    N60 = ictrl[ii++];
    N61 = dctrl[dd++];
    

    P10 = ictrl[ii++];
    P11 = ictrl[ii++];
    P12 = ictrl[ii++];
    P15 = ictrl[ii++];
    P16 = ictrl[ii++];
    P20 = ictrl[ii++];
    P21 = ictrl[ii++];
    P22 = dctrl[dd++];
    P23 = ictrl[ii++];
    P24 = ictrl[ii++];
    P25 = ictrl[ii++];
    P26 = ictrl[ii++];
    P27 = ictrl[ii++];
    P28 = ictrl[ii++];
    P29 = ictrl[ii++];
    P30 = dctrl[dd++];
    P34 = dctrl[dd++];
    P35 = ictrl[ii++];
    P40 = ictrl[ii++];
    P41 = ictrl[ii++];
    P42 = dctrl[dd++];
    P43 = ictrl[ii++];
    P43_xs = dctrl[dd++];
    P43_xe = dctrl[dd++];
    P43_ys = dctrl[dd++];
    P43_ye = dctrl[dd++];
    P44 = ictrl[ii++];
    P45 = ictrl[ii++];
    P46 = ictrl[ii++];
    P46_is = ictrl[ii++];
    P46_ie = ictrl[ii++];
    P47 = ictrl[ii++];
    P47_ts = dctrl[dd++];
    P47_te = dctrl[dd++];
    P50 = ictrl[ii++];
    P51 = ictrl[ii++];
    P52 = ictrl[ii++];
    P53 = ictrl[ii++];
    P54 = ictrl[ii++];
    P55 = dctrl[dd++];
    P56 = ictrl[ii++];
    P57 = ictrl[ii++];
    P58 = ictrl[ii++];
    P59 = ictrl[ii++];
    P61 = ictrl[ii++];
    P62 = ictrl[ii++];
    P63 = ictrl[ii++];
    P64 = ictrl[ii++];
    P65 = ictrl[ii++];
    P66 = ictrl[ii++];
    P71 = ictrl[ii++];
    P72 = ictrl[ii++];
    P73 = ictrl[ii++];
    P74 = ictrl[ii++];
    P75 = ictrl[ii++];
    P76 = ictrl[ii++];
    P77 = ictrl[ii++];
    P78 = ictrl[ii++];
    P79 = ictrl[ii++];
    P80 = ictrl[ii++];
    P81 = ictrl[ii++];
    P82 = ictrl[ii++];
    P85 = ictrl[ii++];
    P88 = ictrl[ii++];
    P91 = dctrl[dd++];
    P92 = ictrl[ii++];
    P101 = ictrl[ii++];
    P101_xm = dctrl[dd++];
    P101_ym = dctrl[dd++];
    P101_zs = dctrl[dd++];
    P101_ze = dctrl[dd++];
    P101_r1 = dctrl[dd++];
    P101_r2 = dctrl[dd++];
    P110 = ictrl[ii++];
    P111 = dctrl[dd++];
    P120 = ictrl[ii++];
    P121 = ictrl[ii++];
    P122 = ictrl[ii++];
    P123 = ictrl[ii++];
    P124 = ictrl[ii++];
    P125 = ictrl[ii++];
    P126 = ictrl[ii++];
    P131 = ictrl[ii++];
    P132 = ictrl[ii++];
    P133 = ictrl[ii++];
    P134 = ictrl[ii++];
    P140 = ictrl[ii++];
    P141 = dctrl[dd++];
    P151 = ictrl[ii++];
    P152 = ictrl[ii++];
    P166 = ictrl[ii++];
    P167 = ictrl[ii++];
    P168 = ictrl[ii++];
    P180 = ictrl[ii++];
    P181 = ictrl[ii++];
    P182 = dctrl[dd++];
    P184 = ictrl[ii++];
    P185 = ictrl[ii++];
    P190 = ictrl[ii++];
    P191 = ictrl[ii++];
    P192 = dctrl[dd++];
    P194 = ictrl[ii++];
    P195 = ictrl[ii++];
    P230 = ictrl[ii++];
    P240 = ictrl[ii++];
    P351 = ictrl[ii++];
    P352 = ictrl[ii++];

    Q10 = ictrl[ii++];
    Q11 = ictrl[ii++];
    Q12 = ictrl[ii++];
    Q13 = ictrl[ii++];
    Q14 = dctrl[dd++];
    Q15 = dctrl[dd++];
    Q16 = dctrl[dd++];
    Q17 = dctrl[dd++];
    Q20 = ictrl[ii++];
    Q22 = dctrl[dd++];
    Q23 = dctrl[dd++];
    Q24 = ictrl[ii++];
    Q25 = dctrl[dd++];
    Q29 = ictrl[ii++];
    Q30 = dctrl[dd++];
    Q41 = dctrl[dd++];
    Q43 = ictrl[ii++];
    Q61 = ictrl[ii++];
    Q73 = ictrl[ii++];
    Q101 = ictrl[ii++];
    Q102 = dctrl[dd++];
    Q110 = ictrl[ii++];
    Q111 = ictrl[ii++];
    Q120 = ictrl[ii++];
    Q121 = ictrl[ii++];
    Q122 = ictrl[ii++];
    Q180 = ictrl[ii++];
    Q181 = ictrl[ii++];
    Q182 = dctrl[dd++];
    Q183 = ictrl[ii++];
    Q200 = ictrl[ii++];
    Q201 = ictrl[ii++];
    Q202 = ictrl[ii++];

    S10 = ictrl[ii++];
    S11 = ictrl[ii++];
    S12 = ictrl[ii++];
    S13 = dctrl[dd++];
    S14 = dctrl[dd++];
    S15 = ictrl[ii++];
    S16 = ictrl[ii++];
    S17 = ictrl[ii++];
    S19 = dctrl[dd++];
    S20 = dctrl[dd++];
    S21 = dctrl[dd++];
    S22 = dctrl[dd++];
    S23 = dctrl[dd++];
    S24 = dctrl[dd++];
    S25 = ictrl[ii++];
    S26_a = dctrl[dd++];
    S26_b = dctrl[dd++];
    S27 = ictrl[ii++];
    S30 = dctrl[dd++];
    S31 = ictrl[ii++];
    S32 = ictrl[ii++];
    S33 = ictrl[ii++];
    S34 = ictrl[ii++];
    S37 = ictrl[ii++];
    S41 = ictrl[ii++];
    S42 = ictrl[ii++];
    S43 = ictrl[ii++];
    S44 = ictrl[ii++];
    S45 = dctrl[dd++];
    S46 = dctrl[dd++];
    S47 = dctrl[dd++];
    S48 = dctrl[dd++];
    S50 = ictrl[ii++];
    S57 = dctrl[dd++];
    S60 = dctrl[dd++];
    S71 = dctrl[dd++];
    S72 = dctrl[dd++];
    S73 = ictrl[ii++];
    S77 = ictrl[ii++];
    S77_xs = dctrl[dd++];
    S77_xe = dctrl[dd++];
    S78 = ictrl[ii++];
    S79 = ictrl[ii++];
    S80 = ictrl[ii++];
    S81 = dctrl[dd++];
    S82 = dctrl[dd++];
    S83 = ictrl[ii++];
    S84 = ictrl[ii++];
    S85 = ictrl[ii++];
    S90 = ictrl[ii++];
    S91 = ictrl[ii++];
    S92 = dctrl[dd++];
    S93 = dctrl[dd++];
    S100 = ictrl[ii++];
    S101 = ictrl[ii++];
    
    T10 = ictrl[ii++];
    T12 = ictrl[ii++];
    T21 = ictrl[ii++];
    T31 = dctrl[dd++];
    T32 = dctrl[dd++];
    T33 = ictrl[ii++];
    T35 = dctrl[dd++];
    T36 = ictrl[ii++];
    T37 = dctrl[dd++];
    T38 = dctrl[dd++];
    T39 = ictrl[ii++];
    T41 = ictrl[ii++];
    T42 = dctrl[dd++];
    T43 = dctrl[dd++];
    T44 = dctrl[dd++];
    T45 = ictrl[ii++];

    W1  = dctrl[dd++];
    W2  = dctrl[dd++];
    W3  = dctrl[dd++];
    W4  = dctrl[dd++];
    W5  = dctrl[dd++];
    W6  = dctrl[dd++];
    W7  = dctrl[dd++];
    W10 = dctrl[dd++];
    W11 = ictrl[ii++];
    W11_u = dctrl[dd++];
    W11_v = dctrl[dd++];
    W11_w = dctrl[dd++];
    W12 = ictrl[ii++];
    W12_u = dctrl[dd++];
    W12_v = dctrl[dd++];
    W12_w = dctrl[dd++];
    W13 = ictrl[ii++];
    W13_u = dctrl[dd++];
    W13_v = dctrl[dd++];
    W13_w = dctrl[dd++];
    W14 = ictrl[ii++];
    W14_u = dctrl[dd++];
    W14_v = dctrl[dd++];
    W14_w = dctrl[dd++];
    W15 = ictrl[ii++];
    W15_u = dctrl[dd++];
    W15_v = dctrl[dd++];
    W15_w = dctrl[dd++];
    W16 = ictrl[ii++];
    W16_u = dctrl[dd++];
    W16_v = dctrl[dd++];
    W16_w = dctrl[dd++];
    W20 = dctrl[dd++];
    W21 = dctrl[dd++];
    W22 = dctrl[dd++];
    W29_x = dctrl[dd++];
    W29_y = dctrl[dd++];
    W29_z = dctrl[dd++];
    W30 = ictrl[ii++];
    W31 = dctrl[dd++];
    W41 = ictrl[ii++];
    W50 = dctrl[dd++];
    W50_air = ictrl[ii++];
    W90 = ictrl[ii++];
    W95 = dctrl[dd++];
    W96 = dctrl[dd++];
    W97 = dctrl[dd++];
    W98 = dctrl[dd++];
    W101 = ictrl[ii++];
    W102_phi = dctrl[dd++];
    W102_c = dctrl[dd++];
    W103 = dctrl[dd++];
    W104 = dctrl[dd++];
    W110 = ictrl[ii++];
    W111 = ictrl[ii++];
    W112 = dctrl[dd++];
    
    X10 = ictrl[ii++];
    X11_u = ictrl[ii++];
    X11_v = ictrl[ii++];
    X11_w = ictrl[ii++];
    X11_p = ictrl[ii++];
    X11_q = ictrl[ii++];
    X11_r = ictrl[ii++];
    X12 = ictrl[ii++];
    X14 = ictrl[ii++];
    X19 = ictrl[ii++];
    X21 = ictrl[ii++];
    X21_d = dctrl[dd++];
    X22 = ictrl[ii++];
    X22_m = dctrl[dd++];
    X23 = ictrl[ii++];
    X23_x = dctrl[dd++];
    X23_y = dctrl[dd++];
    X23_z = dctrl[dd++];
    X24 = ictrl[ii++];
    X24_Ix = dctrl[dd++];
    X24_Iy = dctrl[dd++];
    X24_Iz = dctrl[dd++];
    X25_Cp = dctrl[dd++];
    X25_Cq = dctrl[dd++];
    X25_Cr = dctrl[dd++];
    X26_Cu = dctrl[dd++];
    X26_Cv = dctrl[dd++];
    X26_Cw = dctrl[dd++];
    X31 = ictrl[ii++];
    X32 = ictrl[ii++];
    X33 = ictrl[ii++];
    X34 = ictrl[ii++];
    X39 = ictrl[ii++];
    X40 = ictrl[ii++];
    X41 = dctrl[dd++];
    X42 = dctrl[dd++];
    X43 = dctrl[dd++];
    X44 = dctrl[dd++];
    X45 = ictrl[ii++];
    X46 = ictrl[ii++];
    X48 = ictrl[ii++];
    X49 = ictrl[ii++];
    X50 = ictrl[ii++];
    X60 = ictrl[ii++];
    X100 = ictrl[ii++];
    X100_x = dctrl[dd++];
    X100_y = dctrl[dd++];
    X100_z = dctrl[dd++];
    X101 = ictrl[ii++];
    X101_phi = dctrl[dd++];
    X101_theta = dctrl[dd++];
    X101_psi = dctrl[dd++];
    X102 = ictrl[ii++];
    X102_u = dctrl[dd++];
    X102_v = dctrl[dd++];
    X102_w = dctrl[dd++];
    X103 = ictrl[ii++];
    X103_p = dctrl[dd++];
    X103_q = dctrl[dd++];
    X103_r = dctrl[dd++];
    X110 = ictrl[ii++];
    X120 = ictrl[ii++];
    X120_rad = dctrl[dd++];
    X120_xc = dctrl[dd++];
    X120_yc = dctrl[dd++];
    X120_zc = dctrl[dd++];
    X131 = ictrl[ii++];
    X131_rad = dctrl[dd++];
    X131_h = dctrl[dd++];
    X131_xc = dctrl[dd++];
    X131_yc = dctrl[dd++];
    X131_zc = dctrl[dd++];
    X132 = ictrl[ii++];
    X132_rad = dctrl[dd++];
    X132_h = dctrl[dd++];
    X132_xc = dctrl[dd++];
    X132_yc = dctrl[dd++];
    X132_zc = dctrl[dd++];
    X133 = ictrl[ii++];
    X133_rad = dctrl[dd++];
    X133_h = dctrl[dd++];
    X133_xc = dctrl[dd++];
    X133_yc = dctrl[dd++];
    X133_zc = dctrl[dd++];
    X153 = ictrl[ii++];
    X153_xs = dctrl[dd++];
    X153_xe = dctrl[dd++];
    X153_ys = dctrl[dd++];
    X153_ye = dctrl[dd++];
    X153_zs = dctrl[dd++];
    X153_ze = dctrl[dd++];
    X163 = ictrl[ii++];
    X164 = ictrl[ii++];
    X180 = ictrl[ii++];
    X181 = ictrl[ii++];
    X181_x = dctrl[dd++];
    X181_y = dctrl[dd++];
    X181_z = dctrl[dd++];
    X182 = ictrl[ii++];
    X182_x = dctrl[dd++];
    X182_y = dctrl[dd++];
    X182_z = dctrl[dd++];
    X183 = ictrl[ii++];
    X183_x = dctrl[dd++];
    X183_y = dctrl[dd++];
    X183_z = dctrl[dd++];
    X183_phi = dctrl[dd++];
    X183_theta = dctrl[dd++];
    X183_psi = dctrl[dd++];
    X185 = ictrl[ii++];
    X186 = dctrl[dd++];
    X188 = ictrl[ii++];
    X205 = ictrl[ii++];
    X206 = ictrl[ii++];
    X206_ts = dctrl[dd++];
    X206_te = dctrl[dd++];
    X207 = ictrl[ii++];
    X207_ts = dctrl[dd++];
    X207_te = dctrl[dd++];
    X210 = ictrl[ii++];
    X210_u = dctrl[dd++];
    X210_v = dctrl[dd++];
    X210_w = dctrl[dd++];
    X211 = ictrl[ii++];
    X211_p = dctrl[dd++];
    X211_q = dctrl[dd++];
    X211_r = dctrl[dd++];    
    X240 = ictrl[ii++];
    X241 = dctrl[dd++];
    X242_x = dctrl[dd++];
    X242_y = dctrl[dd++];
    X242_z = dctrl[dd++];
    X243 = dctrl[dd++];
    X310 = ictrl[ii++];    
    X311 = ictrl[ii++];
    X312 = ictrl[ii++];    
    X313 = ictrl[ii++];    
    X314 = ictrl[ii++];    
    X315 = ictrl[ii++];    
    X320 = ictrl[ii++];
    X321 = ictrl[ii++];
    X323_m = dctrl[dd++];
    X323_d = dctrl[dd++];
    X323_l = dctrl[dd++];
    X324 = ictrl[ii++];
    X325_dt = dctrl[dd++];
    X325_relX = dctrl[dd++];
    X325_relY = dctrl[dd++];
    X325_relZ = dctrl[dd++];
    X400 = ictrl[ii++];
    X401_p0 = dctrl[dd++];
    X401_cl = dctrl[dd++]; 
    X401_cb = dctrl[dd++]; 
    X401_a = dctrl[dd++];
    

    Y1 = ictrl[ii++];
    Y2 = ictrl[ii++];
    Y3 = ictrl[ii++];
    Y4 = ictrl[ii++];
    Y5 = ictrl[ii++];
    Y40 = ictrl[ii++];
    Y50 = ictrl[ii++];
    Y60 = ictrl[ii++];
    Y71 = ictrl[ii++];
    Y72 = ictrl[ii++];
    Y73 = ictrl[ii++];
    Y74 = ictrl[ii++];
    
    Z10 = ictrl[ii++];
    Z11 = ictrl[ii++];
    Z12_cdx = dctrl[dd++];
    Z12_cdy = dctrl[dd++];
    Z12_cdz = dctrl[dd++];
    Z12_ckx = dctrl[dd++];
    Z12_cky = dctrl[dd++];
    Z12_ckz = dctrl[dd++];

// --------------------------    
    
    if(A581>0)
    {
    Darray(A581_xs,A581);
    Darray(A581_xe,A581);
    Darray(A581_ys,A581);
    Darray(A581_ye,A581);
    Darray(A581_zs,A581);
    Darray(A581_ze,A581);
    }
    
    if(A583>0)
    {
    Darray(A583_xc,A583);
    Darray(A583_zc,A583);
    Darray(A583_ys,A583);
    Darray(A583_ye,A583);
    Darray(A583_r,A583);
    }

    if(A584>0)
    {
    Darray(A584_xc,A584);
    Darray(A584_yc,A584);
    Darray(A584_zs,A584);
    Darray(A584_ze,A584);
    Darray(A584_r,A584);
    }
    
    if(A585>0)
    {
    Darray(A585_xm1,A585);
    Darray(A585_ym1,A585);
    Darray(A585_zm1,A585);
    Darray(A585_r1,A585);
    Darray(A585_xm2,A585);
    Darray(A585_ym2,A585);
    Darray(A585_zm2,A585);
    Darray(A585_r2,A585);
    }
    
    if(A586>0)
    {
    Darray(A586_xm,A586);
    Darray(A586_ym,A586);
    Darray(A586_zm,A586);
    Darray(A586_r,A586);
    }
    
    if(A587>0)
    {
    Darray(A587_xs,A587);
    Darray(A587_xe,A587);
    Darray(A587_ys,A587);
    Darray(A587_ye,A587);
    Darray(A587_zs,A587);
    Darray(A587_ze,A587);
    }
    
    if(A588>0)
    {
    Darray(A588_xs,A588);
    Darray(A588_xe,A588);
    Darray(A588_ys,A588);
    Darray(A588_ye,A588);
    Darray(A588_zs,A588);
    Darray(A588_ze,A588);
    }
    
    if(A589>0)
    {
    Darray(A589_xs,A589);
    Darray(A589_xe,A589);
    Darray(A589_ys,A589);
    Darray(A589_ye,A589);
    Darray(A589_zs,A589);
    Darray(A589_ze,A589);
    }

    if(B71>0)
    {
    Darray(B71_val,B71);
    Darray(B71_dist,B71);
    Darray(B71_b,B71);
    Darray(B71_x,B71);
    Darray(B71_y,B71);
    }
    
    if(B106>0)
    {
    Darray(B106_b,B106);
    Darray(B106_x,B106);
    Darray(B106_y,B106);
    }
    
    if(B107>0)
    {
    Darray(B107_xs,B107);
    Darray(B107_xe,B107);
    Darray(B107_ys,B107);
    Darray(B107_ye,B107);
    Darray(B107_d,B107);
    }
    
    if(B108>0)
    {
    Darray(B108_xs,B108);
    Darray(B108_xe,B108);
    Darray(B108_ys,B108);
    Darray(B108_ye,B108);
    Darray(B108_d,B108);
    }
    
    if(B240>0)
    {    
    Darray(B240_C,B240);
    Darray(B240_D,B240);
    Darray(B240_xs,B240);
    Darray(B240_xe,B240);
    Darray(B240_ys,B240);
    Darray(B240_ye,B240);
    Darray(B240_zs,B240);
    Darray(B240_ze,B240);
    }
    
    if(B270>0)
    {    
    Darray(B270_xs,B270);
    Darray(B270_xe,B270);
    Darray(B270_ys,B270);
    Darray(B270_ye,B270);
    Darray(B270_zs,B270);
    Darray(B270_ze,B270);
    Darray(B270_n,B270);
    Darray(B270_d50,B270);
    Darray(B270_alpha,B270);
    Darray(B270_beta,B270);
    }
    
    if(B274>0)
    {    
    Darray(B274_xc,B274);
    Darray(B274_yc,B274);
    Darray(B274_zs,B274);
    Darray(B274_ze,B274);
    Darray(B274_r,B274);
    Darray(B274_n,B274);
    Darray(B274_d50,B274);
    Darray(B274_alpha,B274);
    Darray(B274_beta,B274);
    }
    
    if(B281>0)
    {    
    Darray(B281_xs,B281);
    Darray(B281_xe,B281);
    Darray(B281_ys,B281);
    Darray(B281_ye,B281);
    Darray(B281_zs,B281);
    Darray(B281_ze,B281);
    Darray(B281_n,B281);
    Darray(B281_d50,B281);
    Darray(B281_alpha,B281);
    Darray(B281_beta,B281);
    }

    if(B282>0)
    {    
    Darray(B282_xs,B282);
    Darray(B282_xe,B282);
    Darray(B282_ys,B282);
    Darray(B282_ye,B282);
    Darray(B282_zs,B282);
    Darray(B282_ze,B282);
    Darray(B282_n,B282);
    Darray(B282_d50,B282);
    Darray(B282_alpha,B282);
    Darray(B282_beta,B282);
    }
    
    if(B291>0)
    {    
    Darray(B291_xs,B291);
    Darray(B291_xe,B291);
    Darray(B291_ys,B291);
    Darray(B291_ye,B291);
    Darray(B291_zs,B291);
    Darray(B291_ze,B291);
    Darray(B291_d,B291);
    Darray(B291_n,B291);
    Darray(B291_d50,B291);
    Darray(B291_alpha,B291);
    Darray(B291_beta,B291);
    }
    
    if(B310>0)
    {    
    Darray(B310_xs,B310);
    Darray(B310_xe,B310);
    Darray(B310_ys,B310);
    Darray(B310_ye,B310);
    Darray(B310_zs,B310);
    Darray(B310_ze,B310);
    Darray(B310_N,B310);
    Darray(B310_D,B310);
    Darray(B310_Cd,B310);
    }

    if(B321>0)
    {    
    Darray(B321_xs,B321);
    Darray(B321_xe,B321);
    Darray(B321_ys,B321);
    Darray(B321_ye,B321);
    Darray(B321_zs,B321);
    Darray(B321_ze,B321);
    Darray(B321_N,B321);
    Darray(B321_D,B321);
    Darray(B321_Cd,B321);
    }

    if(B322>0)
    {    
    Darray(B322_xs,B322);
    Darray(B322_xe,B322);
    Darray(B322_ys,B322);
    Darray(B322_ye,B322);
    Darray(B322_zs,B322);
    Darray(B322_ze,B322);
    Darray(B322_N,B322);
    Darray(B322_D,B322);
    Darray(B322_Cd,B322);
    }
    
    if(B411>0)
    {
    Iarray(B411_ID,B411);
    Darray(B411_Q,B411);
    }
    
    if(B412>0)
    {
    Iarray(B412_ID,B412);
    Darray(B412_pressBC,B412);
    }
    
    if(B413>0)
    {
    Iarray(B413_ID,B413);
    Darray(B413_h,B413);
    }
    
    if(B414>0)
    {
    Iarray(B414_ID,B414);
    Darray(B414_Uio,B414);
    }
    
    if(B415>0)
    {
    Iarray(B415_ID,B415);
    Darray(B415_U,B415);
    Darray(B415_V,B415);
    Darray(B415_W,B415);
    }
    
    if(B416>0)
    {
    Iarray(B416_ID,B416);
    Darray(B416_alpha,B416);
    }
    
    if(B417>0)
    {
    Iarray(B417_ID,B417);
    Darray(B417_Nx,B417);
    Darray(B417_Ny,B417);
    Darray(B417_Nz,B417);
    }
    
    if(B418>0)
    {
    Iarray(B418_ID,B418);
    Iarray(B418_pio,B418);
    }
    
    if(B421>0)
    {
    Iarray(B421_ID,B421);
    Iarray(B421_Q,B421);
    }
    
    if(B422>0)
    {
    Iarray(B422_ID,B422);
    Iarray(B422_FSF,B422);
    }
    
    if(B440>0)
    {
    Iarray(B440_ID,B440);
    Iarray(B440_face,B440);
    Darray(B440_xs,B440);
    Darray(B440_xe,B440);
    Darray(B440_ys,B440);
    Darray(B440_ye,B440);
    }
    
    if(B441>0)
    {
    Iarray(B441_ID,B441);
    Iarray(B441_face,B441);
    Darray(B441_xs,B441);
    Darray(B441_xe,B441);
    Darray(B441_ys,B441);
    Darray(B441_ye,B441);
    Darray(B441_zs,B441);
    Darray(B441_ze,B441);
    }
    
    if(B442>0)
    {
    Iarray(B442_ID,B442);
    Iarray(B442_face,B442);
    Darray(B442_xm,B442);
    Darray(B442_ym,B442);
    Darray(B442_zm,B442);
    Darray(B442_r,B442);
    }
    
    
    
    if(C75>0)
    {
    Darray(C75_x,C75);   
    Darray(C75_z,C75);  
    
    Darray(C75_a,C75);      
    Darray(C75_s,C75);  
    Darray(C75_l,C75);  
    Darray(C75_v,C75);  
    }

    if(F70>0)
    {
    Darray(F70_xs,F70);  
    Darray(F70_xe,F70);  
    
    Darray(F70_ys,F70);  
    Darray(F70_ye,F70);  
    
    Darray(F70_zs,F70);  
    Darray(F70_ze,F70);  
    }
    
    if(F71>0)
    {
    Darray(F71_xs,F71);  
    Darray(F71_xe,F71);  
    
    Darray(F71_ys,F71);  
    Darray(F71_ye,F71);  
    
    Darray(F71_zs,F71);  
    Darray(F71_ze,F71);  
    }
    
    if(F72>0)
    {
    Darray(F72_xs,F72);  
    Darray(F72_xe,F72);  
    
    Darray(F72_ys,F72);  
    Darray(F72_ye,F72);  
    
    Darray(F72_h,F72);  
    }

    if(F369>0)
    {
    Darray(F369_x,F369);   
    Darray(F369_z,F369);  
    
    Darray(F369_a,F369);      
    Darray(F369_s,F369);  
    Darray(F369_l,F369);  
    Darray(F369_v,F369);  
    }
    
    if(F370>0)
    {
    Darray(F370_xs,F370);  
    Darray(F370_xe,F370);  
    
    Darray(F370_ys,F370);  
    Darray(F370_ye,F370);  
    
    Darray(F370_zs,F370);  
    Darray(F370_ze,F370);  
    }
    
    if(F371>0)
    {
    Darray(F371_xs,F371);  
    Darray(F371_xe,F371);  
    
    Darray(F371_ys,F371);  
    Darray(F371_ye,F371);  
    
    Darray(F371_zs,F371);  
    Darray(F371_ze,F371);  
    }
    
    if(F374>0)
    {
    Darray(F374_xc,F374);  
    Darray(F374_zc,F374);
    Darray(F374_r,F374);  
    }
    
    if(F375>0)
    {
    Darray(F375_xc,F375);  
    Darray(F375_zc,F375);
    Darray(F375_r,F375);  
    }
    
    if(F378>0)
    {
    Darray(F378_xc,F378);  
    Darray(F378_yc,F378);  
    Darray(F378_zc,F378);
    Darray(F378_r,F378);  
    }
    
    if(F379>0)
    {
    Darray(F379_xc,F379);  
    Darray(F379_yc,F379);  
    Darray(F379_zc,F379);
    Darray(F379_r,F379);  
    }
    
    if(F390>0)
    {
    Darray(F390_xs,F390);  
    Darray(F390_xe,F390);  
    
    Darray(F390_ys,F390);  
    Darray(F390_ye,F390);  
    
    Darray(F390_zs,F390);  
    Darray(F390_ze,F390);  
    }
    
    if(F391>0)
    {
    Darray(F391_xs,F391);  
    Darray(F391_xe,F391);  
    
    Darray(F391_ys,F391);  
    Darray(F391_ye,F391);  
    
    Darray(F391_zs,F391);  
    Darray(F391_ze,F391);  
    }
    
    if(F394>0)
    {
    Darray(F394_xc,F394);  
    Darray(F394_zc,F394);
    Darray(F394_r,F394);  
    }
    
    if(F395>0)
    {
    Darray(F395_xc,F395);  
    Darray(F395_zc,F395);
    Darray(F395_r,F395);  
    }
    
    if(F398>0)
    {
    Darray(F398_xc,F398);  
    Darray(F398_yc,F398);  
    Darray(F398_zc,F398);
    Darray(F398_r,F398);  
    }
    
    if(F399>0)
    {
    Darray(F399_xc,F399);  
    Darray(F399_yc,F399);  
    Darray(F399_zc,F399);
    Darray(F399_r,F399);  
    }
    
    if(P35>0)
    {
    Darray(P35_ts,P35);  
    Darray(P35_te,P35);  
    Darray(P35_dt,P35);  
    }
    
    if(P50>0)
    {
    Darray(P50_x,P50);  
    Darray(P50_y,P50);  
    }
    
    if(P51>0)
    {
    Darray(P51_x,P51);  
    Darray(P51_y,P51);  
    }
    
    if(P52>0)
    Darray(P52_y,P52);

    if(P56>0)
    Darray(P56_x,P56);  

    if(P58>0)
    {
    Darray(P58_x,P58);  
    Darray(P58_y,P58); 
    Darray(P58_T,P58);  
    }
    
    if(P61>0)
    {
    Darray(P61_x,P61);  
    Darray(P61_y,P61); 
    Darray(P61_z,P61);  
    }
    
    if(P62>0)
    {
    Darray(P62_xs,P62); 
    Darray(P62_xe,P62); 
    
    Darray(P62_ys,P62); 
    Darray(P62_ye,P62); 
    
    Darray(P62_zs,P62); 
    Darray(P62_ze,P62); 
    }
    
    if(P63>0)
    {
    Darray(P63_x,P63);  
    Darray(P63_y,P63); 
    }

    if(P64>0)
    {
    Darray(P64_x,P64);  
    Darray(P64_y,P64); 
    Darray(P64_z,P64);  
    }
    
    if(P65>0)
    {
    Darray(P65_x,P65);  
    Darray(P65_y,P65); 
    Darray(P65_z,P65);  
    }

    if(P66>0)
    {
    Darray(P66_x,P66);  
    Darray(P66_y,P66); 
    Darray(P66_z,P66);  
    }
    
    if(P81>0)
    {
    Darray(P81_xs,P81); 
    Darray(P81_xe,P81); 
    
    Darray(P81_ys,P81); 
    Darray(P81_ye,P81); 
    
    Darray(P81_zs,P81); 
    Darray(P81_ze,P81); 
    }
    
    if(P85>0)
    {
    Darray(P85_x,P85); 
    Darray(P85_y,P85); 
    Darray(P85_r,P85);
    Darray(P85_cd,P85);
    Darray(P85_cm,P85);
    }
    
    if(P88>0)
    {
    Darray(P88_x,P88); 
    Darray(P88_y,P88); 
    }
    
    if(P121>0)
    {
    Darray(P121_x,P121);  
    Darray(P121_y,P121);  
    }
    
    if(P123>0)
    Darray(P123_y,P123);
    
    if(P124>0)
    Darray(P124_x,P124);
    
    if(P125>0)
    {
    Darray(P125_x,P125);  
    Darray(P125_y,P125);  
    }
    
    if(P133>0)
    Darray(P133_y,P133);
    
    if(P134>0)
    Darray(P134_y,P134);

    if(P140>0)
    {
    Darray(P140_x,P140);
    Darray(P140_y,P140);
    }

    if(P167>0)
    Darray(P167_x,P167); 

    if(P168>0)
    {
    Darray(P168_x,P168);
    Darray(P168_zs,P168);
    Darray(P168_ze,P168);
    }
    
    if(P184>0)
    {
    Iarray(P184_its,P184);  
    Iarray(P184_ite,P184);  
    Iarray(P184_dit,P184);  
    }
    
    if(P185>0)
    {
    Darray(P185_ts,P185);  
    Darray(P185_te,P185);  
    Darray(P185_dt,P185);  
    }

    if(P194>0)
    {
    Iarray(P194_its,P194);  
    Iarray(P194_ite,P194);  
    Iarray(P194_dit,P194);  
    }
    
    if(P195>0)
    {
    Darray(P195_ts,P195);  
    Darray(P195_te,P195);  
    Darray(P195_dt,P195);  
    }
    
    if(P230>0)
    {
    Darray(P230_x,P230);  
    }
    
    if(P240>0)
    {
    Darray(P240_x,P240);  
    }
    
    if(P351>0)
    {
    Darray(P351_x,P351);  
    Darray(P351_y,P351);  
    }
    
    if(P352>0)
    {
    Darray(P352_x,P352);  
    Darray(P352_y,P352);  
    }

    if(Q61>0)
    {
    Darray(Q61_x,Q61); 
    Darray(Q61_y,Q61);
    Darray(Q61_z,Q61);
    Iarray(Q61_i,Q61);
    }
    
    if(Q73>0)
    {
    Darray(Q73_val,Q73);
    Darray(Q73_dist,Q73);
    Darray(Q73_b,Q73);
    Darray(Q73_x,Q73);
    Darray(Q73_y,Q73);
    }
    
    if(Q110>0)
    {
    Darray(Q110_xs,Q110);  
    Darray(Q110_xe,Q110);  
    
    Darray(Q110_ys,Q110);  
    Darray(Q110_ye,Q110);  
    
    Darray(Q110_zs,Q110);  
    Darray(Q110_ze,Q110);  
    }
    if(Q111>0)
    {
    Darray(Q111_xs,Q111);  
    Darray(Q111_xe,Q111);

    Darray(Q111_ys,Q111);  
    Darray(Q111_ye,Q111);  
    
    Darray(Q111_zs,Q111);  
    Darray(Q111_ze,Q111);  
    }

    if(S73>0)
    {
    Darray(S73_val,S73);
    Darray(S73_dist,S73);
    Darray(S73_b,S73);
    Darray(S73_x,S73);
    Darray(S73_y,S73);
    }
    
    if(W41>0)
    {
    Darray(W41_xc,W41);
    Darray(W41_yc,W41);
    Darray(W41_zs,W41);
    Darray(W41_ze,W41);
    Darray(W41_vel,W41);
    Darray(W41_beta,W41);
    }
    
    if(X110>0)
    {
    Darray(X110_xs,X110);  
    Darray(X110_xe,X110);  
    
    Darray(X110_ys,X110);  
    Darray(X110_ye,X110);  
    
    Darray(X110_zs,X110);  
    Darray(X110_ze,X110);  
    }
    
    if(X163>0)
    {
    Darray(X163_x1,X163);
    Darray(X163_y1,X163);
    Darray(X163_z1,X163);
    Darray(X163_x2,X163);
    Darray(X163_y2,X163);
    Darray(X163_z2,X163);
    Darray(X163_x3,X163);
    Darray(X163_y3,X163);
    Darray(X163_z3,X163);
    Darray(X163_x4,X163);
    Darray(X163_y4,X163);
    Darray(X163_z4,X163);
    Darray(X163_x5,X163);
    Darray(X163_y5,X163);
    Darray(X163_z5,X163);
    Darray(X163_x6,X163);
    Darray(X163_y6,X163);
    Darray(X163_z6,X163);
    }
    
    if(X164>0)
    {
    Darray(X164_x1,X164);
    Darray(X164_y1,X164);
    Darray(X164_z1,X164);
    Darray(X164_x2,X164);
    Darray(X164_y2,X164);
    Darray(X164_z2,X164);
    Darray(X164_x3,X164);
    Darray(X164_y3,X164);
    Darray(X164_z3,X164);
    Darray(X164_x4,X164);
    Darray(X164_y4,X164);
    Darray(X164_z4,X164);
    Darray(X164_x5,X164);
    Darray(X164_y5,X164);
    Darray(X164_z5,X164);
    Darray(X164_x6,X164);
    Darray(X164_y6,X164);
    Darray(X164_z6,X164);
    Darray(X164_x7,X164);
    Darray(X164_y7,X164);
    Darray(X164_z7,X164);
    Darray(X164_x8,X164);
    Darray(X164_y8,X164);
    Darray(X164_z8,X164);
    }
    
    if(X311>0)
    {
        Darray(X311_xs,X311);  
        Darray(X311_xe,X311);  
        Darray(X311_ys,X311);  
        Darray(X311_ye,X311);  
        Darray(X311_zs,X311);  
        Darray(X311_ze,X311);  
    
        Darray(X311_w,X311); 
        Darray(X311_rho_c,X311); 
        Darray(X311_EA,X311); 
        Darray(X311_d,X311);  
        Darray(X311_l,X311); 
        Darray(X311_H,X311); 
        Darray(X311_P,X311); 
        Darray(X311_facT,X311);  
        
        Darray(X314_T,X311);  
        Darray(X315_t,X311);  
    }

    if(X312>0)
    {
        Darray(X311_xs,X312);
        Darray(X311_xe,X312);
        Darray(X311_ys,X312);
        Darray(X311_ye,X312);
        Darray(X311_zs,X312);
        Darray(X311_ze,X312);
        Darray(X312_k,X312);
        Darray(X312_T0,X312);
        
        Darray(X314_T,X312);  
        Darray(X315_t,X312);  
    }

    if(X320>0)
    {
        Iarray(X320_type,X320);
    }

    if(X321>0)
    {
        Darray(X321_Sn,X321);
        Darray(X321_d,X321);
        Darray(X321_lambda,X321);
        Darray(X321_dk,X321);
        Darray(X321_rho,X321);
        Darray(X321_nd,X321);
        Darray(X321_nl,X321);

        Darray(X322_D,X321);
        Darray(X322_L,X321);
        Darray(X322_x0,X321);
        Darray(X322_y0,X321);
        Darray(X322_z0,X321);
        Darray(X322_phi,X321);
        Darray(X322_theta,X321);
        Darray(X322_psi,X321);
    }

    if(X324>0)
    {
    Darray(X324_x,X324);
    Darray(X324_y,X324);
    Darray(X324_z,X324);
    }
    
    if(Z11>0)
    {
        Darray(Z11_x,Z11);  
        Darray(Z11_y,Z11);  
        Darray(Z11_z,Z11);  
        Darray(Z11_l,Z11);  
        Darray(Z11_w,Z11);  
        Darray(Z11_t,Z11);  
        Darray(Z11_rho,Z11);  
        Darray(Z11_e,Z11);  
        Darray(Z11_ix,Z11);  
        Darray(Z11_iy,Z11);  
        Darray(Z11_iz,Z11);  
        Darray(Z11_nu,Z11);  
        Darray(Z11_n,Z11);  
    }

// --------------------------

    for(n=0;n<A581;++n)
    {
    A581_xs[n] = dctrl[dd++];
    A581_xe[n] = dctrl[dd++];
    A581_ys[n] = dctrl[dd++];
    A581_ye[n] = dctrl[dd++];
    A581_zs[n] = dctrl[dd++];
    A581_ze[n] = dctrl[dd++];
    }

    for(n=0;n<A583;++n)
    {
    A583_xc[n] = dctrl[dd++];
    A583_zc[n] = dctrl[dd++];
    A583_ys[n] = dctrl[dd++];
    A583_ye[n] = dctrl[dd++];
    A583_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<A584;++n)
    {
    A584_xc[n] = dctrl[dd++];
    A584_yc[n] = dctrl[dd++];
    A584_zs[n] = dctrl[dd++];
    A584_ze[n] = dctrl[dd++];
    A584_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<A585;++n)
    {
    A585_xm1[n] = dctrl[dd++];
    A585_ym1[n] = dctrl[dd++];
    A585_zm1[n] = dctrl[dd++];
    A585_r1[n] = dctrl[dd++];
    A585_xm2[n] = dctrl[dd++];
    A585_ym2[n] = dctrl[dd++];
    A585_zm2[n] = dctrl[dd++];
    A585_r2[n] = dctrl[dd++];
    }
    
    for(n=0;n<A586;++n)
    {
    A586_xm[n] = dctrl[dd++];
    A586_ym[n] = dctrl[dd++];
    A586_zm[n] = dctrl[dd++];
    A586_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<A587;++n)
    {
    A587_xs[n] = dctrl[dd++];
    A587_xe[n] = dctrl[dd++];
    A587_ys[n] = dctrl[dd++];
    A587_ye[n] = dctrl[dd++];
    A587_zs[n] = dctrl[dd++];
    A587_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<A588;++n)
    {
    A588_xs[n] = dctrl[dd++];
    A588_xe[n] = dctrl[dd++];
    A588_ys[n] = dctrl[dd++];
    A588_ye[n] = dctrl[dd++];
    A588_zs[n] = dctrl[dd++];
    A588_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<A589;++n)
    {
    A589_xs[n] = dctrl[dd++];
    A589_xe[n] = dctrl[dd++];
    A589_ys[n] = dctrl[dd++];
    A589_ye[n] = dctrl[dd++];
    A589_zs[n] = dctrl[dd++];
    A589_ze[n] = dctrl[dd++];
    }

    for(n=0;n<B71;++n)
    {
    B71_val[n]= dctrl[dd++];
    B71_dist[n]= dctrl[dd++];
    B71_b[n]  = dctrl[dd++];
    B71_x[n]  = dctrl[dd++];
    B71_y[n]  = dctrl[dd++];
    }
        
    for(n=0;n<B106;++n)
    {
    B106_b[n]  = dctrl[dd++];
    B106_x[n]  = dctrl[dd++];
    B106_y[n]  = dctrl[dd++];
    }
    
    for(n=0;n<B107;++n)
    {
    B107_xs[n]  = dctrl[dd++];
    B107_xe[n]  = dctrl[dd++];
    B107_ys[n]  = dctrl[dd++];
    B107_ye[n]  = dctrl[dd++];
    B107_d[n]  = dctrl[dd++];
    }
    
    for(n=0;n<B108;++n)
    {
    B108_xs[n]  = dctrl[dd++];
    B108_xe[n]  = dctrl[dd++];
    B108_ys[n]  = dctrl[dd++];
    B108_ye[n]  = dctrl[dd++];
    B108_d[n]  = dctrl[dd++];
    }
    
    for(n=0;n<B240;++n)
    {
    B240_C[n]  = dctrl[dd++];
    B240_D[n]  = dctrl[dd++];
    B240_xs[n] = dctrl[dd++];
    B240_xe[n] = dctrl[dd++];
    B240_ys[n] = dctrl[dd++];
    B240_ye[n] = dctrl[dd++];
    B240_zs[n] = dctrl[dd++];
    B240_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<B270;++n)
    {
    B270_xs[n] = dctrl[dd++];
    B270_xe[n] = dctrl[dd++];
    B270_ys[n] = dctrl[dd++];
    B270_ye[n] = dctrl[dd++];
    B270_zs[n] = dctrl[dd++];
    B270_ze[n] = dctrl[dd++];
    B270_n[n]  = dctrl[dd++];
    B270_d50[n]= dctrl[dd++];
    B270_alpha[n]= dctrl[dd++];
    B270_beta[n]= dctrl[dd++];
    }
    
    for(n=0;n<B274;++n)
    {
    B274_xc[n] = dctrl[dd++];
    B274_yc[n] = dctrl[dd++];
    B274_zs[n] = dctrl[dd++];
    B274_ze[n] = dctrl[dd++];
    B274_r[n] = dctrl[dd++];
    B274_n[n]  = dctrl[dd++];
    B274_d50[n]= dctrl[dd++];
    B274_alpha[n]= dctrl[dd++];
    B274_beta[n]= dctrl[dd++];
    }
    
    for(n=0;n<B281;++n)
    {
    B281_xs[n] = dctrl[dd++];
    B281_xe[n] = dctrl[dd++];
    B281_ys[n] = dctrl[dd++];
    B281_ye[n] = dctrl[dd++];
    B281_zs[n] = dctrl[dd++];
    B281_ze[n] = dctrl[dd++];
    B281_n[n]  = dctrl[dd++];
    B281_d50[n]= dctrl[dd++];
    B281_alpha[n]= dctrl[dd++];
    B281_beta[n]= dctrl[dd++];
    }

    for(n=0;n<B282;++n)
    {
    B282_xs[n] = dctrl[dd++];
    B282_xe[n] = dctrl[dd++];
    B282_ys[n] = dctrl[dd++];
    B282_ye[n] = dctrl[dd++];
    B282_zs[n] = dctrl[dd++];
    B282_ze[n] = dctrl[dd++];
    B282_n[n]  = dctrl[dd++];
    B282_d50[n]= dctrl[dd++];
    B282_alpha[n]= dctrl[dd++];
    B282_beta[n]= dctrl[dd++];
    }
    
    for(n=0;n<B291;++n)
    {
    B291_xs[n] = dctrl[dd++];
    B291_xe[n] = dctrl[dd++];
    B291_ys[n] = dctrl[dd++];
    B291_ye[n] = dctrl[dd++];
    B291_zs[n] = dctrl[dd++];
    B291_ze[n] = dctrl[dd++];
    B291_d[n] = dctrl[dd++];
    B291_n[n]  = dctrl[dd++];
    B291_d50[n]= dctrl[dd++];
    B291_alpha[n]= dctrl[dd++];
    B291_beta[n]= dctrl[dd++];
    }
    
    for(n=0;n<B310;++n)
    {
    B310_xs[n] = dctrl[dd++];
    B310_xe[n] = dctrl[dd++];
    B310_ys[n] = dctrl[dd++];
    B310_ye[n] = dctrl[dd++];
    B310_zs[n] = dctrl[dd++];
    B310_ze[n] = dctrl[dd++];
    B310_N[n]  = dctrl[dd++];
    B310_D[n]  = dctrl[dd++];
    B310_Cd[n] = dctrl[dd++];
    }

    for(n=0;n<B321;++n)
    {
    B321_xs[n] = dctrl[dd++];
    B321_xe[n] = dctrl[dd++];
    B321_ys[n] = dctrl[dd++];
    B321_ye[n] = dctrl[dd++];
    B321_zs[n] = dctrl[dd++];
    B321_ze[n] = dctrl[dd++];
    B321_N[n]  = dctrl[dd++];
    B321_D[n]  = dctrl[dd++];
    B321_Cd[n] = dctrl[dd++];
    }

    for(n=0;n<B322;++n)
    {
    B322_xs[n] = dctrl[dd++];
    B322_xe[n] = dctrl[dd++];
    B322_ys[n] = dctrl[dd++];
    B322_ye[n] = dctrl[dd++];
    B322_zs[n] = dctrl[dd++];
    B322_ze[n] = dctrl[dd++];
    B322_N[n]  = dctrl[dd++];
    B322_D[n]  = dctrl[dd++];
    B322_Cd[n] = dctrl[dd++];
    }
    
    for(n=0;n<B411;++n)
    {
    B411_ID[n] = ictrl[ii++];
    B411_Q[n] = dctrl[dd++];
    }
    
    for(n=0;n<B412;++n)
    {
    B412_ID[n] = ictrl[ii++];
    B412_pressBC[n] = dctrl[dd++];
    }
    
    for(n=0;n<B413;++n)
    {
    B413_ID[n] = ictrl[ii++];
    B413_h[n] = dctrl[dd++];
    }
    
    for(n=0;n<B414;++n)
    {
    B414_ID[n] = ictrl[ii++];
    B414_Uio[n] = dctrl[dd++];
    }
    
    for(n=0;n<B415;++n)
    {
    B415_ID[n] = ictrl[ii++];
    B415_U[n] = dctrl[dd++];
    B415_V[n] = dctrl[dd++];
    B415_W[n] = dctrl[dd++];
    }
    
    for(n=0;n<B416;++n)
    {
    B416_ID[n] = ictrl[ii++];
    B416_alpha[n] = dctrl[dd++];
    }
    
    for(n=0;n<B417;++n)
    {
    B417_ID[n] = ictrl[ii++];
    B417_Nx[n] = dctrl[dd++];
    B417_Ny[n] = dctrl[dd++];
    B417_Nz[n] = dctrl[dd++];
    }
    
    for(n=0;n<B418;++n)
    {
    B418_ID[n] = ictrl[ii++];
    B418_pio[n] = ictrl[ii++];
    }
    
    for(n=0;n<B421;++n)
    {
    B421_ID[n] = ictrl[ii++];
    B421_Q[n] = ictrl[ii++];
    }
    
    for(n=0;n<B422;++n)
    {
    B422_ID[n] = ictrl[ii++];
    B422_FSF[n] = ictrl[ii++];
    }
    
    for(n=0;n<B440;++n)
    {
    B440_ID[n]  = ictrl[ii++];
    B440_face[n]  = ictrl[ii++];
    B440_xs[n] = dctrl[dd++];
    B440_xe[n] = dctrl[dd++];
    B440_ys[n] = dctrl[dd++];
    B440_ye[n] = dctrl[dd++];
    }
    
    for(n=0;n<B441;++n)
    {
    B441_ID[n]  = ictrl[ii++];
    B441_face[n]  = ictrl[ii++];
    B441_xs[n] = dctrl[dd++];
    B441_xe[n] = dctrl[dd++];
    B441_ys[n] = dctrl[dd++];
    B441_ye[n] = dctrl[dd++];
    B441_zs[n] = dctrl[dd++];
    B441_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<B442;++n)
    {
    B442_ID[n]  = ictrl[ii++];
    B442_face[n]  = ictrl[ii++];
    B442_xm[n] = dctrl[dd++];
    B442_ym[n] = dctrl[dd++];
    B442_zm[n] = dctrl[dd++];
    B442_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<C75;++n)
    {
    C75_x[n] = dctrl[dd++];
    C75_z[n] = dctrl[dd++];
    C75_a[n] = dctrl[dd++];
    C75_s[n] = dctrl[dd++];
    C75_l[n] = dctrl[dd++];
    C75_v[n] = dctrl[dd++];
    }
    
    for(n=0;n<F70;++n)
    {
    F70_xs[n] = dctrl[dd++];
    F70_xe[n] = dctrl[dd++];
    F70_ys[n] = dctrl[dd++];
    F70_ye[n] = dctrl[dd++];
    F70_zs[n] = dctrl[dd++];
    F70_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<F71;++n)
    {
    F71_xs[n] = dctrl[dd++];
    F71_xe[n] = dctrl[dd++];
    F71_ys[n] = dctrl[dd++];
    F71_ye[n] = dctrl[dd++];
    F71_zs[n] = dctrl[dd++];
    F71_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<F72;++n)
    {
    F72_xs[n] = dctrl[dd++];
    F72_xe[n] = dctrl[dd++];
    F72_ys[n] = dctrl[dd++];
    F72_ye[n] = dctrl[dd++];
    F72_h[n] = dctrl[dd++];
    }

for(n=0;n<F369;++n)
    {
    F369_x[n] = dctrl[dd++];
    F369_z[n] = dctrl[dd++];
    F369_a[n] = dctrl[dd++];
    F369_s[n] = dctrl[dd++];
    F369_l[n] = dctrl[dd++];
    F369_v[n] = dctrl[dd++];
    }
    
    for(n=0;n<F370;++n)
    {
    F370_xs[n] = dctrl[dd++];
    F370_xe[n] = dctrl[dd++];
    F370_ys[n] = dctrl[dd++];
    F370_ye[n] = dctrl[dd++];
    F370_zs[n] = dctrl[dd++];
    F370_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<F371;++n)
    {
    F371_xs[n] = dctrl[dd++];
    F371_xe[n] = dctrl[dd++];
    F371_ys[n] = dctrl[dd++];
    F371_ye[n] = dctrl[dd++];
    F371_zs[n] = dctrl[dd++];
    F371_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<F374;++n)
    {
    F374_xc[n] = dctrl[dd++];
    F374_zc[n] = dctrl[dd++];
    F374_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<F375;++n)
    {
    F375_xc[n] = dctrl[dd++];
    F375_zc[n] = dctrl[dd++];
    F375_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<F378;++n)
    {
    F378_xc[n] = dctrl[dd++];
    F378_yc[n] = dctrl[dd++];
    F378_zc[n] = dctrl[dd++];
    F378_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<F379;++n)
    {
    F379_xc[n] = dctrl[dd++];
    F379_yc[n] = dctrl[dd++];
    F379_zc[n] = dctrl[dd++];
    F379_r[n] = dctrl[dd++];
    }
    
    
    for(n=0;n<F390;++n)
    {
    F390_xs[n] = dctrl[dd++];
    F390_xe[n] = dctrl[dd++];
    F390_ys[n] = dctrl[dd++];
    F390_ye[n] = dctrl[dd++];
    F390_zs[n] = dctrl[dd++];
    F390_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<F391;++n)
    {
    F391_xs[n] = dctrl[dd++];
    F391_xe[n] = dctrl[dd++];
    F391_ys[n] = dctrl[dd++];
    F391_ye[n] = dctrl[dd++];
    F391_zs[n] = dctrl[dd++];
    F391_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<F394;++n)
    {
    F394_xc[n] = dctrl[dd++];
    F394_zc[n] = dctrl[dd++];
    F394_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<F395;++n)
    {
    F395_xc[n] = dctrl[dd++];
    F395_zc[n] = dctrl[dd++];
    F395_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<F398;++n)
    {
    F398_xc[n] = dctrl[dd++];
    F398_yc[n] = dctrl[dd++];
    F398_zc[n] = dctrl[dd++];
    F398_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<F399;++n)
    {
    F399_xc[n] = dctrl[dd++];
    F399_yc[n] = dctrl[dd++];
    F399_zc[n] = dctrl[dd++];
    F399_r[n] = dctrl[dd++];
    }
    
    for(n=0;n<P35;++n)
    {
    P35_ts[n] = dctrl[dd++];
    P35_te[n] = dctrl[dd++];
    P35_dt[n] = dctrl[dd++];
    }
    
    for(n=0;n<P50;++n)
    {
    P50_x[n] = dctrl[dd++];
    P50_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P51;++n)
    {
    P51_x[n] = dctrl[dd++];
    P51_y[n] = dctrl[dd++];
    }

    for(n=0;n<P52;++n)
    {
    P52_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P56;++n)
    {
    P56_x[n] = dctrl[dd++];
    }

    for(n=0;n<P58;++n)
    {
    P58_x[n] = dctrl[dd++];
    P58_y[n] = dctrl[dd++];
    P58_T[n] = dctrl[dd++];
    }
    
    for(n=0;n<P61;++n)
    {
    P61_x[n] = dctrl[dd++];
    P61_y[n] = dctrl[dd++];
    P61_z[n] = dctrl[dd++];
    }
        
    for(n=0;n<P62;++n)
    {
    P62_xs[n] = dctrl[dd++];
    P62_xe[n] = dctrl[dd++];
    P62_ys[n] = dctrl[dd++];
    P62_ye[n] = dctrl[dd++];
    P62_zs[n] = dctrl[dd++];
    P62_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<P63;++n)
    {
    P63_x[n] = dctrl[dd++];
    P63_y[n] = dctrl[dd++];
    }

    for(n=0;n<P64;++n)
    {
    P64_x[n] = dctrl[dd++];
    P64_y[n] = dctrl[dd++];
    P64_z[n] = dctrl[dd++];
    }
    
    for(n=0;n<P65;++n)
    {
    P65_x[n] = dctrl[dd++];
    P65_y[n] = dctrl[dd++];
    P65_z[n] = dctrl[dd++];
    }

    for(n=0;n<P66;++n)
    {
    P66_x[n] = dctrl[dd++];
    P66_y[n] = dctrl[dd++];
    P66_z[n] = dctrl[dd++];
    }
    
    for(n=0;n<P81;++n)
    {
    P81_xs[n] = dctrl[dd++];
    P81_xe[n] = dctrl[dd++];
    P81_ys[n] = dctrl[dd++];
    P81_ye[n] = dctrl[dd++];
    P81_zs[n] = dctrl[dd++];
    P81_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<P85;++n)
    {
    P85_x[n] = dctrl[dd++];
    P85_y[n] = dctrl[dd++];
     P85_r[n] = dctrl[dd++];
     P85_cd[n] = dctrl[dd++];
     P85_cm[n] = dctrl[dd++];
    }
    
    for(n=0;n<P88;++n)
    {
    P88_x[n] = dctrl[dd++];
    P88_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P121;++n)
    {
    P121_x[n] = dctrl[dd++];
    P121_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P123;++n)
    {
    P123_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P124;++n)
    {
    P124_x[n] = dctrl[dd++];
    }
    
    for(n=0;n<P125;++n)
    {
    P125_x[n] = dctrl[dd++];
    P125_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P133;++n)
    {
    P133_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P134;++n)
    {
    P134_y[n] = dctrl[dd++];
    }

    for(n=0;n<P140;++n)
    {
    P140_x[n] = dctrl[dd++];
    P140_y[n] = dctrl[dd++];
    }

    for(n=0;n<P167;++n)
    {
    P167_x[n] = dctrl[dd++];
    }

    for(n=0;n<P168;++n)
    {
    P168_x[n] = dctrl[dd++];
    P168_zs[n] = dctrl[dd++];
    P168_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<P184;++n)
    {
    P184_its[n] = ictrl[ii++];
    P184_ite[n] = ictrl[ii++];
    P184_dit[n] = ictrl[ii++];
    }
    
    for(n=0;n<P185;++n)
    {
    P185_ts[n] = dctrl[dd++];
    P185_te[n] = dctrl[dd++];
    P185_dt[n] = dctrl[dd++];
    }

    for(n=0;n<P194;++n)
    {
    P194_its[n] = ictrl[ii++];
    P194_ite[n] = ictrl[ii++];
    P194_dit[n] = ictrl[ii++];
    }
    
    for(n=0;n<P195;++n)
    {
    P195_ts[n] = dctrl[dd++];
    P195_te[n] = dctrl[dd++];
    P195_dt[n] = dctrl[dd++];
    }
    
    for(n=0;n<P230;++n)
    {
    P230_x[n] = dctrl[dd++];
    }
    
    for(n=0;n<P240;++n)
    {
    P240_x[n] = dctrl[dd++];
    }
    
    for(n=0;n<P351;++n)
    {
    P351_x[n] = dctrl[dd++];
    P351_y[n] = dctrl[dd++];
    }
    
    for(n=0;n<P352;++n)
    {
    P352_x[n] = dctrl[dd++];
    P352_y[n] = dctrl[dd++];
    }

    for(n=0;n<Q61;++n)
    {
    Q61_x[n] = dctrl[dd++];
    Q61_y[n] = dctrl[dd++];
    Q61_z[n] = dctrl[dd++];
    Q61_i[n] = ictrl[ii++];
    }
    
    for(n=0;n<Q73;++n)
    {
    Q73_val[n]= dctrl[dd++];
    Q73_dist[n]= dctrl[dd++];
    Q73_b[n]  = dctrl[dd++];
    Q73_x[n]  = dctrl[dd++];
    Q73_y[n]  = dctrl[dd++];
    }

    for(n=0;n<Q110;++n)
    {
    Q110_xs[n] = dctrl[dd++];
    Q110_xe[n] = dctrl[dd++];
    Q110_ys[n] = dctrl[dd++];
    Q110_ye[n] = dctrl[dd++];
    Q110_zs[n] = dctrl[dd++];
    Q110_ze[n] = dctrl[dd++];
    }
    for(n=0;n<Q111;++n)
    {
    Q111_xs[n] = dctrl[dd++];
    Q111_xe[n] = dctrl[dd++];
    Q111_ys[n] = dctrl[dd++];
    Q111_ye[n] = dctrl[dd++];
    Q111_zs[n] = dctrl[dd++];
    Q111_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<S73;++n)
    {
    S73_val[n]= dctrl[dd++];
    S73_dist[n]= dctrl[dd++];
    S73_b[n]  = dctrl[dd++];
    S73_x[n]  = dctrl[dd++];
    S73_y[n]  = dctrl[dd++];
    }
    
    for(n=0;n<W41;++n)
    {
    W41_xc[n] = dctrl[dd++];
    W41_yc[n] = dctrl[dd++];
    W41_zs[n] = dctrl[dd++];
    W41_ze[n] = dctrl[dd++];
    W41_vel[n] = dctrl[dd++];
    W41_beta[n] = dctrl[dd++];
    }
    
    for(n=0;n<X110;++n)
    {
    X110_xs[n] = dctrl[dd++];
    X110_xe[n] = dctrl[dd++];
    X110_ys[n] = dctrl[dd++];
    X110_ye[n] = dctrl[dd++];
    X110_zs[n] = dctrl[dd++];
    X110_ze[n] = dctrl[dd++];
    }
    
    for(n=0;n<X163;++n)
    {
    X163_x1[n] = dctrl[dd++];
    X163_y1[n] = dctrl[dd++];
    X163_z1[n] = dctrl[dd++];
    X163_x2[n] = dctrl[dd++];
    X163_y2[n] = dctrl[dd++];
    X163_z2[n] = dctrl[dd++];
    X163_x3[n] = dctrl[dd++];
    X163_y3[n] = dctrl[dd++];
    X163_z3[n] = dctrl[dd++];
    X163_x4[n] = dctrl[dd++];
    X163_y4[n] = dctrl[dd++];
    X163_z4[n] = dctrl[dd++];
    X163_x5[n] = dctrl[dd++];
    X163_y5[n] = dctrl[dd++];
    X163_z5[n] = dctrl[dd++];
    X163_x6[n] = dctrl[dd++];
    X163_y6[n] = dctrl[dd++];
    X163_z6[n] = dctrl[dd++];
    }
    
    for(n=0;n<X164;++n)
    {
    X164_x1[n] = dctrl[dd++];
    X164_y1[n] = dctrl[dd++];
    X164_z1[n] = dctrl[dd++];
    X164_x2[n] = dctrl[dd++];
    X164_y2[n] = dctrl[dd++];
    X164_z2[n] = dctrl[dd++];
    X164_x3[n] = dctrl[dd++];
    X164_y3[n] = dctrl[dd++];
    X164_z3[n] = dctrl[dd++];
    X164_x4[n] = dctrl[dd++];
    X164_y4[n] = dctrl[dd++];
    X164_z4[n] = dctrl[dd++];
    X164_x5[n] = dctrl[dd++];
    X164_y5[n] = dctrl[dd++];
    X164_z5[n] = dctrl[dd++];
    X164_x6[n] = dctrl[dd++];
    X164_y6[n] = dctrl[dd++];
    X164_z6[n] = dctrl[dd++];
    X164_x7[n] = dctrl[dd++];
    X164_y7[n] = dctrl[dd++];
    X164_z7[n] = dctrl[dd++];
    X164_x8[n] = dctrl[dd++];
    X164_y8[n] = dctrl[dd++];
    X164_z8[n] = dctrl[dd++];
    }
    
    for(n=0;n<X311;++n)
    {
    X311_xs[n] = dctrl[dd++];
    X311_xe[n] = dctrl[dd++];
    X311_ys[n] = dctrl[dd++];
    X311_ye[n] = dctrl[dd++];
    X311_zs[n] = dctrl[dd++];
    X311_ze[n] = dctrl[dd++];
    X311_w[n] = dctrl[dd++];
    X311_rho_c[n] = dctrl[dd++];
    X311_EA[n] = dctrl[dd++];
    X311_d[n] = dctrl[dd++];    
    X311_l[n] = dctrl[dd++]; 
    X311_H[n] = dctrl[dd++];
    X311_P[n] = dctrl[dd++];
    X311_facT[n] = dctrl[dd++];
    }

    for(n=0;n<X312;++n)
    {
    X311_xs[n] = dctrl[dd++];
    X311_xe[n] = dctrl[dd++];
    X311_ys[n] = dctrl[dd++];
    X311_ye[n] = dctrl[dd++];
    X311_zs[n] = dctrl[dd++];
    X311_ze[n] = dctrl[dd++];
    X312_k[n] = dctrl[dd++];
    X312_T0[n] = dctrl[dd++];
    }

    if (X314 > 0)
    {
        for(n=0;n<X311;++n)
        {
            X314_T[n] = dctrl[dd++]; 
        }
        for(n=0;n<X312;++n)
        {
            X314_T[n] = dctrl[dd++]; 
        }
    }
    if (X315 > 0)
    {
        for(n=0;n<X311;++n)
        {
            X315_t[n] = dctrl[dd++];
        }
        for(n=0;n<X312;++n)
        {
            X315_t[n] = dctrl[dd++]; 
        }
    }

    for(n=0;n<X320;++n)
    {
    X320_type[n] = ictrl[ii++];
    }

    for(n=0;n<X321;++n)
    {
    X321_Sn[n] = dctrl[dd++];
    X321_d[n] = dctrl[dd++];
    X321_lambda[n] = dctrl[dd++];
    X321_dk[n] = dctrl[dd++];
    X321_rho[n] = dctrl[dd++];
    X321_nd[n] = dctrl[dd++];
    X321_nl[n] = dctrl[dd++];

    X322_D[n] = dctrl[dd++];
    X322_L[n] = dctrl[dd++];
    X322_x0[n] = dctrl[dd++];
    X322_y0[n] = dctrl[dd++];
    X322_z0[n] = dctrl[dd++];
    X322_phi[n] = dctrl[dd++];
    X322_theta[n] = dctrl[dd++];
    X322_psi[n] = dctrl[dd++];
    }

    for(n=0;n<X324;++n)
    {
    X324_x[n] = dctrl[dd++];
    X324_y[n] = dctrl[dd++];
    X324_z[n] = dctrl[dd++];
    }

    for(n=0;n<Z11;++n)
    {
    Z11_x[n] = dctrl[dd++];
    Z11_y[n] = dctrl[dd++];
    Z11_z[n] = dctrl[dd++];
    Z11_l[n] = dctrl[dd++];
    Z11_w[n] = dctrl[dd++];
    Z11_t[n] = dctrl[dd++];
    Z11_rho[n] = dctrl[dd++];
    Z11_e[n] = dctrl[dd++];
    Z11_ix[n] = dctrl[dd++];
    Z11_iy[n] = dctrl[dd++];
    Z11_iz[n] = dctrl[dd++];
    Z11_nu[n] = dctrl[dd++];
    Z11_n[n] = dctrl[dd++];
    }
}
