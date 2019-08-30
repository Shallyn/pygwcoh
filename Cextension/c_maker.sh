#! /bin/bash

# ./g_maker.sh Type Object GObject G_TYPE_OBJECT

mylocation=`pwd`

name=$1
name_C=$(echo $name | tr [a-z] [A-Z])

fname_h=${name}.h
fname_c=${name}.c

echo -e "/**\n\
* Writer: Xiaolin.liu\n\
* xiaolin.liu@mail.bnu.edu.cn\n\
*\n\
* This module contains basic functions for  calculation.\n\
* Functions list:\n\
* Kernel: \n\
* 20xx.xx.xx, LOC\n\
**/\n\
\n\
#ifndef __INCLUDE_${name_C}__\n\
#define __INCLUDE_${name_C}__\n\
\n\
#endif\n\
" > ${fname_h}

echo -e "/**\n\
* Writer: Xiaolin.liu\n\
* xiaolin.liu@mail.bnu.edu.cn\n\
*\n\
* This module contains basic functions for  calculation.\n\
* Functions list:\n\
* Kernel: \n\
* 20xx.xx.xx, LOC\n\
**/\n\
\n\
#include \"${fname_h}\"\n\
" > ${fname_c}


