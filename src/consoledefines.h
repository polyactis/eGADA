/*=================================================================
 * consoledefines.h 
 *
 *=================================================================*/
/*
 This File is part of eGADA

 eGADA: enhanced Genome Alteration Detection Algorithm

 eGADA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 eGADA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with eGADA.  If not, see <http://www.gnu.org/licenses/>.

 Authors:
 Yu S. Huang  polyactis@gmail.com
 Roger Pique-Regi    piquereg@usc.edu

 */
#ifndef _CONSOLEDEFINES_H_
#define _CONSOLEDEFINES_H_

//Matlab printing and memory allocation 
#include <stdio.h>
#include <malloc.h>

#define myPrintf printf //Syntax of myPrintf similar to that of fprintf('asfd',var,var)
//memory allocation in Matlab
//#define myDoubleMAlloc(mysize) mxGetPr(mxCreateDoubleMatrix(1,mysize,mxREAL))
//#define myIntMalloc(mysize)    mxGetPr(mxCreateNumericMatrix(1,mysize,mxUINT32_CLASS,mxREAL))
#define myCalloc(n,size) calloc(n, size)
#define myFree(ptr) free(ptr)
#define log2(x) log(x)/log(2)


//#define _DEBUG_BEwTandMinLen_
//#define _DebugBEwTscore_

#endif //_CONSOLEDEFINES_H_
