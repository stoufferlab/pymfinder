/************************************************************************
*
*  File name: bits.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __BITS_H
#define __BITS_H

#include <stdio.h>
#include <stdlib.h>

#ifdef UNIX
typedef long long           int64;
typedef unsigned long long  uint64;
#else
typedef __int64           int64;
typedef unsigned __int64 uint64;
#endif
typedef int               int32;
typedef unsigned int     uint32;
typedef short             int16;
typedef unsigned short   uint16;
typedef char              int8;
typedef unsigned char    uint8;

typedef uint32           uint;

#endif
