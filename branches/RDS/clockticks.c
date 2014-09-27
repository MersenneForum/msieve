#include <windows.h>
#include <stdio.h>
typedef __int64 BIGINT;

/* I believe this code is straight out of the Intel manual. I first
 * saw it in a comp.lang.c thread. It gives the number of processor
 * ticks since startup (for example, if you're on a 300MHz machine
 * like me, you'll get 300,000,000 ticks per second).
 *
 * (There's similar code available for Linux; look on the Web.)
 */

BIGINT prectime(void)
   {
   BIGINT t;
   unsigned int a,b;
   unsigned int *c = (unsigned int *)&t;
   _asm
      {
      _emit 0x0f;
      _emit 0x31;
      mov a,eax;
      mov b,edx;
      }

   c[0]=a;
   c[1]=b;
   return t;
   }

