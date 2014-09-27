PRDPPA(a,b,e)
unsigned long a,b,e[2];
{
unsigned long f,g;
_asm
   {
   mov eax,a		
   mov ebx,b	
   mul ebx		
   mov f,edx	
   mov g,eax	
   }
e[0] += g;
e[1] += (f + (e[0] < g));
}
