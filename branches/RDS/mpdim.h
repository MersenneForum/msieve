#define RADIX		1073741824

#define	BASE		(double)1073741824.0
#define SQRT_RAD	46340
#define MPDIM		150
#define MPCDIM		300
#define MPBITS		30
#define LOGTWO		50
#define EOS			'\0'
#define BLANK		' '
#define	SINGLE		2
#define	DOUBLE		3
#define	TRIPLE		4
#define	QUAD		5


#define COPY1(a,b,ind)	for (ind=0; ind<SIZE(a); ind++) b[ind] = a[ind]
#define SIZE(x)			x[0] 
#define LASTDIGIT(x)	x[1]
#define FIRST(x)		x[SIZE(x)-1]

#define INIT(x)		for (ini=0; ini<MPDIM; ini++) x[ini] = 0;
 
unsigned char numstr[MPCDIM];

#define small_mpcopy(a,b)				\
	{									\
	b[0] = a[0];						\
	b[1] = a[1];						\
	b[2] = a[2];						\
	b[3] = a[3];						\
	b[4] = a[4];						\
	b[5] = a[5];						\
    b[6] = a[6];                        \
    b[7] = a[7];                        \
	}
