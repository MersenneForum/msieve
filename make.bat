cl /Ox newset.c asmmpp.c 64aux.c int64.c basic.c pentium.c
cl /Ox /I..\mpir\include extend30.c cofactorize.c cofactorize_siqs.c cofactorize_squfof.c asmmpp.c 64aux.c int64.c basic.c pentium.c clockticks.c ws2_32.lib ..\mpir\lib32\mpir.lib
del *.obj
	
