#include "hsp3_64.as"
#regcmd "hsp3cmdinit","HSPFLUID64.dll", 1
//#cmd _HSPFLUID64_FLUIDPROCESS $01
//#cmd _HSPFLUID64_GETSPEEDPOINTS $02

#cmd FluidProcess $01
#cmd GetSpeedPoints $02
#cmd SetPoissonsLoopNum $03

/*
#module

#deffunc FluidProcess array vx,array vy,array p,array wall
	xl=length(vx)
	yl=length2(vx)
	_HSPFLUID64_FLUIDPROCESS vx,vy,p,wall,xl,yl
	return

#deffunc GetSpeedPoints array vx,array vy,array posx,array posy,array outspeedx,array outspeedy
	xl=length(vx)
	yl=length2(vx)
	lc=length(posx)
	_HSPFLUID64_GETSPEEDPOINTS vx,vy,posx,posy,outspeedx,outspeedy,xl,yl,lc
	return


#global
*/