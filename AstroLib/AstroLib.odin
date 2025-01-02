package AstroLib

import "core:fmt"
import "core:math"

eOpt :: enum {
	e80,
	e96,
	e00a,
	e00b,
	e00cio,
}

jpldesize :: 60000

jpldedata :: struct {
	rsun, rmoon:       [3]f64,
	year, mon, day:    int,
	rsmag, rmmag, mjd: f64,
}


equintype :: enum {
	ea_m,
	en_m,
	ep_m,
	ea_nu,
	en_nu,
	ep_nu,
}


mum :: 3.986004415e14
mu :: 398600.4415
re :: 6378.1363
velkmps :: 7.905366149846074
earthrot :: 7.292115e-05
speedoflight :: 2.99792458e8
