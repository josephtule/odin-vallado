package MathTimeLib

import "core:fmt"
import "core:math"
import la "core:math/linalg"
import "base:intrinsics"

pi :: 3.14159265358979323846
twopi :: 2.0 * pi
infinite :: 999999.9
undefined :: 999999.1

edirection :: enum {
	eTo,
	eFrom,
}

round :: proc(x: f64) -> f64 {
	return math.floor(x + 0.5)
}

cot :: proc(x: f64) -> f64 {
	temp := math.tan(x)
	if (math.abs(temp) < 0.00000001) {
		return infinite
	} else {
		return 1.0 / temp
	}
}

acosh :: proc(x: f64) -> f64 {
	temp: f64
	if (x * x - 1. < 0.) {
		temp = undefined
		fmt.eprintln("ERROR: arccosh function error")
	} else {
		temp = math.ln(x + math.sqrt(x * x - 1.))
	}
	return temp
}

asinh :: proc(x: f64) -> f64 {
	return math.ln(x + math.sqrt(x * x + 1.))
}

dot :: proc(x, y: [3]f64) -> f64 {
	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
}

mag :: proc(x: [3]f64) -> f64 {
	return math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
}

cross :: proc(vec1, vec2: [3]f64) -> [3]f64 {
	out: [3]f64
	out[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]
	out[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]
	out[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0]
	return out
}

norm :: proc(vec: [3]f64) -> [3]f64 {
	small: f64 : 0.00001
	magv: f64 = mag(vec)
	out: [3]f64
	if magv > small {
		out /= magv
	} else {
		out = {0., 0., 0.}
	}
	return out
}

angle :: proc(vec1, vec2: [3]f64) -> f64 {
	small := 0.00000001
	magv1 := mag(vec1)
	magv2 := mag(vec2)
	magv1v2 := magv1 * magv2
	temp: f64
	if (magv1v2 > small * small) {
		temp = dot(vec1, vec2) / (magv1v2)
		if (math.abs(temp) > 1.0) {temp = sgn(temp) * 1.0}
		return math.acos(temp)
	} else {
		return undefined
	}
}


sgn :: proc(x: f64) -> f64 {
	if x < 0.0 {
		return -1.0
	} else {
		return 1.0
	}
}

rotmat :: proc(angle: f64, axis: $I) -> matrix[3, 3]f64 where intrinsics.type_is_integer(I) {
	// replaces rotimat
	R: matrix[3, 3]f64
	c := math.cos(angle)
	s := math.sin(angle)
	switch axis {
	// odinfmt: disable
    case 1:
        R = {1., 0., 0.,
             0., c, s,
             0., -s, c}
    case 2:
        R = {c, 0., -s,
             0., 1., 0.,
             s, 0., c}
	case 3:
        R = {c, s, 0.,
             -s, c, 0.,
             0., 0., 1.}
	// odinfmt: enable
	case:
		fmt.eprintln("ERROR: invalid rotation axis")
	}
	return R
}

rotvec :: proc(vec: [3]f64, angle: f64, axis: $I) -> [3]f64 where intrinsics.type_is_integer(I) {
	// replaces roti
	return rotmat(angle, axis) * v
}

// addvec, addvec3: adds 2 or 3 vectors with coefficients

// matvecmult: post-multiplies a 3x3 matrix with a 3x1 matrix

vecouter :: proc(vec1, vec2: [3]f64) -> matrix[3, 3]f64 {
	return la.outer_product(vec1, vec2)
}


