package MathTimeLib

import "base:intrinsics"
import "core:fmt"
import "core:math"
import la "core:math/linalg"

is_float :: intrinsics.type_is_float
is_int :: intrinsics.type_is_integer

pi :: 3.14159265358979323846
twopi :: 2.0 * pi
infinite :: 999999.9
undefined :: 999999.1

edirection :: enum {
	eTo,
	eFrom,
}

round :: proc(x: $T) -> T where is_float(T) {
	return math.floor(x + 0.5)
}

cot :: proc(x: $T) -> T where is_float(T) {
	temp := math.tan(x)
	if (math.abs(temp) < 0.00000001) {
		return infinite
	} else {
		return 1.0 / temp
	}
}

acosh :: proc(x: $T) -> T where is_float(T) {
	temp: f64
	if (x * x - 1. < 0.) {
		temp = undefined
		fmt.eprintln("ERROR: arccosh function error")
	} else {
		temp = math.ln(x + math.sqrt(x * x - 1.))
	}
	return temp
}

asinh :: proc(x: $T) -> T where is_float(T) {
	return math.ln(x + math.sqrt(x * x + 1.))
}

dot :: proc(x, y: [$N]$T) -> T where is_float(T) {
	return la.dot(x, y)
}

mag :: proc(x: [$N]$T) -> T where is_float(T) {
	return la.vector_length(x)
}

cross :: proc(vec1, vec2: [3]$T) -> [3]T where is_float(T) {
	out: [3]f64
	out[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]
	out[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]
	out[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0]
	return out
	// return la.cross(vec1,vec2)
}

norm :: proc(vec: [$N]$T) -> [N]T where is_float(T) {
	small: f64 : 0.00001
	magv: f64 = mag(vec)
	out: [3]f64
	if magv > small {
		out /= magv
	} else {
		out = {0., 0., 0.}
	}
	return out
	// return norm
}

angle :: proc(vec1, vec2: [$N]$T) -> T where is_float(T) {
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


sgn :: proc(x: $T) -> T where is_float(T) {
	if x < 0.0 {
		return -1.0
	} else {
		return 1.0
	}
}

rotmat :: proc(angle: $T, axis: $I) -> matrix[3, 3]f64 where is_int(I),
	is_float(T) {
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

rotvec :: proc(vec: [3]$T, angle: T, axis: $I) -> [3]f64 where is_int(I),
	is_float(T) {
	// replaces roti
	return rotmat(angle, axis) * v
}

// addvec, addvec3: adds 2 or 3 vectors with coefficients

// matvecmult: post-multiplies a 3x3 matrix with a 3x1 matrix

vecouter :: proc(vec1, vec2: [3]f64) -> matrix[3, 3]f64 {
	return la.outer_product(vec1, vec2)
}

repeat_vec :: proc(vec: ^[$N]$T, val: T) where intrinsics.type_is_numeric(T) {
	for i in 0 ..< len(vec) {
		vec[i] = val
	}
}


ludecomp :: proc(mat: matrix[$N, N]$T) -> (matrix[N, N]T, [N]int, int) where is_float(T) {
	small: T = cast(T)(0.000001)
	order: int = N

	scale: [N]T
	index: [N]int
	lu: matrix[N, N]T = mat

	imax: int = 0
	for i in 0 ..< order {
		amax: T = 0.0
		for j in 0 ..< order {
			if math.abs(lu[i, j]) > amax {
				amax = math.abs(lu[i, j])
			}
		}
		if math.abs(amax) < small {
			fmt.eprintln("ERROR: Singular matrix while finding LU decomposition")
		}
		scale[i] = 1.0 / amax
	}

	for j in 0 ..< order {
		for i in 0 ..< j {
			sum: T = lu[i, j]
			for k in 0 ..< i {
				sum -= lu[i, k] * lu[k, j]
			}
			lu[i, j] = sum
		}

		amax: T = 0.0
		for i in j ..< order {
			sum: T = lu[i, j]
			for k in 0 ..< j {
				sum -= lu[i, k] * lu[k, j]
			}
			lu[i, j] = sum
			dum: T = scale[i] * math.abs(sum)
			if dum >= amax {
				imax = i
				amax = dum
			}
		}

		if j != imax {
			for k in 0 ..< order {
				dum: T = lu[imax, k]
				lu[imax, k] = lu[j, k]
				lu[j, k] = dum
			}
			scale[imax] = scale[j]
		}
		index[j] = imax

		if math.abs(lu[j, j]) < small {
			fmt.eprintln("ERROR: Singular matrix while finding LU decomposition")
		}

		if j != order - 1 {
			dum: T = 1.0 / lu[j, j]
			for i in (j + 1) ..< order {
				lu[i, j] *= dum
			}
		}
	}

	return lu, index, order
}

lubksub :: proc(lu: matrix[$N, N]$T, b: [N]T, index: [N]int) -> [N]T where is_float(T) {
	i, j, iptr, ii: int
	sum: T
	sol: [N]T = b
	order := N
	ii = -1
	for i = 0; i < order; i += 1 {
		iptr = index[i]
		sum = sol[iptr]
		sol[iptr] = sol[i]
		if (ii != -1) {
			for j = ii; j <= i - 1; j += 1 {
				sum -= lu[i, j] * sol[j]
			}
		} else {
			if sum != 0. {
				ii = i
			}
		}
		sol[i] = sum
	}

	sol[order - 1] = sol[order - 1] / lu[order - 1, order - 1]

	for i = order - 2; i >= 0; i -= 1 {
		sum = sol[i]
		for j = i + 1; j < order; j += 1 {
			sum -= lu[i, j] * sol[j]
		}
		sol[i] = sum / lu[i, i]
	}

	return sol
}


matinverse :: proc(mat: matrix[$N, N]$T) -> matrix[N, N]T where is_float(T) {
	i, j: int
	if N <= 4 {
		return la.inverse(mat)
	}
	matinv := mat
	lu, index, order := ludecomp(mat)
	b: [N]T;repeat_vec(&b, T(order))
	for j = 0; j < order; j += 1 {
		for i = 0; i < order; i += 1 {
			if i == j {
				b[i] = 1.
			} else {
				b[i] = 0.
			}
		}
		b = lubksub(lu, b, index)
		for i = 0; i < order; i += 1 {
			matinv[i, j] = b[i]
		}
	}
	return matinv
}


determinant :: proc(mat: matrix[$N, N]$T) -> T where is_float(T) {
	return la.determinant(mat)
}

cholesky :: proc(mat: matrix[$N, N]$T) -> matrix[N, N]T where is_float(T) {
	// gives lower chol, transpose to get upper (matlab)
	res: matrix[N, N]T
	n := N
	for r := 0; r < n; r += 1 {
		for c := 0; c <= r; c += 1 {
			if c == r {
				sum: T = 0.
				for j := 0; j < c; j += 1 {
					sum += res[c, j] * res[c, j]
				}
				res[c, c] = math.sqrt(mat[c, c] - sum)
			} else {
				sum: T = 0
				for j := 0; j < c; j += 1 {
					sum += res[r, j] * res[c, j]
				}
				res[r, c] = 1. / res[c, c] * (mat[r, c] - sum)
			}
		}
	}
	return res
}

RootsOptions :: enum {
	All, // all roots
	Real, // real roots only
	Unique, // unique real root only
}

// TODO: check
quadratic :: proc(
	a, b, c: $T,
	opt: RootsOptions,
) -> (
	real1, imag1, real2, imag2: T,
) where is_float(T) {
	small := 0.0000001
	real1, imag1, real2, imag2: T
	discrim: T = b * b - 4. * a * c

	if math.abs(discrim) < small {
		real1 = -b / (2. * a)
		real2 = real1
		if opt == RootsOptions.Unique {
			real2 = 99999.9
		}

	} else {
		if math.abs(a) < small {
			real1 = -c / b
		} else {
			if opt == RootsOptions.All {
				real1 = -b / (2. * a)
				real2 = real1
				imag1 = sqrt(-discrim) / (2. * a)
				imag2 = -imag1
			} else {
				real1 = 99999.9
				real2 = 99999.9
			}
		}
	}
	return real1, imag1, real2, imag2
}

// TODO: check
cubicspl :: proc(p1, p2, p3, p4: $T) -> (T, T, T, T) where is_float(T) {
	acu0: T = p2
	acu1: T = -p1 / 3. - 0.5 * p2 + p3 - p4 / 6.
	acu2: T = 0.5 * p1 - p2 + 0.5 * p3
	acu3: T = -p1 / 6. + 0.5 * p2 - 0.5 * p3 + p4 / 6.
}

// TODO: check 
cubic :: proc(
	a, b, c, d: $T,
	opt: RootsOptions,
) -> (
	real1, imag1, real2, imag2, real3, imag3: T,
) where is_float(T) {
	rad := 57.29577951308230
	onethird := 1. / 3.
	small = 0.00000001
	temp1, temp2, p, q, r, delta, e0, cosphi, sinphi, phi: T
	real1, real2, real3: T
	imag1, imag2, imag3: T

	if math.abs(a) > small {
		// force coefs into std form
		p = b / a
		q = c / a
		r = d / a
		a = onethird * (3. * q - p * p)
		b = (1. / 27.) * (2. * p * p * p - 9. * p * q + 27 * r)

		delta = (a * a * a / 27.) + (b * b * 0.25)

		// cardan's formula
		if (delta > small) {
			temp1 = (-b * 0.5) + sqrt(delta)
			temp2 = (-b * 0.5) - sqrt(delta)
			temp1 = sgn(temp1) * math.pow(math.abs(temp1), onethird)
			temp2 = sgn(temp2) * math.pow(math.abs(temp2), onethird)
			real1 = temp1 + temp2 - p * onethird

			if opt == RootsOptions.All {
				real2 = -0.5 * (temp1 + temp2) - p * onethird
				imag2 = -0.5 * sqrt(3.) * (temp1 - temp2)
				real3 = -0.5 * (temp1 + temp2) - p * onethird
				imag3 = -imag2
			} else {
				real2 = 99999.9
				real3 = 99999.9
			}
		} else {
			if math.abs(delta) < small {
				real1 = -2. * sgn(b) * math.pow(math.abs(b * 0.5), onethird) - p * onethird
				real2 = sgn(b) * math.pow(math.abs(b * 0.5), onethird) - p * onethird
				if opt == RootsOptions.Unique {
					real3 = 99999.9
				} else {
					real3 = real2
				}
			} else {
				e0 = 2. * sqrt(-a * onethird)
				cosphi = (-b / (2. * math.sqrt(-a * a * a / 27)))
				sinphi = math.sqrt(1 - cosphi * cosphi)
				phi = math.atan2(sinphi, cosphi)
				if phi < 0. {
					phi += 2. * math.pi
				}
				real1 = e0 * math.cos(phi * onethird) - p * onethird
				real2 = e0 * math.cos(phi * onethird + 120. / rad) - p * onethird
				real3 = e0 * math.cos(phi * onethird + 240. / rad) - p * onethird
			}

		}
	} else {
		real1, imag1, real2, imag2 = quadratic(b, c, d, opt)
		real3 = 99999.9
		imag3 = 99999.9
	}
	return real1, imag1, real2, imag2, real3, imag3
}


cubicinterp :: proc()