package MathTimeLib

import "base:intrinsics"
import "core:fmt"
import "core:math"
import la "core:math/linalg"

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

repeat_vec :: proc(vec: ^[$N]$T, val: T) where intrinsics.type_is_numeric(T) {
	for i in 0 ..< len(vec) {
		vec[i] = val
	}
}


ludecomp :: proc(
	mat: matrix[$N, N]$T,
) -> (
	matrix[N, N]T,
	[N]int,
	int,
) where intrinsics.type_is_float(T) {
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

lubksub :: proc(
	lu: matrix[$N, N]$T,
	b: [N]T,
	index: [N]int,
) -> [N]T where intrinsics.type_is_float(T) {
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


matinverse :: proc(mat: matrix[$N, N]$T) -> matrix[N, N]T where intrinsics.type_is_float(T) {
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


determinant :: proc(mat: matrix[$N, N]$T) -> T where intrinsics.type_is_float(T) {
	return la.determinant(mat)
}

cholesky :: proc(mat: matrix[$N, N]$T) -> matrix[N, N]T where intrinsics.type_is_float(T) {
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

quadratic :: proc(
	a, b, c: $T,
) -> (
	real1, comp1, real2, comp2: T,
) where intrinsics.type_is_float(T) {
	small := 0.0000001
	real1, comp1, real2, comp2: T
	discrim: T = b * b - 4. * a * c

    if math.abs(discrim) < small {
        real1 = -b/(2.*a);
        real2 =real1;
        // TODO: complete this
    }
}
