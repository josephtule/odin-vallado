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
// :Math Functions -------------------------------------------------------------
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


ludecomp :: proc(
	mat: matrix[$N, N]$T,
) -> (
	matrix[N, N]T,
	[N]int,
	int,
) where is_float(T) {
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
) -> [N]T where is_float(T) {
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
					phi += 2. * pi
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


cubicinterp :: proc(
	a1, b1, c1, d1, a2, b2, c2, d2, vin: $T,
) where is_float(T) {
	kc0, kc1, kc2, kc3: T
	ac0, ac1, ac2, ac3: T
	real1, imag1: T
	real2, imag2: T
	real3, imag3: T
	vout: T

	ac0, ac1, ac2, ac3 = cubicspl(a1, b1, c1, d1)
	kc0, kc1, kc2, kc3 = cubicspl(a2, b2, c2, d2)
	real1, imag1, real2, imag2, real3, imag3 = cubic(
		kc3,
		kc2,
		kc1,
		kc0 - vin,
		RootsOptions.Real,
	)

	if real1 >= -0.000001 && real1 <= 1.001 {
		vout = real1
	} else {
		if real2 >= -0.000001 && real2 <= 1.001 {
			vout = real2
		} else {
			if real3 >= -0.000001 && real3 <= 1.001 {
				vout = real3
			} else {
				vout = 0.
				panic("ERROR: cubicinterp could not find root")
			}
		}
	}

	return ac3 * vout * vout * vout + ac2 * vout * vout + ac1 * vout + ac0
}

percentile :: proc(sequence: [$N]$T, pcnt: T) -> T where is_float(T) {
	// NOTE: this function does not do anything in Vallado's code 
	b: [100]T
	if N > 0 {
		if N == 1 {
			return b[0]
		} else {
			n := T((N - 1) * pcnt + 1.)
			if n == N {
				return b[N]
			} else {
				k: int = int(n)
				d: T = T(n - k)
				return b[k - 1] + d * (b[k] - b[k - 1])
			}
		}
	} else {return 0.}
}


factorial :: proc(n: int) -> f64 {
	if n == 0 {
		return 1.
	} else {
		return f64(n) * factorial(n - 1)
	}
}


// :Time Functions -------------------------------------------------------------

Month :: enum u8 {
	Janurary  = 1,
	Jan       = 1,
	February  = 2,
	Feb       = 2,
	March     = 3,
	Mar       = 3,
	April     = 4,
	Apr       = 4,
	May       = 5,
	June      = 6,
	Jun       = 6,
	July      = 7,
	Jul       = 7,
	August    = 8,
	Aug       = 8,
	September = 9,
	Sep       = 9,
	October   = 10,
	Oct       = 10,
	November  = 11,
	Nov       = 11,
	December  = 12,
	Dec       = 12,
}

getIntMonth :: proc(month: Month) -> uint {
	return uint(month)
}

Day :: enum u8 {
	Sunday    = 1,
	Sun       = 1,
	Monday    = 2,
	Mon       = 2,
	Tuesday   = 3,
	Tue       = 3,
	Wednesday = 4,
	Wed       = 4,
	Thursday  = 5,
	Thr       = 5,
	Thu       = 5,
	Friday    = 6,
	Fri       = 6,
	Saturday  = 7,
	Sat       = 7,
}

getIntDay :: proc(day: Day) -> uint {
	return uint(day)
}

getIntDayofweek :: proc(jd: f64) -> uint {
	jdtemp := math.floor(jd + 0.5)
	return uint(math.floor(jdtemp - 7. * math.floor((jdtemp + 1.) / 7) + 2))
}

jday :: proc(
	year, month, day, hr, minute: int,
	sec: f64,
) -> (
	jd, jd_frac: f64,
) {
	// convert because likely the other functions will read the following in as integers
	yearf := f64(year)
	monf := f64(month)
	dayf := f64(day)
	hrf := f64(hr)
	minf := f64(minute)

	jd =
		367.0 * yearf -
		math.floor((7. * (yearf + math.floor((monf + 9.) / 12.0))) * 0.25) +
		math.floor(275 * monf / 9.0) +
		dayf +
		1721013.5
	jd_frac = (sec + minf * 60. + hrf * 3600.) / 86400.

	if math.abs(jd_frac) >= 1. {
		dt: f64 = math.floor(jd_frac)
		jd = jd + dt
		jd_frac = jd_frac - dt
	}
	return jd, jd_frac
}

DateType :: enum {
	Juilian,
	Gregorian,
}
jdayall :: proc(
	year, month, day, hr, minute: int,
	sec: f64,
	type: DateType,
) -> (
	jd: f64,
) {
	yearf := f64(year)
	monf := f64(month)
	dayf := f64(day)
	hrf := f64(hr)
	minf := f64(minute)
	b: f64

	if monf <= 2. {
		yearf -= 1
		monf += 12
	}


	if type == .Juilian {
		b = 0.
	} else {
		b =
			2. - math.floor(yearf * 0.01) + math.floor(math.floor(yearf * 0.01) * 0.25)
		jd =
			math.floor(365.25 * (yearf + 4716)) +
			math.floor(30.6001 * (monf + 1.)) +
			dayf +
			b -
			1524.5 +
			((sec / 60. + minf) / 60. + hrf) / 24.
	}
	return jd
}

DaysInMonth :: enum u8 {
	Jan = 31,
	Feb = 28,
	Mar = 31,
	Apr = 30,
	May = 31,
	Jun = 30,
	Jul = 31,
	Aug = 31,
	Sep = 30,
	Oct = 31,
	Nov = 30,
	Dec = 31,
}

days2mdhms :: proc(
	year: int,
	days: f64,
) -> (
	month, day, hr, minute: int,
	sec: f64,
) {
	i, inttemp, dayofyr: int
	temp: f64
	lmonth: [13]int = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
	dayofyr = int(math.floor(days))

	if (year % 4) == 0 {
		lmonth[2] = 29
	}
	i = 1
	inttemp = 0
	for (dayofyr > inttemp + lmonth[i]) && (i < 12) {
		inttemp += lmonth[i]
		i += 1
	}
	month = i
	day = dayofyr - inttemp

	temp = (days - f64(dayofyr)) * 24.
	hr = int(temp)
	temp = (temp - f64(hr)) * 60.
	minute = int(temp)
	sec = (temp - f64(minute)) * 60.

	return month, day, hr, minute, sec
}


invjday :: proc(
	jd, jd_frac: f64,
) -> (
	year, month, day, hr, minute: int,
	sec: f64,
) {
	leapyears: int
	dt, days, tu, temp: f64
	jd := jd
	jd_frac := jd_frac
	if math.abs(jd_frac) >= 1. {
		jd += math.floor(jd_frac)
		jd_frac -= math.floor(jd_frac)
	}

	dt = jd - math.floor(jd) - 0.5
	if math.abs(dt) > 0.00000001 {
		jd -= dt
		jd_frac += dt
	}

	temp = jd - 2415019.5
	tu = temp / 365.25
	year = 1900 + int(tu)
	leapyears = int(f64(year - 1901) * 0.25)

	days = math.floor(temp - (f64(year - 1900) * 365. + f64(leapyears)))

	if days + jd_frac < 1. {
		year -= 1
		leapyears = int(f64(year - 1901) * 0.25)
		days = math.floor(temp - (f64(year - 1900) * 365. + f64(leapyears)))
	}

	month, day, hr, minute, sec = days2mdhms(year, days + jd_frac)
	return year, month, day, hr, minute, sec
}

finddays :: proc(year, month, day, hr, minute: int, sec: f64) -> (days: f64) {
	lmonth: [13]int = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
	i: int

	if (year - 1900) % 4 == 0 {
		lmonth[2] = 29
		if year % 100 == 0 && !(year % 400 == 0) {
			lmonth[2] = 28
		}
	}

	i = 1
	days = 0.
	for i < month && i < 12 {
		days += f64(lmonth[i])
		i += 1
	}
	days += f64(day) + f64(hr) / 24. + f64(minute) / 1440. + sec / 86400.
	return days
}


hms_sec :: proc(
	hour, minute: ^int,
	second, utsec: ^f64,
	direction: edirection,
) {
	temp: f64

	if direction == .eTo {
		utsec^ = f64(hour^) * 3600. + f64(minute^) * 60. + second^

	} else {
		temp = utsec^ / 3600.
		hour^ = int(temp)
		minute^ = int((temp - f64(hour^)) * 60.)
		second^ = (temp - f64(hour^) - f64(minute^) / 60.) * 3600.
	}
}

hms_ut :: proc(hr, min: ^int, sec, ut: ^f64, direction: edirection) {
	if direction == .eTo {
		ut^ = f64(hr^) * 100.0 + f64(min^) + sec^ * 0.01
	} else {
		hr^ = int(ut^ * 0.01)
		min^ = int(ut^ - f64(hr^) * 100.)
		sec^ = (ut^ - f64(hr^) * 100. - f64(min^)) * 100.
	}
}

hms_rad :: proc(hr, min: ^int, sec, hms: ^f64, direction: edirection) {
	rad2deg :: 57.29577951308230
	temp: f64

	temp = 15. / rad2deg
	if direction == .eTo {
		hms^ = f64(hr^) + f64(min^) / 60. + sec^ / 3600.
	} else {
		temp = hms^ / temp
		hr^ = int(temp)
		min^ = int((temp - f64(hr^)) * 60)
		sec^ = (temp - f64(hr^) - f64(min^) / 60.) * 3600.
	}
}

dms_rad :: proc(deg, min: ^int, sec, dms: ^f64, direction: edirection) {
	rad2deg :: 57.29577951308230
	temp: f64

	if direction == .eTo {
		dms^ = (f64(deg^) + f64(min^) / 60. + sec^ / 3600.) / rad2deg
	} else {
		temp = dms^ * rad2deg
		deg^ = int(temp)
		min^ = int((temp - f64(deg^)) * 60.)
		sec^ = (temp - f64(deg^) - f64(min^) / 60.) * 3600.
	}
}

jd2sse :: proc(jd: f64) -> f64 {
	return (jd - 2451544.5) * 86400.
}


convtime :: proc(
	year, mon, day, hr, min: int,
	sec: f64,
	timezone: int,
	dut1: f64,
	dat: int,
	ut1: ^f64,
	tut1: ^f64,
	jdut1: ^f64,
	jdut1Frac: ^f64,
	utc: ^f64,
	tai: ^f64,
	tt: ^f64,
	ttt: ^f64,
	jdtt: ^f64,
	jdttFrac: ^f64,
	tcg: ^f64,
	tdb: ^f64,
	ttdb: ^f64,
	jdtdb: ^f64,
	jdtdbFrac: ^f64,
	tcb: ^f64,
) {

	year := year
	mon := mon
	day := day
	hr := hr
	min := min
	sec := sec

	deg2rad, jd, jdFrac, sectemp, me: f64
	localhr, hrtemp, mintemp: int

	deg2rad = pi / 180.0

	// ------------------------  implementation   ------------------
	jd, jdFrac = jday(year, mon, day, 0, 0, 0.0)
	//     mjd  = jd - 2400000.5;
	//     mfme = hr*60.0 + min + sec/60.0;

	// ------------------ start if ( ut1 is known ------------------
	localhr = timezone + hr

	hms_sec(&localhr, &min, &sec, utc, .eTo)
	ut1^ = utc^ + dut1
	hms_sec(&hrtemp, &mintemp, &sectemp, ut1, .eFrom)
	jdut1^, jdut1Frac^ = jday(year, mon, day, hrtemp, mintemp, sectemp)
	tut1^ = (jdut1^ + jdut1Frac^ - 2451545.0) / 36525.0

	tai^ = utc^ + f64(dat)

	tt^ = tai^ + 32.184 // sec
	hms_sec(&hrtemp, &mintemp, &sectemp, tt, .eFrom)
	jdtt^, jdttFrac^ = jday(year, mon, day, hrtemp, mintemp, sectemp)
	ttt^ = (jdtt^ + jdttFrac^ - 2451545.0) / 36525.0

	tcg^ = tt^ + 6.969290134e-10 * (jdut1^ - 2443144.5) * 86400.0 // sec

	me = 357.5277233 + 35999.05034 * ttt^
	me = math.mod(me, 360.0)
	me = me * deg2rad
	tdb^ = tt^ + 0.001657 * math.sin(me) + 0.00001385 * math.sin(2.0 * me)
	hms_sec(&hrtemp, &mintemp, &sectemp, tdb, .eFrom)
	jdtdb^, jdtdbFrac^ = jday(year, mon, day, hrtemp, mintemp, sectemp)
	ttdb^ = (jdtdb^ + jdtdbFrac^ - 2451545.0) / 36525.0

	tcb^ = tdb^ + 1.55051976772e-8 * (jdtt^ - 2443144.5) * 86400.0 // sec
}
