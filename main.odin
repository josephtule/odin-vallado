package main

import "base:intrinsics"
import "core:fmt"
import "core:math"
import la "core:math/linalg"
import "core:strconv"
import "core:strings"
import "core:unicode/utf8"

import asl "AstroLib"
import mtl "MathTimeLib"
import sgp "SGP4"


main :: proc() {
	mat: matrix[3, 3]f64 = {1, 2, 1, 1, 2, 3, 3, 3, 3}
	ans, index, order := mtl.ludecomp(mat)
	fmt.println(ans)
	b: [3]f64 = {6, 14, 18}
	sol := mtl.lubksub(ans, b, index)
	fmt.println(sol)
	minv := mtl.matinverse(mat)
	fmt.println(minv)
	mat2: matrix[3, 3]f64 = {
		1.9042,
		1.0286,
		1.2999,
		1.0286,
		0.8214,
		0.9853,
		1.2999,
		0.9853,
		1.1919,
	}
	fmt.println(mtl.cholesky(mat2))

	fmt.println(mtl.factorial(4))
	fmt.println(u8(mtl.Month.April) == 4)
	fmt.println(uint(mtl.Month.March))
	r: [3]f64 = {0., 0., 0.}
	fmt.println(testproc(3., &r))
	fmt.println(r)


	strin := "323'42"
	builder := strings.builder_make()
	strings.write_string(&builder, strin)
	builder.buf[2] = '_'
	strout := strings.to_string(builder)
	fmt.println(strout)
	fmt.println(strconv.parse_int(strout[0:]))
	// strfin := transmute(string)str
	// fmt.println(strconv.parse_int(strfin[0:2]))


	longstr1_in := "1 25544  98067A   08264.51782528 -.00002182 -23000-4 -11606-4 0  2927"
	// longstr1_in := "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
	longstr2_in := "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"

	longstr1_mut := strings.builder_make()
	strings.write_string(&longstr1_mut, longstr1_in)
	longstr2_mut := strings.builder_make()
	strings.write_string(&longstr2_mut, longstr2_in)

	satrec: sgp.elsetrec
	ok: bool

	ibexp, nexp, j: int
	nndotind: int = 44
	bstarind: int = 53
	for j = 10; j <= 15; j += 1 {
		if longstr1_mut.buf[j] == ' ' {
			longstr1_mut.buf[j] = '_'
		}
	}
	if longstr1_mut.buf[44] != ' ' {
		longstr1_mut.buf[43] = longstr1_mut.buf[44]
		nndotind = 43
	}
	longstr1_mut.buf[44] = '.'
	if (longstr1_mut.buf[7] == ' ') {
		longstr1_mut.buf[7] = 'U'
	}
	if (longstr1_mut.buf[9] == ' ') {
		longstr1_mut.buf[9] = '.'
	}
	for j = 45; j <= 49; j += 1 {
		if (longstr1_mut.buf[j] == ' ') {
			longstr1_mut.buf[j] = '0'
		}
	}
	if (longstr1_mut.buf[51] == ' ') {
		longstr1_mut.buf[51] = '0'
	}
	if (longstr1_mut.buf[53] != ' ') {
		longstr1_mut.buf[52] = longstr1_mut.buf[53]
		bstarind = 52
	}
	longstr1_mut.buf[53] = '.'
	longstr2_mut.buf[25] = '.'
	for j = 26; j <= 32; j += 1 {
		if (longstr2_mut.buf[j] == ' ') {
			longstr2_mut.buf[j] = '0'
		}
	}
	if (longstr1_mut.buf[62] == ' ') {
		longstr1_mut.buf[62] = '0'
	}
	if (longstr1_mut.buf[68] == ' ') {
		longstr1_mut.buf[68] = '0'
	}

	longstr1 := strings.to_string(longstr1_mut)
	fmt.println("original string")
	fmt.println(longstr1_in)
	fmt.println("modified string")
	fmt.println(longstr1)

	cardnumb, _ := strconv.parse_int(longstr1[:2])
	satrec.satnum, _ = strconv.parse_int(longstr1[2:7])
	if longstr1[7] == 'U' {
		satrec.classification = .unclassified
	} else {
		satrec.classification = .classified
	}
	satrec.intldesg, _ = strconv.parse_int(longstr1[9:14])
	satrec.intldesg_piece = strings.clone(longstr1[14:17])
	satrec.epochyr, _ = strconv.parse_int(longstr1[18:20])
	satrec.epochdays, _ = strconv.parse_f64(longstr1[20:32])
	satrec.ndot, _ = strconv.parse_f64(longstr1[33:43])
	satrec.nddot, _ = strconv.parse_f64(longstr1[nndotind:50])
	nexp, _ = strconv.parse_int(longstr1[50:52])
	satrec.bstar, _ = strconv.parse_f64(longstr1[bstarind:59])
	ibexp , _ = strconv.parse_int(longstr1[59:61])
	satrec.ephtype, _ = strconv.parse_int(longstr1[61:64])
	elnumind := 64
	if longstr1[64] == ' ' {
		elnumind += 1
	}
	satrec.elnum, ok = strconv.parse_i64(longstr1[elnumind:69])

	



	fmt.println("This is the end")
	fmt.println(longstr1[64:])
	fmt.println(satrec.elnum)
	fmt.println(ok)

}

testproc :: proc(a: f64, r: ^[3]f64) -> f64 {
	a := 4.
	a += a
	r[2] = 5.3
	return a
}
