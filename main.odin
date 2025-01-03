package main

import "base:intrinsics"
import "core:fmt"
import "core:math"
import la "core:math/linalg"

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
}
