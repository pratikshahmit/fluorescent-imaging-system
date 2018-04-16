function sub2ind(rows, cols, row, col) {
	return (row * cols + col) * 4;
}

function canvas2matrix(data, cols, rows) {
	var im = new Object();
	im.rows = rows;
	im.cols = cols;
	im.pixels = [[],[],[]]
	for (var row = 0; row < im.rows; row += 1) {
		im.pixels[0].push([]);
		im.pixels[1].push([]);
		im.pixels[2].push([]);
		for (var col = 0; col < im.cols; col += 1) {
			idx = sub2ind(im.rows, im.cols, row, col);
			// ignore alpha
			im.pixels[0][row].push(data[idx]);
			im.pixels[1][row].push(data[idx + 1]);
			im.pixels[2][row].push(data[idx + 2]);
		}
	}
	return im;
}

function matrix2imgData(im, canvas) {
	imgData = canvas.createImageData(im.cols, im.rows);
	for (var row = 0; row < im.rows; row += 1) {
		for (var col = 0; col < im.cols; col += 1) {
			idx = sub2ind(im.rows, im.cols, row, col);
			// alpha defaults to full 255
			imgData.data[idx]     = im.pixels[0][row][col];
			imgData.data[idx + 1] = im.pixels[1][row][col];
			imgData.data[idx + 2] = im.pixels[2][row][col];
			imgData.data[idx + 3] = 255;
		}
	}
	return imgData;
}

function matrixNormalized2imgData(im, canvas) {
	imgData = canvas.createImageData(im.cols, im.rows);
	for (var row = 0; row < im.rows; row += 1) {
		for (var col = 0; col < im.cols; col += 1) {
			idx = sub2ind(im.rows, im.cols, row, col);
			// alpha defaults to full 255
			imgData.data[idx]     = im.pixels[0][row][col] * 255;
			imgData.data[idx + 1] = im.pixels[1][row][col] * 255;
			imgData.data[idx + 2] = im.pixels[2][row][col] * 255;
			imgData.data[idx + 3] = 255;
		}
	}
	return imgData;
}

function labf(t) {
	if (t > Math.pow(6/29, 3)) {
		return Math.pow(t, 1/3);
	} else {
		return (1/3) * (Math.pow(29/6, 2) * t) + (4/29);
	}
}

function rgb2lab(im) {
	var nIm = new Object();
	nIm.rows = im.rows;
	nIm.cols = im.cols;
	nIm.pixels = [[],[],[]]
	for (var row = 0; row < im.rows; row += 1) {
		nIm.pixels[0].push([]);
		nIm.pixels[1].push([]);
		nIm.pixels[2].push([]);
		for (var col = 0; col < im.cols; col += 1) {
			r = im.pixels[0][row][col] / 255;
			g = im.pixels[1][row][col] / 255;
			b = im.pixels[2][row][col] / 255;
			// rgb to sRGB to CIEXYZ
			one_over_b21 = 1 / .17697;
			x = one_over_b21 * (.49000 * r + .31000 * g + .20000 * b);
			y = one_over_b21 * (.17697 * r + .81240 * g + .01063 * b);
			z = one_over_b21 * (     0 * r + .01000 * g + .99000 * b);
			// CIEXYZ to CIELAB
			Xn = 95.047;
			Yn = 100;
			Zn = 108.883;
			nIm.pixels[0][row].push(116 * labf(y/Yn) - 16);
			nIm.pixels[1][row].push(500 * (labf(x/Xn) - labf(y/Yn)));
			nIm.pixels[2][row].push(200 * (labf(y/Yn) - labf(z/Zn)));
		}
	}
	return nIm;
}


function rgb2hsv(im) {
	var nIm = new Object();
	nIm.rows = im.rows;
	nIm.cols = im.cols;
	nIm.pixels = [[],[],[]]
	for (var row = 0; row < im.rows; row += 1) {
		nIm.pixels[0].push([]);
		nIm.pixels[1].push([]);
		nIm.pixels[2].push([]);
		for (var col = 0; col < im.cols; col += 1) {
			r = im.pixels[0][row][col] / 255;
			g = im.pixels[1][row][col] / 255;
			b = im.pixels[2][row][col] / 255;
			cmax = Math.max(r, Math.max(g, b));
			cmin = Math.min(r, Math.min(g, b));
			delta = cmax - cmin;
			h = 0;
			if (cmax == r) {
				h = 60 * (((g-b)/delta) % 6);
			} else if (cmax == g) {
				h = 60 * (((b-r)/delta) + 2);
			} else {
				h = 60 * (((r-g)/delta) + 4);
			}
			s = 0;
			if (cmax == 0) {
				s = 0;
			} else {
				s = delta / cmax;
			}
			v = cmax
			nIm.pixels[0][row].push(h);
			nIm.pixels[1][row].push(s);
			nIm.pixels[2][row].push(v);
		}
	}
	return nIm;
}


function rgb2hsv_retry(im) {
	var nIm = new Object();
	nIm.rows = im.rows;
	nIm.cols = im.cols;
	nIm.pixels = [[],[],[]]
	for (var row = 0; row < im.rows; row += 1) {
		nIm.pixels[0].push([]);
		nIm.pixels[1].push([]);
		nIm.pixels[2].push([]);
		for (var col = 0; col < im.cols; col += 1) {
			var answer = rgb2hsv_internet(im.pixels[0][row][col], im.pixels[1][row][col], im.pixels[2][row][col])
			nIm.pixels[0][row].push(answer.h);
			nIm.pixels[1][row].push(answer.s);
			nIm.pixels[2][row].push(answer.v);
		}
	}
	return nIm;
}

function rgb2hsv_internet () {
    var rr, gg, bb,
        r = arguments[0] / 255,
        g = arguments[1] / 255,
        b = arguments[2] / 255,
        h, s,
        v = Math.max(r, g, b),
        diff = v - Math.min(r, g, b),
        diffc = function(c){
            return (v - c) / 6 / diff + 1 / 2;
        };

    if (diff == 0) {
        h = s = 0;
    } else {
        s = diff / v;
        rr = diffc(r);
        gg = diffc(g);
        bb = diffc(b);

        if (r === v) {
            h = bb - gg;
        }else if (g === v) {
            h = (1 / 3) + rr - bb;
        }else if (b === v) {
            h = (2 / 3) + gg - rr;
        }
        if (h < 0) {
            h += 1;
        }else if (h > 1) {
            h -= 1;
        }
    }
    return {
        h: h,
        s: s,
        v: v
    };
}

function histeq(im, pln) {
	// construct cdf
	//    we need to know: cumulative pixels for each value, which pixels have each value
	//    construct normal histogram, then cumulate it
	// initialize normal histogram
	var hist = new Object();
	var values = new Array();
	// populate normal histogram
	for (var row = 0; row < im.rows; row += 1) {
		for (var col = 0; col < im.cols; col += 1) {
			v = im.pixels[pln][row][col];
			if (!(v in hist)) {
				hist[v] = 0;
				values.push(v)
			}
			hist[v] = hist[v] + 1;
		}
	}
	// get the values in the correct order
	values.sort(function(a,b){return a - b});
	// create cumulative histogram
	var cdf = new Object();
	for (var i = 0; i < values.length; i++) {
		var prev = 0;
		if (i == 0) {
			prev = 0;
		} else {
			prev = cdf[values[i-1]];
		}
		cdf[values[i]] = prev + hist[values[i]];
	}
	// equalize!
	min = hist[values[0]];
	numpixels = im.rows * im.cols;
	var equalized = new Object();
	for (var i = 0; i < values.length; i++) {
		equalized[values[i]] = ((cdf[values[i]] - min) / (numpixels - min));
	}

	// substitute equalized values in the image
	for (var row = 0; row < im.rows; row += 1) {
		for (var col = 0; col < im.cols; col += 1) {
			v = im.pixels[pln][row][col];
			im.pixels[pln][row][col] = equalized[v];
		}
	}
	return im;
}

function inverselabf(t) {
	if (t > 6/29) {
		return Math.pow(t, 3);
	} else {
		return 3 * Math.pow(6/29, 2) * (t - (4/29));
	}
}

function lab2rgb(im) {
	for (var row = 0; row < im.rows; row += 1) {
		for (var col = 0; col < im.cols; col += 1) {
			l = im.pixels[0][row][col];
			a = im.pixels[1][row][col];
			b = im.pixels[2][row][col];
			// CIELAB to CIEXYZ
			Xn = 95.047;
			Yn = 100;
			Zn = 108.883;
			x = Xn * inverselabf(1/116*(l + 16) + 1/500*a)
			y = Yn * inverselabf(1/116*(l + 16))
			z = Zn * inverselabf(1/116*(l + 16) - 1/200*b)
			// CIEXYZ to sRGB to rgb
			im.pixels[0][row][col] = ( .4184700 * x - .1586600 * y - .082835 * z) * 255;
			im.pixels[1][row][col] = (-.0911690 * x + .2524300 * y + .015708 * z) * 255;
			im.pixels[2][row][col] = ( .0009209 * x - .0025498 * y + .178600 * z) * 255;
		}
	}
	return im;
}

function normalize(im) {

	var nIm = new Object();
	nIm.rows = im.rows;
	nIm.cols = im.cols;
	nIm.pixels = [[],[],[]]

	lMin = Number.POSITIVE_INFINITY;
	lMax = Number.NEGATIVE_INFINITY;
	aMin = Number.POSITIVE_INFINITY;
	aMax = Number.NEGATIVE_INFINITY;
	bMin = Number.POSITIVE_INFINITY;
	bMax = Number.NEGATIVE_INFINITY;

	for (var row = 0; row < im.rows; row += 1) {
		for (var col = 0; col < im.cols; col += 1) {
			if (im.pixels[0][row][col] < lMin) lMin = im.pixels[0][row][col];
			if (im.pixels[0][row][col] > lMax) lMax = im.pixels[0][row][col];
			if (im.pixels[1][row][col] < aMin) aMin = im.pixels[1][row][col];
			if (im.pixels[1][row][col] > aMax) aMax = im.pixels[1][row][col];
			if (im.pixels[2][row][col] < bMin) bMin = im.pixels[2][row][col];
			if (im.pixels[2][row][col] > bMax) bMax = im.pixels[2][row][col];
		}
	}

	for (var row = 0; row < im.rows; row += 1) {
		nIm.pixels[0].push([]);
		nIm.pixels[1].push([]);
		nIm.pixels[2].push([]);
		for (var col = 0; col < im.cols; col += 1) {
			nIm.pixels[0][row].push((im.pixels[0][row][col] - lMin) / (lMax - lMin));
			nIm.pixels[1][row].push((im.pixels[1][row][col] - aMin) / (aMax - aMin));
			nIm.pixels[2][row].push((im.pixels[2][row][col] - bMin) / (bMax - bMin));
		}
	}

	return nIm;
}

function segment(im, lMin, lMax, aMin, aMax, bMin, bMax) {
	labelMatrix = [];
	numTruePixels = 0;
	for (var row = 0; row < im.rows; row += 1) {
		labelMatrix.push([]);
		for (var col = 0; col < im.cols; col += 1) {
			if   ((im.pixels[0][row][col] >= lMin || im.pixels[0][row][col] <= lMax)
			   && (im.pixels[1][row][col] >= aMin && im.pixels[1][row][col] <= aMax)
			   && (im.pixels[2][row][col] >= bMin && im.pixels[2][row][col] <= bMax)) 
			{
				labelMatrix[row].push(true);
				numTruePixels++;
			} else {
				labelMatrix[row].push(false);
			}
		}
	}
	return {
		bw: labelMatrix,
		n: numTruePixels
	}
}

function segmentTeeth(im, lMin, lMax, aMin, aMax, bMin, bMax) {
	labelMatrix = [];
	numTruePixels = 0;
	for (var row = 0; row < im.rows; row += 1) {
		labelMatrix.push([]);
		for (var col = 0; col < im.cols; col += 1) {
			if   ((im.pixels[0][row][col] >= lMin && im.pixels[0][row][col] <= lMax)
			   && (im.pixels[1][row][col] >= aMin && im.pixels[1][row][col] <= aMax)
			   && (im.pixels[2][row][col] >= bMin && im.pixels[2][row][col] <= bMax)) 
			{
				labelMatrix[row].push(true);
				numTruePixels++;
			} else {
				labelMatrix[row].push(false);
			}
		}
	}
	return {
		bw: labelMatrix,
		n: numTruePixels
	}
}

function mask(im, bw) {
	for (var row = 0; row < im.rows; row += 1) {
		for (var col = 0; col < im.cols; col += 1) {
			if (bw[row][col] == true) {
				im.pixels[0][row][col] = 0;
				im.pixels[1][row][col] = 255;
				im.pixels[2][row][col] = 255;
			}
		}
	}
	return im;
}