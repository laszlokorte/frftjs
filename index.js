function frft(f, a) {
    const N = f.length;
    a = a % 4;
    
    // for integer values of a the fractional fourier transform is the same as the discrete fourier transform
    if (a === 0) return f.slice();
    if (a === 2) return f.slice().reverse();
    if (a === 1) return fft(f.slice());
    if (a === 3) return ifft(f.slice());
    
    let f0 = fftShift(f.slice());

    // reduce a into the range (0.5, 1.5)
    if (a > 2.0) {
        a = a - 2;
        f0 = f0.reverse();
    }
    if (a > 1.5) {
        a = a - 1;
        f0 = fftShift(fft(fftShift(f0)));
    }
    if (a < 0.5) {
        a = a + 1;
        f0 = fftShift(ifft(fftShift(f0)));
    }

    const alpha = a * Math.PI / 2;
    const s = Math.PI / (N + 1) / Math.sin(alpha) / 4;
    const t = Math.PI / (N + 1) * Math.tan(alpha / 2) / 4;
    const Cs = complScale(complExp(compl(0, -1 * (1 - a) * Math.PI / 4)), Math.sqrt(s / Math.PI));

    const sincArry = Array.from({length: 2 * N - 2}, (_, i) => -(2 * N - 3) + 2 * i).map(sinc).map(x => complScale(x, 1/2))

    const f1 = ifft(fconv(f0, sincArry)).reverse().slice(N, 2 * N - 1);
    
    const chrpA = Array.from({length: 2 * N - 1}, (_, i) => complExp(compl(0, -1 * t * Math.pow(-N + 1 + i, 2))));
    const l0 = arrayMod(chrpA, 2, 0)
    const l1 = arrayMod(chrpA, 2, 1)
    const f0m = complZipArrays(f0, l0, complMul);
    const f1m = complZipArrays(f1, l1, complMul);
    
    const chrpB = Array.from({length: 4 * N - 1}, (_, i) => complExp(compl(0, 1 * s * Math.pow(-(2 * N - 1) + i, 2))));
    const e1 = arrayMod(chrpB, 2, 0);
    const e0 = arrayMod(chrpB, 2, 1);
    const f0c = fconv(f0m, e0);
    const f1c = fconv(f1m, e1);
    const h0 = ifft(complZipArrays(f0c, f1c, complAdd));

    return fftShift(l0.map((l, i) => complMul(Cs, l, h0[N + i])).map(x => complScale(x, Math.sqrt(N))))
}

function complZipArrays(a,b, op) {
  return a.map((value, i) => op(value, b[i]))
}

function arrayMod(a, m, c) {
  return a.filter((_, i) => i % m === c)
}

function fconv(x, y, c) {
    const N = x.length + y.length - 1;
    const P = Math.pow(2, Math.ceil(Math.log2(N)));
    const fft_x = fft(padd(x, P))
    const fft_y = fft(padd(y, P))
    const z = complZipArrays(fft_x, fft_y, complMul)

    return z
}

function compl(re, im) {
  return [re, im]
}

function complExp(z) {
  return [
    Math.exp(z[0]) * Math.cos(z[1]),
    Math.exp(z[0]) * Math.sin(z[1]),
  ]
}

function complAdd(...args) {
  return args.reduce(function (a, b) {
    return [
      a[0] + b[0],
      a[1] + b[1],
    ]
  }, compl(0,0))
}

function complMul(...args) {
  return args.reduce(function(a, b) {
    return [
      a[0] * b[0] - a[1] * b[1],
      a[0] * b[1] + a[1] * b[0],
    ]
  }, compl(1, 0)) 
}

function complScale(z, s) {
  return [
    z[0]*s,
    z[1]*s,
  ]
}


function padd(a, p) {
  return [...a, ...Array(p - a.length).fill(null).map(_ => compl(0,0))]
}

function sinc(x) {
  return [Math.sin(Math.PI*x) / (Math.PI*x), 0]
}

const NEG_TWO_PI = -2*Math.PI;


function _fft(amplitudes) {
  const N = amplitudes.length

  if(N <= 1){
    return amplitudes
  }

  const hN = N / 2;

  const evenT = _fft(arrayMod(amplitudes, 2, 0))
  const oddT = _fft(arrayMod(amplitudes, 2, 1))


  for(let k = 0; k < hN; ++k){
    const p = k/N;

    const t = complMul(complExp(compl(0, NEG_TWO_PI*p)), oddT[k]);

    amplitudes[k] = [evenT[k][0] + t[0], evenT[k][1] + t[1]];
    amplitudes[k + hN] = [evenT[k][0] - t[0], evenT[k][1] - t[1]];
  }


  return amplitudes;
}

function fftShift(values) {
    return values.map((_,i) => values[(i + values.length/2) % values.length])
}


function fft(amplitudes) {
  const N = amplitudes.length;
  const scaleFactor = 1 / Math.sqrt(N);

  amplitudes = _fft(amplitudes)

  complScaleArray(amplitudes, scaleFactor)

  return amplitudes;
}

function ifft(amplitudes) {
  const N = amplitudes.length;
  const scaleFactor = 1 / Math.sqrt(N);

  complConjugateArray(amplitudes)

  amplitudes = _fft(amplitudes)

  complConjugateArray(amplitudes)
  complScaleArray(amplitudes, scaleFactor)

  return amplitudes;
}

function complConjugateArray(amplitudes) {
  const N = amplitudes.length

  for(let i = 0; i < N; ++i){
    amplitudes[i][1] = -amplitudes[i][1]
  }
}

function complScaleArray(amplitudes, scale) {
  const N = amplitudes.length

  for(let i = 0; i < N; ++i){
    amplitudes[i][0] = amplitudes[i][0] * scale
    amplitudes[i][1] = amplitudes[i][1] * scale
  }
}