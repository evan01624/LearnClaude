'use strict';

const canvas = document.getElementById('canvas');
const ctx    = canvas.getContext('2d');

// ── Fluid params ──────────────────────────────────────────────────────────────
const N     = 128;       // grid cells per axis
const ITER  = 4;         // Gauss-Seidel iterations (pressure projection)
const VISC  = 1e-6;      // near-zero viscosity (water)
const DT    = 0.016;     // simulation Δt (one 60fps frame)
const DAMP  = 0.998;     // velocity decay per frame
const FSCALE = 50;       // mouse force multiplier
const FRAD  = 5;         // mouse force radius (grid cells)

// ── Particle params ───────────────────────────────────────────────────────────
const NP      = 1500;
const MAXSPD  = 15;      // pixel speed mapped to full red
const SPRING_K = 0.008;  // restoring force strength (0=none, 0.02=snappy)

// ── Grid ──────────────────────────────────────────────────────────────────────
const SZ = (N + 2) * (N + 2);
const IX = (i, j) => i + (N + 2) * j;

let u  = new Float32Array(SZ), u0 = new Float32Array(SZ);
let v  = new Float32Array(SZ), v0 = new Float32Array(SZ);
let p  = new Float32Array(SZ), dv = new Float32Array(SZ);

// ── Fluid solver ──────────────────────────────────────────────────────────────
function setBnd(b, x) {
  for (let i = 1; i <= N; i++) {
    x[IX(0,   i)] = b === 1 ? -x[IX(1, i)] : x[IX(1, i)];
    x[IX(N+1, i)] = b === 1 ? -x[IX(N, i)] : x[IX(N, i)];
    x[IX(i,   0)] = b === 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, N+1)] = b === 2 ? -x[IX(i, N)] : x[IX(i, N)];
  }
  x[IX(0,   0)]   = 0.5 * (x[IX(1,   0)] + x[IX(0,   1)]);
  x[IX(N+1, 0)]   = 0.5 * (x[IX(N,   0)] + x[IX(N+1, 1)]);
  x[IX(0,   N+1)] = 0.5 * (x[IX(1, N+1)] + x[IX(0,   N)]);
  x[IX(N+1, N+1)] = 0.5 * (x[IX(N, N+1)] + x[IX(N+1, N)]);
}

function linSolve(b, x, x0, a, c) {
  const r = 1 / c;
  for (let k = 0; k < ITER; k++) {
    for (let j = 1; j <= N; j++) {
      for (let i = 1; i <= N; i++) {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (
          x[IX(i-1, j)] + x[IX(i+1, j)] +
          x[IX(i, j-1)] + x[IX(i, j+1)]
        )) * r;
      }
    }
    setBnd(b, x);
  }
}

function diffuse(b, x, x0, diff, dt) {
  const a = dt * diff * N * N;
  linSolve(b, x, x0, a, 1 + 4 * a);
}

function project(ux, vy, pArr, dvArr) {
  const h = 1 / N;
  for (let j = 1; j <= N; j++) {
    for (let i = 1; i <= N; i++) {
      dvArr[IX(i, j)] = -0.5 * h * (
        ux[IX(i+1, j)] - ux[IX(i-1, j)] +
        vy[IX(i, j+1)] - vy[IX(i, j-1)]
      );
      pArr[IX(i, j)] = 0;
    }
  }
  setBnd(0, dvArr); setBnd(0, pArr);
  linSolve(0, pArr, dvArr, 1, 4);
  for (let j = 1; j <= N; j++) {
    for (let i = 1; i <= N; i++) {
      ux[IX(i, j)] -= 0.5 * (pArr[IX(i+1, j)] - pArr[IX(i-1, j)]) / h;
      vy[IX(i, j)] -= 0.5 * (pArr[IX(i, j+1)] - pArr[IX(i, j-1)]) / h;
    }
  }
  setBnd(1, ux); setBnd(2, vy);
}

function advect(b, d, d0, ux, vy, dt) {
  const dt0 = dt * N;
  for (let j = 1; j <= N; j++) {
    for (let i = 1; i <= N; i++) {
      const x = Math.max(0.5, Math.min(N + 0.5, i - dt0 * ux[IX(i, j)]));
      const y = Math.max(0.5, Math.min(N + 0.5, j - dt0 * vy[IX(i, j)]));
      const i0 = Math.floor(x), i1 = i0 + 1;
      const j0 = Math.floor(y), j1 = j0 + 1;
      const s1 = x - i0, s0 = 1 - s1;
      const t1 = y - j0, t0 = 1 - t1;
      d[IX(i, j)] =
        s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
        s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    }
  }
  setBnd(b, d);
}

function fluidStep() {
  // Diffuse
  [u, u0] = [u0, u]; diffuse(1, u, u0, VISC, DT);
  [v, v0] = [v0, v]; diffuse(2, v, v0, VISC, DT);
  project(u, v, p, dv);

  // Advect
  [u, u0] = [u0, u]; [v, v0] = [v0, v];
  advect(1, u, u0, u0, v0, DT);
  advect(2, v, v0, u0, v0, DT);
  project(u, v, p, dv);

  for (let i = 0; i < SZ; i++) { u[i] *= DAMP; v[i] *= DAMP; }
}

// ── Particles ─────────────────────────────────────────────────────────────────
let W = 0, H = 0;
const px    = new Float32Array(NP);
const py    = new Float32Array(NP);
const restX = new Float32Array(NP);
const restY = new Float32Array(NP);

function initParticles() {
  for (let i = 0; i < NP; i++) {
    px[i] = restX[i] = Math.random() * W;
    py[i] = restY[i] = Math.random() * H;
  }
}

function sampleVel(sx, sy) {
  const gx = Math.max(0.5, Math.min(N + 0.5, 1 + sx / W * N));
  const gy = Math.max(0.5, Math.min(N + 0.5, 1 + sy / H * N));
  const i0 = Math.floor(gx), i1 = i0 + 1;
  const j0 = Math.floor(gy), j1 = j0 + 1;
  const s1 = gx - i0, s0 = 1 - s1;
  const t1 = gy - j0, t0 = 1 - t1;
  return [
    s0 * (t0 * u[IX(i0, j0)] + t1 * u[IX(i0, j1)]) + s1 * (t0 * u[IX(i1, j0)] + t1 * u[IX(i1, j1)]),
    s0 * (t0 * v[IX(i0, j0)] + t1 * v[IX(i0, j1)]) + s1 * (t0 * v[IX(i1, j0)] + t1 * v[IX(i1, j1)])
  ];
}

// ── Color ─────────────────────────────────────────────────────────────────────
function speedColor(spd) {
  const t   = Math.min(spd / MAXSPD, 1);
  const hue = (1 - t) * 240;       // 240° blue → 0° red
  const lit = 35 + t * 25;
  return `hsl(${hue | 0},100%,${lit | 0}%)`;
}

// ── Input ─────────────────────────────────────────────────────────────────────
let mx = 0, my = 0, lastMx = 0, lastMy = 0, pressing = false;

canvas.addEventListener('mousedown', e => { pressing = true;  mx = e.clientX; my = e.clientY; lastMx = mx; lastMy = my; });
canvas.addEventListener('mousemove', e => { mx = e.clientX;  my = e.clientY; });
canvas.addEventListener('mouseup',   ()  => { pressing = false; });
canvas.addEventListener('mouseleave',()  => { pressing = false; });

canvas.addEventListener('touchstart', e => {
  e.preventDefault();
  pressing = true;
  mx = e.touches[0].clientX; my = e.touches[0].clientY;
  lastMx = mx; lastMy = my;
}, { passive: false });
canvas.addEventListener('touchmove', e => {
  e.preventDefault();
  mx = e.touches[0].clientX; my = e.touches[0].clientY;
}, { passive: false });
canvas.addEventListener('touchend', e => { e.preventDefault(); pressing = false; }, { passive: false });

// ── Main loop ─────────────────────────────────────────────────────────────────
function loop() {
  if (pressing) {
    const dx = mx - lastMx;
    const dy = my - lastMy;
    const gi = Math.round(1 + mx / W * N);
    const gj = Math.round(1 + my / H * N);
    for (let dj = -FRAD; dj <= FRAD; dj++) {
      for (let di = -FRAD; di <= FRAD; di++) {
        const d = Math.sqrt(di * di + dj * dj);
        if (d > FRAD) continue;
        const ii = gi + di, jj = gj + dj;
        if (ii < 1 || ii > N || jj < 1 || jj > N) continue;
        const f = 1 - d / FRAD;
        u[IX(ii, jj)] += dx / W * FSCALE * f;
        v[IX(ii, jj)] += dy / H * FSCALE * f;
      }
    }
  }
  lastMx = mx; lastMy = my;

  fluidStep();

  ctx.fillStyle = 'rgba(2, 8, 30, 0.18)';
  ctx.fillRect(0, 0, W, H);

  for (let i = 0; i < NP; i++) {
    const [ux, vy] = sampleVel(px[i], py[i]);
    // Grid velocity → pixel displacement: vel * DT * domain_size
    const vx  = ux * DT * W + (restX[i] - px[i]) * SPRING_K;
    const vy2 = vy * DT * H + (restY[i] - py[i]) * SPRING_K;
    px[i] += vx;
    py[i] += vy2;

    // Clamp to screen bounds (fluid setBnd handles velocity reflection)
    if (px[i] < 0) px[i] = 0;
    else if (px[i] > W) px[i] = W;
    if (py[i] < 0) py[i] = 0;
    else if (py[i] > H) py[i] = H;

    const spd = Math.sqrt(vx * vx + vy2 * vy2);
    ctx.fillStyle = speedColor(spd);
    ctx.fillRect(px[i] - 1, py[i] - 1, 3, 3);
  }

  requestAnimationFrame(loop);
}

function resize() {
  W = canvas.width  = window.innerWidth;
  H = canvas.height = window.innerHeight;
  initParticles();
}

window.addEventListener('resize', resize);
resize();
loop();
