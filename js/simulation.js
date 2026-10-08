(() => {
'use strict';

// ---------- grid and constants ----------
const NX = 352, NY = 162, N = NX * NY;
const S = 13.2;               // cells per metre
const DX = 1 / S;
const GR = 4;                 // ground rows
const U0 = 0.1;               // inlet velocity (lattice units)
const RHO_AIR = 1.225, G0 = 9.80665, NU_AIR = 1.5e-5;
const WEIGHT_KGF = 300, AREA = 52;
const SPAN = [0, 10.4, 3.0, 0.6];   // effective span per body: wings, canard, fuselage
const CHORD = 2.5, CHORD_CELLS = CHORD * S;
const CX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const CY = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];
const WT = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
const OFF = CX.map((c, i) => c + CY[i] * NX);
const GW = CX.map((c, i) => 6 * WT[i] * c * U0);     // moving-wall term (ground)
const CSM = 18 * Math.SQRT2 * 0.16 * 0.16;           // Smagorinsky
const GROUND = 4;

// geometry in the aircraft frame (m): x aft, y up, origin at the middle of the wing cell
const WING_LE = -1.25, WING_Y = 0.75, WING_CAMBER = 0.045;
const CAN_LE = -8.0, CAN_C = 2.0, CPX = -7.0, CPY = 0.15, CAN_GAP = 0.75, CAN_CAMBER = 0.02;
const WHEEL_X = 0.0, WHEEL_Y = -1.3, WHEEL_R = 0.35;
const PROP_X = 1.6, PROP_Y = 0.2, PROP_R = 1.0;
const FUSE = [[-7.3, 0.02], [-7.3, 0.28], [-1.6, 0.5], [1.2, 0.5], [1.2, -0.5], [-1.6, -0.5]];
const OX = 11.5;

let cur = [], nxt = [];
for (let i = 0; i < 9; i++) { cur.push(new Float32Array(N)); nxt.push(new Float32Array(N)); }
const ux = new Float32Array(N), uy = new Float32Array(N), rho = new Float32Array(N);
const mask = new Uint8Array(N), newMask = new Uint8Array(N), flag = new Uint8Array(N);
const bfx = new Float32Array(N), bfy = new Float32Array(N);
const sponge = new Float32Array(N);
const E0 = new Float32Array(9);
const FA = new Float64Array(10);       // force accumulator per body (x,y)
const tmp = new Float32Array(9);

function feq(i, r, u, v) {
  const cu = CX[i] * u + CY[i] * v;
  return WT[i] * r * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * (u * u + v * v));
}
for (let i = 0; i < 9; i++) E0[i] = feq(i, 1, U0, 0);

for (let y = 0; y < NY; y++) for (let x = 0; x < NX; x++) {
  let t = 0.5;
  const dxo = x - (NX - 34); if (dxo > 0) t = Math.max(t, 0.5 + 0.7 * (dxo / 33) * (dxo / 33));
  const dyo = y - (NY - 12); if (dyo > 0) t = Math.max(t, 0.5 + 0.5 * (dyo / 11) * (dyo / 11));
  sponge[y * NX + x] = t;
}

// ---------- parameters ----------
const P = { v: 37, aoa: 6, canard: 0, alt: 3, prop: 60, re: 600, center: false, view: 'vel', smoke: true };
let tau0 = 0.52;
let ca = 1, sa = 0, cd = 1, sd = 0, OY = 5;
let running = true, dirty = true;
let simT = 0, simSteps = 0, settle = 0, groundShift = 0, propPhase = 0;
let stepsPF = 4, msAvg = 2;
const EMA = new Float64Array(10);
let hist = [], histAcc = 0;

// ---------- geometry ----------
function bodyW(x, y) { return [OX + x * ca + y * sa, OY - x * sa + y * ca]; }
function canW(x, y) {
  const dx = x - CPX, dy = y - CPY;
  return bodyW(CPX + dx * cd + dy * sd, CPY - dx * sd + dy * cd);
}
function section(xle, c, y0, camber, n) {
  const pts = [];
  for (let k = 0; k <= n; k++) { const s = k / n; pts.push([xle + s * c, y0 + camber * 4 * s * (1 - s) * c]); }
  return pts;
}
const WING_LO = section(WING_LE, CHORD, -WING_Y, WING_CAMBER, 14);
const WING_UP = section(WING_LE, CHORD, WING_Y, WING_CAMBER, 14);
const CAN_LO = section(CAN_LE, CAN_C, CPY - CAN_GAP, CAN_CAMBER, 10);
const CAN_UP = section(CAN_LE, CAN_C, CPY + CAN_GAP, CAN_CAMBER, 10);

function toGrid(w) { return [w[0] * S, GR + w[1] * S]; }
function segCells(a, b, r, cb) {
  const x0 = a[0], y0 = a[1], x1 = b[0], y1 = b[1];
  const xa = Math.max(1, Math.floor(Math.min(x0, x1) - r - 1)), xb = Math.min(NX - 2, Math.ceil(Math.max(x0, x1) + r + 1));
  const ya = Math.max(GR, Math.floor(Math.min(y0, y1) - r - 1)), yb = Math.min(NY - 2, Math.ceil(Math.max(y0, y1) + r + 1));
  const dx = x1 - x0, dy = y1 - y0, l2 = dx * dx + dy * dy || 1e-9, r2 = r * r;
  for (let y = ya; y <= yb; y++) for (let x = xa; x <= xb; x++) {
    const px = x + 0.5 - x0, py = y + 0.5 - y0;
    let t = (px * dx + py * dy) / l2; t = t < 0 ? 0 : t > 1 ? 1 : t;
    const ex = px - t * dx, ey = py - t * dy;
    if (ex * ex + ey * ey <= r2) cb(y * NX + x);
  }
}
function polyline(pts, tf, r, id) {
  let prev = toGrid(tf(pts[0][0], pts[0][1]));
  for (let k = 1; k < pts.length; k++) {
    const p = toGrid(tf(pts[k][0], pts[k][1]));
    segCells(prev, p, r, idx => { newMask[idx] = id; });
    prev = p;
  }
}
function polyFill(pts, tf, id) {
  const g = pts.map(p => toGrid(tf(p[0], p[1])));
  let xa = NX, xb = 0, ya = NY, yb = 0;
  for (const p of g) { xa = Math.min(xa, p[0]); xb = Math.max(xb, p[0]); ya = Math.min(ya, p[1]); yb = Math.max(yb, p[1]); }
  xa = Math.max(1, Math.floor(xa)); xb = Math.min(NX - 2, Math.ceil(xb));
  ya = Math.max(GR, Math.floor(ya)); yb = Math.min(NY - 2, Math.ceil(yb));
  for (let y = ya; y <= yb; y++) for (let x = xa; x <= xb; x++) {
    const px = x + 0.5, py = y + 0.5; let inside = false;
    for (let i = 0, j = g.length - 1; i < g.length; j = i++) {
      const xi = g[i][0], yi = g[i][1], xj = g[j][0], yj = g[j][1];
      if ((yi > py) !== (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi) + xi) inside = !inside;
    }
    if (inside) newMask[y * NX + x] = id;
  }
}

function rebuild() {
  const al = P.aoa * Math.PI / 180, de = P.canard * Math.PI / 180;
  ca = Math.cos(al); sa = Math.sin(al); cd = Math.cos(de); sd = Math.sin(de);
  OY = P.alt + WHEEL_R - (-WHEEL_X * sa + WHEEL_Y * ca);

  newMask.fill(0);
  for (let y = 0; y < GR; y++) newMask.fill(GROUND, y * NX, (y + 1) * NX);
  polyline(WING_LO, bodyW, 0.95, 1);
  polyline(WING_UP, bodyW, 0.95, 1);
  polyline(CAN_LO, canW, 0.95, 2);
  polyline(CAN_UP, canW, 0.95, 2);
  if (P.center) {
    polyFill(FUSE, bodyW, 3);
    polyFill([[-2.2, -0.62], [-2.2, 0.1], [-1.5, 0.1], [-1.5, -0.62]], bodyW, 3);
    polyline([[-1.85, 0.1], [-1.85, 1.1]], bodyW, 2.0, 3);
    polyFill([[0.3, 0.5], [0.3, 0.95], [1.1, 0.95], [1.1, 0.5]], bodyW, 3);
    const wc = toGrid(bodyW(WHEEL_X, WHEEL_Y));
    segCells(wc, wc, WHEEL_R * S, idx => { newMask[idx] = 3; });
  }
  for (let i = 0; i < N; i++) {
    const o = mask[i], n = newMask[i];
    if (n) { ux[i] = 0; uy[i] = 0; }
    else if (o) { for (let k = 0; k < 9; k++) cur[k][i] = WT[k]; ux[i] = 0; uy[i] = 0; rho[i] = 1; }
    mask[i] = n;
  }
  // cells next to a wall
  flag.fill(0);
  for (let y = 1; y < NY - 1; y++) for (let x = 1; x < NX - 1; x++) {
    const i = y * NX + x;
    if (mask[i]) continue;
    if (mask[i - 1] | mask[i + 1] | mask[i - NX] | mask[i + NX] | mask[i - NX - 1] | mask[i - NX + 1] | mask[i + NX - 1] | mask[i + NX + 1]) flag[i] = 1;
  }
  // propeller: a strip that pushes the air rearwards
  bfx.fill(0); bfy.fill(0);
  const k = P.prop / 100;
  if (k > 0) {
    const cells = [];
    const a = toGrid(bodyW(PROP_X, PROP_Y - PROP_R)), b = toGrid(bodyW(PROP_X, PROP_Y + PROP_R));
    segCells(a, b, 1.6, idx => { if (!mask[idx]) cells.push(idx); });
    const thrust = (2 * PROP_R * S) * U0 * U0 * (1 + k / 2) * k;
    const fc = cells.length ? thrust / cells.length : 0;
    for (const idx of cells) { bfx[idx] = fc * ca; bfy[idx] = -fc * sa; flag[idx] |= 2; }
  }
  settle = 2300;
  dirty = true;
}

// ---------- smoke ----------
const RAKES = Math.floor((NY - GR - 6) / 4), RK = 400, EMIT = 12;
const SX = new Float32Array(RAKES * RK), SY = new Float32Array(RAKES * RK);
const head = new Int32Array(RAKES);
let emitAcc = 0;
function rakeY(r) { return GR + 2.5 + r * 4; }
function seedSmoke() {
  SX.fill(NaN);
  for (let r = 0; r < RAKES; r++) {
    const y = rakeY(r); let h = 0;
    for (let x = NX - 3; x >= 1.5 && h < RK; x -= EMIT * U0) {
      const o = r * RK + h;
      SX[o] = mask[Math.floor(y) * NX + Math.floor(x)] ? NaN : x; SY[o] = y; h++;
    }
    head[r] = h % RK;
  }
  emitAcc = 0;
}
function advectSmoke(n) {
  for (let o = 0; o < SX.length; o++) {
    let x = SX[o];
    if (x !== x) continue;
    let y = SY[o];
    const fx = x - 0.5, fy = y - 0.5;
    const ix = Math.floor(fx), iy = Math.floor(fy);
    if (ix < 0 || ix > NX - 2 || iy < 0 || iy > NY - 2) { SX[o] = NaN; continue; }
    const tx = fx - ix, ty = fy - iy, i = iy * NX + ix;
    const u = (ux[i] * (1 - tx) + ux[i + 1] * tx) * (1 - ty) + (ux[i + NX] * (1 - tx) + ux[i + NX + 1] * tx) * ty;
    const v = (uy[i] * (1 - tx) + uy[i + 1] * tx) * (1 - ty) + (uy[i + NX] * (1 - tx) + uy[i + NX + 1] * tx) * ty;
    x += u * n; y += v * n;
    if (x < 0.5 || x > NX - 1.5 || y < GR || y > NY - 1.5 || mask[Math.floor(y) * NX + Math.floor(x)]) { SX[o] = NaN; continue; }
    SX[o] = x; SY[o] = y;
  }
  emitAcc += n;
  while (emitAcc >= EMIT) {
    emitAcc -= EMIT;
    for (let r = 0; r < RAKES; r++) {
      const o = r * RK + head[r];
      SX[o] = 1.2 + emitAcc * U0; SY[o] = rakeY(r);
      head[r] = (head[r] + 1) % RK;
    }
  }
}

// ---------- solver ----------
function initFlow() {
  for (let i = 0; i < N; i++) {
    const solid = mask[i];
    for (let k = 0; k < 9; k++) { cur[k][i] = solid ? WT[k] : E0[k]; nxt[k][i] = cur[k][i]; }
    ux[i] = solid ? 0 : U0; uy[i] = 0; rho[i] = 1;
  }
  simT = 0; simSteps = 0; settle = 2300; EMA.fill(0); hist = []; histAcc = 0;
  seedSmoke();
  dirty = true;
}

function step() {
  const a = cur, b = nxt;
  const a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7], a8 = a[8];
  const b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7], b8 = b[8];
  const t0 = tau0, t0sq = tau0 * tau0;
  for (let y = GR; y < NY - 1; y++) {
    let idx = y * NX + 1;
    for (let x = 1; x < NX - 1; x++, idx++) {
      if (mask[idx]) continue;
      const fl = flag[idx];
      let f0 = a0[idx], f1, f2, f3, f4, f5, f6, f7, f8;
      if (fl & 1) {
        for (let i = 1; i < 9; i++) {
          const n = idx - OFF[i], m = mask[n];
          if (m) {
            let v = a[OPP[i]][idx];
            if (m === GROUND) v += GW[i];
            else { FA[m * 2] -= 2 * v * CX[i]; FA[m * 2 + 1] -= 2 * v * CY[i]; }
            tmp[i] = v;
          } else tmp[i] = a[i][n];
        }
        f1 = tmp[1]; f2 = tmp[2]; f3 = tmp[3]; f4 = tmp[4]; f5 = tmp[5]; f6 = tmp[6]; f7 = tmp[7]; f8 = tmp[8];
      } else {
        f1 = a1[idx - 1]; f2 = a2[idx - NX]; f3 = a3[idx + 1]; f4 = a4[idx + NX];
        f5 = a5[idx - 1 - NX]; f6 = a6[idx + 1 - NX]; f7 = a7[idx + 1 + NX]; f8 = a8[idx - 1 + NX];
      }
      let r = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
      let u, v;
      if (r > 0.6 && r < 1.6) {
        u = (f1 - f3 + f5 - f6 - f7 + f8) / r;
        v = (f2 - f4 + f5 + f6 - f7 - f8) / r;
        const sp = u * u + v * v;
        if (sp > 0.1225) { const sc = 0.35 / Math.sqrt(sp); u *= sc; v *= sc; }
      } else {
        // degenerate cell: back to inlet equilibrium
        r = 1; u = U0; v = 0;
        f0 = E0[0]; f1 = E0[1]; f2 = E0[2]; f3 = E0[3]; f4 = E0[4]; f5 = E0[5]; f6 = E0[6]; f7 = E0[7]; f8 = E0[8];
      }
      rho[idx] = r; ux[idx] = u; uy[idx] = v;
      const usq = 1.5 * (u * u + v * v), r1 = r / 9, r2 = r / 36;
      const u3 = 3 * u, v3 = 3 * v, uu = 4.5 * u * u, vv = 4.5 * v * v;
      const p = u + v, q = v - u, pp = 4.5 * p * p, qq = 4.5 * q * q, p3 = 3 * p, q3 = 3 * q;
      const n0 = f0 - (4 / 9) * r * (1 - usq);
      const n1 = f1 - r1 * (1 + u3 + uu - usq);
      const n3 = f3 - r1 * (1 - u3 + uu - usq);
      const n2 = f2 - r1 * (1 + v3 + vv - usq);
      const n4 = f4 - r1 * (1 - v3 + vv - usq);
      const n5 = f5 - r2 * (1 + p3 + pp - usq);
      const n7 = f7 - r2 * (1 - p3 + pp - usq);
      const n6 = f6 - r2 * (1 + q3 + qq - usq);
      const n8 = f8 - r2 * (1 - q3 + qq - usq);
      const dg = n5 + n6 + n7 + n8;
      const pxx = n1 + n3 + dg, pyy = n2 + n4 + dg, pxy = n5 - n6 + n7 - n8;
      const Q = Math.sqrt(pxx * pxx + 2 * pxy * pxy + pyy * pyy);
      let tau = 0.5 * (t0 + Math.sqrt(t0sq + CSM * Q / r));
      const s = sponge[idx]; if (tau < s) tau = s;
      const om = 1 / tau;
      f0 -= om * n0; f1 -= om * n1; f2 -= om * n2; f3 -= om * n3; f4 -= om * n4;
      f5 -= om * n5; f6 -= om * n6; f7 -= om * n7; f8 -= om * n8;
      if (fl & 2) {
        const gx = bfx[idx], gy = bfy[idx];
        f1 += gx / 3; f3 -= gx / 3; f2 += gy / 3; f4 -= gy / 3;
        f5 += (gx + gy) / 12; f7 -= (gx + gy) / 12; f6 += (gy - gx) / 12; f8 -= (gy - gx) / 12;
      }
      b0[idx] = f0; b1[idx] = f1; b2[idx] = f2; b3[idx] = f3; b4[idx] = f4;
      b5[idx] = f5; b6[idx] = f6; b7[idx] = f7; b8[idx] = f8;
    }
  }
  // inlet (left), outlet (right, zero gradient), top (free air)
  for (let y = GR; y < NY; y++) {
    const i0 = y * NX, i1 = i0 + NX - 1;
    for (let i = 0; i < 9; i++) { b[i][i0] = E0[i]; b[i][i1] = b[i][i1 - 1]; }
    ux[i1] = ux[i1 - 1]; uy[i1] = uy[i1 - 1]; rho[i1] = rho[i1 - 1];
  }
  const top = (NY - 1) * NX;
  for (let i = 0; i < 9; i++) b[i].fill(E0[i], top, top + NX);
  cur = b; nxt = a;
}

// ---------- colours ----------
function lut(stops) {
  const out = new Uint32Array(256);
  for (let i = 0; i < 256; i++) {
    const t = i / 255; let k = 0;
    while (k < stops.length - 2 && t > stops[k + 1][0]) k++;
    const s0 = stops[k], s1 = stops[k + 1];
    const f = Math.min(1, Math.max(0, (t - s0[0]) / (s1[0] - s0[0])));
    const c0 = parseInt(s0[1].slice(1), 16), c1 = parseInt(s1[1].slice(1), 16);
    const R = Math.round(((c0 >> 16) & 255) * (1 - f) + ((c1 >> 16) & 255) * f);
    const G = Math.round(((c0 >> 8) & 255) * (1 - f) + ((c1 >> 8) & 255) * f);
    const B = Math.round((c0 & 255) * (1 - f) + (c1 & 255) * f);
    out[i] = (255 << 24) | (B << 16) | (G << 8) | R;
  }
  return out;
}
const STOPS_SEQ = [[0, '#030a13'], [0.3, '#0f3554'], [0.556, '#2b7a92'], [0.74, '#97cdbf'], [0.88, '#f2e7b4'], [1, '#f6ad5c']];
const STOPS_DIV = [[0, '#c9f0f8'], [0.22, '#37a9cf'], [0.5, '#0a1d2e'], [0.78, '#d68527'], [1, '#ffe6ae']];
const LUT_SEQ = lut(STOPS_SEQ), LUT_DIV = lut(STOPS_DIV);
const SOLID_C = (255 << 24) | (0x0c << 16) | (0x07 << 8) | 0x03;
const GROUND_C = (255 << 24) | (0x22 << 16) | (0x1c << 8) | 0x14;
function gradCss(stops) { return 'linear-gradient(90deg,' + stops.map(s => s[1] + ' ' + (s[0] * 100).toFixed(1) + '%').join(',') + ')'; }

// ---------- drawing ----------
const cv = document.getElementById('tunnel'), ctx = cv.getContext('2d');
const off = document.createElement('canvas'); off.width = NX; off.height = NY;
const octx = off.getContext('2d');
const img = octx.createImageData(NX, NY), pix = new Uint32Array(img.data.buffer);
let cw = 800, ch = 368, dpr = 1;

function fitCanvas() {
  dpr = Math.min(2, window.devicePixelRatio || 1);
  const w = Math.max(200, Math.round(cv.clientWidth * dpr)), h = Math.round(w * NY / NX);
  if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
  cw = w; ch = h; dirty = true;
}

function paintField() {
  const view = P.view;
  let rref = 0, cnt = 0;
  if (view === 'pres') { for (let y = GR + 1; y < NY - 1; y++) { rref += rho[y * NX + 3]; cnt++; } rref /= cnt; }
  const vs = 255 / (1.8 * U0), ws = 127 / 0.035, ps = 127 / (1.5 * 1.5 * U0 * U0);
  for (let y = 0; y < NY; y++) {
    const row = (NY - 1 - y) * NX, base = y * NX;
    for (let x = 0; x < NX; x++) {
      const i = base + x, m = mask[i];
      if (m) { pix[row + x] = m === GROUND ? GROUND_C : SOLID_C; continue; }
      let c;
      if (view === 'vel') {
        let k = Math.sqrt(ux[i] * ux[i] + uy[i] * uy[i]) * vs;
        c = LUT_SEQ[k > 255 ? 255 : k | 0];
      } else if (view === 'vort') {
        let w = 0;
        if (x > 0 && x < NX - 1 && y > GR && y < NY - 1) w = 0.5 * (uy[i + 1] - uy[i - 1] - ux[i + NX] + ux[i - NX]);
        let k = 128 + w * ws;
        c = LUT_DIV[k < 0 ? 0 : k > 255 ? 255 : k | 0];
      } else {
        let k = 128 + (rho[i] - rref) * ps;
        c = LUT_DIV[k < 0 ? 0 : k > 255 ? 255 : k | 0];
      }
      pix[row + x] = c;
    }
  }
  octx.putImageData(img, 0, 0);
}

function render() {
  const kx = cw / NX, ky = ch / NY;
  paintField();
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = 'high';
  ctx.drawImage(off, 0, 0, cw, ch);
  const px = dpr;
  const G = (gx, gy) => [gx * kx, ch - gy * ky];
  const W = w => G(w[0] * S, GR + w[1] * S);
  const B = (x, y) => W(bodyW(x, y));
  const C = (x, y) => W(canW(x, y));

  // smoke
  if (P.smoke) {
    ctx.beginPath();
    for (let r = 0; r < RAKES; r++) {
      let pen = false, lx = 0, ly = 0;
      for (let k = 0; k < RK; k++) {
        const o = r * RK + (head[r] + k) % RK, x = SX[o];
        if (x !== x) { pen = false; continue; }
        const y = SY[o];
        if (pen && Math.abs(x - lx) + Math.abs(y - ly) < 7) ctx.lineTo(x * kx, ch - y * ky);
        else ctx.moveTo(x * kx, ch - y * ky);
        pen = true; lx = x; ly = y;
      }
    }
    ctx.strokeStyle = P.view === 'vel' ? 'rgba(241,233,214,0.42)' : 'rgba(241,233,214,0.5)';
    ctx.lineWidth = 1 * px; ctx.lineJoin = 'round';
    ctx.stroke();
  }

  // moving ground
  const gy0 = ch - GR * ky;
  ctx.strokeStyle = 'rgba(241,233,214,0.55)'; ctx.lineWidth = 1 * px;
  ctx.beginPath(); ctx.moveTo(0, gy0); ctx.lineTo(cw, gy0);
  const hs = 0.6 * S * kx, shift = (groundShift * kx) % hs;
  for (let x = shift - hs; x < cw + hs; x += hs) { ctx.moveTo(x, gy0); ctx.lineTo(x - GR * ky * 0.9, ch); }
  ctx.stroke();

  const line = (T, a, b) => { const p = T(a[0], a[1]), q = T(b[0], b[1]); ctx.moveTo(p[0], p[1]); ctx.lineTo(q[0], q[1]); };
  const path = (T, pts, close) => { pts.forEach((p, i) => { const q = T(p[0], p[1]); if (i) ctx.lineTo(q[0], q[1]); else ctx.moveTo(q[0], q[1]); }); if (close) ctx.closePath(); };
  const SILK = '#f1e9d6', BAMBOO = '#d9a441';

  // height dimension
  const wb = bodyW(WHEEL_X, WHEEL_Y), wp = W([wb[0], wb[1] - WHEEL_R]), wg = W([wb[0], 0]);
  ctx.strokeStyle = 'rgba(241,233,214,0.7)'; ctx.lineWidth = 1 * px;
  ctx.beginPath(); ctx.moveTo(wp[0], wp[1]); ctx.lineTo(wg[0], wg[1]);
  ctx.moveTo(wp[0] - 5 * px, wp[1]); ctx.lineTo(wp[0] + 5 * px, wp[1]); ctx.stroke();
  ctx.font = `${11 * px}px 'IBM Plex Mono', ui-monospace, Menlo, monospace`;
  ctx.fillStyle = 'rgba(241,233,214,0.9)'; ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
  if (wg[1] - wp[1] > 16 * px) ctx.fillText(fmt(P.alt, 1) + ' m', wp[0] + 7 * px, (wp[1] + wg[1]) / 2);

  // fuselage
  ctx.beginPath(); path(B, FUSE, true);
  if (P.center) { ctx.fillStyle = 'rgba(241,233,214,0.16)'; ctx.fill(); }
  ctx.strokeStyle = 'rgba(241,233,214,0.62)'; ctx.lineWidth = 1 * px; ctx.stroke();
  ctx.beginPath();
  for (let k = 0; k < 9; k++) {
    const xa = -7.3 + k * 0.95, xb = Math.min(1.2, xa + 0.95);
    const top = x => 0.28 + (x + 7.3) * (0.22 / 5.7), bot = x => 0.02 - (x + 7.3) * (0.52 / 5.7);
    const cl = x => Math.min(-1.6, x);
    line(B, [xa, k % 2 ? top(cl(xa)) : bot(cl(xa))], [xb, k % 2 ? bot(cl(xb)) : top(cl(xb))]);
  }
  ctx.strokeStyle = 'rgba(217,164,65,0.5)'; ctx.stroke();

  // struts and wires
  ctx.beginPath();
  for (const x of [-0.95, 0.95]) line(B, [x, -WING_Y + 0.07], [x, WING_Y + 0.07]);
  for (const x of [-0.8, 0.8]) line(C, [CPX + x, CPY - CAN_GAP], [CPX + x, CPY + CAN_GAP]);
  line(B, [-0.5, -WING_Y], [WHEEL_X, WHEEL_Y]); line(B, [0.5, -WING_Y], [WHEEL_X, WHEEL_Y]);
  ctx.strokeStyle = BAMBOO; ctx.lineWidth = 1.6 * px; ctx.stroke();
  ctx.beginPath();
  line(B, [-0.95, -WING_Y], [0.95, WING_Y]); line(B, [-0.95, WING_Y], [0.95, -WING_Y]);
  line(C, [CPX - 0.8, CPY - CAN_GAP], [CPX + 0.8, CPY + CAN_GAP]); line(C, [CPX - 0.8, CPY + CAN_GAP], [CPX + 0.8, CPY - CAN_GAP]);
  ctx.strokeStyle = 'rgba(241,233,214,0.35)'; ctx.lineWidth = 0.8 * px; ctx.stroke();

  // wings and canard
  ctx.beginPath(); path(B, WING_LO); path(B, WING_UP); path(C, CAN_LO); path(C, CAN_UP);
  ctx.strokeStyle = SILK; ctx.lineWidth = Math.max(2 * px, 1.9 * kx); ctx.lineCap = 'round'; ctx.stroke(); ctx.lineCap = 'butt';

  // wheel
  const wc = B(WHEEL_X, WHEEL_Y), wr = WHEEL_R * S * kx;
  ctx.beginPath(); ctx.arc(wc[0], wc[1], wr, 0, 6.2832);
  if (P.center) { ctx.fillStyle = 'rgba(241,233,214,0.16)'; ctx.fill(); }
  for (let k = 0; k < 4; k++) { const an = k * 0.785 + groundShift * 0.19; ctx.moveTo(wc[0] - wr * Math.cos(an), wc[1] - wr * Math.sin(an)); ctx.lineTo(wc[0] + wr * Math.cos(an), wc[1] + wr * Math.sin(an)); }
  ctx.strokeStyle = 'rgba(241,233,214,0.8)'; ctx.lineWidth = 1 * px; ctx.stroke();

  // basket, pilot and engine
  ctx.beginPath();
  path(B, [[-2.2, -0.62], [-2.2, 0.1], [-1.5, 0.1], [-1.5, -0.62]], true);
  path(B, [[0.3, 0.5], [0.3, 0.95], [1.1, 0.95], [1.1, 0.5]], true);
  if (P.center) { ctx.fillStyle = 'rgba(241,233,214,0.16)'; ctx.fill(); }
  ctx.strokeStyle = 'rgba(241,233,214,0.75)'; ctx.stroke();
  ctx.beginPath(); line(B, [-1.85, 0.1], [-1.85, 0.92]); line(B, [-1.85, 0.7], [-1.55, 0.35]);
  ctx.strokeStyle = SILK; ctx.lineWidth = 2 * px; ctx.stroke();
  const hd = B(-1.85, 1.06);
  ctx.beginPath(); ctx.arc(hd[0], hd[1], 0.13 * S * kx, 0, 6.2832); ctx.fillStyle = SILK; ctx.fill();
  ctx.beginPath(); line(B, [-2.07, 1.19], [-1.63, 1.19]); ctx.lineWidth = 1.6 * px; ctx.stroke();

  // propeller
  const bl = PROP_R * Math.abs(Math.cos(propPhase));
  ctx.beginPath(); line(B, [PROP_X, PROP_Y - PROP_R], [PROP_X, PROP_Y + PROP_R]);
  ctx.strokeStyle = 'rgba(217,164,65,0.4)'; ctx.lineWidth = 1.2 * px; ctx.stroke();
  ctx.beginPath(); line(B, [PROP_X, PROP_Y - bl], [PROP_X, PROP_Y + bl]);
  ctx.strokeStyle = BAMBOO; ctx.lineWidth = 3 * px; ctx.lineCap = 'round'; ctx.stroke(); ctx.lineCap = 'butt';

  // scale bar and wind
  ctx.strokeStyle = 'rgba(241,233,214,0.9)'; ctx.lineWidth = 1 * px;
  const s0 = W([0.7, 0.5]), s1 = W([2.7, 0.5]);
  ctx.beginPath(); ctx.moveTo(s0[0], s0[1] - 4 * px); ctx.lineTo(s0[0], s0[1]); ctx.lineTo(s1[0], s1[1]); ctx.lineTo(s1[0], s1[1] - 4 * px); ctx.stroke();
  ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText('2 m', (s0[0] + s1[0]) / 2, s0[1] - 5 * px);
  ctx.textAlign = 'left'; ctx.textBaseline = 'top';
  ctx.fillText('wind ' + fmt(P.v, 0) + ' km/h →', 10 * px, 9 * px);
}

// ---------- readouts ----------
const $ = id => document.getElementById(id);
function fmt(v, d) { return v.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d }); }
function forces() {
  const V = P.v / 3.6;
  const conv = RHO_AIR * V * V * DX / (U0 * U0) / G0;      // kgf per metre of span
  const lw = EMA[3] * conv * SPAN[1], lc = EMA[5] * conv * SPAN[2], lf = EMA[7] * conv * SPAN[3];
  const dw = EMA[2] * conv * SPAN[1], dc = EMA[4] * conv * SPAN[2], df = EMA[6] * conv * SPAN[3];
  const lift = lw + lc + lf, drag = dw + dc + df;
  return { V, lw, lc, lift, drag, ratio: lift / WEIGHT_KGF };
}
function readouts() {
  const F = forces(), V = F.V;
  $('r-ratio').textContent = fmt(F.ratio, 2);
  $('r-lift').textContent = fmt(F.lift, 0) + ' kgf';
  $('r-wing').textContent = fmt(F.lw, 0) + ' kgf';
  $('r-can').textContent = fmt(F.lc, 0) + ' kgf';
  $('r-drag').textContent = fmt(F.drag, 0) + ' kgf';
  $('r-ld').textContent = F.drag > 0.5 ? fmt(F.lift / F.drag, 1) : '–';
  $('r-cl').textContent = fmt(F.lift * G0 / (0.5 * RHO_AIR * V * V * AREA), 2);
  $('r-thrust').textContent = '+' + fmt(P.prop / 100 * P.v, 0) + ' km/h';
  $('r-re').textContent = '≈ ' + fmt(V * CHORD / NU_AIR / 1e6, 1) + ' million';
  let state = 'wait', text = 'Flow settling';
  if (settle <= 0) {
    if (F.ratio < 0.9) { state = 'warn'; text = 'Does not carry the 300 kgf'; }
    else if (F.ratio <= 1.1) { state = 'good'; text = 'Level flight'; }
    else { state = 'good'; text = 'Excess lift: climbs'; }
  }
  $('r-pill').dataset.state = state; $('r-pill').textContent = text;
  $('r-beam').dataset.state = state;
  $('r-fill').style.width = Math.max(0, Math.min(100, F.ratio * 50)) + '%';
  $('clock').textContent = 't = ' + fmt(simT, 2) + ' s · ' + fmt(stepsPS, 0) + ' steps/s';
  drawHist();
}

const hc = $('hist'), hx = hc.getContext('2d');
function drawHist() {
  const st = getComputedStyle(document.documentElement);
  const ink = st.getPropertyValue('--ink').trim(), muted = st.getPropertyValue('--muted').trim();
  const lineC = st.getPropertyValue('--line').trim(), acc = st.getPropertyValue('--accent').trim();
  const d = Math.min(2, window.devicePixelRatio || 1);
  const w = Math.round(hc.clientWidth * d), h = Math.round(hc.clientHeight * d);
  if (hc.width !== w || hc.height !== h) { hc.width = w; hc.height = h; }
  hx.clearRect(0, 0, w, h);
  const padL = 34 * d, padR = 8 * d, padT = 8 * d, padB = 8 * d;
  let top = 2;
  for (const v of hist) if (v > top) top = v;
  top = Math.ceil(top * 2) / 2;
  let bot = 0;
  for (const v of hist) if (v < bot) bot = v;
  bot = Math.floor(bot * 2) / 2;
  const Y = v => padT + (h - padT - padB) * (1 - (v - bot) / (top - bot));
  hx.font = `${10.5 * d}px 'IBM Plex Mono', ui-monospace, Menlo, monospace`;
  hx.textAlign = 'right'; hx.textBaseline = 'middle'; hx.lineWidth = 1 * d;
  for (const v of [bot, 1, top]) {
    hx.strokeStyle = v === 1 ? muted : lineC;
    hx.setLineDash(v === 1 ? [4 * d, 4 * d] : []);
    hx.beginPath(); hx.moveTo(padL, Y(v)); hx.lineTo(w - padR, Y(v)); hx.stroke();
    hx.fillStyle = v === 1 ? ink : muted;
    hx.fillText(fmt(v, v % 1 ? 1 : 0), padL - 6 * d, Y(v));
  }
  hx.setLineDash([]);
  if (hist.length > 1) {
    const HMAX = 520, x0 = padL, xw = w - padR - padL;
    hx.beginPath();
    hist.forEach((v, i) => { const x = x0 + xw * (i + HMAX - hist.length) / (HMAX - 1), y = Y(v); if (i) hx.lineTo(x, y); else hx.moveTo(x, y); });
    hx.strokeStyle = acc; hx.lineWidth = 1.8 * d; hx.lineJoin = 'round'; hx.stroke();
    const lx = x0 + xw, ly = Y(hist[hist.length - 1]);
    hx.beginPath(); hx.arc(lx, ly, 3 * d, 0, 6.2832); hx.fillStyle = acc; hx.fill();
  }
  const span = 520 * 30 * U0 * DX / (P.v / 3.6);
  $('hist-cap').textContent = 'Lift ÷ weight over the last ' + fmt(span, 1) + ' simulated seconds. The dashed line is the weight.';
}

// ---------- main loop ----------
let stepsPS = 0, lastRead = 0, lastT = performance.now(), stepCount = 0;
function advance(n) {
  FA.fill(0);
  for (let s = 0; s < n; s++) step();
  const al = 1 - Math.exp(-n / 500);
  for (let i = 2; i < 8; i++) EMA[i] += al * (FA[i] / n - EMA[i]);
  if (P.smoke) advectSmoke(n);
  simSteps += n; settle -= n; stepCount += n;
  simT += n * U0 * DX / (P.v / 3.6);
  groundShift += n * U0;
  propPhase += n * 0.012 * (0.2 + P.prop / 100);
  histAcc += n;
  while (histAcc >= 30) { histAcc -= 30; hist.push(forces().ratio); if (hist.length > 520) hist.shift(); }
}
function frame(now) {
  requestAnimationFrame(frame);
  if (running) {
    const t0 = performance.now(), n = stepsPF;
    advance(n);
    const ms = (performance.now() - t0) / n;
    msAvg = msAvg * 0.9 + ms * 0.1;
    stepsPF = Math.max(1, Math.min(12, Math.floor(11 / msAvg)));
    dirty = true;
  }
  if (dirty) { render(); dirty = false; }
  if (now - lastRead > 220) {
    stepsPS = stepCount * 1000 / Math.max(1, now - lastT); stepCount = 0; lastT = now; lastRead = now;
    readouts();
  }
}

// ---------- controls ----------
function reFromSlider(t) { return Math.round(100 * Math.pow(30, t / 100) / 10) * 10; }
function setTau() { tau0 = 0.5 + 3 * U0 * CHORD_CELLS / P.re; }
function syncOutputs() {
  $('o-v').textContent = fmt(P.v, 0) + ' km/h';
  $('o-aoa').textContent = fmt(P.aoa, 1) + '°';
  $('o-canard').textContent = (P.canard > 0 ? '+' : '') + fmt(P.canard, 1) + '°';
  $('o-alt').textContent = fmt(P.alt, 1) + ' m';
  $('o-prop').textContent = fmt(P.prop, 0) + '%';
  $('o-re').textContent = fmt(P.re, 0);
}
function readControls() {
  P.v = +$('k-v').value; P.aoa = +$('k-aoa').value; P.canard = +$('k-canard').value;
  P.alt = +$('k-alt').value; P.prop = +$('k-prop').value; P.re = reFromSlider(+$('k-re').value);
  P.center = $('center').checked; P.smoke = $('smoke').checked;
  P.view = document.querySelector('input[name="view"]:checked').value;
}
function legend() {
  const strip = $('leg-strip');
  if (P.view === 'vel') {
    strip.style.background = gradCss(STOPS_SEQ);
    $('leg-lo').textContent = 'still'; $('leg-hi').textContent = '1.8 × flight (' + fmt(1.8 * P.v, 0) + ' km/h)';
  } else if (P.view === 'vort') {
    strip.style.background = gradCss(STOPS_DIV);
    $('leg-lo').textContent = 'clockwise spin'; $('leg-hi').textContent = 'counter-clockwise spin';
  } else {
    strip.style.background = gradCss(STOPS_DIV);
    $('leg-lo').textContent = 'suction (Cp −1.5)'; $('leg-hi').textContent = 'overpressure (Cp +1.5)';
  }
}
for (const id of ['k-aoa', 'k-canard', 'k-alt', 'k-prop']) $(id).addEventListener('input', () => { readControls(); syncOutputs(); rebuild(); });
$('k-v').addEventListener('input', () => { readControls(); syncOutputs(); legend(); dirty = true; });
$('k-re').addEventListener('input', () => { readControls(); setTau(); syncOutputs(); });
$('center').addEventListener('change', () => { readControls(); rebuild(); });
$('smoke').addEventListener('change', () => { readControls(); if (P.smoke) seedSmoke(); dirty = true; });
document.querySelectorAll('input[name="view"]').forEach(el => el.addEventListener('change', () => { readControls(); legend(); dirty = true; }));
$('pause').addEventListener('click', () => { running = !running; $('pause').textContent = running ? 'Pause' : 'Resume'; });
$('reset').addEventListener('click', () => { initFlow(); });
function preset(v, alt, aoa) {
  $('k-v').value = v; $('k-alt').value = alt; $('k-aoa').value = aoa; $('k-canard').value = 0;
  readControls(); syncOutputs(); legend(); rebuild();
}
$('preset-out').addEventListener('click', () => preset(31, 3, 10));
$('preset-nov').addEventListener('click', () => preset(37, 6, 6));

// ---------- start ----------
document.documentElement.lang = 'en';
readControls(); setTau(); syncOutputs(); legend();
rebuild(); initFlow();
new ResizeObserver(fitCanvas).observe(cv);
fitCanvas();
const calm = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
advance(calm ? 600 : 120);
if (calm) { running = false; $('pause').textContent = 'Resume'; }
readouts();
requestAnimationFrame(frame);
})();
