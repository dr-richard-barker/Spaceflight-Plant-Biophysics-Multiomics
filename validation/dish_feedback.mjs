// validation/dish_feedback.ts
import { writeFileSync as writeFileSync2 } from "node:fs";

// src/solver/lbm/d2q9.ts
var Q = 9;
var CX = new Int8Array([0, 1, 0, -1, 0, 1, -1, -1, 1]);
var CY = new Int8Array([0, 0, 1, 0, -1, 1, 1, -1, -1]);
var W = new Float64Array([
  4 / 9,
  1 / 9,
  1 / 9,
  1 / 9,
  1 / 9,
  1 / 36,
  1 / 36,
  1 / 36,
  1 / 36
]);
var OPP = new Uint8Array([0, 3, 4, 1, 2, 7, 8, 5, 6]);
var CS2 = 1 / 3;
function tauFromViscosity(nu) {
  return nu / CS2 + 0.5;
}
function viscosityFromRe(u, length, Re) {
  return u * length / Re;
}
function feq(i, rho, ux, uy) {
  const eu = CX[i] * ux + CY[i] * uy;
  const usq = ux * ux + uy * uy;
  return W[i] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * usq);
}

// src/solver/lbm/LBMFluid.ts
var LBMFluid = class {
  nx;
  ny;
  size;
  tau;
  omega;
  f;
  ftmp;
  rho;
  ux;
  uy;
  /** 1 = solid (bounce-back wall), 0 = fluid. */
  solid;
  /** Wall velocity for solid cells (used by moving walls; 0 for static). */
  uwx;
  uwy;
  /** Optional body force per cell; null disables the forcing branch entirely. */
  fx = null;
  fy = null;
  constructor(nx, ny, tau) {
    this.nx = nx;
    this.ny = ny;
    this.size = nx * ny;
    this.tau = tau;
    this.omega = 1 / tau;
    this.f = new Float32Array(Q * this.size);
    this.ftmp = new Float32Array(Q * this.size);
    this.rho = new Float32Array(this.size);
    this.ux = new Float32Array(this.size);
    this.uy = new Float32Array(this.size);
    this.solid = new Uint8Array(this.size);
    this.uwx = new Float32Array(this.size);
    this.uwy = new Float32Array(this.size);
  }
  setTau(tau) {
    this.tau = tau;
    this.omega = 1 / tau;
  }
  index(x, y) {
    return x + y * this.nx;
  }
  /** Enable the body-force arrays (allocates fx/fy). Idempotent. */
  enableForcing() {
    if (!this.fx) this.fx = new Float32Array(this.size);
    if (!this.fy) this.fy = new Float32Array(this.size);
  }
  /** Initialise every fluid cell to equilibrium at a uniform density and velocity. */
  setEquilibrium(rho0, ux0, uy0) {
    for (let c = 0; c < this.size; c++) {
      this.rho[c] = rho0;
      this.ux[c] = ux0;
      this.uy[c] = uy0;
      const base = c * Q;
      for (let i = 0; i < Q; i++) this.f[base + i] = feq(i, rho0, ux0, uy0);
    }
  }
  /**
   * One timestep: fused macroscopic + BGK collision (with optional Guo forcing),
   * then push-streaming with bounce-back, then buffer swap.
   */
  collideAndStream() {
    const { nx, ny, f, omega, solid, fx, fy } = this;
    for (let c = 0; c < this.size; c++) {
      if (solid[c]) continue;
      const base = c * Q;
      let rho = 0;
      let mx = 0;
      let my = 0;
      for (let i = 0; i < Q; i++) {
        const fi = f[base + i];
        rho += fi;
        mx += CX[i] * fi;
        my += CY[i] * fi;
      }
      const Fx = fx ? fx[c] : 0;
      const Fy = fy ? fy[c] : 0;
      const ux = (mx + 0.5 * Fx) / rho;
      const uy = (my + 0.5 * Fy) / rho;
      this.rho[c] = rho;
      this.ux[c] = ux;
      this.uy[c] = uy;
      const forceCoeff = 1 - 0.5 * omega;
      for (let i = 0; i < Q; i++) {
        const cix = CX[i];
        const ciy = CY[i];
        const eu = cix * ux + ciy * uy;
        const usq = ux * ux + uy * uy;
        const eqi = W[i] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * usq);
        let fi = f[base + i] - omega * (f[base + i] - eqi);
        if (fx || fy) {
          const gi = forceCoeff * W[i] * ((3 * (cix - ux) + 9 * eu * cix) * Fx + (3 * (ciy - uy) + 9 * eu * ciy) * Fy);
          fi += gi;
        }
        f[base + i] = fi;
      }
    }
    const ftmp = this.ftmp;
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const c = x + y * nx;
        if (solid[c]) continue;
        const base = c * Q;
        const rhoC = this.rho[c];
        for (let i = 0; i < Q; i++) {
          const xn = x + CX[i];
          const yn = y + CY[i];
          const fi = f[base + i];
          if (xn < 0 || xn >= nx || yn < 0 || yn >= ny) {
            ftmp[base + OPP[i]] = fi;
            continue;
          }
          const n = xn + yn * nx;
          if (solid[n]) {
            const uwEu = CX[i] * this.uwx[n] + CY[i] * this.uwy[n];
            ftmp[base + OPP[i]] = fi - 2 * W[i] * rhoC * uwEu / CS2;
          } else {
            ftmp[n * Q + i] = fi;
          }
        }
      }
    }
    const tmp = this.f;
    this.f = this.ftmp;
    this.ftmp = tmp;
  }
  /** Recompute macroscopic fields from the current populations (for diagnostics/render). */
  computeMacroscopic() {
    const { f, solid, fx, fy } = this;
    for (let c = 0; c < this.size; c++) {
      if (solid[c]) {
        this.rho[c] = 1;
        this.ux[c] = 0;
        this.uy[c] = 0;
        continue;
      }
      const base = c * Q;
      let rho = 0;
      let mx = 0;
      let my = 0;
      for (let i = 0; i < Q; i++) {
        const fi = f[base + i];
        rho += fi;
        mx += CX[i] * fi;
        my += CY[i] * fi;
      }
      const Fx = fx ? fx[c] : 0;
      const Fy = fy ? fy[c] : 0;
      this.rho[c] = rho;
      this.ux[c] = (mx + 0.5 * Fx) / rho;
      this.uy[c] = (my + 0.5 * Fy) / rho;
    }
  }
  /** Speed |u| at a cell (call computeMacroscopic first for solid-safe values). */
  speed(c) {
    return Math.hypot(this.ux[c], this.uy[c]);
  }
  /** Vorticity (∂uy/∂x − ∂ux/∂y) via central differences; 0 on the domain border. */
  vorticity(x, y) {
    if (x <= 0 || x >= this.nx - 1 || y <= 0 || y >= this.ny - 1) return 0;
    const duy_dx = this.uy[this.index(x + 1, y)] - this.uy[this.index(x - 1, y)];
    const dux_dy = this.ux[this.index(x, y + 1)] - this.ux[this.index(x, y - 1)];
    return 0.5 * (duy_dx - dux_dy);
  }
};

// src/solver/lbm/ScalarField.ts
var CX5 = new Int8Array([0, 1, 0, -1, 0]);
var CY5 = new Int8Array([0, 0, 1, 0, -1]);
var W5 = new Float64Array([1 / 3, 1 / 6, 1 / 6, 1 / 6, 1 / 6]);
var OPP5 = new Uint8Array([0, 3, 4, 1, 2]);
var CS2_5 = 1 / 3;
var Q5 = 5;
function tauFromDiffusivity(D) {
  return D / CS2_5 + 0.5;
}
var ScalarField = class {
  fluid;
  size;
  tauC;
  omegaC;
  g;
  gtmp;
  /** Concentration field (macroscopic). */
  C;
  /** Optional per-cell source/sink (e.g. stomatal flux). null disables the branch. */
  source = null;
  constructor(fluid, diffusivity, initialC = 0) {
    this.fluid = fluid;
    this.size = fluid.size;
    this.tauC = tauFromDiffusivity(diffusivity);
    this.omegaC = 1 / this.tauC;
    this.g = new Float32Array(Q5 * this.size);
    this.gtmp = new Float32Array(Q5 * this.size);
    this.C = new Float32Array(this.size).fill(initialC);
    this.setEquilibrium(initialC);
  }
  setDiffusivity(D) {
    this.tauC = tauFromDiffusivity(D);
    this.omegaC = 1 / this.tauC;
  }
  enableSource() {
    if (!this.source) this.source = new Float32Array(this.size);
  }
  /**
   * Per-solid-cell Dirichlet value. NaN → the wall stays zero-flux (Neumann /
   * adiabatic). A finite value pins the concentration at the wall midpoint via
   * anti-bounce-back (e.g. hot/cold walls, or a fixed ambient boundary).
   */
  dirichlet = null;
  setDirichlet(cell, value) {
    if (!this.dirichlet) this.dirichlet = new Float32Array(this.size).fill(NaN);
    this.dirichlet[cell] = value;
  }
  feq5(i, C, ux, uy) {
    const eu = CX5[i] * ux + CY5[i] * uy;
    return W5[i] * C * (1 + eu / CS2_5);
  }
  setEquilibrium(C0) {
    for (let c = 0; c < this.size; c++) {
      this.C[c] = C0;
      const base = c * Q5;
      for (let i = 0; i < Q5; i++) this.g[base + i] = this.feq5(i, C0, 0, 0);
    }
  }
  /** Initialise the concentration from a spatial function C0(x, y) (u = 0). */
  initField(fn) {
    const nx = this.fluid.nx;
    for (let c = 0; c < this.size; c++) {
      const x = c % nx;
      const y = c / nx | 0;
      const C0 = fn(x, y);
      this.C[c] = C0;
      const base = c * Q5;
      for (let i = 0; i < Q5; i++) this.g[base + i] = this.feq5(i, C0, 0, 0);
    }
  }
  /** One advection–diffusion step using the fluid's current velocity field. */
  step() {
    const { fluid, size, omegaC, source } = this;
    const { nx, ny, solid, ux, uy } = fluid;
    const g = this.g;
    for (let c = 0; c < size; c++) {
      if (solid[c]) continue;
      const base = c * Q5;
      let C = 0;
      for (let i = 0; i < Q5; i++) C += g[base + i];
      const S = source ? source[c] : 0;
      C += 0.5 * S;
      this.C[c] = C;
      const vx = ux[c];
      const vy = uy[c];
      for (let i = 0; i < Q5; i++) {
        const eqi = this.feq5(i, C, vx, vy);
        let gi = g[base + i] - omegaC * (g[base + i] - eqi);
        if (source) gi += (1 - 0.5 * omegaC) * W5[i] * S;
        g[base + i] = gi;
      }
    }
    const gtmp = this.gtmp;
    const dir = this.dirichlet;
    for (let y = 0; y < ny; y++) {
      for (let x = 0; x < nx; x++) {
        const c = x + y * nx;
        if (solid[c]) continue;
        const base = c * Q5;
        for (let i = 0; i < Q5; i++) {
          const xn = x + CX5[i];
          const yn = y + CY5[i];
          const gi = g[base + i];
          if (xn < 0 || xn >= nx || yn < 0 || yn >= ny) {
            gtmp[base + OPP5[i]] = gi;
            continue;
          }
          const n = xn + yn * nx;
          if (solid[n]) {
            const cw = dir ? dir[n] : NaN;
            if (!Number.isNaN(cw)) {
              gtmp[base + OPP5[i]] = -gi + 2 * W5[i] * cw;
            } else {
              gtmp[base + OPP5[i]] = gi;
            }
          } else {
            gtmp[n * Q5 + i] = gi;
          }
        }
      }
    }
    const tmp = this.g;
    this.g = this.gtmp;
    this.gtmp = tmp;
  }
};

// src/solver/lbm/buoyancy.ts
var GRAVITY = {
  microgravity: 0,
  moon: 1.62,
  mars: 3.71,
  earth: 9.81
};
function applyBoussinesqForce(fluid, gx, gy, contributors) {
  fluid.enableForcing();
  const fx = fluid.fx;
  const fy = fluid.fy;
  const { size, solid, rho } = fluid;
  for (let c = 0; c < size; c++) {
    if (solid[c]) {
      fx[c] = 0;
      fy[c] = 0;
      continue;
    }
    let expansion = 0;
    for (let k = 0; k < contributors.length; k++) {
      const ct = contributors[k];
      expansion += ct.beta * (ct.field.C[c] - ct.ref);
    }
    const factor = rho[c] * expansion;
    fx[c] = gx * factor;
    fy[c] = gy * factor;
  }
}

// src/solver/diagnostics/dimensionless.ts
function reynolds(u, length, nu) {
  return u * length / nu;
}
function grashof(g, beta, deltaC, length, nu) {
  return g * beta * deltaC * length ** 3 / (nu * nu);
}
function peclet(u, length, D) {
  return u * length / D;
}
function computeDimensionless(params) {
  const { u, length, nu } = params;
  const D = params.D ?? nu;
  const g = params.g ?? 0;
  const beta = params.beta ?? 0;
  const deltaC = params.deltaC ?? 0;
  const Sc = nu / D;
  const Gr = grashof(g, beta, deltaC, length, nu);
  return {
    Re: reynolds(u, length, nu),
    Gr,
    Ra: Gr * Sc,
    Pe: peclet(u, length, D),
    Sc
  };
}

// src/solver/diagnostics/ghia.ts
var GHIA_RE100 = {
  Re: 100,
  // vertical centreline: y (0..1) and u/U_lid
  yU: [
    0,
    0.0547,
    0.0625,
    0.0703,
    0.1016,
    0.1719,
    0.2813,
    0.4531,
    0.5,
    0.6172,
    0.7344,
    0.8516,
    0.9531,
    0.9609,
    0.9688,
    0.9766,
    1
  ],
  u: [
    0,
    -0.03717,
    -0.04192,
    -0.04775,
    -0.06434,
    -0.1015,
    -0.15662,
    -0.2109,
    -0.20581,
    -0.13641,
    332e-5,
    0.23151,
    0.68717,
    0.73722,
    0.78871,
    0.84123,
    1
  ],
  // horizontal centreline: x (0..1) and v/U_lid
  xV: [
    0,
    0.0625,
    0.0703,
    0.0781,
    0.0938,
    0.1563,
    0.2266,
    0.2344,
    0.5,
    0.8047,
    0.8594,
    0.9063,
    0.9453,
    0.9531,
    0.9609,
    0.9688,
    1
  ],
  v: [
    0,
    0.09233,
    0.10091,
    0.1089,
    0.12317,
    0.16077,
    0.17507,
    0.17527,
    0.05454,
    -0.24533,
    -0.22445,
    -0.16914,
    -0.10313,
    -0.08864,
    -0.07391,
    -0.05906,
    0
  ]
};
function ghiaL2Error(sampleU, sampleV) {
  let sq = 0;
  let n = 0;
  for (let i = 0; i < GHIA_RE100.yU.length; i++) {
    const d = sampleU(GHIA_RE100.yU[i]) - GHIA_RE100.u[i];
    sq += d * d;
    n++;
  }
  for (let i = 0; i < GHIA_RE100.xV.length; i++) {
    const d = sampleV(GHIA_RE100.xV[i]) - GHIA_RE100.v[i];
    sq += d * d;
    n++;
  }
  return Math.sqrt(sq / n);
}

// src/solver/diagnostics/strouhal.ts
var StrouhalProbe = class {
  samples = [];
  crossings = [];
  prev = 0;
  step = 0;
  warmup;
  constructor(warmupSteps = 2e3) {
    this.warmup = warmupSteps;
  }
  /** Feed the current transverse velocity at the probe each timestep. */
  record(uy) {
    this.step++;
    if (this.step < this.warmup) {
      this.prev = uy;
      return;
    }
    this.samples.push(uy);
    if (this.prev < 0 && uy >= 0) this.crossings.push(this.step);
    this.prev = uy;
  }
  /** Estimated shedding period (steps), or null if not enough cycles yet. */
  period() {
    if (this.crossings.length < 3) return null;
    let sum = 0;
    for (let i = 1; i < this.crossings.length; i++) {
      sum += this.crossings[i] - this.crossings[i - 1];
    }
    return sum / (this.crossings.length - 1);
  }
  /** Strouhal number St = f D / U = D / (period · U), or null if not ready. */
  strouhal(diameter, u) {
    const T = this.period();
    if (!T) return null;
    return diameter / (T * u);
  }
  cycles() {
    return Math.max(0, this.crossings.length - 1);
  }
  reset() {
    this.samples.length = 0;
    this.crossings.length = 0;
    this.prev = 0;
    this.step = 0;
  }
};

// src/solver/diagnostics/analytic.ts
function erf(x) {
  const sign = x < 0 ? -1 : 1;
  const ax = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * ax);
  const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-ax * ax);
  return sign * y;
}
function erfc(x) {
  return 1 - erf(x);
}
function diffusionStep(x, t, D, x0, cLeft, cRight) {
  if (t <= 0) return x < x0 ? cLeft : cRight;
  const arg = (x - x0) / (2 * Math.sqrt(D * t));
  return cRight + (cLeft - cRight) * 0.5 * erfc(arg);
}

// src/leaf/types.ts
var NEUTRAL_DEFORM = {
  curl: 0,
  cup: 1,
  twist: 0,
  droop: 0,
  thicknessScale: 1
};

// src/leaf/sdf.ts
function polygonSDF(poly, px, py) {
  const n = poly.length;
  if (n < 3) return Infinity;
  let minSq = Infinity;
  let inside = false;
  for (let i = 0, j = n - 1; i < n; j = i++) {
    const [ax, ay] = poly[j];
    const [bx, by] = poly[i];
    const ex = bx - ax;
    const ey = by - ay;
    const wx = px - ax;
    const wy = py - ay;
    const len2 = ex * ex + ey * ey;
    const t = len2 > 0 ? Math.max(0, Math.min(1, (wx * ex + wy * ey) / len2)) : 0;
    const dx = wx - t * ex;
    const dy = wy - t * ey;
    const sq = dx * dx + dy * dy;
    if (sq < minSq) minSq = sq;
    if (ay > py !== by > py) {
      const xCross = ax + (py - ay) / (by - ay) * ex;
      if (px < xCross) inside = !inside;
    }
  }
  const d = Math.sqrt(minSq);
  return inside ? -d : d;
}
function signedArea(poly) {
  const n = poly.length;
  let a = 0;
  for (let i = 0, j = n - 1; i < n; j = i++) {
    a += poly[j][0] * poly[i][1] - poly[i][0] * poly[j][1];
  }
  return a / 2;
}
function polygonArea(poly) {
  return Math.abs(signedArea(poly));
}
function polygonPerimeter(poly) {
  const n = poly.length;
  let p = 0;
  for (let i = 0, j = n - 1; i < n; j = i++) {
    p += Math.hypot(poly[i][0] - poly[j][0], poly[i][1] - poly[j][1]);
  }
  return p;
}
function ensureCCW(poly) {
  return signedArea(poly) < 0 ? [...poly].reverse() : [...poly];
}
function bounds(poly) {
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  for (const [x, y] of poly) {
    if (x < minX) minX = x;
    if (y < minY) minY = y;
    if (x > maxX) maxX = x;
    if (y > maxY) maxY = y;
  }
  return [minX, minY, maxX, maxY];
}
function crossingsAtY(poly, yq) {
  const n = poly.length;
  const xs = [];
  for (let i = 0, j = n - 1; i < n; j = i++) {
    const [ax, ay] = poly[j];
    const [bx, by] = poly[i];
    if (ay > yq !== by > yq) {
      xs.push(ax + (yq - ay) / (by - ay) * (bx - ax));
    }
  }
  xs.sort((a, b) => a - b);
  return xs;
}

// src/leaf/resample.ts
function cumulative(points, closed) {
  const n = points.length;
  const acc = new Array(closed ? n + 1 : n);
  acc[0] = 0;
  for (let i = 1; i < n; i++) {
    acc[i] = acc[i - 1] + Math.hypot(points[i][0] - points[i - 1][0], points[i][1] - points[i - 1][1]);
  }
  if (closed) {
    acc[n] = acc[n - 1] + Math.hypot(points[0][0] - points[n - 1][0], points[0][1] - points[n - 1][1]);
  }
  return acc;
}
function polylineLength(line) {
  const acc = cumulative(line, false);
  return acc[line.length - 1] ?? 0;
}
function sampleTable(table, key) {
  const n = table.length;
  if (n === 0) return 0;
  if (key <= table[0][0]) return table[0][1];
  if (key >= table[n - 1][0]) return table[n - 1][1];
  for (let i = 1; i < n; i++) {
    const [k1, v1] = table[i];
    if (key <= k1) {
      const [k0, v0] = table[i - 1];
      const span = k1 - k0;
      return span > 0 ? v0 + (key - k0) / span * (v1 - v0) : v1;
    }
  }
  return table[n - 1][1];
}

// src/leaf/loft.ts
var DEFAULT_LOFT = { nu: 25, nv: 61 };
function midribXAt(midrib, y) {
  const n = midrib.length;
  if (n === 0) return 0;
  if (n === 1) return midrib[0][0];
  if (y <= midrib[0][1]) return midrib[0][0];
  if (y >= midrib[n - 1][1]) return midrib[n - 1][0];
  for (let i = 1; i < n; i++) {
    if (y <= midrib[i][1]) {
      const span = midrib[i][1] - midrib[i - 1][1];
      const t = span > 0 ? (y - midrib[i - 1][1]) / span : 0;
      return midrib[i - 1][0] + t * (midrib[i][0] - midrib[i - 1][0]);
    }
  }
  return midrib[n - 1][0];
}
function loft(asset, deform = NEUTRAL_DEFORM, options = DEFAULT_LOFT) {
  const { nu, nv } = options;
  const outline = ensureCCW(asset.outline);
  const [, minY, , maxY] = bounds(outline);
  const lengthMm = polylineLength(asset.midrib) || maxY - minY;
  const span = maxY - minY;
  const eps = span * 1e-4;
  const y = new Float64Array(nv);
  const xMid = new Float64Array(nv);
  const halfWidth = new Float64Array(nv);
  const xL = new Float64Array(nv);
  const xR = new Float64Array(nv);
  let halfWidthMax = 0;
  for (let j = 0; j < nv; j++) {
    const yj = minY + eps + (span - 2 * eps) * j / (nv - 1);
    y[j] = yj;
    const xs = crossingsAtY(outline, yj);
    let left;
    let right;
    if (xs.length >= 2) {
      left = xs[0];
      right = xs[xs.length - 1];
    } else {
      left = right = midribXAt(asset.midrib, yj);
    }
    xL[j] = left;
    xR[j] = right;
    xMid[j] = Math.min(right, Math.max(left, midribXAt(asset.midrib, yj)));
    halfWidth[j] = (right - left) / 2;
    if (halfWidth[j] > halfWidthMax) halfWidthMax = halfWidth[j];
  }
  const thickness = asset.thicknessMm * deform.thicknessScale;
  const half = thickness / 2;
  const surf = new Float64Array(nu * nv * 3);
  for (let j = 0; j < nv; j++) {
    const s = j / (nv - 1);
    const widthScale = halfWidthMax > 0 ? halfWidth[j] / halfWidthMax : 0;
    const zMid = sampleTable(asset.profile, s) + deform.droop * s * s * lengthMm;
    const theta = deform.twist * s;
    const ct = Math.cos(theta);
    const st = Math.sin(theta);
    for (let i = 0; i < nu; i++) {
      const t = nu > 1 ? -1 + 2 * i / (nu - 1) : 0;
      const x = t < 0 ? xMid[j] + t * (xMid[j] - xL[j]) : xMid[j] + t * (xR[j] - xMid[j]);
      const z = zMid + deform.cup * sampleTable(asset.camber, t) * widthScale + deform.curl * t * t * t * t * halfWidth[j];
      const dx = x - xMid[j];
      const dz = z - zMid;
      const o = (j * nu + i) * 3;
      surf[o] = xMid[j] + dx * ct - dz * st;
      surf[o + 1] = y[j];
      surf[o + 2] = zMid + dx * st + dz * ct;
    }
  }
  const count = nu * nv;
  const normals = new Float64Array(count * 3);
  const top = new Float64Array(count * 3);
  const bottom = new Float64Array(count * 3);
  const at = (j, i, k) => surf[(j * nu + i) * 3 + k];
  for (let j = 0; j < nv; j++) {
    for (let i = 0; i < nu; i++) {
      const i0 = Math.max(0, i - 1);
      const i1 = Math.min(nu - 1, i + 1);
      const j0 = Math.max(0, j - 1);
      const j1 = Math.min(nv - 1, j + 1);
      const ux = at(j, i1, 0) - at(j, i0, 0);
      const uy = at(j, i1, 1) - at(j, i0, 1);
      const uz = at(j, i1, 2) - at(j, i0, 2);
      const vx = at(j1, i, 0) - at(j0, i, 0);
      const vy = at(j1, i, 1) - at(j0, i, 1);
      const vz = at(j1, i, 2) - at(j0, i, 2);
      let nx = uy * vz - uz * vy;
      let ny = uz * vx - ux * vz;
      let nz = ux * vy - uy * vx;
      const len = Math.hypot(nx, ny, nz);
      if (len > 1e-12) {
        nx /= len;
        ny /= len;
        nz /= len;
      } else {
        nx = 0;
        ny = 0;
        nz = 1;
      }
      const o = (j * nu + i) * 3;
      normals[o] = nx;
      normals[o + 1] = ny;
      normals[o + 2] = nz;
      top[o] = surf[o] + nx * half;
      top[o + 1] = surf[o + 1] + ny * half;
      top[o + 2] = surf[o + 2] + nz * half;
      bottom[o] = surf[o] - nx * half;
      bottom[o + 1] = surf[o + 1] - ny * half;
      bottom[o + 2] = surf[o + 2] - nz * half;
    }
  }
  const vertices = new Float64Array(count * 6);
  const vNormals = new Float64Array(count * 6);
  vertices.set(top, 0);
  vertices.set(bottom, count * 3);
  vNormals.set(normals, 0);
  for (let k = 0; k < count * 3; k++) vNormals[count * 3 + k] = -normals[k];
  const T = (j, i) => j * nu + i;
  const B = (j, i) => count + j * nu + i;
  const faces = [];
  for (let j = 0; j < nv - 1; j++) {
    for (let i = 0; i < nu - 1; i++) {
      faces.push(T(j, i), T(j, i + 1), T(j + 1, i + 1));
      faces.push(T(j, i), T(j + 1, i + 1), T(j + 1, i));
      faces.push(B(j, i), B(j + 1, i + 1), B(j, i + 1));
      faces.push(B(j, i), B(j + 1, i), B(j + 1, i + 1));
    }
  }
  const loop = [];
  for (let i = 0; i < nu; i++) loop.push([0, i]);
  for (let j = 1; j < nv; j++) loop.push([j, nu - 1]);
  for (let i = nu - 2; i >= 0; i--) loop.push([nv - 1, i]);
  for (let j = nv - 2; j >= 1; j--) loop.push([j, 0]);
  for (let k = 0; k < loop.length; k++) {
    const [jp, ip] = loop[k];
    const [jq, iq] = loop[(k + 1) % loop.length];
    faces.push(T(jp, ip), B(jp, ip), B(jq, iq));
    faces.push(T(jp, ip), B(jq, iq), T(jq, iq));
  }
  return {
    grid: { nu, nv, y, xMid, halfWidth, top, bottom },
    mesh: { vertices, normals: vNormals, faces: Uint32Array.from(faces) },
    lengthMm,
    halfWidthMaxMm: halfWidthMax
  };
}

// src/leaf/section.ts
function measure(polygon) {
  let minA = Infinity;
  let maxA = -Infinity;
  let minB = Infinity;
  let maxB = -Infinity;
  for (const [a, b] of polygon) {
    if (a < minA) minA = a;
    if (a > maxA) maxA = a;
    if (b < minB) minB = b;
    if (b > maxB) maxB = b;
  }
  return {
    polygon,
    chordMm: polygon.length ? maxA - minA : 0,
    heightMm: polygon.length ? maxB - minB : 0
  };
}
function transverseSection(result, y0) {
  const { grid } = result;
  const { nu, nv, y, top, bottom, halfWidth } = grid;
  let station = y0;
  if (station === void 0) {
    let jWidest = 0;
    for (let j = 1; j < nv; j++) if (halfWidth[j] > halfWidth[jWidest]) jWidest = j;
    station = y[jWidest];
  }
  let j0 = 0;
  while (j0 < nv - 2 && y[j0 + 1] < station) j0++;
  const j1 = Math.min(nv - 1, j0 + 1);
  const span = y[j1] - y[j0];
  const w = span > 0 ? Math.min(1, Math.max(0, (station - y[j0]) / span)) : 0;
  const blend = (sheet, i, k) => {
    const a = sheet[(j0 * nu + i) * 3 + k];
    const b = sheet[(j1 * nu + i) * 3 + k];
    return a + (b - a) * w;
  };
  const upper = [];
  const lower = [];
  for (let i = 0; i < nu; i++) {
    upper.push([blend(top, i, 0), blend(top, i, 2)]);
    lower.push([blend(bottom, i, 0), blend(bottom, i, 2)]);
  }
  return measure([...upper, ...lower.reverse()]);
}
function longitudinalSection(result, x0 = 0) {
  const { grid } = result;
  const { nu, nv, top, bottom } = grid;
  const crossRow = (sheet, j) => {
    for (let i = 0; i < nu - 1; i++) {
      const xa = sheet[(j * nu + i) * 3];
      const xb = sheet[(j * nu + i + 1) * 3];
      if (xa === xb) continue;
      if (xa > x0 !== xb > x0) {
        const t = (x0 - xa) / (xb - xa);
        const ya = sheet[(j * nu + i) * 3 + 1];
        const yb = sheet[(j * nu + i + 1) * 3 + 1];
        const za = sheet[(j * nu + i) * 3 + 2];
        const zb = sheet[(j * nu + i + 1) * 3 + 2];
        return [ya + (yb - ya) * t, za + (zb - za) * t];
      }
    }
    return null;
  };
  const upper = [];
  const lower = [];
  for (let j = 0; j < nv; j++) {
    const t = crossRow(top, j);
    const b = crossRow(bottom, j);
    if (t && b) {
      upper.push(t);
      lower.push(b);
    }
  }
  if (upper.length < 2) return { polygon: [], chordMm: 0, heightMm: 0 };
  return measure([...upper, ...lower.reverse()]);
}

// src/leaf/assets/provisional-cotyledon.json
var provisional_cotyledon_default = {
  id: "provisional-cotyledon",
  label: "Cotyledon (provisional, 2.5 mm)",
  provisional: true,
  source: {
    series: "provisional",
    note: "Near-circular cotyledon blade, flat. Analytic stand-in. Generated by tools/leaf-assets/make_provisional.mjs \u2014 NOT a fit to imagery."
  },
  outline: [
    [
      0,
      0
    ],
    [
      0.0279,
      0.026
    ],
    [
      0.0592,
      0.0521
    ],
    [
      0.0916,
      0.0781
    ],
    [
      0.1243,
      0.1042
    ],
    [
      0.1571,
      0.1302
    ],
    [
      0.1899,
      0.1563
    ],
    [
      0.2225,
      0.1823
    ],
    [
      0.2548,
      0.2083
    ],
    [
      0.2868,
      0.2344
    ],
    [
      0.3183,
      0.2604
    ],
    [
      0.3494,
      0.2865
    ],
    [
      0.38,
      0.3125
    ],
    [
      0.41,
      0.3385
    ],
    [
      0.4395,
      0.3646
    ],
    [
      0.4683,
      0.3906
    ],
    [
      0.4966,
      0.4167
    ],
    [
      0.5242,
      0.4427
    ],
    [
      0.5511,
      0.4688
    ],
    [
      0.5774,
      0.4948
    ],
    [
      0.603,
      0.5208
    ],
    [
      0.6279,
      0.5469
    ],
    [
      0.652,
      0.5729
    ],
    [
      0.6754,
      0.599
    ],
    [
      0.6981,
      0.625
    ],
    [
      0.72,
      0.651
    ],
    [
      0.7412,
      0.6771
    ],
    [
      0.7616,
      0.7031
    ],
    [
      0.7812,
      0.7292
    ],
    [
      0.8,
      0.7552
    ],
    [
      0.818,
      0.7813
    ],
    [
      0.8352,
      0.8073
    ],
    [
      0.8515,
      0.8333
    ],
    [
      0.8671,
      0.8594
    ],
    [
      0.8818,
      0.8854
    ],
    [
      0.8957,
      0.9115
    ],
    [
      0.9087,
      0.9375
    ],
    [
      0.9209,
      0.9635
    ],
    [
      0.9323,
      0.9896
    ],
    [
      0.9428,
      1.0156
    ],
    [
      0.9524,
      1.0417
    ],
    [
      0.9611,
      1.0677
    ],
    [
      0.969,
      1.0938
    ],
    [
      0.976,
      1.1198
    ],
    [
      0.9821,
      1.1458
    ],
    [
      0.9873,
      1.1719
    ],
    [
      0.9917,
      1.1979
    ],
    [
      0.9951,
      1.224
    ],
    [
      0.9976,
      1.25
    ],
    [
      0.9992,
      1.276
    ],
    [
      1,
      1.3021
    ],
    [
      0.9998,
      1.3281
    ],
    [
      0.9987,
      1.3542
    ],
    [
      0.9966,
      1.3802
    ],
    [
      0.9937,
      1.4063
    ],
    [
      0.9898,
      1.4323
    ],
    [
      0.985,
      1.4583
    ],
    [
      0.9792,
      1.4844
    ],
    [
      0.9726,
      1.5104
    ],
    [
      0.9649,
      1.5365
    ],
    [
      0.9564,
      1.5625
    ],
    [
      0.9469,
      1.5885
    ],
    [
      0.9364,
      1.6146
    ],
    [
      0.925,
      1.6406
    ],
    [
      0.9127,
      1.6667
    ],
    [
      0.8993,
      1.6927
    ],
    [
      0.8851,
      1.7188
    ],
    [
      0.8698,
      1.7448
    ],
    [
      0.8536,
      1.7708
    ],
    [
      0.8365,
      1.7969
    ],
    [
      0.8184,
      1.8229
    ],
    [
      0.7993,
      1.849
    ],
    [
      0.7792,
      1.875
    ],
    [
      0.7581,
      1.901
    ],
    [
      0.7361,
      1.9271
    ],
    [
      0.7131,
      1.9531
    ],
    [
      0.6891,
      1.9792
    ],
    [
      0.6641,
      2.0052
    ],
    [
      0.6382,
      2.0313
    ],
    [
      0.6112,
      2.0573
    ],
    [
      0.5833,
      2.0833
    ],
    [
      0.5543,
      2.1094
    ],
    [
      0.5244,
      2.1354
    ],
    [
      0.4935,
      2.1615
    ],
    [
      0.4616,
      2.1875
    ],
    [
      0.4287,
      2.2135
    ],
    [
      0.3947,
      2.2396
    ],
    [
      0.3598,
      2.2656
    ],
    [
      0.3239,
      2.2917
    ],
    [
      0.2869,
      2.3177
    ],
    [
      0.249,
      2.3438
    ],
    [
      0.21,
      2.3698
    ],
    [
      0.1701,
      2.3958
    ],
    [
      0.1291,
      2.4219
    ],
    [
      0.0871,
      2.4479
    ],
    [
      0.044,
      2.474
    ],
    [
      0,
      2.5
    ],
    [
      -0.044,
      2.474
    ],
    [
      -0.0871,
      2.4479
    ],
    [
      -0.1291,
      2.4219
    ],
    [
      -0.1701,
      2.3958
    ],
    [
      -0.21,
      2.3698
    ],
    [
      -0.249,
      2.3438
    ],
    [
      -0.2869,
      2.3177
    ],
    [
      -0.3239,
      2.2917
    ],
    [
      -0.3598,
      2.2656
    ],
    [
      -0.3947,
      2.2396
    ],
    [
      -0.4287,
      2.2135
    ],
    [
      -0.4616,
      2.1875
    ],
    [
      -0.4935,
      2.1615
    ],
    [
      -0.5244,
      2.1354
    ],
    [
      -0.5543,
      2.1094
    ],
    [
      -0.5833,
      2.0833
    ],
    [
      -0.6112,
      2.0573
    ],
    [
      -0.6382,
      2.0313
    ],
    [
      -0.6641,
      2.0052
    ],
    [
      -0.6891,
      1.9792
    ],
    [
      -0.7131,
      1.9531
    ],
    [
      -0.7361,
      1.9271
    ],
    [
      -0.7581,
      1.901
    ],
    [
      -0.7792,
      1.875
    ],
    [
      -0.7993,
      1.849
    ],
    [
      -0.8184,
      1.8229
    ],
    [
      -0.8365,
      1.7969
    ],
    [
      -0.8536,
      1.7708
    ],
    [
      -0.8698,
      1.7448
    ],
    [
      -0.8851,
      1.7188
    ],
    [
      -0.8993,
      1.6927
    ],
    [
      -0.9127,
      1.6667
    ],
    [
      -0.925,
      1.6406
    ],
    [
      -0.9364,
      1.6146
    ],
    [
      -0.9469,
      1.5885
    ],
    [
      -0.9564,
      1.5625
    ],
    [
      -0.9649,
      1.5365
    ],
    [
      -0.9726,
      1.5104
    ],
    [
      -0.9792,
      1.4844
    ],
    [
      -0.985,
      1.4583
    ],
    [
      -0.9898,
      1.4323
    ],
    [
      -0.9937,
      1.4063
    ],
    [
      -0.9966,
      1.3802
    ],
    [
      -0.9987,
      1.3542
    ],
    [
      -0.9998,
      1.3281
    ],
    [
      -1,
      1.3021
    ],
    [
      -0.9992,
      1.276
    ],
    [
      -0.9976,
      1.25
    ],
    [
      -0.9951,
      1.224
    ],
    [
      -0.9917,
      1.1979
    ],
    [
      -0.9873,
      1.1719
    ],
    [
      -0.9821,
      1.1458
    ],
    [
      -0.976,
      1.1198
    ],
    [
      -0.969,
      1.0938
    ],
    [
      -0.9611,
      1.0677
    ],
    [
      -0.9524,
      1.0417
    ],
    [
      -0.9428,
      1.0156
    ],
    [
      -0.9323,
      0.9896
    ],
    [
      -0.9209,
      0.9635
    ],
    [
      -0.9087,
      0.9375
    ],
    [
      -0.8957,
      0.9115
    ],
    [
      -0.8818,
      0.8854
    ],
    [
      -0.8671,
      0.8594
    ],
    [
      -0.8515,
      0.8333
    ],
    [
      -0.8352,
      0.8073
    ],
    [
      -0.818,
      0.7813
    ],
    [
      -0.8,
      0.7552
    ],
    [
      -0.7812,
      0.7292
    ],
    [
      -0.7616,
      0.7031
    ],
    [
      -0.7412,
      0.6771
    ],
    [
      -0.72,
      0.651
    ],
    [
      -0.6981,
      0.625
    ],
    [
      -0.6754,
      0.599
    ],
    [
      -0.652,
      0.5729
    ],
    [
      -0.6279,
      0.5469
    ],
    [
      -0.603,
      0.5208
    ],
    [
      -0.5774,
      0.4948
    ],
    [
      -0.5511,
      0.4688
    ],
    [
      -0.5242,
      0.4427
    ],
    [
      -0.4966,
      0.4167
    ],
    [
      -0.4683,
      0.3906
    ],
    [
      -0.4395,
      0.3646
    ],
    [
      -0.41,
      0.3385
    ],
    [
      -0.38,
      0.3125
    ],
    [
      -0.3494,
      0.2865
    ],
    [
      -0.3183,
      0.2604
    ],
    [
      -0.2868,
      0.2344
    ],
    [
      -0.2548,
      0.2083
    ],
    [
      -0.2225,
      0.1823
    ],
    [
      -0.1899,
      0.1563
    ],
    [
      -0.1571,
      0.1302
    ],
    [
      -0.1243,
      0.1042
    ],
    [
      -0.0916,
      0.0781
    ],
    [
      -0.0592,
      0.0521
    ],
    [
      -0.0279,
      0.026
    ]
  ],
  midrib: [
    [
      0,
      0
    ],
    [
      0,
      0.2778
    ],
    [
      0,
      0.5556
    ],
    [
      0,
      0.8333
    ],
    [
      0,
      1.1111
    ],
    [
      0,
      1.3889
    ],
    [
      0,
      1.6667
    ],
    [
      0,
      1.9444
    ],
    [
      0,
      2.2222
    ],
    [
      0,
      2.5
    ]
  ],
  profile: [
    [
      0,
      0
    ],
    [
      1,
      0
    ]
  ],
  camber: [
    [
      -1,
      0
    ],
    [
      0,
      0
    ],
    [
      1,
      0
    ]
  ],
  thicknessMm: 0.15
};

// src/leaf/assets/provisional-juvenile.json
var provisional_juvenile_default = {
  id: "provisional-juvenile",
  label: "Juvenile leaf (provisional, 6 mm)",
  provisional: true,
  source: {
    series: "provisional",
    note: "First true leaf, faint upward cup. Analytic stand-in. Generated by tools/leaf-assets/make_provisional.mjs \u2014 NOT a fit to imagery."
  },
  outline: [
    [
      0,
      0
    ],
    [
      0.0124,
      0.05
    ],
    [
      0.0325,
      0.1
    ],
    [
      0.0569,
      0.15
    ],
    [
      0.0844,
      0.2
    ],
    [
      0.1143,
      0.25
    ],
    [
      0.1463,
      0.3
    ],
    [
      0.1799,
      0.35
    ],
    [
      0.215,
      0.4
    ],
    [
      0.2513,
      0.45
    ],
    [
      0.2886,
      0.5
    ],
    [
      0.3268,
      0.55
    ],
    [
      0.3658,
      0.6
    ],
    [
      0.4053,
      0.65
    ],
    [
      0.4455,
      0.7
    ],
    [
      0.486,
      0.75
    ],
    [
      0.5269,
      0.8
    ],
    [
      0.5681,
      0.85
    ],
    [
      0.6094,
      0.9
    ],
    [
      0.6509,
      0.95
    ],
    [
      0.6924,
      1
    ],
    [
      0.7339,
      1.05
    ],
    [
      0.7754,
      1.1
    ],
    [
      0.8168,
      1.15
    ],
    [
      0.858,
      1.2
    ],
    [
      0.899,
      1.25
    ],
    [
      0.9398,
      1.3
    ],
    [
      0.9802,
      1.35
    ],
    [
      1.0203,
      1.4
    ],
    [
      1.06,
      1.45
    ],
    [
      1.0993,
      1.5
    ],
    [
      1.1382,
      1.55
    ],
    [
      1.1766,
      1.6
    ],
    [
      1.2144,
      1.65
    ],
    [
      1.2517,
      1.7
    ],
    [
      1.2884,
      1.75
    ],
    [
      1.3244,
      1.8
    ],
    [
      1.3598,
      1.85
    ],
    [
      1.3945,
      1.9
    ],
    [
      1.4286,
      1.95
    ],
    [
      1.4618,
      2
    ],
    [
      1.4943,
      2.05
    ],
    [
      1.526,
      2.1
    ],
    [
      1.5569,
      2.15
    ],
    [
      1.587,
      2.2
    ],
    [
      1.6161,
      2.25
    ],
    [
      1.6444,
      2.3
    ],
    [
      1.6718,
      2.35
    ],
    [
      1.6982,
      2.4
    ],
    [
      1.7237,
      2.45
    ],
    [
      1.7482,
      2.5
    ],
    [
      1.7716,
      2.55
    ],
    [
      1.7941,
      2.6
    ],
    [
      1.8154,
      2.65
    ],
    [
      1.8358,
      2.7
    ],
    [
      1.855,
      2.75
    ],
    [
      1.8731,
      2.8
    ],
    [
      1.8901,
      2.85
    ],
    [
      1.906,
      2.9
    ],
    [
      1.9206,
      2.95
    ],
    [
      1.9341,
      3
    ],
    [
      1.9464,
      3.05
    ],
    [
      1.9575,
      3.1
    ],
    [
      1.9673,
      3.15
    ],
    [
      1.9759,
      3.2
    ],
    [
      1.9832,
      3.25
    ],
    [
      1.9892,
      3.3
    ],
    [
      1.9939,
      3.35
    ],
    [
      1.9973,
      3.4
    ],
    [
      1.9993,
      3.45
    ],
    [
      2,
      3.5
    ],
    [
      1.9993,
      3.55
    ],
    [
      1.9972,
      3.6
    ],
    [
      1.9938,
      3.65
    ],
    [
      1.9889,
      3.7
    ],
    [
      1.9825,
      3.75
    ],
    [
      1.9748,
      3.8
    ],
    [
      1.9655,
      3.85
    ],
    [
      1.9548,
      3.9
    ],
    [
      1.9426,
      3.95
    ],
    [
      1.9289,
      4
    ],
    [
      1.9137,
      4.05
    ],
    [
      1.8969,
      4.1
    ],
    [
      1.8786,
      4.15
    ],
    [
      1.8587,
      4.2
    ],
    [
      1.8373,
      4.25
    ],
    [
      1.8143,
      4.3
    ],
    [
      1.7896,
      4.35
    ],
    [
      1.7634,
      4.4
    ],
    [
      1.7355,
      4.45
    ],
    [
      1.706,
      4.5
    ],
    [
      1.6749,
      4.55
    ],
    [
      1.642,
      4.6
    ],
    [
      1.6075,
      4.65
    ],
    [
      1.5714,
      4.7
    ],
    [
      1.5335,
      4.75
    ],
    [
      1.4939,
      4.8
    ],
    [
      1.4525,
      4.85
    ],
    [
      1.4095,
      4.9
    ],
    [
      1.3647,
      4.95
    ],
    [
      1.3181,
      5
    ],
    [
      1.2698,
      5.05
    ],
    [
      1.2197,
      5.1
    ],
    [
      1.1677,
      5.15
    ],
    [
      1.114,
      5.2
    ],
    [
      1.0585,
      5.25
    ],
    [
      1.0011,
      5.3
    ],
    [
      0.9419,
      5.35
    ],
    [
      0.8808,
      5.4
    ],
    [
      0.8179,
      5.45
    ],
    [
      0.7531,
      5.5
    ],
    [
      0.6865,
      5.55
    ],
    [
      0.6179,
      5.6
    ],
    [
      0.5474,
      5.65
    ],
    [
      0.4751,
      5.7
    ],
    [
      0.4007,
      5.75
    ],
    [
      0.3245,
      5.8
    ],
    [
      0.2463,
      5.85
    ],
    [
      0.1662,
      5.9
    ],
    [
      0.0841,
      5.95
    ],
    [
      0,
      6
    ],
    [
      -0.0841,
      5.95
    ],
    [
      -0.1662,
      5.9
    ],
    [
      -0.2463,
      5.85
    ],
    [
      -0.3245,
      5.8
    ],
    [
      -0.4007,
      5.75
    ],
    [
      -0.4751,
      5.7
    ],
    [
      -0.5474,
      5.65
    ],
    [
      -0.6179,
      5.6
    ],
    [
      -0.6865,
      5.55
    ],
    [
      -0.7531,
      5.5
    ],
    [
      -0.8179,
      5.45
    ],
    [
      -0.8808,
      5.4
    ],
    [
      -0.9419,
      5.35
    ],
    [
      -1.0011,
      5.3
    ],
    [
      -1.0585,
      5.25
    ],
    [
      -1.114,
      5.2
    ],
    [
      -1.1677,
      5.15
    ],
    [
      -1.2197,
      5.1
    ],
    [
      -1.2698,
      5.05
    ],
    [
      -1.3181,
      5
    ],
    [
      -1.3647,
      4.95
    ],
    [
      -1.4095,
      4.9
    ],
    [
      -1.4525,
      4.85
    ],
    [
      -1.4939,
      4.8
    ],
    [
      -1.5335,
      4.75
    ],
    [
      -1.5714,
      4.7
    ],
    [
      -1.6075,
      4.65
    ],
    [
      -1.642,
      4.6
    ],
    [
      -1.6749,
      4.55
    ],
    [
      -1.706,
      4.5
    ],
    [
      -1.7355,
      4.45
    ],
    [
      -1.7634,
      4.4
    ],
    [
      -1.7896,
      4.35
    ],
    [
      -1.8143,
      4.3
    ],
    [
      -1.8373,
      4.25
    ],
    [
      -1.8587,
      4.2
    ],
    [
      -1.8786,
      4.15
    ],
    [
      -1.8969,
      4.1
    ],
    [
      -1.9137,
      4.05
    ],
    [
      -1.9289,
      4
    ],
    [
      -1.9426,
      3.95
    ],
    [
      -1.9548,
      3.9
    ],
    [
      -1.9655,
      3.85
    ],
    [
      -1.9748,
      3.8
    ],
    [
      -1.9825,
      3.75
    ],
    [
      -1.9889,
      3.7
    ],
    [
      -1.9938,
      3.65
    ],
    [
      -1.9972,
      3.6
    ],
    [
      -1.9993,
      3.55
    ],
    [
      -2,
      3.5
    ],
    [
      -1.9993,
      3.45
    ],
    [
      -1.9973,
      3.4
    ],
    [
      -1.9939,
      3.35
    ],
    [
      -1.9892,
      3.3
    ],
    [
      -1.9832,
      3.25
    ],
    [
      -1.9759,
      3.2
    ],
    [
      -1.9673,
      3.15
    ],
    [
      -1.9575,
      3.1
    ],
    [
      -1.9464,
      3.05
    ],
    [
      -1.9341,
      3
    ],
    [
      -1.9206,
      2.95
    ],
    [
      -1.906,
      2.9
    ],
    [
      -1.8901,
      2.85
    ],
    [
      -1.8731,
      2.8
    ],
    [
      -1.855,
      2.75
    ],
    [
      -1.8358,
      2.7
    ],
    [
      -1.8154,
      2.65
    ],
    [
      -1.7941,
      2.6
    ],
    [
      -1.7716,
      2.55
    ],
    [
      -1.7482,
      2.5
    ],
    [
      -1.7237,
      2.45
    ],
    [
      -1.6982,
      2.4
    ],
    [
      -1.6718,
      2.35
    ],
    [
      -1.6444,
      2.3
    ],
    [
      -1.6161,
      2.25
    ],
    [
      -1.587,
      2.2
    ],
    [
      -1.5569,
      2.15
    ],
    [
      -1.526,
      2.1
    ],
    [
      -1.4943,
      2.05
    ],
    [
      -1.4618,
      2
    ],
    [
      -1.4286,
      1.95
    ],
    [
      -1.3945,
      1.9
    ],
    [
      -1.3598,
      1.85
    ],
    [
      -1.3244,
      1.8
    ],
    [
      -1.2884,
      1.75
    ],
    [
      -1.2517,
      1.7
    ],
    [
      -1.2144,
      1.65
    ],
    [
      -1.1766,
      1.6
    ],
    [
      -1.1382,
      1.55
    ],
    [
      -1.0993,
      1.5
    ],
    [
      -1.06,
      1.45
    ],
    [
      -1.0203,
      1.4
    ],
    [
      -0.9802,
      1.35
    ],
    [
      -0.9398,
      1.3
    ],
    [
      -0.899,
      1.25
    ],
    [
      -0.858,
      1.2
    ],
    [
      -0.8168,
      1.15
    ],
    [
      -0.7754,
      1.1
    ],
    [
      -0.7339,
      1.05
    ],
    [
      -0.6924,
      1
    ],
    [
      -0.6509,
      0.95
    ],
    [
      -0.6094,
      0.9
    ],
    [
      -0.5681,
      0.85
    ],
    [
      -0.5269,
      0.8
    ],
    [
      -0.486,
      0.75
    ],
    [
      -0.4455,
      0.7
    ],
    [
      -0.4053,
      0.65
    ],
    [
      -0.3658,
      0.6
    ],
    [
      -0.3268,
      0.55
    ],
    [
      -0.2886,
      0.5
    ],
    [
      -0.2513,
      0.45
    ],
    [
      -0.215,
      0.4
    ],
    [
      -0.1799,
      0.35
    ],
    [
      -0.1463,
      0.3
    ],
    [
      -0.1143,
      0.25
    ],
    [
      -0.0844,
      0.2
    ],
    [
      -0.0569,
      0.15
    ],
    [
      -0.0325,
      0.1
    ],
    [
      -0.0124,
      0.05
    ]
  ],
  midrib: [
    [
      0,
      0
    ],
    [
      0,
      0.6667
    ],
    [
      0,
      1.3333
    ],
    [
      0,
      2
    ],
    [
      0,
      2.6667
    ],
    [
      0,
      3.3333
    ],
    [
      0,
      4
    ],
    [
      0,
      4.6667
    ],
    [
      0,
      5.3333
    ],
    [
      0,
      6
    ]
  ],
  profile: [
    [
      0,
      0
    ],
    [
      0.5,
      0.08
    ],
    [
      1,
      0.1
    ]
  ],
  camber: [
    [
      -1,
      0.06
    ],
    [
      -0.5,
      0.01
    ],
    [
      0,
      0
    ],
    [
      0.5,
      0.01
    ],
    [
      1,
      0.06
    ]
  ],
  thicknessMm: 0.18
};

// src/leaf/assets/provisional-mature-cupped.json
var provisional_mature_cupped_default = {
  id: "provisional-mature-cupped",
  label: "Mature rosette leaf, cupped (provisional, 15 mm)",
  provisional: true,
  source: {
    series: "provisional",
    note: "Same blade as provisional-mature-flat, arched along the midrib and cupped across it. The pair isolates deformation at fixed outline. Generated by tools/leaf-assets/make_provisional.mjs \u2014 NOT a fit to imagery."
  },
  outline: [
    [
      0,
      0
    ],
    [
      67e-4,
      0.0938
    ],
    [
      0.0203,
      0.1875
    ],
    [
      0.0388,
      0.2813
    ],
    [
      0.0614,
      0.375
    ],
    [
      0.0875,
      0.4688
    ],
    [
      0.1168,
      0.5625
    ],
    [
      0.149,
      0.6563
    ],
    [
      0.1838,
      0.75
    ],
    [
      0.221,
      0.8438
    ],
    [
      0.2603,
      0.9375
    ],
    [
      0.3017,
      1.0313
    ],
    [
      0.3447,
      1.125
    ],
    [
      0.3893,
      1.2188
    ],
    [
      0.4353,
      1.3125
    ],
    [
      0.4825,
      1.4063
    ],
    [
      0.5307,
      1.5
    ],
    [
      0.5797,
      1.5938
    ],
    [
      0.6295,
      1.6875
    ],
    [
      0.6797,
      1.7813
    ],
    [
      0.7304,
      1.875
    ],
    [
      0.7814,
      1.9688
    ],
    [
      0.8326,
      2.0625
    ],
    [
      0.8839,
      2.1563
    ],
    [
      0.9352,
      2.25
    ],
    [
      0.9865,
      2.3438
    ],
    [
      1.0378,
      2.4375
    ],
    [
      1.089,
      2.5313
    ],
    [
      1.1402,
      2.625
    ],
    [
      1.1914,
      2.7188
    ],
    [
      1.2426,
      2.8125
    ],
    [
      1.2938,
      2.9063
    ],
    [
      1.3453,
      3
    ],
    [
      1.3969,
      3.0938
    ],
    [
      1.4488,
      3.1875
    ],
    [
      1.5011,
      3.2813
    ],
    [
      1.5539,
      3.375
    ],
    [
      1.6073,
      3.4688
    ],
    [
      1.6614,
      3.5625
    ],
    [
      1.7162,
      3.6563
    ],
    [
      1.7719,
      3.75
    ],
    [
      1.8285,
      3.8437
    ],
    [
      1.886,
      3.9375
    ],
    [
      1.9445,
      4.0313
    ],
    [
      2.004,
      4.125
    ],
    [
      2.0646,
      4.2188
    ],
    [
      2.1261,
      4.3125
    ],
    [
      2.1886,
      4.4063
    ],
    [
      2.2519,
      4.5
    ],
    [
      2.3161,
      4.5938
    ],
    [
      2.3809,
      4.6875
    ],
    [
      2.4462,
      4.7813
    ],
    [
      2.5118,
      4.875
    ],
    [
      2.5777,
      4.9688
    ],
    [
      2.6435,
      5.0625
    ],
    [
      2.7091,
      5.1563
    ],
    [
      2.7743,
      5.25
    ],
    [
      2.8388,
      5.3438
    ],
    [
      2.9024,
      5.4375
    ],
    [
      2.9648,
      5.5313
    ],
    [
      3.0258,
      5.625
    ],
    [
      3.0852,
      5.7188
    ],
    [
      3.1427,
      5.8125
    ],
    [
      3.1982,
      5.9063
    ],
    [
      3.2513,
      6
    ],
    [
      3.3021,
      6.0938
    ],
    [
      3.3502,
      6.1875
    ],
    [
      3.3955,
      6.2813
    ],
    [
      3.4381,
      6.375
    ],
    [
      3.4777,
      6.4688
    ],
    [
      3.5144,
      6.5625
    ],
    [
      3.548,
      6.6563
    ],
    [
      3.5788,
      6.75
    ],
    [
      3.6066,
      6.8438
    ],
    [
      3.6316,
      6.9375
    ],
    [
      3.6538,
      7.0313
    ],
    [
      3.6735,
      7.125
    ],
    [
      3.6908,
      7.2188
    ],
    [
      3.7059,
      7.3125
    ],
    [
      3.7189,
      7.4063
    ],
    [
      3.7301,
      7.5
    ],
    [
      3.7398,
      7.5938
    ],
    [
      3.7481,
      7.6875
    ],
    [
      3.7553,
      7.7813
    ],
    [
      3.7617,
      7.875
    ],
    [
      3.7675,
      7.9688
    ],
    [
      3.773,
      8.0625
    ],
    [
      3.7783,
      8.1563
    ],
    [
      3.7836,
      8.25
    ],
    [
      3.7891,
      8.3438
    ],
    [
      3.7951,
      8.4375
    ],
    [
      3.8015,
      8.5313
    ],
    [
      3.8085,
      8.625
    ],
    [
      3.8161,
      8.7188
    ],
    [
      3.8244,
      8.8125
    ],
    [
      3.8333,
      8.9063
    ],
    [
      3.8429,
      9
    ],
    [
      3.8529,
      9.0938
    ],
    [
      3.8634,
      9.1875
    ],
    [
      3.8741,
      9.2813
    ],
    [
      3.8849,
      9.375
    ],
    [
      3.8955,
      9.4688
    ],
    [
      3.9057,
      9.5625
    ],
    [
      3.9153,
      9.6563
    ],
    [
      3.924,
      9.75
    ],
    [
      3.9314,
      9.8438
    ],
    [
      3.9374,
      9.9375
    ],
    [
      3.9415,
      10.0313
    ],
    [
      3.9435,
      10.125
    ],
    [
      3.943,
      10.2188
    ],
    [
      3.9399,
      10.3125
    ],
    [
      3.9337,
      10.4063
    ],
    [
      3.9244,
      10.5
    ],
    [
      3.9115,
      10.5938
    ],
    [
      3.895,
      10.6875
    ],
    [
      3.8746,
      10.7813
    ],
    [
      3.8502,
      10.875
    ],
    [
      3.8217,
      10.9688
    ],
    [
      3.7889,
      11.0625
    ],
    [
      3.752,
      11.1563
    ],
    [
      3.7108,
      11.25
    ],
    [
      3.6654,
      11.3438
    ],
    [
      3.6159,
      11.4375
    ],
    [
      3.5623,
      11.5313
    ],
    [
      3.5049,
      11.625
    ],
    [
      3.4437,
      11.7188
    ],
    [
      3.3789,
      11.8125
    ],
    [
      3.3107,
      11.9063
    ],
    [
      3.2394,
      12
    ],
    [
      3.1652,
      12.0938
    ],
    [
      3.0883,
      12.1875
    ],
    [
      3.0089,
      12.2813
    ],
    [
      2.9273,
      12.375
    ],
    [
      2.8438,
      12.4688
    ],
    [
      2.7584,
      12.5625
    ],
    [
      2.6715,
      12.6563
    ],
    [
      2.5832,
      12.75
    ],
    [
      2.4936,
      12.8438
    ],
    [
      2.4029,
      12.9375
    ],
    [
      2.3112,
      13.0313
    ],
    [
      2.2185,
      13.125
    ],
    [
      2.1249,
      13.2188
    ],
    [
      2.0303,
      13.3125
    ],
    [
      1.9347,
      13.4063
    ],
    [
      1.838,
      13.5
    ],
    [
      1.7401,
      13.5938
    ],
    [
      1.641,
      13.6875
    ],
    [
      1.5404,
      13.7813
    ],
    [
      1.4381,
      13.875
    ],
    [
      1.3339,
      13.9688
    ],
    [
      1.2277,
      14.0625
    ],
    [
      1.1191,
      14.1563
    ],
    [
      1.008,
      14.25
    ],
    [
      0.894,
      14.3438
    ],
    [
      0.777,
      14.4375
    ],
    [
      0.6568,
      14.5313
    ],
    [
      0.5331,
      14.625
    ],
    [
      0.4057,
      14.7188
    ],
    [
      0.2745,
      14.8125
    ],
    [
      0.1393,
      14.9063
    ],
    [
      0,
      15
    ],
    [
      -0.1393,
      14.9063
    ],
    [
      -0.2745,
      14.8125
    ],
    [
      -0.4057,
      14.7188
    ],
    [
      -0.5331,
      14.625
    ],
    [
      -0.6568,
      14.5313
    ],
    [
      -0.777,
      14.4375
    ],
    [
      -0.894,
      14.3438
    ],
    [
      -1.008,
      14.25
    ],
    [
      -1.1191,
      14.1563
    ],
    [
      -1.2277,
      14.0625
    ],
    [
      -1.3339,
      13.9688
    ],
    [
      -1.4381,
      13.875
    ],
    [
      -1.5404,
      13.7813
    ],
    [
      -1.641,
      13.6875
    ],
    [
      -1.7401,
      13.5938
    ],
    [
      -1.838,
      13.5
    ],
    [
      -1.9347,
      13.4063
    ],
    [
      -2.0303,
      13.3125
    ],
    [
      -2.1249,
      13.2188
    ],
    [
      -2.2185,
      13.125
    ],
    [
      -2.3112,
      13.0313
    ],
    [
      -2.4029,
      12.9375
    ],
    [
      -2.4936,
      12.8438
    ],
    [
      -2.5832,
      12.75
    ],
    [
      -2.6715,
      12.6563
    ],
    [
      -2.7584,
      12.5625
    ],
    [
      -2.8438,
      12.4688
    ],
    [
      -2.9273,
      12.375
    ],
    [
      -3.0089,
      12.2813
    ],
    [
      -3.0883,
      12.1875
    ],
    [
      -3.1652,
      12.0938
    ],
    [
      -3.2394,
      12
    ],
    [
      -3.3107,
      11.9063
    ],
    [
      -3.3789,
      11.8125
    ],
    [
      -3.4437,
      11.7188
    ],
    [
      -3.5049,
      11.625
    ],
    [
      -3.5623,
      11.5313
    ],
    [
      -3.6159,
      11.4375
    ],
    [
      -3.6654,
      11.3438
    ],
    [
      -3.7108,
      11.25
    ],
    [
      -3.752,
      11.1563
    ],
    [
      -3.7889,
      11.0625
    ],
    [
      -3.8217,
      10.9688
    ],
    [
      -3.8502,
      10.875
    ],
    [
      -3.8746,
      10.7813
    ],
    [
      -3.895,
      10.6875
    ],
    [
      -3.9115,
      10.5938
    ],
    [
      -3.9244,
      10.5
    ],
    [
      -3.9337,
      10.4063
    ],
    [
      -3.9399,
      10.3125
    ],
    [
      -3.943,
      10.2188
    ],
    [
      -3.9435,
      10.125
    ],
    [
      -3.9415,
      10.0313
    ],
    [
      -3.9374,
      9.9375
    ],
    [
      -3.9314,
      9.8438
    ],
    [
      -3.924,
      9.75
    ],
    [
      -3.9153,
      9.6563
    ],
    [
      -3.9057,
      9.5625
    ],
    [
      -3.8955,
      9.4688
    ],
    [
      -3.8849,
      9.375
    ],
    [
      -3.8741,
      9.2813
    ],
    [
      -3.8634,
      9.1875
    ],
    [
      -3.8529,
      9.0938
    ],
    [
      -3.8429,
      9
    ],
    [
      -3.8333,
      8.9063
    ],
    [
      -3.8244,
      8.8125
    ],
    [
      -3.8161,
      8.7188
    ],
    [
      -3.8085,
      8.625
    ],
    [
      -3.8015,
      8.5313
    ],
    [
      -3.7951,
      8.4375
    ],
    [
      -3.7891,
      8.3438
    ],
    [
      -3.7836,
      8.25
    ],
    [
      -3.7783,
      8.1563
    ],
    [
      -3.773,
      8.0625
    ],
    [
      -3.7675,
      7.9688
    ],
    [
      -3.7617,
      7.875
    ],
    [
      -3.7553,
      7.7813
    ],
    [
      -3.7481,
      7.6875
    ],
    [
      -3.7398,
      7.5938
    ],
    [
      -3.7301,
      7.5
    ],
    [
      -3.7189,
      7.4063
    ],
    [
      -3.7059,
      7.3125
    ],
    [
      -3.6908,
      7.2188
    ],
    [
      -3.6735,
      7.125
    ],
    [
      -3.6538,
      7.0313
    ],
    [
      -3.6316,
      6.9375
    ],
    [
      -3.6066,
      6.8438
    ],
    [
      -3.5788,
      6.75
    ],
    [
      -3.548,
      6.6563
    ],
    [
      -3.5144,
      6.5625
    ],
    [
      -3.4777,
      6.4688
    ],
    [
      -3.4381,
      6.375
    ],
    [
      -3.3955,
      6.2813
    ],
    [
      -3.3502,
      6.1875
    ],
    [
      -3.3021,
      6.0938
    ],
    [
      -3.2513,
      6
    ],
    [
      -3.1982,
      5.9063
    ],
    [
      -3.1427,
      5.8125
    ],
    [
      -3.0852,
      5.7188
    ],
    [
      -3.0258,
      5.625
    ],
    [
      -2.9648,
      5.5313
    ],
    [
      -2.9024,
      5.4375
    ],
    [
      -2.8388,
      5.3438
    ],
    [
      -2.7743,
      5.25
    ],
    [
      -2.7091,
      5.1563
    ],
    [
      -2.6435,
      5.0625
    ],
    [
      -2.5777,
      4.9688
    ],
    [
      -2.5118,
      4.875
    ],
    [
      -2.4462,
      4.7813
    ],
    [
      -2.3809,
      4.6875
    ],
    [
      -2.3161,
      4.5938
    ],
    [
      -2.2519,
      4.5
    ],
    [
      -2.1886,
      4.4063
    ],
    [
      -2.1261,
      4.3125
    ],
    [
      -2.0646,
      4.2188
    ],
    [
      -2.004,
      4.125
    ],
    [
      -1.9445,
      4.0313
    ],
    [
      -1.886,
      3.9375
    ],
    [
      -1.8285,
      3.8437
    ],
    [
      -1.7719,
      3.75
    ],
    [
      -1.7162,
      3.6563
    ],
    [
      -1.6614,
      3.5625
    ],
    [
      -1.6073,
      3.4688
    ],
    [
      -1.5539,
      3.375
    ],
    [
      -1.5011,
      3.2813
    ],
    [
      -1.4488,
      3.1875
    ],
    [
      -1.3969,
      3.0938
    ],
    [
      -1.3453,
      3
    ],
    [
      -1.2938,
      2.9063
    ],
    [
      -1.2426,
      2.8125
    ],
    [
      -1.1914,
      2.7188
    ],
    [
      -1.1402,
      2.625
    ],
    [
      -1.089,
      2.5313
    ],
    [
      -1.0378,
      2.4375
    ],
    [
      -0.9865,
      2.3438
    ],
    [
      -0.9352,
      2.25
    ],
    [
      -0.8839,
      2.1563
    ],
    [
      -0.8326,
      2.0625
    ],
    [
      -0.7814,
      1.9688
    ],
    [
      -0.7304,
      1.875
    ],
    [
      -0.6797,
      1.7813
    ],
    [
      -0.6295,
      1.6875
    ],
    [
      -0.5797,
      1.5938
    ],
    [
      -0.5307,
      1.5
    ],
    [
      -0.4825,
      1.4063
    ],
    [
      -0.4353,
      1.3125
    ],
    [
      -0.3893,
      1.2188
    ],
    [
      -0.3447,
      1.125
    ],
    [
      -0.3017,
      1.0313
    ],
    [
      -0.2603,
      0.9375
    ],
    [
      -0.221,
      0.8438
    ],
    [
      -0.1838,
      0.75
    ],
    [
      -0.149,
      0.6563
    ],
    [
      -0.1168,
      0.5625
    ],
    [
      -0.0875,
      0.4688
    ],
    [
      -0.0614,
      0.375
    ],
    [
      -0.0388,
      0.2813
    ],
    [
      -0.0203,
      0.1875
    ],
    [
      -67e-4,
      0.0938
    ]
  ],
  midrib: [
    [
      0,
      0
    ],
    [
      0,
      1.6667
    ],
    [
      0,
      3.3333
    ],
    [
      0,
      5
    ],
    [
      0,
      6.6667
    ],
    [
      0,
      8.3333
    ],
    [
      0,
      10
    ],
    [
      0,
      11.6667
    ],
    [
      0,
      13.3333
    ],
    [
      0,
      15
    ]
  ],
  profile: [
    [
      0,
      0
    ],
    [
      0.3,
      0.35
    ],
    [
      0.7,
      0.6
    ],
    [
      1,
      0.55
    ]
  ],
  camber: [
    [
      -1,
      0.9
    ],
    [
      -0.6,
      0.35
    ],
    [
      -0.25,
      0.06
    ],
    [
      0,
      0
    ],
    [
      0.25,
      0.06
    ],
    [
      0.6,
      0.35
    ],
    [
      1,
      0.9
    ]
  ],
  thicknessMm: 0.21
};

// src/leaf/assets/provisional-mature-flat.json
var provisional_mature_flat_default = {
  id: "provisional-mature-flat",
  label: "Mature rosette leaf, flat (provisional, 15 mm)",
  provisional: true,
  source: {
    series: "provisional",
    note: "Flat 15 mm blade sized to the T3 calibration (1.5 cm \u2194 52 cells). This is the like-for-like replacement for the stampEllipse reference geometry. Generated by tools/leaf-assets/make_provisional.mjs \u2014 NOT a fit to imagery."
  },
  outline: [
    [
      0,
      0
    ],
    [
      67e-4,
      0.0938
    ],
    [
      0.0203,
      0.1875
    ],
    [
      0.0388,
      0.2813
    ],
    [
      0.0614,
      0.375
    ],
    [
      0.0875,
      0.4688
    ],
    [
      0.1168,
      0.5625
    ],
    [
      0.149,
      0.6563
    ],
    [
      0.1838,
      0.75
    ],
    [
      0.221,
      0.8438
    ],
    [
      0.2603,
      0.9375
    ],
    [
      0.3017,
      1.0313
    ],
    [
      0.3447,
      1.125
    ],
    [
      0.3893,
      1.2188
    ],
    [
      0.4353,
      1.3125
    ],
    [
      0.4825,
      1.4063
    ],
    [
      0.5307,
      1.5
    ],
    [
      0.5797,
      1.5938
    ],
    [
      0.6295,
      1.6875
    ],
    [
      0.6797,
      1.7813
    ],
    [
      0.7304,
      1.875
    ],
    [
      0.7814,
      1.9688
    ],
    [
      0.8326,
      2.0625
    ],
    [
      0.8839,
      2.1563
    ],
    [
      0.9352,
      2.25
    ],
    [
      0.9865,
      2.3438
    ],
    [
      1.0378,
      2.4375
    ],
    [
      1.089,
      2.5313
    ],
    [
      1.1402,
      2.625
    ],
    [
      1.1914,
      2.7188
    ],
    [
      1.2426,
      2.8125
    ],
    [
      1.2938,
      2.9063
    ],
    [
      1.3453,
      3
    ],
    [
      1.3969,
      3.0938
    ],
    [
      1.4488,
      3.1875
    ],
    [
      1.5011,
      3.2813
    ],
    [
      1.5539,
      3.375
    ],
    [
      1.6073,
      3.4688
    ],
    [
      1.6614,
      3.5625
    ],
    [
      1.7162,
      3.6563
    ],
    [
      1.7719,
      3.75
    ],
    [
      1.8285,
      3.8437
    ],
    [
      1.886,
      3.9375
    ],
    [
      1.9445,
      4.0313
    ],
    [
      2.004,
      4.125
    ],
    [
      2.0646,
      4.2188
    ],
    [
      2.1261,
      4.3125
    ],
    [
      2.1886,
      4.4063
    ],
    [
      2.2519,
      4.5
    ],
    [
      2.3161,
      4.5938
    ],
    [
      2.3809,
      4.6875
    ],
    [
      2.4462,
      4.7813
    ],
    [
      2.5118,
      4.875
    ],
    [
      2.5777,
      4.9688
    ],
    [
      2.6435,
      5.0625
    ],
    [
      2.7091,
      5.1563
    ],
    [
      2.7743,
      5.25
    ],
    [
      2.8388,
      5.3438
    ],
    [
      2.9024,
      5.4375
    ],
    [
      2.9648,
      5.5313
    ],
    [
      3.0258,
      5.625
    ],
    [
      3.0852,
      5.7188
    ],
    [
      3.1427,
      5.8125
    ],
    [
      3.1982,
      5.9063
    ],
    [
      3.2513,
      6
    ],
    [
      3.3021,
      6.0938
    ],
    [
      3.3502,
      6.1875
    ],
    [
      3.3955,
      6.2813
    ],
    [
      3.4381,
      6.375
    ],
    [
      3.4777,
      6.4688
    ],
    [
      3.5144,
      6.5625
    ],
    [
      3.548,
      6.6563
    ],
    [
      3.5788,
      6.75
    ],
    [
      3.6066,
      6.8438
    ],
    [
      3.6316,
      6.9375
    ],
    [
      3.6538,
      7.0313
    ],
    [
      3.6735,
      7.125
    ],
    [
      3.6908,
      7.2188
    ],
    [
      3.7059,
      7.3125
    ],
    [
      3.7189,
      7.4063
    ],
    [
      3.7301,
      7.5
    ],
    [
      3.7398,
      7.5938
    ],
    [
      3.7481,
      7.6875
    ],
    [
      3.7553,
      7.7813
    ],
    [
      3.7617,
      7.875
    ],
    [
      3.7675,
      7.9688
    ],
    [
      3.773,
      8.0625
    ],
    [
      3.7783,
      8.1563
    ],
    [
      3.7836,
      8.25
    ],
    [
      3.7891,
      8.3438
    ],
    [
      3.7951,
      8.4375
    ],
    [
      3.8015,
      8.5313
    ],
    [
      3.8085,
      8.625
    ],
    [
      3.8161,
      8.7188
    ],
    [
      3.8244,
      8.8125
    ],
    [
      3.8333,
      8.9063
    ],
    [
      3.8429,
      9
    ],
    [
      3.8529,
      9.0938
    ],
    [
      3.8634,
      9.1875
    ],
    [
      3.8741,
      9.2813
    ],
    [
      3.8849,
      9.375
    ],
    [
      3.8955,
      9.4688
    ],
    [
      3.9057,
      9.5625
    ],
    [
      3.9153,
      9.6563
    ],
    [
      3.924,
      9.75
    ],
    [
      3.9314,
      9.8438
    ],
    [
      3.9374,
      9.9375
    ],
    [
      3.9415,
      10.0313
    ],
    [
      3.9435,
      10.125
    ],
    [
      3.943,
      10.2188
    ],
    [
      3.9399,
      10.3125
    ],
    [
      3.9337,
      10.4063
    ],
    [
      3.9244,
      10.5
    ],
    [
      3.9115,
      10.5938
    ],
    [
      3.895,
      10.6875
    ],
    [
      3.8746,
      10.7813
    ],
    [
      3.8502,
      10.875
    ],
    [
      3.8217,
      10.9688
    ],
    [
      3.7889,
      11.0625
    ],
    [
      3.752,
      11.1563
    ],
    [
      3.7108,
      11.25
    ],
    [
      3.6654,
      11.3438
    ],
    [
      3.6159,
      11.4375
    ],
    [
      3.5623,
      11.5313
    ],
    [
      3.5049,
      11.625
    ],
    [
      3.4437,
      11.7188
    ],
    [
      3.3789,
      11.8125
    ],
    [
      3.3107,
      11.9063
    ],
    [
      3.2394,
      12
    ],
    [
      3.1652,
      12.0938
    ],
    [
      3.0883,
      12.1875
    ],
    [
      3.0089,
      12.2813
    ],
    [
      2.9273,
      12.375
    ],
    [
      2.8438,
      12.4688
    ],
    [
      2.7584,
      12.5625
    ],
    [
      2.6715,
      12.6563
    ],
    [
      2.5832,
      12.75
    ],
    [
      2.4936,
      12.8438
    ],
    [
      2.4029,
      12.9375
    ],
    [
      2.3112,
      13.0313
    ],
    [
      2.2185,
      13.125
    ],
    [
      2.1249,
      13.2188
    ],
    [
      2.0303,
      13.3125
    ],
    [
      1.9347,
      13.4063
    ],
    [
      1.838,
      13.5
    ],
    [
      1.7401,
      13.5938
    ],
    [
      1.641,
      13.6875
    ],
    [
      1.5404,
      13.7813
    ],
    [
      1.4381,
      13.875
    ],
    [
      1.3339,
      13.9688
    ],
    [
      1.2277,
      14.0625
    ],
    [
      1.1191,
      14.1563
    ],
    [
      1.008,
      14.25
    ],
    [
      0.894,
      14.3438
    ],
    [
      0.777,
      14.4375
    ],
    [
      0.6568,
      14.5313
    ],
    [
      0.5331,
      14.625
    ],
    [
      0.4057,
      14.7188
    ],
    [
      0.2745,
      14.8125
    ],
    [
      0.1393,
      14.9063
    ],
    [
      0,
      15
    ],
    [
      -0.1393,
      14.9063
    ],
    [
      -0.2745,
      14.8125
    ],
    [
      -0.4057,
      14.7188
    ],
    [
      -0.5331,
      14.625
    ],
    [
      -0.6568,
      14.5313
    ],
    [
      -0.777,
      14.4375
    ],
    [
      -0.894,
      14.3438
    ],
    [
      -1.008,
      14.25
    ],
    [
      -1.1191,
      14.1563
    ],
    [
      -1.2277,
      14.0625
    ],
    [
      -1.3339,
      13.9688
    ],
    [
      -1.4381,
      13.875
    ],
    [
      -1.5404,
      13.7813
    ],
    [
      -1.641,
      13.6875
    ],
    [
      -1.7401,
      13.5938
    ],
    [
      -1.838,
      13.5
    ],
    [
      -1.9347,
      13.4063
    ],
    [
      -2.0303,
      13.3125
    ],
    [
      -2.1249,
      13.2188
    ],
    [
      -2.2185,
      13.125
    ],
    [
      -2.3112,
      13.0313
    ],
    [
      -2.4029,
      12.9375
    ],
    [
      -2.4936,
      12.8438
    ],
    [
      -2.5832,
      12.75
    ],
    [
      -2.6715,
      12.6563
    ],
    [
      -2.7584,
      12.5625
    ],
    [
      -2.8438,
      12.4688
    ],
    [
      -2.9273,
      12.375
    ],
    [
      -3.0089,
      12.2813
    ],
    [
      -3.0883,
      12.1875
    ],
    [
      -3.1652,
      12.0938
    ],
    [
      -3.2394,
      12
    ],
    [
      -3.3107,
      11.9063
    ],
    [
      -3.3789,
      11.8125
    ],
    [
      -3.4437,
      11.7188
    ],
    [
      -3.5049,
      11.625
    ],
    [
      -3.5623,
      11.5313
    ],
    [
      -3.6159,
      11.4375
    ],
    [
      -3.6654,
      11.3438
    ],
    [
      -3.7108,
      11.25
    ],
    [
      -3.752,
      11.1563
    ],
    [
      -3.7889,
      11.0625
    ],
    [
      -3.8217,
      10.9688
    ],
    [
      -3.8502,
      10.875
    ],
    [
      -3.8746,
      10.7813
    ],
    [
      -3.895,
      10.6875
    ],
    [
      -3.9115,
      10.5938
    ],
    [
      -3.9244,
      10.5
    ],
    [
      -3.9337,
      10.4063
    ],
    [
      -3.9399,
      10.3125
    ],
    [
      -3.943,
      10.2188
    ],
    [
      -3.9435,
      10.125
    ],
    [
      -3.9415,
      10.0313
    ],
    [
      -3.9374,
      9.9375
    ],
    [
      -3.9314,
      9.8438
    ],
    [
      -3.924,
      9.75
    ],
    [
      -3.9153,
      9.6563
    ],
    [
      -3.9057,
      9.5625
    ],
    [
      -3.8955,
      9.4688
    ],
    [
      -3.8849,
      9.375
    ],
    [
      -3.8741,
      9.2813
    ],
    [
      -3.8634,
      9.1875
    ],
    [
      -3.8529,
      9.0938
    ],
    [
      -3.8429,
      9
    ],
    [
      -3.8333,
      8.9063
    ],
    [
      -3.8244,
      8.8125
    ],
    [
      -3.8161,
      8.7188
    ],
    [
      -3.8085,
      8.625
    ],
    [
      -3.8015,
      8.5313
    ],
    [
      -3.7951,
      8.4375
    ],
    [
      -3.7891,
      8.3438
    ],
    [
      -3.7836,
      8.25
    ],
    [
      -3.7783,
      8.1563
    ],
    [
      -3.773,
      8.0625
    ],
    [
      -3.7675,
      7.9688
    ],
    [
      -3.7617,
      7.875
    ],
    [
      -3.7553,
      7.7813
    ],
    [
      -3.7481,
      7.6875
    ],
    [
      -3.7398,
      7.5938
    ],
    [
      -3.7301,
      7.5
    ],
    [
      -3.7189,
      7.4063
    ],
    [
      -3.7059,
      7.3125
    ],
    [
      -3.6908,
      7.2188
    ],
    [
      -3.6735,
      7.125
    ],
    [
      -3.6538,
      7.0313
    ],
    [
      -3.6316,
      6.9375
    ],
    [
      -3.6066,
      6.8438
    ],
    [
      -3.5788,
      6.75
    ],
    [
      -3.548,
      6.6563
    ],
    [
      -3.5144,
      6.5625
    ],
    [
      -3.4777,
      6.4688
    ],
    [
      -3.4381,
      6.375
    ],
    [
      -3.3955,
      6.2813
    ],
    [
      -3.3502,
      6.1875
    ],
    [
      -3.3021,
      6.0938
    ],
    [
      -3.2513,
      6
    ],
    [
      -3.1982,
      5.9063
    ],
    [
      -3.1427,
      5.8125
    ],
    [
      -3.0852,
      5.7188
    ],
    [
      -3.0258,
      5.625
    ],
    [
      -2.9648,
      5.5313
    ],
    [
      -2.9024,
      5.4375
    ],
    [
      -2.8388,
      5.3438
    ],
    [
      -2.7743,
      5.25
    ],
    [
      -2.7091,
      5.1563
    ],
    [
      -2.6435,
      5.0625
    ],
    [
      -2.5777,
      4.9688
    ],
    [
      -2.5118,
      4.875
    ],
    [
      -2.4462,
      4.7813
    ],
    [
      -2.3809,
      4.6875
    ],
    [
      -2.3161,
      4.5938
    ],
    [
      -2.2519,
      4.5
    ],
    [
      -2.1886,
      4.4063
    ],
    [
      -2.1261,
      4.3125
    ],
    [
      -2.0646,
      4.2188
    ],
    [
      -2.004,
      4.125
    ],
    [
      -1.9445,
      4.0313
    ],
    [
      -1.886,
      3.9375
    ],
    [
      -1.8285,
      3.8437
    ],
    [
      -1.7719,
      3.75
    ],
    [
      -1.7162,
      3.6563
    ],
    [
      -1.6614,
      3.5625
    ],
    [
      -1.6073,
      3.4688
    ],
    [
      -1.5539,
      3.375
    ],
    [
      -1.5011,
      3.2813
    ],
    [
      -1.4488,
      3.1875
    ],
    [
      -1.3969,
      3.0938
    ],
    [
      -1.3453,
      3
    ],
    [
      -1.2938,
      2.9063
    ],
    [
      -1.2426,
      2.8125
    ],
    [
      -1.1914,
      2.7188
    ],
    [
      -1.1402,
      2.625
    ],
    [
      -1.089,
      2.5313
    ],
    [
      -1.0378,
      2.4375
    ],
    [
      -0.9865,
      2.3438
    ],
    [
      -0.9352,
      2.25
    ],
    [
      -0.8839,
      2.1563
    ],
    [
      -0.8326,
      2.0625
    ],
    [
      -0.7814,
      1.9688
    ],
    [
      -0.7304,
      1.875
    ],
    [
      -0.6797,
      1.7813
    ],
    [
      -0.6295,
      1.6875
    ],
    [
      -0.5797,
      1.5938
    ],
    [
      -0.5307,
      1.5
    ],
    [
      -0.4825,
      1.4063
    ],
    [
      -0.4353,
      1.3125
    ],
    [
      -0.3893,
      1.2188
    ],
    [
      -0.3447,
      1.125
    ],
    [
      -0.3017,
      1.0313
    ],
    [
      -0.2603,
      0.9375
    ],
    [
      -0.221,
      0.8438
    ],
    [
      -0.1838,
      0.75
    ],
    [
      -0.149,
      0.6563
    ],
    [
      -0.1168,
      0.5625
    ],
    [
      -0.0875,
      0.4688
    ],
    [
      -0.0614,
      0.375
    ],
    [
      -0.0388,
      0.2813
    ],
    [
      -0.0203,
      0.1875
    ],
    [
      -67e-4,
      0.0938
    ]
  ],
  midrib: [
    [
      0,
      0
    ],
    [
      0,
      1.6667
    ],
    [
      0,
      3.3333
    ],
    [
      0,
      5
    ],
    [
      0,
      6.6667
    ],
    [
      0,
      8.3333
    ],
    [
      0,
      10
    ],
    [
      0,
      11.6667
    ],
    [
      0,
      13.3333
    ],
    [
      0,
      15
    ]
  ],
  profile: [
    [
      0,
      0
    ],
    [
      1,
      0
    ]
  ],
  camber: [
    [
      -1,
      0
    ],
    [
      0,
      0
    ],
    [
      1,
      0
    ]
  ],
  thicknessMm: 0.21
};

// src/leaf/registry.ts
var RAW = [provisional_mature_flat_default, provisional_mature_cupped_default, provisional_juvenile_default, provisional_cotyledon_default];
function widestChord(outline, samples = 200) {
  const [, minY, , maxY] = bounds(outline);
  const span = maxY - minY;
  if (!(span > 0)) return 0;
  let widest = 0;
  for (let k = 1; k < samples; k++) {
    const xs = crossingsAtY(outline, minY + span * k / samples);
    if (xs.length >= 2) {
      const w = xs[xs.length - 1] - xs[0];
      if (w > widest) widest = w;
    }
  }
  return widest;
}
function computeTraits(asset) {
  const outline = ensureCCW(asset.outline);
  const [, minY, , maxY] = bounds(outline);
  const lengthMm = polylineLength(asset.midrib) || maxY - minY;
  const widthMm = widestChord(outline);
  return {
    areaMm2: polygonArea(outline),
    lengthMm,
    widthMm,
    aspectRatio: widthMm > 0 ? lengthMm / widthMm : 0,
    perimeterMm: polygonPerimeter(outline)
  };
}
function validate(asset) {
  const where = `leaf asset '${asset.id ?? "(no id)"}'`;
  if (!asset.id) throw new Error(`${where}: missing id`);
  if (!Array.isArray(asset.outline) || asset.outline.length < 3) {
    throw new Error(`${where}: outline needs at least 3 points`);
  }
  if (!Array.isArray(asset.midrib) || asset.midrib.length < 2) {
    throw new Error(`${where}: midrib needs at least 2 points`);
  }
  if (!(asset.thicknessMm > 0)) throw new Error(`${where}: thicknessMm must be positive`);
  const normalised = { ...asset, outline: ensureCCW(asset.outline) };
  return { ...normalised, traits: computeTraits(normalised) };
}
var LEAF_ASSETS = RAW.map(validate);
function getLeafAsset(id) {
  const found = LEAF_ASSETS.find((a) => a.id === id);
  if (!found) {
    throw new Error(`unknown leaf asset '${id}'; have: ${LEAF_ASSETS.map((a) => a.id).join(", ")}`);
  }
  return found;
}
var ALL_PROVISIONAL = LEAF_ASSETS.every((a) => a.provisional);

// src/scenarios/scenarios.ts
function buildCavity() {
  const N = 130;
  const U = 0.1;
  const L = N - 2;
  const Re = 100;
  const nu = viscosityFromRe(U, L, Re);
  const tau = tauFromViscosity(nu);
  const fluid = new LBMFluid(N, N, tau);
  fluid.setEquilibrium(1, 0, 0);
  for (let x = 0; x < N; x++) {
    for (let y = 0; y < N; y++) {
      const c = fluid.index(x, y);
      const border = x === 0 || x === N - 1 || y === 0 || y === N - 1;
      if (!border) continue;
      fluid.solid[c] = 1;
      if (y === N - 1) {
        fluid.uwx[c] = U;
        fluid.uwy[c] = 0;
      }
    }
  }
  const sampleU = (yNorm) => {
    const xc = Math.floor(N / 2);
    const y = Math.min(N - 2, Math.max(1, Math.round(1 + yNorm * (N - 3))));
    return fluid.ux[fluid.index(xc, y)] / U;
  };
  const sampleV = (xNorm) => {
    const yc = Math.floor(N / 2);
    const x = Math.min(N - 2, Math.max(1, Math.round(1 + xNorm * (N - 3))));
    return fluid.uy[fluid.index(x, yc)] / U;
  };
  return {
    id: "cavity",
    label: "Lid-driven cavity (Re 100)",
    fluid,
    renderMode: "speed",
    renderScale: U,
    postStream() {
    },
    onAfterStep() {
    },
    diagnostics(step) {
      fluid.computeMacroscopic();
      const dl = computeDimensionless({ u: U, length: L, nu });
      const err = ghiaL2Error(sampleU, sampleV);
      const converged = step > 6e3;
      const gateStatus = !converged ? "pending" : err < 0.05 ? "pass" : "warn";
      return [
        { label: "Reynolds", value: dl.Re.toFixed(0) },
        { label: "tau", value: tau.toFixed(3) },
        { label: "steps", value: step.toLocaleString() },
        {
          label: "Ghia L2 error",
          value: err.toFixed(4),
          status: gateStatus
        },
        {
          label: "Gate 1 (cavity)",
          value: converged ? err < 0.05 ? "PASS (<0.05)" : "off target" : "converging\u2026",
          status: gateStatus
        }
      ];
    }
  };
}
function buildCylinder() {
  const nx = 420;
  const ny = 120;
  const U = 0.05;
  const D = 20;
  const Re = 100;
  const nu = viscosityFromRe(U, D, Re);
  const tau = tauFromViscosity(nu);
  const fluid = new LBMFluid(nx, ny, tau);
  fluid.setEquilibrium(1, U, 0);
  for (let x = 0; x < nx; x++) {
    fluid.solid[fluid.index(x, 0)] = 1;
    fluid.solid[fluid.index(x, ny - 1)] = 1;
  }
  const cx = 90;
  const cy = ny / 2 - 2;
  const r = D / 2;
  for (let x = 0; x < nx; x++) {
    for (let y = 0; y < ny; y++) {
      if ((x - cx) ** 2 + (y - cy) ** 2 <= r * r) {
        fluid.solid[fluid.index(x, y)] = 1;
      }
    }
  }
  const probe = new StrouhalProbe(4e3);
  const probeX = Math.min(nx - 2, cx + 3 * D);
  const probeY = Math.round(cy);
  const probeCell = fluid.index(probeX, probeY);
  const applyInletOutlet = () => {
    for (let y = 1; y < ny - 1; y++) {
      const c = fluid.index(0, y);
      if (fluid.solid[c]) continue;
      const base = c * 9;
      for (let i = 0; i < 9; i++) fluid.f[base + i] = feq(i, 1, U, 0);
    }
    for (let y = 1; y < ny - 1; y++) {
      const c = fluid.index(nx - 1, y);
      if (fluid.solid[c]) continue;
      const src = fluid.index(nx - 2, y) * 9;
      const dst = c * 9;
      for (let i = 0; i < 9; i++) fluid.f[dst + i] = fluid.f[src + i];
    }
  };
  return {
    id: "cylinder",
    label: "Flow past cylinder (Re 100)",
    fluid,
    renderMode: "vorticity",
    renderScale: 0.02,
    postStream() {
      applyInletOutlet();
    },
    onAfterStep(step) {
      const CXa = [0, 1, 0, -1, 0, 1, -1, -1, 1];
      const CYa = [0, 0, 1, 0, -1, 1, 1, -1, -1];
      if (step < 600) {
        const py = Math.round(cy) + 3;
        const x0 = cx + Math.round(r) + 2;
        for (let x = x0; x < x0 + 4; x++) {
          const cc = fluid.index(x, py);
          if (fluid.solid[cc]) continue;
          const b = cc * 9;
          let rho2 = 0;
          let mx = 0;
          let my2 = 0;
          for (let i = 0; i < 9; i++) {
            const fi = fluid.f[b + i];
            rho2 += fi;
            mx += CXa[i] * fi;
            my2 += CYa[i] * fi;
          }
          const ux = mx / rho2;
          const uy = my2 / rho2 + 0.2 * U;
          for (let i = 0; i < 9; i++) fluid.f[b + i] = feq(i, rho2, ux, uy);
        }
      }
      const base = probeCell * 9;
      let rho = 0;
      let my = 0;
      for (let i = 0; i < 9; i++) {
        const fi = fluid.f[base + i];
        rho += fi;
        my += CYa[i] * fi;
      }
      probe.record(my / rho);
    },
    diagnostics(step) {
      const dl = computeDimensionless({ u: U, length: D, nu });
      const st = probe.strouhal(D, U);
      const cycles = probe.cycles();
      const inRange = st !== null && st > 0.16 && st < 0.22;
      const status = st === null ? "pending" : inRange ? "pass" : "warn";
      return [
        { label: "Reynolds", value: dl.Re.toFixed(0) },
        { label: "tau", value: tau.toFixed(3) },
        { label: "steps", value: step.toLocaleString() },
        { label: "shedding cycles", value: cycles.toString() },
        {
          label: "Strouhal",
          value: st === null ? "measuring\u2026" : st.toFixed(3),
          status
        },
        {
          label: "Gate 2 (cylinder)",
          value: st === null ? "waiting for wake\u2026" : inRange ? "PASS (0.16\u20130.22, confined)" : "off target",
          status
        }
      ];
    }
  };
}
function buildNaturalConvection() {
  const N = 66;
  const H = N - 2;
  const nu = 0.02;
  const Pr = 0.71;
  const kappa = nu / Pr;
  const Ra = 1e4;
  const dT = 1;
  const A = Ra * nu * kappa / (dT * H ** 3);
  const tauF = tauFromViscosity(nu);
  const fluid = new LBMFluid(N, N, tauF);
  fluid.setEquilibrium(1, 0, 0);
  fluid.enableForcing();
  const temp = new ScalarField(fluid, kappa, 0.5);
  for (let x = 0; x < N; x++) {
    for (let y = 0; y < N; y++) {
      const c = fluid.index(x, y);
      const border = x === 0 || x === N - 1 || y === 0 || y === N - 1;
      if (!border) continue;
      fluid.solid[c] = 1;
      if (x === 0) temp.setDirichlet(c, 1);
      else if (x === N - 1) temp.setDirichlet(c, 0);
    }
  }
  const nusseltHot = () => {
    let sum = 0;
    let n = 0;
    for (let y = 1; y < N - 1; y++) {
      const t1 = temp.C[fluid.index(1, y)];
      sum += (1 - t1) * H / (0.5 * dT);
      n++;
    }
    return sum / n;
  };
  return {
    id: "natconv",
    label: "Natural-convection cavity (Ra 1e4)",
    fluid,
    scalarField: temp,
    renderMode: "scalar",
    renderScale: 1,
    postStream() {
    },
    onAfterStep() {
      temp.step();
      applyBoussinesqForce(fluid, 0, A, [{ field: temp, beta: 1, ref: 0.5 }]);
    },
    diagnostics(step) {
      const Nu = nusseltHot();
      const converged = step > 4e4;
      const inRange = Nu > 1.9 && Nu < 2.5;
      const status = !converged ? "pending" : inRange ? "pass" : "warn";
      return [
        { label: "Rayleigh", value: Ra.toExponential(0) },
        { label: "Prandtl", value: Pr.toFixed(2) },
        { label: "steps", value: step.toLocaleString() },
        { label: "Nu (hot wall)", value: Nu.toFixed(3), status },
        {
          label: "Gate 3 (buoyancy)",
          value: !converged ? "converging\u2026" : inRange ? "PASS (Nu\u22482.24)" : "off target",
          status
        }
      ];
    }
  };
}
function buildDiffusion() {
  const nx = 200;
  const ny = 8;
  const D = 1 / 6;
  const x0 = nx / 2;
  const iface = x0 - 0.5;
  const fluid = new LBMFluid(nx, ny, tauFromViscosity(1 / 6));
  fluid.setEquilibrium(1, 0, 0);
  const conc = new ScalarField(fluid, D, 0);
  conc.initField((x) => x < x0 ? 1 : 0);
  let latchedPass = false;
  const l2Error = (step) => {
    const yrow = ny >> 1;
    let sq = 0;
    let n = 0;
    for (let x = 1; x < nx - 1; x++) {
      const analytic = diffusionStep(x, step, D, iface, 1, 0);
      const d = conc.C[fluid.index(x, yrow)] - analytic;
      sq += d * d;
      n++;
    }
    return Math.sqrt(sq / n);
  };
  return {
    id: "diffusion",
    label: "Pure diffusion (erfc limit)",
    fluid,
    scalarField: conc,
    renderMode: "scalar",
    renderScale: 1,
    postStream() {
    },
    onAfterStep() {
      conc.step();
    },
    diagnostics(step) {
      const err = l2Error(step);
      const inWindow = step >= 400 && step <= 1600;
      if (inWindow && err < 0.02) latchedPass = true;
      const status = step < 400 ? "pending" : latchedPass ? "pass" : inWindow ? "warn" : "warn";
      return [
        { label: "diffusivity D", value: D.toFixed(3) },
        { label: "steps", value: step.toLocaleString() },
        { label: "erfc L2 error", value: err.toFixed(4), status },
        {
          label: "Gate 4 (diffusion)",
          value: step < 400 ? "diffusing\u2026" : latchedPass ? "PASS (<0.02)" : step > 1600 ? "window closed" : "measuring\u2026",
          status
        }
      ];
    }
  };
}
function chamberWalls(fluid, sp) {
  const { nx, ny } = fluid;
  for (let x = 0; x < nx; x++) {
    for (let y = 0; y < ny; y++) {
      if (x !== 0 && x !== nx - 1 && y !== 0 && y !== ny - 1) continue;
      const c = fluid.index(x, y);
      fluid.solid[c] = 1;
      sp.co2.setDirichlet(c, 0);
      sp.o2.setDirichlet(c, 0);
      sp.h2o.setDirichlet(c, 0);
    }
  }
}
function stampEllipse(fluid, isLeaf, cx, cy, a, b, angle) {
  const ca = Math.cos(angle);
  const sa = Math.sin(angle);
  const R = Math.ceil(Math.max(a, b)) + 1;
  const x0 = Math.max(0, Math.floor(cx - R));
  const x1 = Math.min(fluid.nx - 1, Math.ceil(cx + R));
  const y0 = Math.max(0, Math.floor(cy - R));
  const y1 = Math.min(fluid.ny - 1, Math.ceil(cy + R));
  for (let x = x0; x <= x1; x++) {
    for (let y = y0; y <= y1; y++) {
      const dx = x - cx;
      const dy = y - cy;
      const lx = dx * ca + dy * sa;
      const ly = -dx * sa + dy * ca;
      if ((lx / a) ** 2 + (ly / b) ** 2 <= 1) {
        const c = fluid.index(x, y);
        fluid.solid[c] = 1;
        isLeaf[c] = 1;
      }
    }
  }
}
function stampPolygon(fluid, isLeaf, poly, cx, cy, cellsPerMm) {
  let minA = Infinity;
  let maxA = -Infinity;
  let minB = Infinity;
  let maxB = -Infinity;
  for (const [a, b] of poly) {
    if (a < minA) minA = a;
    if (a > maxA) maxA = a;
    if (b < minB) minB = b;
    if (b > maxB) maxB = b;
  }
  const midA = (minA + maxA) / 2;
  const midB = (minB + maxB) / 2;
  const x0 = Math.max(0, Math.floor(cx + (minA - midA) * cellsPerMm) - 1);
  const x1 = Math.min(fluid.nx - 1, Math.ceil(cx + (maxA - midA) * cellsPerMm) + 1);
  const y0 = Math.max(0, Math.floor(cy + (minB - midB) * cellsPerMm) - 1);
  const y1 = Math.min(fluid.ny - 1, Math.ceil(cy + (maxB - midB) * cellsPerMm) + 1);
  let filled = 0;
  for (let x = x0; x <= x1; x++) {
    for (let y = y0; y <= y1; y++) {
      const a = midA + (x - cx) / cellsPerMm;
      const b = midB + (y - cy) / cellsPerMm;
      if (polygonSDF(poly, a, b) <= 0) {
        const c = fluid.index(x, y);
        fluid.solid[c] = 1;
        isLeaf[c] = 1;
        filled++;
      }
    }
  }
  return filled;
}
function leafGeometry(fluid, sp) {
  chamberWalls(fluid, sp);
  const isLeaf = new Uint8Array(fluid.size);
  stampEllipse(fluid, isLeaf, fluid.nx / 2, fluid.ny / 2, 26, 4, 0);
  return { isLeaf, charLen: 52 };
}
function dishLeafGeometry(cfg) {
  return (fluid, sp) => {
    chamberWalls(fluid, sp);
    for (let y = 1; y <= cfg.agarCells; y++) {
      for (let x = 1; x < fluid.nx - 1; x++) fluid.solid[fluid.index(x, y)] = 1;
    }
    const isLeaf = new Uint8Array(fluid.size);
    const b = 4;
    const cy = cfg.agarCells + 1 + cfg.gapCells + b - 0.5;
    if (cy + b >= fluid.ny - 1) throw new Error(`dishLeafGeometry: leaf (top at ${cy + b}) does not fit under the lid (ny ${fluid.ny})`);
    stampEllipse(fluid, isLeaf, fluid.nx / 2, cy, 26, b, 0);
    return { isLeaf, charLen: 52 };
  };
}
function rosetteGeometry(fluid, sp) {
  chamberWalls(fluid, sp);
  const isLeaf = new Uint8Array(fluid.size);
  const cx = fluid.nx / 2;
  const cy = fluid.ny / 2;
  const a = 20;
  const b = 3.4;
  const off = 15;
  const anglesDeg = [22, 8, -8, -22, 158, 172, 188, 202];
  for (const deg of anglesDeg) {
    const ang = deg * Math.PI / 180;
    stampEllipse(fluid, isLeaf, cx + off * Math.cos(ang), cy + off * Math.sin(ang), a, b, ang);
  }
  return { isLeaf, charLen: 70 };
}
function canopyGeometry(fluid, sp) {
  const { nx, ny } = fluid;
  const isLeaf = new Uint8Array(fluid.size);
  for (let x = 0; x < nx; x++) {
    for (let y = 0; y < ny; y++) {
      if (x !== 0 && x !== nx - 1 && y !== ny - 1) continue;
      const c = fluid.index(x, y);
      fluid.solid[c] = 1;
      sp.co2.setDirichlet(c, 0);
      sp.o2.setDirichlet(c, 0);
      sp.h2o.setDirichlet(c, 0);
    }
  }
  for (let x = 0; x < nx; x++) fluid.solid[fluid.index(x, 0)] = 1;
  const spacing = 12;
  const halfW = 1.7;
  const height = 13;
  for (let bx = spacing; bx < nx - spacing + 1; bx += spacing) {
    stampEllipse(fluid, isLeaf, bx, 1 + height, halfW, height, 0);
  }
  return { isLeaf, charLen: 26 };
}
var U_STAR_CM_S = 166;
var DX_M = 288e-6;
var CELLS_PER_MM = 1 / (DX_M * 1e3);
var REF_THICKNESS_CELLS = 8;
function assetLeafGeometry(cfg) {
  return (fluid, sp) => {
    chamberWalls(fluid, sp);
    const isLeaf = new Uint8Array(fluid.size);
    const asset = getLeafAsset(cfg.assetId);
    const cut = cfg.cut ?? "longitudinal";
    const thicknessCells = cfg.thicknessCells ?? REF_THICKNESS_CELLS;
    const thicknessScale = thicknessCells / CELLS_PER_MM / asset.thicknessMm;
    const deform = { ...NEUTRAL_DEFORM, ...cfg.deform, thicknessScale };
    const posed = loft(asset, deform);
    const section = cut === "transverse" ? transverseSection(posed) : longitudinalSection(posed);
    if (section.polygon.length < 3) {
      throw new Error(`leaf asset '${cfg.assetId}': ${cut} cut produced no section`);
    }
    const filled = stampPolygon(
      fluid,
      isLeaf,
      section.polygon,
      fluid.nx / 2,
      fluid.ny / 2,
      CELLS_PER_MM
    );
    if (filled < 16) {
      throw new Error(
        `leaf asset '${cfg.assetId}': ${cut} cut filled only ${filled} cells \u2014 too thin to resolve`
      );
    }
    const chordCells = section.chordMm * CELLS_PER_MM;
    const heightCells = section.heightMm * CELLS_PER_MM;
    return {
      isLeaf,
      charLen: chordCells,
      notes: [
        { label: "leaf geometry", value: `${asset.id}${asset.provisional ? " (provisional)" : ""}` },
        { label: "section", value: `${cut} \xB7 ${chordCells.toFixed(0)} \xD7 ${heightCells.toFixed(1)} cells` },
        {
          label: "blade thickness",
          value: `${thicknessCells} cells (${(thicknessCells / CELLS_PER_MM).toFixed(2)} mm; ${thicknessScale.toFixed(0)}\xD7 life size)`
        }
      ]
    };
  };
}
var S_CO2_BASE = 5e-4;
var GBL_REF = 1;
var DC_CO2_REF = 0.161;
var K_GBL = GBL_REF * DC_CO2_REF / S_CO2_BASE;
var D_CO2_PHYS = 16e-6;
var MOLAR_DENSITY = 41.6;
function makeLeafScene(cfg) {
  return () => {
    const { nx, ny } = cfg;
    const nu = 0.02;
    const r = cfg.refine ?? 1;
    const U = (cfg.forcedU ?? 0) / r;
    const fluid = new LBMFluid(nx, ny, tauFromViscosity(nu));
    fluid.setEquilibrium(1, U, 0);
    fluid.enableForcing();
    const co2 = new ScalarField(fluid, 0.033, 0);
    const o2 = new ScalarField(fluid, 0.042, 0);
    const h2o = new ScalarField(fluid, 0.05, 0);
    co2.enableSource();
    o2.enableSource();
    h2o.enableSource();
    const species = [co2, o2, h2o];
    const { isLeaf, charLen, notes: geometryNotes } = cfg.geometry(fluid, { co2, o2, h2o });
    if (cfg.forcedU !== void 0) {
      for (let y = 1; y < ny - 1; y++) {
        for (const x of [0, nx - 1]) {
          const c = fluid.index(x, y);
          fluid.solid[c] = 0;
          for (const sp of species) if (sp.dirichlet) sp.dirichlet[c] = NaN;
        }
      }
    }
    const fluidInletOutlet = () => {
      for (let y = 1; y < ny - 1; y++) {
        const ci = fluid.index(0, y) * 9;
        for (let i = 0; i < 9; i++) fluid.f[ci + i] = feq(i, 1, U, 0);
        const co = fluid.index(nx - 1, y) * 9;
        const cs = fluid.index(nx - 2, y) * 9;
        if ((cfg.outletBC ?? "pressure") === "pressure") {
          let rho = 0;
          let mx = 0;
          let my = 0;
          for (let i = 0; i < 9; i++) {
            const fi = fluid.f[cs + i];
            rho += fi;
            mx += CX[i] * fi;
            my += CY[i] * fi;
          }
          const vx = mx / rho;
          const vy = my / rho;
          for (let i = 0; i < 9; i++) {
            fluid.f[co + i] = feq(i, 1, vx, vy) + (fluid.f[cs + i] - feq(i, rho, vx, vy));
          }
        } else {
          for (let i = 0; i < 9; i++) fluid.f[co + i] = fluid.f[cs + i];
        }
      }
    };
    const age = cfg.ageOfAir && cfg.forcedU !== void 0 ? new ScalarField(fluid, 0.033, 0) : null;
    if (age) {
      age.enableSource();
      for (let c = 0; c < fluid.size; c++) if (!fluid.solid[c]) age.source[c] = 1;
    }
    const speciesInletOutlet = () => {
      for (const sp of age ? [...species, age] : species) {
        for (let y = 1; y < ny - 1; y++) {
          const ci = fluid.index(0, y) * 5;
          for (let i = 0; i < 5; i++) sp.g[ci + i] = 0;
          const co = fluid.index(nx - 1, y) * 5;
          const cs = fluid.index(nx - 2, y) * 5;
          for (let i = 0; i < 5; i++) sp.g[co + i] = sp.g[cs + i];
        }
      }
    };
    const sScale = (cfg.sourceScale ?? 1) * (cfg.dark ? -0.31 : 1) / r;
    const S_CO2 = S_CO2_BASE * sScale;
    const S_O2 = 4e-4 * sScale;
    const S_H2O = 5e-4 * (cfg.sourceScale ?? 1) * (cfg.dark ? 0.31 : 1) / r;
    const surfaceCells = [];
    for (let x = 1; x < nx - 1; x++) {
      for (let y = 1; y < ny - 1; y++) {
        const c = fluid.index(x, y);
        if (fluid.solid[c]) continue;
        const touches = isLeaf[fluid.index(x + 1, y)] || isLeaf[fluid.index(x - 1, y)] || isLeaf[fluid.index(x, y + 1)] || isLeaf[fluid.index(x, y - 1)];
        if (touches) {
          surfaceCells.push(c);
          co2.source[c] = -S_CO2;
          o2.source[c] = S_O2;
          h2o.source[c] = S_H2O;
        }
      }
    }
    const Ca = cfg.co2Ambient ?? 0;
    const GAMMA = 0.125 * Ca;
    const KC = 0.75 * Ca;
    const mAmbient = Ca > 0 ? (Ca - GAMMA) / (Ca - GAMMA + KC) : 1;
    let aEff = 1;
    let eEff = 1;
    const updateTranspiration = () => {
      const sat = cfg.h2oSatExcess;
      if (!(sat !== void 0 && sat > 0)) return;
      let sum = 0;
      for (const c of surfaceCells) {
        const m = Math.max(0, 1 - h2o.C[c] / sat);
        h2o.source[c] = S_H2O * m;
        sum += m;
      }
      eEff = sum / (surfaceCells.length || 1);
    };
    const updatePhotosynthesis = () => {
      if (!(Ca > 0) || cfg.dark) return;
      let sum = 0;
      for (const c of surfaceCells) {
        const cAbs = Ca + co2.C[c];
        const m = cAbs <= GAMMA ? 0 : (cAbs - GAMMA) / (cAbs - GAMMA + KC) / mAmbient;
        co2.source[c] = -S_CO2 * m;
        o2.source[c] = S_O2 * m;
        sum += m;
      }
      aEff = sum / surfaceCells.length;
    };
    const membraneCells = [];
    if (cfg.membraneK !== void 0) {
      for (let x = 0; x < nx; x++) {
        for (let y = 0; y < ny; y++) {
          if (x !== 0 && x !== nx - 1 && y !== 0 && y !== ny - 1) continue;
          const c = fluid.index(x, y);
          if (co2.dirichlet) co2.dirichlet[c] = NaN;
          if (o2.dirichlet) o2.dirichlet[c] = NaN;
          if (h2o.dirichlet) h2o.dirichlet[c] = NaN;
        }
      }
      if (cfg.membraneK > 0) {
        const onBorder = (xx, yy) => xx === 0 || xx === nx - 1 || yy === ny - 1;
        for (let x = 1; x < nx - 1; x++) {
          for (let y = 1; y < ny - 1; y++) {
            const c = fluid.index(x, y);
            if (fluid.solid[c]) continue;
            if (onBorder(x + 1, y) || onBorder(x - 1, y) || onBorder(x, y + 1) || onBorder(x, y - 1))
              membraneCells.push(c);
          }
        }
      }
    }
    const leak = cfg.membraneK ?? 0;
    const applyMembrane = () => {
      const f = 1 - leak;
      for (const c of membraneCells) {
        for (const sp of species) {
          const b = c * 5;
          for (let i = 0; i < 5; i++) sp.g[b + i] *= f;
        }
      }
    };
    const B = 8e-4;
    const gLat = B * cfg.gRatio / r ** 3;
    const include = new Set(cfg.buoyancySpecies ?? ["h2o", "co2", "o2"]);
    const contributors = [
      { name: "h2o", field: h2o, beta: -1, ref: 0 },
      { name: "co2", field: co2, beta: 0.7, ref: 0 },
      { name: "o2", field: o2, beta: 0.5, ref: 0 }
    ].filter((c) => include.has(c.name));
    const surfaceStats = (field) => {
      let sum = 0;
      let peak = 0;
      for (const c of surfaceCells) {
        const v = field.C[c];
        sum += v;
        if (Math.abs(v) > Math.abs(peak)) peak = v;
      }
      return { mean: sum / surfaceCells.length, peak };
    };
    const maxSpeed = () => {
      let m = 0;
      for (let c = 0; c < fluid.size; c++) {
        if (fluid.solid[c]) continue;
        const s = fluid.speed(c);
        if (s > m) m = s;
      }
      return m;
    };
    const domainMean = (field) => {
      let sum = 0;
      let n = 0;
      for (let c = 0; c < fluid.size; c++) {
        if (fluid.solid[c]) continue;
        sum += field.C[c];
        n++;
      }
      return sum / n;
    };
    const boundaryLabel = cfg.forcedU !== void 0 ? "vented (forced airflow)" : cfg.membraneK === void 0 ? "ambient walls" : cfg.membraneK === 0 ? "sealed (BRIC)" : "micropore tape (CARA)";
    return {
      id: cfg.id,
      label: cfg.label,
      fluid,
      scalarField: h2o,
      renderMode: "scalar",
      renderScale: cfg.renderScale,
      postStream() {
        if (cfg.forcedU !== void 0) fluidInletOutlet();
      },
      onAfterStep() {
        updatePhotosynthesis();
        updateTranspiration();
        co2.step();
        o2.step();
        h2o.step();
        if (age) age.step();
        if (cfg.forcedU !== void 0) speciesInletOutlet();
        if (membraneCells.length) applyMembrane();
        if (cfg.buoyancyRef === "enclosureMean" && gLat !== 0) {
          for (const ct of contributors) ct.ref = domainMean(ct.field);
        }
        applyBoussinesqForce(fluid, 0, -gLat, contributors);
      },
      diagnostics(step) {
        fluid.computeMacroscopic();
        const umax = maxSpeed();
        const w = surfaceStats(h2o);
        const cc = surfaceStats(co2);
        const oo = surfaceStats(o2);
        const Sc = nu / 0.05;
        const Ra = gLat * Math.abs(w.mean) * charLen ** 3 / (nu * nu) * Sc;
        const gAbs = cfg.gRatio * GRAVITY.earth;
        const out = [
          ...geometryNotes ?? [],
          { label: "enclosure", value: boundaryLabel },
          { label: "gravity", value: `${gAbs.toFixed(2)} m/s\xB2  (${cfg.gRatio.toFixed(3)} g)` }
        ];
        if (cfg.forcedU !== void 0) {
          out.push({ label: "forced airflow", value: `${(U * r * U_STAR_CM_S).toFixed(1)} cm/s` });
        }
        out.push(
          { label: "steps", value: step.toLocaleString() },
          { label: "u_max (convection)", value: umax.toExponential(2) }
        );
        if (cfg.membraneK !== void 0 || cfg.dark) {
          out.push({ label: "dish-mean CO\u2082 excess", value: domainMean(co2).toFixed(3) });
          if (cfg.h2oSatExcess !== void 0) {
            out.push(
              { label: "dish-mean H\u2082O excess", value: domainMean(h2o).toFixed(4) },
              { label: "transpiration (% potential)", value: (eEff * 100).toFixed(1) }
            );
          }
        } else {
          out.push({ label: "Rayleigh (H\u2082O)", value: Ra.toExponential(1) });
        }
        out.push(
          { label: "\u0394C H\u2082O  mean/peak", value: `${w.mean.toFixed(3)} / ${w.peak.toFixed(3)}` },
          { label: "\u0394C CO\u2082  mean/peak", value: `${cc.mean.toFixed(3)} / ${cc.peak.toFixed(3)}` },
          { label: "\u0394C O\u2082  mean/peak", value: `${oo.mean.toFixed(3)} / ${oo.peak.toFixed(3)}` }
        );
        {
          const v = surfaceCells.map((c) => Math.abs(co2.C[c])).sort((a, b) => a - b);
          const n = v.length || 1;
          const mean = v.reduce((s, x) => s + x, 0) / n;
          const sd = Math.sqrt(v.reduce((s, x) => s + (x - mean) ** 2, 0) / n);
          const q = (p) => v[Math.min(v.length - 1, Math.floor(p * (v.length - 1)))] ?? NaN;
          const worse = v.filter((x) => x > DC_CO2_REF).length / n;
          out.push(
            { label: "surface |\u0394C CO\u2082| CoV", value: (mean > 0 ? sd / mean : NaN).toFixed(4) },
            { label: "surface |\u0394C CO\u2082| p10/p50/p90", value: `${q(0.1).toFixed(4)} / ${q(0.5).toFixed(4)} / ${q(0.9).toFixed(4)}` },
            { label: "surface frac > Earth-leaf mean", value: worse.toFixed(4) }
          );
        }
        let sAbs = 0;
        for (const c of surfaceCells) sAbs += Math.abs(co2.source[c]);
        sAbs /= surfaceCells.length || 1;
        const dCco2 = Math.abs(cc.mean);
        if (dCco2 > 5e-3 && sAbs > 0) {
          const gbl = K_GBL * r * sAbs / dCco2;
          const gblMS = gbl / MOLAR_DENSITY;
          const deltaMM = D_CO2_PHYS / gblMS * 1e3;
          const Sh = charLen * (DX_M / r) / (deltaMM * 1e-3);
          out.push(
            { label: "g_bl CO\u2082", value: `${gbl.toFixed(3)} mol m\u207B\xB2 s\u207B\xB9` },
            { label: "\u03B4 film (CO\u2082)", value: `${deltaMM.toFixed(2)} mm` },
            { label: "Sherwood Sh", value: Sh.toFixed(1) }
          );
        }
        {
          const dCrel = Math.abs(cc.mean - domainMean(co2));
          if (dCrel > 5e-3 && sAbs > 0) {
            out.push({ label: "g_bl vs enclosure bulk", value: `${(K_GBL * r * sAbs / dCrel).toFixed(3)} mol m\u207B\xB2 s\u207B\xB9` });
          }
        }
        if (age) {
          let nFluid = 0;
          let ageSum = 0;
          for (let c = 0; c < fluid.size; c++) {
            if (fluid.solid[c]) continue;
            nFluid++;
            ageSum += age.C[c];
          }
          let qIn = 0;
          let qOut = 0;
          let aOut = 0;
          for (let y = 1; y < ny - 1; y++) {
            const ci = fluid.index(1, y);
            const co = fluid.index(nx - 2, y);
            if (!fluid.solid[ci]) qIn += fluid.ux[ci];
            if (!fluid.solid[co]) {
              qOut += fluid.ux[co];
              aOut += fluid.ux[co] * age.C[co];
            }
          }
          const tau = nFluid / qIn;
          const meanAge = ageSum / nFluid;
          out.push(
            { label: "age \u03C4 (nominal)", value: `${tau.toFixed(1)} steps` },
            { label: "age outlet / \u03C4", value: (aOut / qOut / tau).toFixed(4) },
            { label: "age mean / \u03C4", value: (meanAge / tau).toFixed(4) },
            { label: "\u03B5_a ventilation", value: (tau / (2 * meanAge)).toFixed(4) },
            { label: "outlet / inlet flux", value: (qOut / qIn).toFixed(4) }
          );
        }
        if (cfg.co2Ambient !== void 0 && !cfg.dark) {
          out.push({ label: "net assimilation (% potential)", value: (aEff * 100).toFixed(1) });
        }
        return out;
      }
    };
  };
}
var R2E = GRAVITY.mars / GRAVITY.earth;
var M2E = GRAVITY.moon / GRAVITY.earth;
var LEAF_DESC = "Single Arabidopsis leaf in a closed chamber \u2014 stomatal CO\u2082 uptake + O\u2082/H\u2082O release with solutal buoyancy. Sweep gravity (select each preset, Run ~30k steps) and watch u_max fall and the surface gaps \u0394C steepen.";
var ROSETTE_DESC = "A rosette of overlapping leaves. Leaves shade each other\u2019s airflow and trap air near the crown, so the interior surface gaps (\u0394C peak) run steeper than an isolated leaf \u2014 an effect that worsens as gravity falls.";
var FAN_DESC = "A leaf in microgravity with forced ventilation \u2014 the engineered substitute for the buoyant convection \xB5g removes. As fan speed rises, the inlet flow thins the boundary layer and the surface gaps \u0394C fall back toward Earth-1 g levels (compare with \u201CLeaf \xB7 Earth\u201D).";
var BRIC_DESC = "BRIC spaceflight hardware \u2014 a hermetically SEALED Petri dish in \xB5g. No gas exchange with the cabin, so the enclosure atmosphere drifts without bound: in light the leaf depletes its own CO\u2082 (carbon starvation); in the dark respiration consumes O\u2082 and builds CO\u2082 (hypoxia). Watch \u201Cdish-mean CO\u2082 excess\u201D run away \u2014 it never reaches steady state.";
var CARA_DESC = "CARA spaceflight hardware \u2014 a Petri dish sealed with gas-permeable micropore surgical tape in \xB5g. The tape vents the dish toward the cabin, so the dish-mean stays bounded near ambient \u2014 but with no convection inside, a diffusive boundary layer still steepens the gases at the leaf surface. Compare its dish-mean (bounded) with BRIC (runaway).";
var VEGGIE_DESC = "VEGGIE spaceflight hardware \u2014 light + forced airflow actively ventilating the growth volume in \xB5g. The circulated air thins the leaf boundary layer as well as venting the enclosure, restoring near-Earth surface gas gradients. (If the dish stayed taped, internal gradients would resemble CARA; VEGGIE\u2019s design intent is to ventilate the plants directly.)";
var FB_DESC = "Closed-loop photosynthesis: the leaf\u2019s CO\u2082 uptake is now CO\u2082-limited (compensation point + saturation), so the boundary-layer CO\u2082 depletion the model builds up feeds back to suppress assimilation. Watch \u201Cnet assimilation (% potential)\u201D \u2014 it falls as the boundary layer thickens (\xB5g) and where air is trapped (rosette crown), and is highest under ventilation.";
var CANOPY_DESC = "A microgreen \u201Clawn\u201D \u2014 a dense row of upright shoots on soil, ambient above. Convection ventilates only the canopy top; the within-canopy air stagnates, and in microgravity the whole stand\u2019s gas gaps blow out.";
var SHAPE_DESC = "The same single-leaf chamber, but the blade is a real leaf outline cut along its midrib instead of an ellipse. Chord (52 cells) and thickness (8 cells) match \u201CLeaf \xB7 Earth/\xB5g\u201D exactly, so the difference is down to the outline alone \u2014 a tapered blade with a narrow base and a rounded tip, not a symmetric ellipse. Measured (T14): the real outline conducts LESS than the ellipse \u2014 g_bl 0.941 vs 0.997 at 1 g (\u22125.6%) and 0.426 vs 0.443 in \xB5g (\u22123.8%), at steady state. The ellipse has been flattering the leaf by a few percent. NOTE: the shipped outline is a provisional analytic stand-in until the NeuraLeaf fitting pipeline has been run on real imagery; the readout says so.";
var CURL_DESC = "Curl \xD7 gravity. The cut runs ACROSS the blade \u2014 the plane cupping actually lives in \u2014 so a curled leaf appears as a U that traps air underneath. Measured (T14): curl costs about the same in blade-MEAN conductance at either gravity (\xD70.81 at 1 g, \xD70.83 in \xB5g), so curl and gravity do NOT compound on the mean. They do compound on the WORST SPOT: the peak surface gap runs \u22120.208 \u2192 \u22120.326 at 1 g but \u22120.338 \u2192 \u22120.540 in \xB5g, so a curled blade in \xB5g sits at 2.6\xD7 the peak CO\u2082 gap of a flat blade at 1 g. CAVEAT: the curled section is 35 cells wide against the flat one's 27, because the blade thickness has to be inflated ~11\xD7 to be resolvable at all and that fattens a tilted margin sideways. Curl is therefore not the only variable here \u2014 treat these as indicative until the blade is resolved on a finer lattice. Not comparable with the longitudinal-cut scenarios either, which have a different characteristic length.";
var DOM = { nx: 128, ny: 96 };
var SCENARIOS = [
  {
    id: "leaf-earth",
    label: "Leaf \xB7 Earth (1 g)",
    description: LEAF_DESC,
    build: makeLeafScene({ ...DOM, id: "leaf-earth", label: "Single leaf \u2014 Earth (1 g)", gRatio: 1, renderScale: 0.12, geometry: leafGeometry })
  },
  {
    id: "leaf-mars",
    label: "Leaf \xB7 Mars (0.38 g)",
    description: LEAF_DESC,
    build: makeLeafScene({ ...DOM, id: "leaf-mars", label: "Single leaf \u2014 Mars (0.38 g)", gRatio: R2E, renderScale: 0.12, geometry: leafGeometry })
  },
  {
    id: "leaf-moon",
    label: "Leaf \xB7 Moon (0.17 g)",
    description: LEAF_DESC,
    build: makeLeafScene({ ...DOM, id: "leaf-moon", label: "Single leaf \u2014 Moon (0.17 g)", gRatio: M2E, renderScale: 0.12, geometry: leafGeometry })
  },
  {
    id: "leaf-ug",
    label: "Leaf \xB7 microgravity (0 g)",
    description: LEAF_DESC,
    build: makeLeafScene({ ...DOM, id: "leaf-ug", label: "Single leaf \u2014 microgravity (0 g)", gRatio: 0, renderScale: 0.12, geometry: leafGeometry })
  },
  {
    id: "rosette-earth",
    label: "Rosette \xB7 Earth (1 g)",
    description: ROSETTE_DESC,
    build: makeLeafScene({ ...DOM, id: "rosette-earth", label: "Rosette \u2014 Earth (1 g)", gRatio: 1, renderScale: 0.18, geometry: rosetteGeometry })
  },
  {
    id: "rosette-ug",
    label: "Rosette \xB7 microgravity (0 g)",
    description: ROSETTE_DESC,
    build: makeLeafScene({ ...DOM, id: "rosette-ug", label: "Rosette \u2014 microgravity (0 g)", gRatio: 0, renderScale: 0.18, geometry: rosetteGeometry })
  },
  {
    id: "canopy-earth",
    label: "Microgreen canopy \xB7 Earth (1 g)",
    description: CANOPY_DESC,
    build: makeLeafScene({ ...DOM, id: "canopy-earth", label: "Microgreen canopy \u2014 Earth (1 g)", gRatio: 1, renderScale: 0.3, sourceScale: 0.3, geometry: canopyGeometry })
  },
  {
    id: "canopy-ug",
    label: "Microgreen canopy \xB7 microgravity (0 g)",
    description: CANOPY_DESC,
    build: makeLeafScene({ ...DOM, id: "canopy-ug", label: "Microgreen canopy \u2014 microgravity (0 g)", gRatio: 0, renderScale: 0.3, sourceScale: 0.3, geometry: canopyGeometry })
  },
  {
    id: "leaf-ug-fan-lo",
    label: "Leaf \xB7 \xB5g + fan 3 cm/s",
    description: FAN_DESC,
    build: makeLeafScene({ ...DOM, id: "leaf-ug-fan-lo", label: "Leaf \u2014 \xB5g + fan ~3 cm/s", gRatio: 0, renderScale: 0.12, geometry: leafGeometry, forcedU: 0.02 })
  },
  {
    id: "leaf-ug-fan-mid",
    label: "Leaf \xB7 \xB5g + fan 8 cm/s",
    description: FAN_DESC,
    build: makeLeafScene({ ...DOM, id: "leaf-ug-fan-mid", label: "Leaf \u2014 \xB5g + fan ~8 cm/s", gRatio: 0, renderScale: 0.12, geometry: leafGeometry, forcedU: 0.05 })
  },
  {
    id: "leaf-ug-fan-hi",
    label: "Leaf \xB7 \xB5g + fan 17 cm/s",
    description: FAN_DESC,
    build: makeLeafScene({ ...DOM, id: "leaf-ug-fan-hi", label: "Leaf \u2014 \xB5g + fan ~17 cm/s", gRatio: 0, renderScale: 0.12, geometry: leafGeometry, forcedU: 0.1 })
  },
  {
    id: "hw-bric-light",
    label: "Hardware \xB7 BRIC sealed (light)",
    description: BRIC_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-bric-light", label: "BRIC \u2014 sealed dish, \xB5g, light", gRatio: 0, renderScale: 0.5, geometry: leafGeometry, membraneK: 0 })
  },
  {
    id: "hw-bric-dark",
    label: "Hardware \xB7 BRIC sealed (dark)",
    description: BRIC_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-bric-dark", label: "BRIC \u2014 sealed dish, \xB5g, dark (respiration)", gRatio: 0, renderScale: 0.5, geometry: leafGeometry, membraneK: 0, dark: true })
  },
  {
    id: "hw-cara",
    label: "Hardware \xB7 CARA micropore tape",
    description: CARA_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-cara", label: "CARA \u2014 micropore-taped dish, \xB5g, light", gRatio: 0, renderScale: 0.3, geometry: leafGeometry, membraneK: 0.01 })
  },
  {
    id: "hw-cara-dark",
    label: "Hardware \xB7 CARA micropore tape (dark)",
    description: CARA_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-cara-dark", label: "CARA \u2014 micropore-taped dish, \xB5g, dark", gRatio: 0, renderScale: 0.3, geometry: leafGeometry, membraneK: 0.01, dark: true })
  },
  {
    id: "hw-veggie",
    label: "Hardware \xB7 VEGGIE vented",
    description: VEGGIE_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-veggie", label: "VEGGIE \u2014 vented + light, \xB5g", gRatio: 0, renderScale: 0.12, geometry: leafGeometry, forcedU: 0.05 })
  },
  {
    id: "hw-rosette-bric",
    label: "Rosette \xB7 BRIC sealed",
    description: BRIC_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-rosette-bric", label: "Rosette \u2014 BRIC sealed, \xB5g", gRatio: 0, renderScale: 0.6, geometry: rosetteGeometry, membraneK: 0 })
  },
  {
    id: "hw-rosette-veggie",
    label: "Rosette \xB7 VEGGIE vented",
    description: VEGGIE_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-rosette-veggie", label: "Rosette \u2014 vented + light, \xB5g", gRatio: 0, renderScale: 0.18, geometry: rosetteGeometry, forcedU: 0.05 })
  },
  {
    id: "hw-canopy-bric",
    label: "Canopy \xB7 BRIC sealed",
    description: BRIC_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-canopy-bric", label: "Canopy \u2014 BRIC sealed, \xB5g", gRatio: 0, renderScale: 0.7, geometry: canopyGeometry, sourceScale: 0.3, membraneK: 0 })
  },
  {
    id: "hw-canopy-veggie",
    label: "Canopy \xB7 VEGGIE vented",
    description: VEGGIE_DESC,
    build: makeLeafScene({ ...DOM, id: "hw-canopy-veggie", label: "Canopy \u2014 vented + light, \xB5g", gRatio: 0, renderScale: 0.3, geometry: canopyGeometry, sourceScale: 0.3, forcedU: 0.05 })
  },
  {
    id: "fb-leaf-earth",
    label: "Feedback \xB7 Leaf Earth (1 g)",
    description: FB_DESC,
    build: makeLeafScene({ ...DOM, id: "fb-leaf-earth", label: "Leaf \u2014 Earth, CO\u2082-limited photosynthesis", gRatio: 1, renderScale: 0.12, geometry: leafGeometry, co2Ambient: 8 })
  },
  {
    id: "fb-leaf-ug",
    label: "Feedback \xB7 Leaf \xB5g (0 g)",
    description: FB_DESC,
    build: makeLeafScene({ ...DOM, id: "fb-leaf-ug", label: "Leaf \u2014 \xB5g, CO\u2082-limited photosynthesis", gRatio: 0, renderScale: 0.12, geometry: leafGeometry, co2Ambient: 8 })
  },
  {
    id: "fb-rosette-ug",
    label: "Feedback \xB7 Rosette \xB5g (0 g)",
    description: FB_DESC,
    build: makeLeafScene({ ...DOM, id: "fb-rosette-ug", label: "Rosette \u2014 \xB5g, CO\u2082-limited photosynthesis", gRatio: 0, renderScale: 0.18, geometry: rosetteGeometry, co2Ambient: 8 })
  },
  {
    id: "shape-flat-earth",
    label: "Shape \xB7 real outline, Earth (1 g)",
    description: SHAPE_DESC,
    build: makeLeafScene({ ...DOM, id: "shape-flat-earth", label: "Leaf outline \u2014 Earth (1 g)", gRatio: 1, renderScale: 0.12, geometry: assetLeafGeometry({ assetId: "provisional-mature-flat" }) })
  },
  {
    id: "shape-flat-ug",
    label: "Shape \xB7 real outline, \xB5g (0 g)",
    description: SHAPE_DESC,
    build: makeLeafScene({ ...DOM, id: "shape-flat-ug", label: "Leaf outline \u2014 microgravity (0 g)", gRatio: 0, renderScale: 0.12, geometry: assetLeafGeometry({ assetId: "provisional-mature-flat" }) })
  },
  {
    id: "curl-none-earth",
    label: "Curl \xB7 flat blade, Earth (1 g)",
    description: CURL_DESC,
    build: makeLeafScene({ ...DOM, id: "curl-none-earth", label: "Transverse cut, flat \u2014 Earth (1 g)", gRatio: 1, renderScale: 0.12, geometry: assetLeafGeometry({ assetId: "provisional-mature-flat", cut: "transverse" }) })
  },
  {
    id: "curl-strong-earth",
    label: "Curl \xB7 curled blade, Earth (1 g)",
    description: CURL_DESC,
    build: makeLeafScene({ ...DOM, id: "curl-strong-earth", label: "Transverse cut, curl 0.8 \u2014 Earth (1 g)", gRatio: 1, renderScale: 0.12, geometry: assetLeafGeometry({ assetId: "provisional-mature-flat", cut: "transverse", deform: { curl: 0.8 } }) })
  },
  {
    id: "curl-none-ug",
    label: "Curl \xB7 flat blade, \xB5g (0 g)",
    description: CURL_DESC,
    build: makeLeafScene({ ...DOM, id: "curl-none-ug", label: "Transverse cut, flat \u2014 microgravity (0 g)", gRatio: 0, renderScale: 0.12, geometry: assetLeafGeometry({ assetId: "provisional-mature-flat", cut: "transverse" }) })
  },
  {
    id: "curl-strong-ug",
    label: "Curl \xB7 curled blade, \xB5g (0 g)",
    description: CURL_DESC,
    build: makeLeafScene({ ...DOM, id: "curl-strong-ug", label: "Transverse cut, curl 0.8 \u2014 microgravity (0 g)", gRatio: 0, renderScale: 0.12, geometry: assetLeafGeometry({ assetId: "provisional-mature-flat", cut: "transverse", deform: { curl: 0.8 } }) })
  },
  {
    id: "cavity",
    label: "Lid-driven cavity",
    description: "Validation gate 1 \u2014 compares the vertical/horizontal centreline velocity profiles to Ghia et al. (1982) at Re 100.",
    build: buildCavity
  },
  {
    id: "cylinder",
    label: "Flow past cylinder",
    description: "Validation gate 2 \u2014 measures the K\xE1rm\xE1n vortex-shedding Strouhal number. Unconfined St \u2248 0.164 at Re 100; this confined channel (\u224817% blockage) raises it to \u2248 0.19\u20130.20 (measured 0.196).",
    build: buildCylinder
  },
  {
    id: "natconv",
    label: "Natural convection (buoyancy)",
    description: "Validation gate 3 \u2014 differentially heated cavity at Ra 1e4. Buoyancy from the temperature field drives a convection roll; the hot-wall Nusselt number must match de Vahl Davis (1983), Nu \u2248 2.24. This is the proof that the gravity/buoyancy coupling is physically correct.",
    build: buildNaturalConvection
  },
  {
    id: "diffusion",
    label: "Pure diffusion (microgravity limit)",
    description: "Validation gate 4 \u2014 with no flow, an initial concentration step must relax to the error-function profile C(x,t)=\xBD erfc((x\u2212x0)/2\u221A(Dt)). This molecular-diffusion limit is exactly what dominates gas transport once microgravity removes convection.",
    build: buildDiffusion
  }
];

// validation/jobs.ts
import { mkdirSync, readFileSync, readdirSync, writeFileSync } from "node:fs";
var val = (rows, key) => rows.find((r) => r.label.startsWith(key))?.value ?? "";
var num = (rows, key) => {
  const s = val(rows, key);
  return s ? parseFloat(s.replace(/[^0-9eE.+-].*$/, "")) : NaN;
};
function cli(name, jobs, runJob, merge) {
  const dir = `results/tables/.partials/${name}`;
  const [cmd, arg] = process.argv.slice(2);
  const doJob = (j) => {
    const t0 = Date.now();
    const out = runJob(j);
    mkdirSync(dir, { recursive: true });
    writeFileSync(`${dir}/${j.id}.json`, JSON.stringify(out));
    console.log(`${name}/${j.id} done in ${((Date.now() - t0) / 1e3).toFixed(0)} s`);
  };
  const doMerge = () => {
    const results = /* @__PURE__ */ new Map();
    const have = new Set(readdirSync(dir));
    const missing = jobs.filter((j) => !have.has(`${j.id}.json`)).map((j) => j.id);
    if (missing.length) throw new Error(`${name}: missing partials for ${missing.join(", ")}`);
    for (const j of jobs) results.set(j.id, JSON.parse(readFileSync(`${dir}/${j.id}.json`, "utf8")));
    merge(results);
  };
  if (cmd === "list") {
    for (const j of jobs) console.log(j.id);
  } else if (cmd === "job") {
    const j = jobs.find((x) => x.id === arg);
    if (!j) throw new Error(`${name}: unknown job ${arg}`);
    doJob(j);
  } else if (cmd === "merge") {
    doMerge();
  } else {
    jobs.forEach(doJob);
    doMerge();
  }
}

// validation/dishes.ts
var DX = 0.288;
var WALL = 1;
var GAP_MM = 1;
var cells = (mm) => Math.round(mm / DX);
function dishGeom(d) {
  const agarMm = d.agarMl / d.floorCm2 * 10;
  const headMm = d.heightMm - WALL - agarMm;
  const agarCells = cells(agarMm);
  const headCells = cells(headMm);
  return { agarMm, headMm, nx: cells(d.widthMm) + 2, ny: agarCells + headCells + 2, agarCells, headCells };
}
var bricW = (54.8 + 55.6) / 2 - 2 * WALL;
var bricFloor = Math.PI * (bricW / 20) ** 2;
var DISHES = {
  "BRIC-mid": { name: "BRIC 60 mm", widthMm: bricW, floorCm2: bricFloor, heightMm: 14, agarMl: 12.5 },
  "BRIC-min": { name: "BRIC 60 mm", widthMm: bricW, floorCm2: bricFloor, heightMm: 13, agarMl: 15 },
  "BRIC-max": { name: "BRIC 60 mm", widthMm: bricW, floorCm2: bricFloor, heightMm: 15, agarMl: 10 },
  "CARA-mid": { name: "CARA 120 mm sq", widthMm: 120, floorCm2: 144, heightMm: 16, agarMl: 75 },
  "CARA-min": { name: "CARA 120 mm sq", widthMm: 120, floorCm2: 144, heightMm: 15, agarMl: 80 },
  "CARA-max": { name: "CARA 120 mm sq", widthMm: 120, floorCm2: 144, heightMm: 17, agarMl: 70 }
};

// validation/dish_feedback.ts
var DT_S = 173e-6;
var CA = 8;
var GAMMA_EXCESS = 0.125 * CA - CA;
var H2O_SAT = 1.3;
var SAMPLE = 2890;
var SCALE = Number(process.env.DISH_FB_SCALE ?? 1);
var JOBS = [
  { id: "BRIC_sealed_fb_g1", dish: "BRIC-mid", bc: "sealed", g: 1, feedback: true, steps: 24e5 },
  { id: "BRIC_sealed_fb_g0", dish: "BRIC-mid", bc: "sealed", g: 0, feedback: true, steps: 24e5 },
  { id: "CARA_taped_fb_g1", dish: "CARA-mid", bc: "taped", g: 1, feedback: true, steps: 3e5 },
  { id: "CARA_taped_fb_g0", dish: "CARA-mid", bc: "taped", g: 0, feedback: true, steps: 3e5 },
  // Buoyancy-reference fix alone (no feedback), same length as T19, to isolate its effect.
  { id: "BRIC_sealed_nofb_g1", dish: "BRIC-mid", bc: "sealed", g: 1, feedback: false, steps: 15e4 }
];
var rd = (d, k) => num(d, k);
cli(
  "dish_feedback",
  JOBS,
  (j) => {
    const geo = dishGeom(DISHES[j.dish]);
    const inst = makeLeafScene({
      id: j.id,
      label: j.id,
      gRatio: j.g,
      nx: geo.nx,
      ny: geo.ny,
      renderScale: 0.12,
      geometry: dishLeafGeometry({ agarCells: geo.agarCells, gapCells: cells(GAP_MM) }),
      membraneK: j.bc === "sealed" ? 0 : 0.01,
      co2Ambient: j.feedback ? CA : void 0,
      buoyancyRef: "enclosureMean",
      h2oSatExcess: H2O_SAT
    })();
    const f = inst.fluid;
    const series = [];
    const steps = Math.max(SAMPLE, Math.round(j.steps * SCALE));
    for (let s = 1; s <= steps; s++) {
      f.collideAndStream();
      inst.postStream();
      inst.onAfterStep(s);
      if (s % SAMPLE === 0 || s === steps) {
        const d = inst.diagnostics(s);
        series.push([s * DT_S, rd(d, "dish-mean"), j.feedback ? rd(d, "net assimilation") : NaN, rd(d, "g_bl vs enclosure bulk"), rd(d, "g_bl CO\u2082")]);
      }
    }
    return { ...j, series };
  },
  (res) => {
    const rows = ["job,t_s,dish_mean_CO2,net_assimilation_pct,g_bl_vs_dish_air,g_bl_vs_ambient"];
    const summary = ["job,gravity_g,boundary,feedback,t_end_s,dish_mean_CO2_end,A_end_pct,t_A_below_10pct_s,t_dish_to_compensation_s,g_bl_dish_air_at_5s,g_bl_dish_air_at_26s,g_bl_dish_air_end"];
    const fx = (v, dp) => v === null || !Number.isFinite(v) ? "NaN" : v.toFixed(dp);
    const at = (s, t, k) => {
      let best = s[0];
      for (const r of s) if (Math.abs(r[0] - t) < Math.abs(best[0] - t)) best = r;
      return best[k];
    };
    for (const j of JOBS) {
      const s = res.get(j.id).series;
      for (const r of s) rows.push([j.id, fx(r[0], 2), fx(r[1], 4), fx(r[2], 1), fx(r[3], 3), fx(r[4], 3)].join(","));
      const last = s[s.length - 1];
      const tA = s.find((r) => r[2] !== null && Number.isFinite(r[2]) && r[2] < 10);
      const tC = s.find((r) => r[1] <= GAMMA_EXCESS * 0.95);
      summary.push([
        j.id,
        j.g,
        j.bc,
        j.feedback,
        fx(last[0], 1),
        fx(last[1], 3),
        fx(last[2], 1),
        tA ? tA[0].toFixed(1) : "not reached",
        tC ? tC[0].toFixed(1) : "not reached",
        fx(at(s, 5.2, 3), 3),
        fx(at(s, 26, 3), 3),
        fx(last[3], 3)
      ].join(","));
    }
    writeFileSync2("results/tables/T20_dish_feedback.csv", rows.join("\n") + "\n");
    writeFileSync2("results/tables/T20_dish_feedback_summary.csv", summary.join("\n") + "\n");
    console.log(summary.join("\n"));
    console.log("wrote T20_dish_feedback.csv, T20_dish_feedback_summary.csv");
  }
);
