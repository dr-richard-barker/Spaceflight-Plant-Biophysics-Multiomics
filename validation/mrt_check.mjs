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
function feq(i, rho, ux, uy) {
  const eu = CX[i] * ux + CY[i] * uy;
  const usq = ux * ux + uy * uy;
  return W[i] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * usq);
}

// src/solver/lbm/LBMFluid.ts
var MRT_M = [
  [1, 1, 1, 1, 1, 1, 1, 1, 1],
  [-4, -1, -1, -1, -1, 2, 2, 2, 2],
  [4, -2, -2, -2, -2, 1, 1, 1, 1],
  [0, 1, 0, -1, 0, 1, -1, -1, 1],
  [0, -2, 0, 2, 0, 1, -1, -1, 1],
  [0, 0, 1, 0, -1, 1, 1, -1, -1],
  [0, 0, -2, 0, 2, 1, 1, -1, -1],
  [0, 1, -1, 1, -1, 0, 0, 0, 0],
  [0, 0, 0, 0, 0, 1, -1, 1, -1]
];
var MRT_M_FLAT = new Float64Array(MRT_M.flat());
var MRT_MINV = new Float64Array(81);
for (let k = 0; k < 9; k++) {
  const norm = MRT_M[k].reduce((s, v) => s + v * v, 0);
  for (let i = 0; i < 9; i++) MRT_MINV[i * 9 + k] = MRT_M[k][i] / norm;
}
var DEFAULT_MRT_RATES = { se: 1.64, seps: 1.54, sq: null };
var LBMFluid = class _LBMFluid {
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
  /** Collision scheme for fluids created from now on (validation scripts set this
   * once to run every scene with MRT). BGK keeps all published results unchanged. */
  static defaultCollision = "bgk";
  collision = _LBMFluid.defaultCollision;
  mrtRates = { ...DEFAULT_MRT_RATES };
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
    if (this.collision === "mrt") this.collideMrt();
    else for (let c = 0; c < this.size; c++) {
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
  /**
   * MRT collision (Lallemand & Luo 2000) with Guo forcing in moment space:
   * m* = m - S (m - m_eq) + (I - S/2) M F. Conserved moments (rho, jx, jy) have s = 0,
   * so momentum gains exactly F per step, as with the BGK branch; the stresses relax
   * at s_nu = 1/tau, so the shear viscosity matches BGK at the same tau.
   */
  collideMrt() {
    const { f, solid, fx, fy, omega } = this;
    const { se, seps } = this.mrtRates;
    const sq = this.mrtRates.sq ?? 8 * (2 - omega) / (8 - omega);
    const S = [0, se, seps, 0, sq, 0, sq, omega, omega];
    const feqv = new Float64Array(Q);
    const src = new Float64Array(Q);
    const mpost = new Float64Array(Q);
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
      const usq = ux * ux + uy * uy;
      for (let i = 0; i < Q; i++) {
        const eu = CX[i] * ux + CY[i] * uy;
        feqv[i] = W[i] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * usq);
        src[i] = fx || fy ? W[i] * ((3 * (CX[i] - ux) + 9 * eu * CX[i]) * Fx + (3 * (CY[i] - uy) + 9 * eu * CY[i]) * Fy) : 0;
      }
      for (let k = 0; k < Q; k++) {
        let m = 0;
        let meq = 0;
        let fm = 0;
        for (let i = 0; i < Q; i++) {
          const mk = MRT_M_FLAT[k * 9 + i];
          if (mk === 0) continue;
          m += mk * f[base + i];
          meq += mk * feqv[i];
          fm += mk * src[i];
        }
        mpost[k] = m - S[k] * (m - meq) + (1 - 0.5 * S[k]) * fm;
      }
      for (let i = 0; i < Q; i++) {
        let fi = 0;
        for (let k = 0; k < Q; k++) fi += MRT_MINV[i * 9 + k] * mpost[k];
        f[base + i] = fi;
      }
    }
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

// validation/mrt_check.ts
function makeFluid(collision, forced, tau) {
  const nx = 24, ny = 18;
  const f = new LBMFluid(nx, ny, tau);
  f.collision = collision;
  const omega = 1 / tau;
  f.mrtRates = { se: omega, seps: omega, sq: omega };
  f.setEquilibrium(1, 0, 0);
  for (let c = 0; c < f.size; c++) {
    for (let i = 0; i < 9; i++) f.f[c * 9 + i] *= 1 + 0.02 * Math.sin(0.7 * c + 1.3 * i);
  }
  for (let x = 0; x < nx; x++) {
    f.solid[f.index(x, 0)] = 1;
    f.solid[f.index(x, ny - 1)] = 1;
  }
  if (forced) {
    f.enableForcing();
    for (let c = 0; c < f.size; c++) {
      f.fx[c] = 1e-5 * Math.cos(0.3 * c);
      f.fy[c] = -2e-5 * Math.sin(0.2 * c);
    }
  }
  return f;
}
var failures = 0;
for (const forced of [false, true]) {
  for (const tau of [0.55, 0.8, 1.2]) {
    const a = makeFluid("bgk", forced, tau);
    const b = makeFluid("mrt", forced, tau);
    for (let s = 0; s < 200; s++) {
      a.collideAndStream();
      b.collideAndStream();
    }
    let maxd = 0;
    for (let k = 0; k < a.f.length; k++) maxd = Math.max(maxd, Math.abs(a.f[k] - b.f[k]));
    const ok = maxd < 1e-5;
    if (!ok) failures++;
    console.log(`${forced ? "forced  " : "unforced"} tau=${tau}: max |f_BGK - f_MRT| = ${maxd.toExponential(2)} ${ok ? "PASS" : "FAIL"}`);
  }
}
{
  const a = makeFluid("bgk", true, 0.8);
  const b = makeFluid("mrt", true, 0.8);
  b.mrtRates = { se: 1.64, seps: 1.54, sq: null };
  for (let s = 0; s < 200; s++) {
    a.collideAndStream();
    b.collideAndStream();
  }
  let maxd = 0;
  for (let k = 0; k < a.f.length; k++) maxd = Math.max(maxd, Math.abs(a.f[k] - b.f[k]));
  const ok = maxd > 1e-6;
  if (!ok) failures++;
  console.log(`control (default MRT rates differ from BGK): max diff ${maxd.toExponential(2)} ${ok ? "PASS" : "FAIL"}`);
}
process.exit(failures ? 1 : 0);
