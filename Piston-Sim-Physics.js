/* Piston-Sim-Physics.js
 * Single cylinder, four stroke, spark ignition simulator.
 *
 * The chamber is an axisymmetric (r, y) section: y points down from the deck
 * plane, r out from the bore axis. The mesh has three blocks of rows:
 *   dome rows   recess in the head (fixed),
 *   gap rows    between deck plane and piston crown (stretch with the piston),
 *   bowl rows   recess in the piston (move with the piston).
 * Cells outside the recess profile are masked, so the chamber walls are steps.
 *
 * Flow model: low Mach number. Pressure has a uniform thermodynamic part p0(t)
 * plus a hydrodynamic part that only steers the velocity. Each step:
 *   1 sources: valves, injection, droplets, wall film, spark, flame, knock
 *   2 every cell is compressed isentropically to the new p0, chosen so the
 *     gas fills the new chamber volume exactly
 *   3 projection: a Poisson solve makes the face fluxes move each cell's
 *     surplus volume to its neighbours (this is what drives squish, the
 *     intake jet and the flame expansion)
 *   4 conservative upwind transport (MUSCL, minmod) of every quantity
 *   5 implicit diffusion along lines (turbulent viscosity, wall friction,
 *     wall heat transfer), then the k equation
 * Combustion: a front tracked by the burned fraction c that advances into
 * fresh gas at the turbulent flame speed (Godunov upwind |grad c|), burning
 * the reactants it overruns. Rich leftovers oxidise later by mixing.
 * Fuel liquid: Lagrangian droplet parcels with drag, heating, evaporation and
 * a wall film. Knock: Livengood-Wu integral with the Douaud-Eyzat delay.
 */
'use strict';

const PS_RU = 8314.462618;
const PS_ATM = 101325;
const PS_DEG = Math.PI / 180;

// Constant cp per species, effective values over the range each one lives in.
const PS_AIR = { R: 287.05, cp: 1100 };
const PS_PROD = { R: 290.5, cp: 1400 };

// gulder: [W m/s, eta, xi] for SL0 = W phi^eta exp(-xi (phi - 1.075)^2) at 298 K and 1 atm.
const PS_FUELS = {
  gasolina: {
    label: 'Gasolina', LHV: 43.4e6, AFR: 14.7, M: 100, cpV: 1700, gulder: [0.38, -0.326, 4.48], octane: 95,
    phiMin: 0.5, phiMax: 1.8, liquid: true, rhoL: 720, cL: 2200, Lb: 3.3e5, Tb: 370, Tc: 560,
  },
  etanol: {
    label: 'Etanol', LHV: 26.8e6, AFR: 9.0, M: 46.07, cpV: 1900, gulder: [0.465, 0.25, 6.34], octane: 108,
    phiMin: 0.45, phiMax: 2.2, liquid: true, rhoL: 789, cL: 2440, Lb: 8.4e5, Tb: 351, Tc: 514,
  },
  metano: {
    label: 'Metano (GNC)', LHV: 50.0e6, AFR: 17.2, M: 16.04, cpV: 2600, gulder: [0.422, 0.15, 5.18], octane: 120,
    phiMin: 0.55, phiMax: 1.6, liquid: false,
  },
};

const PS_DEFAULTS = {
  bore: 0.086, stroke: 0.086, rod: 0.143, cr: 10.5,
  chamber: 'cupula',            // 'cupula' (head recess), 'bowl' (piston recess), 'plana'
  squishClear: 0.0011, squishFrac: 0.4,
  rpm: 2000, pIntake: 0.95e5, tIntake: 310, pExhaust: 1.05e5,
  lambda: 1.0, fuel: 'gasolina', octane: null, spark: 18,
  injection: 'puerto',          // 'puerto' or 'directa'
  pfiVapor: 0.8, pfiDrop: 50e-6,
  diSoi: 300, diPressure: 150e5, diDrop: 18e-6, diCone: 70,
  ivo: 12, ivc: 48, evo: 50, evc: 12,
  ivLift: 0.0095, evLift: 0.009, ivDiam: 0.34, evDiam: 0.29,
  swirl: 0.35, jetAngle: 35,
  tHead: 480, tLiner: 420, tPiston: 540, heatTransfer: 1.0,
  flameFactor: 1.0,
  resolution: 'normal',
};

const PS_RESOLUTIONS = {
  rapida: { nr: 14, ng: 6, nrec: 6 },
  normal: { nr: 20, ng: 8, nrec: 8 },
  fina: { nr: 28, ng: 10, nrec: 11 },
};

// ---------------------------------------------------------------- geometry

function psBuildGeometry(p) {
  const res = PS_RESOLUTIONS[p.resolution] || PS_RESOLUTIONS.normal;
  const R = p.bore / 2, Ac = Math.PI * R * R;
  const Vd = Ac * p.stroke, Vc = Vd / (p.cr - 1);
  let type = p.chamber;
  const nr = res.nr, dr = R / nr;
  const Az = new Float64Array(nr), rc = new Float64Array(nr), rw = new Float64Array(nr + 1);
  for (let i = 0; i <= nr; i++) rw[i] = i * dr;
  for (let i = 0; i < nr; i++) {
    Az[i] = Math.PI * (rw[i + 1] * rw[i + 1] - rw[i] * rw[i]);
    rc[i] = (2 / 3) * (rw[i + 1] ** 3 - rw[i] ** 3) / (rw[i + 1] ** 2 - rw[i] ** 2);
  }
  let hTdc = Vc / Ac, rRec = 0, nRec = 0, dzRec = 0;
  const recN = new Int32Array(nr);
  if (type !== 'plana') {
    rRec = R * Math.sqrt(Math.min(0.95, Math.max(0.05, 1 - p.squishFrac)));
    hTdc = Math.min(Math.max(p.squishClear, 0.0003), 0.85 * Vc / Ac);
    nRec = res.nrec;
    let sumA = 0;
    for (let i = 0; i < nr; i++) {
      const s = rc[i] < rRec ? Math.sqrt(1 - (rc[i] / rRec) ** 2) : 0;
      recN[i] = Math.min(nRec, Math.round(s * nRec + 0.25));
      sumA += Az[i] * recN[i];
    }
    if (sumA <= 0) { type = 'plana'; hTdc = Vc / Ac; nRec = 0; recN.fill(0); }
    else dzRec = (Vc - Ac * hTdc) / sumA;
  }
  const nd = type === 'cupula' ? nRec : 0, nb = type === 'bowl' ? nRec : 0, ng = res.ng;
  const NZ = nd + ng + nb, N = nr * NZ;
  const rowType = new Int8Array(NZ);   // 0 dome, 1 gap, 2 bowl
  for (let j = 0; j < NZ; j++) rowType[j] = j < nd ? 0 : j < nd + ng ? 1 : 2;
  const act = new Uint8Array(N);
  for (let i = 0; i < nr; i++) for (let j = 0; j < NZ; j++) {
    let a = 1;
    if (rowType[j] === 0) a = recN[i] > (nd - 1 - j) ? 1 : 0;
    else if (rowType[j] === 2) a = recN[i] > (j - nd - ng) ? 1 : 0;
    act[i * NZ + j] = a;
  }
  // Recess depth per column, used to draw the chamber walls.
  const recDepth = new Float64Array(nr);
  for (let i = 0; i < nr; i++) recDepth[i] = recN[i] * dzRec;
  const valveBand = (c0, c1) => {
    const cells = [];
    for (let i = 0; i < nr; i++) if (rc[i] >= c0 * R && rc[i] <= c1 * R) cells.push(i);
    if (!cells.length) cells.push(Math.min(nr - 1, Math.round(0.5 * (c0 + c1) * nr)));
    return cells;
  };
  return {
    type, R, B: p.bore, Ac, Vd, Vc, hTdc, rRec, nRec, dzRec, recN, recDepth,
    nr, dr, rw, rc, Az, nd, ng, nb, NZ, N, rowType, act,
    exCols: valveBand(0.16, 0.40), inCols: valveBand(0.46, 0.72),
    exR: 0.28 * R, inR: 0.59 * R,
  };
}

// Piston distance from TDC and its speed (m, m/s); ca in degrees.
function psPiston(p, ca, omega) {
  const a = p.stroke / 2, l = p.rod, th = ca * PS_DEG;
  const s = Math.sin(th), c = Math.cos(th), root = Math.sqrt(l * l - a * a * s * s);
  const x = a + l - (a * c + root);
  const v = omega * a * s * (1 + a * c / root);
  return { x, v };
}

// ---------------------------------------------------------------- correlations

function psOrifice(Cd, A, p0, T0, R0, g, pd) {
  if (A <= 0 || p0 <= 0) return 0;
  const pr = Math.min(1, Math.max(pd / p0, 0));
  const crit = Math.pow(2 / (g + 1), g / (g - 1));
  let psi;
  if (pr > crit) psi = Math.sqrt(Math.max(0, 2 * g / (g - 1) * (Math.pow(pr, 2 / g) - Math.pow(pr, (g + 1) / g))));
  else psi = Math.sqrt(g) * Math.pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
  return Cd * A * p0 / Math.sqrt(R0 * T0) * psi;
}

function psLaminarSpeed(fuel, phi, Tu, p, xr) {
  if (!(phi > fuel.phiMin && phi < fuel.phiMax)) return 0;
  const [W, eta, xi] = fuel.gulder;
  const s0 = W * Math.pow(phi, eta) * Math.exp(-xi * (phi - 1.075) ** 2);
  const al = 2.18 - 0.8 * (phi - 1), be = -0.16 + 0.22 * (phi - 1);
  const dil = Math.max(0, 1 - 2.06 * Math.pow(Math.max(xr, 0), 0.77));
  return s0 * Math.pow(Tu / 298, al) * Math.pow(Math.max(p, 1e3) / PS_ATM, be) * dil;
}

// Douaud-Eyzat ignition delay (s).
function psIgnitionDelay(octane, p, T) {
  return 17.68e-3 * Math.pow(octane / 100, 3.402) * Math.pow(p / PS_ATM, -1.7) * Math.exp(3800 / T);
}

const psViscosity = T => 1.8e-5 * Math.pow(T / 300, 0.7);
const psConductivity = T => 0.026 * Math.pow(T / 300, 0.8);

// First circumferential acoustic mode of the bore (Draper), Hz.
function psKnockFrequency(bore, T, gamma, R) {
  return 1.841 * Math.sqrt(gamma * R * T) / (Math.PI * bore);
}

// Heywood's squish velocity estimate for a recess of volume Vrec and diameter Drec.
function psSquishEstimate(bore, Drec, Vrec, z, pistonSpeed) {
  const Ac = Math.PI * bore * bore / 4;
  return (Drec / (4 * z)) * ((bore / Drec) ** 2 - 1) * Vrec / (Ac * z + Vrec) * pistonSpeed;
}

// ---------------------------------------------------------------- linear algebra

// In place Cholesky of a symmetric banded matrix stored by rows: a[i*W + d] = A(i, i - d).
function psBandCholesky(a, n, bw) {
  const W = bw + 1;
  for (let j = 0; j < n; j++) {
    const jb = j * W;
    let s = a[jb];
    for (let d = 1; d <= bw && d <= j; d++) { const l = a[jb + d]; s -= l * l; }
    const piv = Math.sqrt(Math.max(s, 1e-300));
    a[jb] = piv;
    const iMax = Math.min(n - 1, j + bw);
    for (let i = j + 1; i <= iMax; i++) {
      const ib = i * W, dij = i - j;
      let t = a[ib + dij];
      if (t === 0) {
        // Only the fill-in can be non zero; still need the full update.
      }
      const kMin = Math.max(0, i - bw);
      for (let k = kMin; k < j; k++) {
        const aik = a[ib + (i - k)];
        if (aik !== 0) t -= aik * a[jb + (j - k)];
      }
      a[ib + dij] = t / piv;
    }
  }
}

function psBandSolve(a, n, bw, b) {
  const W = bw + 1;
  for (let i = 0; i < n; i++) {
    let s = b[i];
    const ib = i * W;
    for (let d = 1; d <= bw && d <= i; d++) s -= a[ib + d] * b[i - d];
    b[i] = s / a[ib];
  }
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let d = 1; d <= bw && i + d < n; d++) s -= a[(i + d) * W + d] * b[i + d];
    b[i] = s / a[i * W];
  }
}

// Thomas algorithm. lo[k] couples k to k-1, up[k] couples k to k+1.
function psTridiag(n, lo, di, up, rhs, x, work) {
  let b = di[0];
  x[0] = rhs[0] / b;
  for (let k = 1; k < n; k++) {
    work[k] = up[k - 1] / b;
    b = di[k] - lo[k] * work[k];
    x[k] = (rhs[k] - lo[k] * x[k - 1]) / b;
  }
  for (let k = n - 2; k >= 0; k--) x[k] -= work[k + 1] * x[k + 1];
}

// ---------------------------------------------------------------- simulator

const PS_FIELDS = ['m', 'mF', 'mP', 'E', 'mu', 'mw', 'mL', 'mk', 'mc', 'mxr', 'mik', 'mphi'];

class PSCylinder {
  constructor(params) {
    this.p = Object.assign({}, PS_DEFAULTS, params || {});
    this.reset();
  }

  setParam(key, value) {
    this.p[key] = value;
    if (key === 'rpm') this.omega = value * 2 * Math.PI / 60;
  }

  // Live values for the readouts.
  diagnostics() {
    const g = this.g, fuel = this.fuel;
    let M = 0, Tm = 0, kk = 0, Lm = 0, I = 0, mF = 0, mP = 0, mb = 0, Tb = 0;
    for (let c = 0; c < g.N; c++) {
      if (!g.act[c]) continue;
      const m = this.m[c], r = g.rc[(c / g.NZ) | 0];
      M += m; Tm += m * this.T[c]; kk += m * this.k[c]; Lm += this.mL[c]; I += m * r * r;
      mF += this.mF[c]; mP += this.mP[c];
      if (this.c[c] > 0.5) { mb += m; Tb += m * this.T[c]; }
    }
    const air = M - mF - mP;
    const liquid = this.liquidTotal || 0, film = this.filmTotal || 0;
    const fuelAll = mF + liquid + film;
    return {
      p: this.p0, V: this.Vtot, M, T: Tm / M, Tu: this.Tu || Tm / M, Tb: mb > 0 ? Tb / mb : null,
      u: Math.sqrt(2 / 3 * kk / M), swirl: I > 0 ? Lm / (I * this.omega) : 0,
      vapor: mF, liquid, film, air, residual: mP / M,
      lambda: fuelAll > 0 ? air / (fuelAll * fuel.AFR) : Infinity, mfb: this.mfb || 0, usq: this.usq || 0,
      parcels: this.parcels.length, liftIn: this.liftIn || 0, liftEx: this.liftEx || 0, dt: this.dt,
    };
  }

  reset() {
    const p = this.p;
    this.fuel = PS_FUELS[p.fuel] || PS_FUELS.gasolina;
    p.octane = p.octane || this.fuel.octane;
    const g = this.g = psBuildGeometry(p);
    const N = g.N, NZ = g.NZ;
    const F = () => new Float64Array(N);
    for (const f of PS_FIELDS) this[f] = F();
    this.u = F(); this.w = F(); this.v = F(); this.k = F(); this.c = F();
    this.T = F(); this.R = F(); this.cp = F(); this.rho = F();
    this.V = F(); this.Vt = F(); this.S = F();
    this.Fe = F(); this.Fd = F();            // volume flux through east and lower faces (m3/s)
    this.phi = F(); this.q = F(); this.flame = F(); this.film = F(); this.filmWall = new Int8Array(N);
    this.Hcol = new Float64Array(g.nr);
    this.dz = new Float64Array(NZ); this.yTop = new Float64Array(NZ + 1); this.gv = new Float64Array(NZ + 1);
    this.band = new Float64Array(N * (NZ + 1));
    const L = Math.max(g.nr, NZ) + 2;
    this.tri = { lo: new Float64Array(L), di: new Float64Array(L), up: new Float64Array(L), rhs: new Float64Array(L), x: new Float64Array(L), work: new Float64Array(L), idx: new Int32Array(L) };
    this.old = {}; for (const f of PS_FIELDS) this.old[f] = F();
    this.parcels = [];
    this.parcelMass = 0;
    this.injAcc = 0;

    this.omega = p.rpm * 2 * Math.PI / 60;
    this.ca = 540 + p.ivc;                   // start at intake valve closing
    this.t = 0;
    this.dt = 1e-6;
    this.geometryAt(this.ca);

    // Fresh charge with some residual gas.
    const fuel = this.fuel, xr = 0.06;
    const f = 1 / (p.lambda * fuel.AFR);
    const YP = xr, YF = (1 - xr) * f / (1 + f), YA = 1 - YP - YF;
    const T0 = p.tIntake + 25;
    const R0 = YA * PS_AIR.R + YP * PS_PROD.R + YF * PS_RU / fuel.M;
    const cp0 = YA * PS_AIR.cp + YP * PS_PROD.cp + YF * fuel.cpV;
    this.p0 = p.pIntake;
    const Sp = 2 * p.stroke * p.rpm / 60, k0 = 1.5 * (0.45 * Sp) ** 2;
    for (let c = 0; c < N; c++) {
      if (!g.act[c]) continue;
      const m = this.p0 * this.V[c] / (R0 * T0);
      this.m[c] = m; this.mF[c] = m * YF; this.mP[c] = m * YP; this.E[c] = m * cp0 * T0;
      this.mk[c] = m * k0;
      const i = (c / NZ) | 0;
      // Weak solid body swirl left from the intake.
      this.mL[c] = m * g.rc[i] * g.rc[i] * this.omega * 0.6 * p.swirl;
    }
    this.derive();

    this.exhaustRes = { YF: 0.0, YP: 1.0, T: 900 };
    this.intakeBuf = { m: 0, mF: 0, mP: 0, E: 0 };
    this.lastTrappedAir = this.airMass();
    this.sparkDone = false; this.ignited = false; this.misfire = false;
    this.knock = null;
    this.stats = this.newCycleStats();
    this.prevCycle = null;
    this.history = [];
    this.cycleCount = 0;
    this.trace = this.newTrace();
  }

  newTrace() { return { ca: [], p: [], V: [], hrr: [], mfb: [], T: [], u: [], usq: [], vap: [], liq: [], film: [], burned: [] }; }

  newCycleStats() {
    return {
      work: 0, workGross: 0, qComb: 0, qWall: 0, fuelIn: 0, fuelBurned: 0, airTrapped: 0, massTrapped: 0, residual: 0,
      lambdaTrapped: 0, pMax: 0, caPMax: 0, uSpark: 0, sqMax: 0, mfb10: null, mfb50: null, mfb90: null, knockCa: null,
      knockFrac: 0, knockFreq: 0, hcEvo: 0, filmMax: 0, liquidMax: 0, misfire: false, ignitionCa: null,
    };
  }

  // Row heights, positions and face speeds at crank angle ca.
  geometryAt(ca) {
    const g = this.g, NZ = g.NZ, kin = psPiston(this.p, ca, this.omega);
    const h = g.hTdc + kin.x, vp = kin.v;
    this.h = h; this.vp = vp;
    const dz = this.dz, yTop = this.yTop, gv = this.gv;
    for (let j = 0; j < NZ; j++) {
      const t = g.rowType[j];
      if (t === 0) { dz[j] = g.dzRec; yTop[j] = -(g.nd - j) * g.dzRec; gv[j] = 0; }
      else if (t === 1) { const m = j - g.nd; dz[j] = h / g.ng; yTop[j] = m * h / g.ng; gv[j] = vp * m / g.ng; }
      else { const m = j - g.nd - g.ng; dz[j] = g.dzRec; yTop[j] = h + m * g.dzRec; gv[j] = vp; }
    }
    yTop[NZ] = yTop[NZ - 1] + dz[NZ - 1]; gv[NZ] = g.nb ? vp : vp;
    let Vtot = 0;
    for (let i = 0; i < g.nr; i++) {
      let H = 0;
      for (let j = 0; j < NZ; j++) {
        const c = i * NZ + j;
        if (g.act[c]) { this.V[c] = g.Az[i] * dz[j]; Vtot += this.V[c]; H += dz[j]; }
        else this.V[c] = 0;
      }
      this.Hcol[i] = H;
    }
    this.Vtot = Vtot;
    return Vtot;
  }

  derive() {
    const g = this.g, fuel = this.fuel, RF = PS_RU / fuel.M;
    for (let c = 0; c < g.N; c++) {
      if (!g.act[c]) continue;
      const m = this.m[c], YF = this.mF[c] / m, YP = this.mP[c] / m, YA = 1 - YF - YP;
      const R = YA * PS_AIR.R + YP * PS_PROD.R + YF * RF;
      const cp = YA * PS_AIR.cp + YP * PS_PROD.cp + YF * fuel.cpV;
      this.R[c] = R; this.cp[c] = cp;
      this.T[c] = this.E[c] / (m * cp);
      this.rho[c] = this.p0 / (R * this.T[c]);
      const i = (c / g.NZ) | 0;
      this.u[c] = this.mu[c] / m; this.w[c] = this.mw[c] / m; this.v[c] = this.mL[c] / (m * g.rc[i]);
      this.k[c] = Math.max(this.mk[c] / m, 1e-4); this.c[c] = Math.min(1, Math.max(0, this.mc[c] / m));
    }
  }

  airMass() {
    let a = 0;
    for (let c = 0; c < this.g.N; c++) if (this.g.act[c]) a += this.m[c] - this.mF[c] - this.mP[c];
    return a;
  }

  cyclePos(ca) { return ((ca % 720) + 720) % 720; }

  valveLift(pos, open, close, lmax) {
    if (pos < open || pos > close) return 0;
    const s = Math.sin(Math.PI * (pos - open) / (close - open));
    return lmax * s * s;
  }

  // Top active cell of column i.
  roofCell(i) {
    const NZ = this.g.NZ;
    for (let j = 0; j < NZ; j++) if (this.g.act[i * NZ + j]) return i * NZ + j;
    return -1;
  }

  // Cell at (r, y), or -1 inside a wall.
  cellAt(r, y) {
    const g = this.g;
    if (r < 0) r = -r;
    if (r >= g.R) return -1;
    const i = Math.min(g.nr - 1, Math.floor(r / g.dr));
    const yTop = this.yTop;
    if (y < yTop[0] || y >= yTop[g.NZ]) return -1;
    let lo = 0, hi = g.NZ - 1;
    while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (yTop[mid] <= y) lo = mid; else hi = mid - 1; }
    const c = i * g.NZ + lo;
    return g.act[c] ? c : -1;
  }

  nearestCell(r, y) {
    const g = this.g, NZ = g.NZ;
    const i = Math.min(g.nr - 1, Math.max(0, Math.floor(Math.abs(r) / g.dr)));
    let best = -1, bd = Infinity;
    for (let j = 0; j < NZ; j++) {
      const c = i * NZ + j;
      if (!g.act[c]) continue;
      const d = Math.abs(this.yTop[j] + 0.5 * this.dz[j] - y);
      if (d < bd) { bd = d; best = c; }
    }
    return best;
  }

  // ------------------------------------------------------------ main step

  step() {
    const p = this.p, g = this.g, N = g.N, NZ = g.NZ, nr = g.nr;
    const degPerS = p.rpm * 6;
    // Time step: transport, front and crank angle limits.
    let dt = Math.min(this.dtLimit || 1e-4, 0.25 / degPerS, this.dt * 1.25);
    dt = Math.max(dt, 2e-8);
    this.dt = dt;
    const caOld = this.ca, caNew = caOld + dt * degPerS;
    const posOld = this.cyclePos(caOld), posNew = this.cyclePos(caNew);
    const st = this.stats;
    this.qStep = 0;

    this.derive();
    this.sources(dt, posOld, posNew);
    this.derive();

    // ---- isentropic compression to the new volume
    const Vold = this.Vtot;
    const Vnew = this.geometryAt(caNew);
    const p0 = this.p0;
    let sumV = 0;
    for (let c = 0; c < N; c++) if (g.act[c]) { this.Vt[c] = this.m[c] * this.R[c] * this.T[c] / p0; sumV += this.Vt[c]; }
    let x = Math.pow(sumV / Vnew, 1.3);
    for (let it = 0; it < 30; it++) {
      let f = -Vnew, df = 0;
      const lx = Math.log(x);
      for (let c = 0; c < N; c++) {
        if (!g.act[c]) continue;
        const gm = this.cp[c] / (this.cp[c] - this.R[c]);
        const vi = this.Vt[c] * Math.exp(-lx / gm);
        f += vi; df -= vi / (gm * x);
      }
      const dx = -f / df;
      x += dx;
      if (x <= 0) x = 1e-3;
      if (Math.abs(dx) < 1e-13 * x) break;
    }
    const pNew = p0 * x, lx = Math.log(x);
    for (let c = 0; c < N; c++) {
      if (!g.act[c]) continue;
      const gm = this.cp[c] / (this.cp[c] - this.R[c]);
      const ratio = Math.exp(-lx / gm);            // volume ratio
      this.Vt[c] *= ratio;
      this.E[c] *= Math.exp(lx * (gm - 1) / gm);
      this.mk[c] *= Math.pow(ratio, -2 / 3);       // rapid distortion of k
      this.S[c] = (this.Vt[c] - this.V[c]) / dt;
    }
    st.work += 0.5 * (p0 + pNew) * (Vnew - Vold);
    if (posNew < 180 || posNew >= 540) st.workGross += 0.5 * (p0 + pNew) * (Vnew - Vold);
    this.p0 = pNew;
    this.derive();

    // ---- centrifugal acceleration from swirl
    for (let c = 0; c < N; c++) {
      if (!g.act[c]) continue;
      const i = (c / NZ) | 0;
      this.u[c] += dt * this.v[c] * this.v[c] / g.rc[i];
    }

    this.project(dt);
    this.advect(dt);
    this.derive();
    this.diffuse(dt);
    this.turbulence(dt);
    this.derive();

    this.ca = caNew; this.t += dt;
    this.bookkeeping(posOld, posNew, dt);
  }

  // ------------------------------------------------------------ projection

  project(dt) {
    const g = this.g, N = g.N, NZ = g.NZ, nr = g.nr, W = NZ + 1;
    const a = this.band, Fe = this.Fe, Fd = this.Fd, phi = this.phi, act = g.act, dz = this.dz;
    a.fill(0);
    for (let c = 0; c < N; c++) phi[c] = act[c] ? this.S[c] : 0;
    let maxDiag = 0;
    for (let i = 0; i < nr; i++) for (let j = 0; j < NZ; j++) {
      const c = i * NZ + j;
      if (!act[c]) { a[c * W] = 1; continue; }
      if (i + 1 < nr && act[c + NZ]) {
        const n = c + NZ, A = 2 * Math.PI * g.rw[i + 1] * dz[j];
        const Fs = A * 0.5 * (this.u[c] + this.u[n]);
        const beta = A * dt / (0.5 * (this.rho[c] + this.rho[n]) * g.dr);
        Fe[c] = Fs; phi[c] -= Fs; phi[n] += Fs;
        a[c * W] += beta; a[n * W] += beta; a[n * W + NZ] -= beta;
        this.bE = this.bE || new Float64Array(N); this.bE[c] = beta;
      } else Fe[c] = 0;
      if (j + 1 < NZ && act[c + 1]) {
        const n = c + 1, A = g.Az[i];
        const Fs = A * (0.5 * (this.w[c] + this.w[n]) - this.gv[j + 1]);
        const beta = A * dt / (0.5 * (this.rho[c] + this.rho[n]) * 0.5 * (dz[j] + dz[j + 1]));
        Fd[c] = Fs; phi[c] -= Fs; phi[n] += Fs;
        a[c * W] += beta; a[n * W] += beta; a[n * W + 1] -= beta;
        this.bD = this.bD || new Float64Array(N); this.bD[c] = beta;
      } else Fd[c] = 0;
    }
    for (let c = 0; c < N; c++) if (act[c]) maxDiag = Math.max(maxDiag, a[c * W]);
    let pin = -1;
    for (let c = 0; c < N && pin < 0; c++) if (act[c]) pin = c;
    a[pin * W] += 1e4 * maxDiag;
    psBandCholesky(a, N, NZ);
    psBandSolve(a, N, NZ, phi);
    // Corrected fluxes and cell velocities.
    const bE = this.bE, bD = this.bD;
    const gr = this.tmpGr || (this.tmpGr = new Float64Array(N)), gz = this.tmpGz || (this.tmpGz = new Float64Array(N));
    const nrc = this.tmpNr || (this.tmpNr = new Float64Array(N)), nzc = this.tmpNz || (this.tmpNz = new Float64Array(N));
    gr.fill(0); gz.fill(0); nrc.fill(0); nzc.fill(0);
    let maxRate = 0;
    const out = this.tmpOut || (this.tmpOut = new Float64Array(N));
    out.fill(0);
    for (let i = 0; i < nr; i++) for (let j = 0; j < NZ; j++) {
      const c = i * NZ + j;
      if (!act[c]) continue;
      if (i + 1 < nr && act[c + NZ]) {
        const n = c + NZ, d = phi[n] - phi[c];
        Fe[c] -= bE[c] * d;
        const du = -dt * d / (0.5 * (this.rho[c] + this.rho[n]) * g.dr);
        gr[c] += du; gr[n] += du; nrc[c]++; nrc[n]++;
        if (Fe[c] > 0) out[c] += Fe[c]; else out[n] -= Fe[c];
      }
      if (j + 1 < NZ && act[c + 1]) {
        const n = c + 1, d = phi[n] - phi[c];
        Fd[c] -= bD[c] * d;
        const dw = -dt * d / (0.5 * (this.rho[c] + this.rho[n]) * 0.5 * (dz[j] + dz[j + 1]));
        gz[c] += dw; gz[n] += dw; nzc[c]++; nzc[n]++;
        if (Fd[c] > 0) out[c] += Fd[c]; else out[n] -= Fd[c];
      }
    }
    for (let c = 0; c < N; c++) {
      if (!act[c]) continue;
      if (nrc[c]) this.u[c] += gr[c] / nrc[c];
      if (nzc[c]) this.w[c] += gz[c] / nzc[c];
      this.mu[c] = this.m[c] * this.u[c]; this.mw[c] = this.m[c] * this.w[c];
      const rate = out[c] / this.Vt[c];
      if (rate > maxRate) maxRate = rate;
    }
    this.outRate = out;
    this.maxFluxRate = maxRate;
  }

  // ------------------------------------------------------------ transport

  advect(dt) {
    const g = this.g, N = g.N, NZ = g.NZ, nr = g.nr, act = g.act;
    const old = this.old, nf = PS_FIELDS.length;
    for (const f of PS_FIELDS) old[f].set(this[f]);
    // Donor limiter so no cell can give away more than 90 % of its gas.
    const lim = this.tmpLim || (this.tmpLim = new Float64Array(N));
    for (let c = 0; c < N; c++) lim[c] = act[c] ? Math.min(1, 0.9 * this.Vt[c] / Math.max(this.outRate[c] * dt, 1e-300)) : 0;
    const fields = PS_FIELDS.map(f => this[f]), olds = PS_FIELDS.map(f => old[f]);
    const om = old.m;
    const spec = (k, c) => olds[k][c] / om[c];
    const minmod = (a, b) => (a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b));
    const face = (P, Nn, U, flux) => {
      // P donor, Nn acceptor, U upstream of donor (or -1).
      const mass = flux * dt * lim[P] * om[P] / this.Vt[P];
      fields[0][P] -= mass; fields[0][Nn] += mass;
      for (let k = 1; k < nf; k++) {
        const sP = spec(k, P);
        let sf = sP;
        if (U >= 0) sf += 0.5 * minmod(sP - spec(k, U), spec(k, Nn) - sP);
        const q = mass * sf;
        fields[k][P] -= q; fields[k][Nn] += q;
      }
    };
    for (let i = 0; i < nr; i++) for (let j = 0; j < NZ; j++) {
      const c = i * NZ + j;
      if (!act[c]) continue;
      if (i + 1 < nr && act[c + NZ]) {
        const n = c + NZ, F = this.Fe[c];
        if (F > 0) face(c, n, i > 0 && act[c - NZ] ? c - NZ : -1, F);
        else if (F < 0) face(n, c, i + 2 < nr && act[n + NZ] ? n + NZ : -1, -F);
      }
      if (j + 1 < NZ && act[c + 1]) {
        const n = c + 1, F = this.Fd[c];
        if (F > 0) face(c, n, j > 0 && act[c - 1] ? c - 1 : -1, F);
        else if (F < 0) face(n, c, j + 2 < NZ && act[n + 1] ? n + 1 : -1, -F);
      }
    }
    for (let c = 0; c < N; c++) {
      if (!act[c]) continue;
      if (this.m[c] < 1e-14) this.m[c] = 1e-14;
      if (this.mF[c] < 0) this.mF[c] = 0;
      if (this.mP[c] < 0) this.mP[c] = 0;
      if (this.mF[c] + this.mP[c] > this.m[c]) { const s = this.m[c] / (this.mF[c] + this.mP[c]); this.mF[c] *= s; this.mP[c] *= s; }
      if (this.mc[c] < 0) this.mc[c] = 0;
      if (this.mk[c] < 0) this.mk[c] = 0;
    }
  }

  // Implicit diffusion along radial lines then axial lines. Momentum, species,
  // temperature and k share the turbulent viscosity; walls add friction and heat.
  diffuse(dt) {
    const g = this.g, N = g.N, NZ = g.NZ, nr = g.nr, act = g.act, p = this.p, dz = this.dz;
    const nut = this.tmpNut || (this.tmpNut = new Float64Array(N));
    const Lm = this.tmpL || (this.tmpL = new Float64Array(N));
    this.lengthScales(Lm);
    for (let c = 0; c < N; c++) if (act[c]) {
      const nu = psViscosity(this.T[c]) / this.rho[c];
      nut[c] = nu + 0.55 * Math.sqrt(this.k[c]) * Lm[c];
    }
    const Y1 = this.tmpY1 || (this.tmpY1 = new Float64Array(N)), Y2 = this.tmpY2 || (this.tmpY2 = new Float64Array(N));
    for (let c = 0; c < N; c++) if (act[c]) { Y1[c] = this.mF[c] / this.m[c]; Y2[c] = this.mP[c] / this.m[c]; }
    const vars = [this.u, this.w, this.v, Y1, Y2, this.T, this.k];
    const wallT = [p.tHead, p.tLiner, p.tPiston];
    const tri = this.tri;
    let qWall = 0;
    const self = this;
    // Solve one line of cells. cells: indices; G[k] conductance between k and k+1 (kg/s);
    // wallLo/wallHi: wall conductance and type at the line ends.
    function line(cells, n, G, walls) {
      for (let vi = 0; vi < vars.length; vi++) {
        const X = vars[vi];
        for (let k = 0; k < n; k++) {
          const c = cells[k];
          const wgt = vi === 5 ? self.m[c] * self.cp[c] : self.m[c];
          const gcoef = vi === 5 ? self.cp[c] / 0.9 : vi >= 3 ? 1 / 0.9 : 1;
          tri.di[k] = wgt / dt; tri.rhs[k] = wgt / dt * X[c];
          tri.lo[k] = 0; tri.up[k] = 0;
          if (k > 0) { const gg = G[k - 1] * gcoef; tri.lo[k] = -gg; tri.di[k] += gg; }
          if (k < n - 1) { const gg = G[k] * gcoef; tri.up[k] = -gg; tri.di[k] += gg; }
          // Walls on this cell
          const wl = walls[k];
          if (wl) {
            if (vi <= 2) {               // friction toward the wall velocity
              tri.di[k] += wl.f;
              if (vi === 1) tri.rhs[k] += wl.f * wl.vw;
            } else if (vi === 5) {        // heat transfer
              tri.di[k] += wl.h; tri.rhs[k] += wl.h * wl.Tw;
            }
          }
        }
        psTridiag(n, tri.lo, tri.di, tri.up, tri.rhs, tri.x, tri.work);
        for (let k = 0; k < n; k++) {
          const c = cells[k];
          if (vi === 5) {
            const wl = walls[k];
            if (wl) qWall += wl.h * (tri.x[k] - wl.Tw) * dt;
          }
          X[c] = tri.x[k];
        }
      }
    }
    const cells = new Int32Array(Math.max(nr, NZ)), G = new Float64Array(Math.max(nr, NZ));
    const walls = new Array(Math.max(nr, NZ));
    const ht = p.heatTransfer;
    const wallFor = (c, A, type, normalDist, vw) => {
      const speed = Math.hypot(this.u[c], this.w[c] - (type === 2 ? this.vp : 0), this.v[c]);
      const sk = Math.sqrt(this.k[c]);
      const h = ht * (0.006 * this.rho[c] * this.cp[c] * (speed + 1.5 * sk) + psConductivity(this.T[c]) / Math.max(normalDist, 1e-4)) * A;
      const f = this.rho[c] * 0.004 * (speed + sk) * A;
      return { h, f, Tw: wallT[type], vw: vw };
    };
    const merge = (a, b) => (a ? (b ? { h: a.h + b.h, f: a.f + b.f, Tw: (a.h * a.Tw + b.h * b.Tw) / (a.h + b.h), vw: (a.f * a.vw + b.f * b.vw) / (a.f + b.f) } : a) : b);
    // Radial lines.
    for (let j = 0; j < NZ; j++) {
      let n = 0;
      const typeSide = g.rowType[j] === 0 ? 0 : g.rowType[j] === 1 ? 1 : 2;
      for (let i = 0; i <= nr; i++) {
        const c = i * NZ + j;
        if (i < nr && act[c]) {
          cells[n] = c;
          if (n > 0) {
            const pc = cells[n - 1];
            const A = 2 * Math.PI * g.rw[i] * dz[j];
            G[n - 1] = 0.5 * (this.rho[pc] * nut[pc] + this.rho[c] * nut[c]) * A / g.dr;
          }
          walls[n] = null;
          n++;
        } else {
          if (n > 0) {
            const last = cells[n - 1], il = (last / NZ) | 0;
            const A = 2 * Math.PI * g.rw[il + 1] * dz[j];
            walls[n - 1] = wallFor(last, A, typeSide, 0.5 * g.dr, typeSide === 2 ? this.vp : 0);
            line(cells, n, G, walls);
          }
          n = 0;
        }
      }
    }
    // Axial lines.
    for (let i = 0; i < nr; i++) {
      let n = 0;
      for (let j = 0; j <= NZ; j++) {
        const c = i * NZ + j;
        if (j < NZ && act[c]) {
          cells[n] = c;
          walls[n] = null;
          if (n > 0) {
            const pc = cells[n - 1];
            G[n - 1] = 0.5 * (this.rho[pc] * nut[pc] + this.rho[c] * nut[c]) * g.Az[i] / (0.5 * (dz[j - 1] + dz[j]));
          } else {
            walls[0] = wallFor(c, g.Az[i], 0, 0.5 * dz[j], 0);   // roof: head
          }
          n++;
        } else {
          if (n > 0) {
            const last = cells[n - 1], jl = last - i * NZ;
            walls[n - 1] = merge(walls[n - 1], wallFor(last, g.Az[i], 2, 0.5 * dz[jl], this.vp));
            line(cells, n, G, walls);
          }
          n = 0;
        }
      }
    }
    // Liner friction and heat already done in the radial pass for gap rows.
    for (let c = 0; c < N; c++) {
      if (!act[c]) continue;
      const m = this.m[c];
      const YF = Math.max(0, Y1[c]), YP = Math.max(0, Y2[c]);
      this.mF[c] = m * YF; this.mP[c] = m * Math.min(YP, 1 - YF);
      const i = (c / NZ) | 0;
      this.mu[c] = m * this.u[c]; this.mw[c] = m * this.w[c]; this.mL[c] = m * this.v[c] * g.rc[i];
      this.mk[c] = m * Math.max(this.k[c], 1e-4);
      const fuel = this.fuel, YA = 1 - YF - YP;
      const cp = YA * PS_AIR.cp + YP * PS_PROD.cp + YF * fuel.cpV;
      this.E[c] = m * cp * this.T[c];
    }
    this.stats.qWall += qWall;
  }

  lengthScales(Lm) {
    const g = this.g, NZ = g.NZ, act = g.act, dz = this.dz;
    for (let i = 0; i < g.nr; i++) {
      let top = -1, bot = -1;
      for (let j = 0; j < NZ; j++) if (act[i * NZ + j]) { if (top < 0) top = j; bot = j; }
      if (top < 0) continue;
      const yT = this.yTop[top], yB = this.yTop[bot] + dz[bot];
      const Lmax = Math.min(0.2 * this.Hcol[i], 0.04 * g.B);
      for (let j = top; j <= bot; j++) {
        const c = i * NZ + j;
        if (!act[c]) continue;
        const yc = this.yTop[j] + 0.5 * dz[j];
        let yw = Math.min(yc - yT, yB - yc, g.R - g.rc[i]);
        // radial walls of the recess
        if (i + 1 < g.nr && !act[c + NZ]) yw = Math.min(yw, g.rw[i + 1] - g.rc[i]);
        Lm[c] = Math.max(5e-5, Math.min(0.41 * yw, Lmax));
      }
    }
  }

  turbulence(dt) {
    const g = this.g, N = g.N, NZ = g.NZ, nr = g.nr, act = g.act, dz = this.dz;
    const Lm = this.tmpL;
    for (let i = 0; i < nr; i++) for (let j = 0; j < NZ; j++) {
      const c = i * NZ + j;
      if (!act[c]) continue;
      let s2 = 0, walls = 0, wallY = 1;
      const add = (n, d) => {
        const du = this.u[n] - this.u[c], dw = this.w[n] - this.w[c], dv = this.v[n] - this.v[c];
        s2 += 0.5 * (du * du + dw * dw + dv * dv) / (d * d);
      };
      const wall = d => { walls++; wallY = Math.min(wallY, 0.5 * d); };
      // Wall faces go through a wall function instead of the velocity jump.
      if (i > 0) add(c - NZ, g.dr);
      if (i + 1 < nr && act[c + NZ]) add(c + NZ, g.dr); else wall(g.dr);
      if (j > 0 && act[c - 1]) add(c - 1, 0.5 * (dz[j] + dz[j - 1])); else wall(dz[j]);
      if (j + 1 < NZ && act[c + 1]) add(c + 1, 0.5 * (dz[j] + dz[j + 1])); else wall(dz[j]);
      const u0 = this.u[c] / g.rc[i];
      s2 += 2 * u0 * u0;
      const k = this.k[c], L = Lm[c];
      const nut = 0.55 * Math.sqrt(k) * L;
      let P = nut * s2;
      if (walls) {
        const ut = 0.063 * Math.hypot(this.u[c], this.w[c] - (g.rowType[j] === 0 ? 0 : this.vp), this.v[c]);
        P += walls * ut * ut * ut / (0.41 * Math.max(wallY, 1e-5)) * 0.5;
      }
      const kn = (k + dt * P) / (1 + dt * 0.1643 * Math.sqrt(k) / L);
      this.mk[c] = this.m[c] * Math.min(Math.max(kn, 1e-4), 5e3);
    }
  }

  // ------------------------------------------------------------ sources

  sources(dt, posOld, posNew) {
    this.valves(dt, posNew);
    this.injection(dt, posOld, posNew);
    this.droplets(dt);
    this.filmEvaporation(dt);
    this.combustion(dt, posOld, posNew);
  }

  addGas(c, dm, YF, YP, T, cp, ux, uy, vt, k) {
    const i = (c / this.g.NZ) | 0;
    this.m[c] += dm; this.mF[c] += dm * YF; this.mP[c] += dm * YP; this.E[c] += dm * cp * T;
    this.mu[c] += dm * ux; this.mw[c] += dm * uy; this.mL[c] += dm * vt * this.g.rc[i]; this.mk[c] += dm * k;
  }

  removeGas(c, dm) {
    const f = Math.min(0.9, dm / this.m[c]);
    const out = { m: this.m[c] * f, mF: this.mF[c] * f, mP: this.mP[c] * f, E: this.E[c] * f };
    for (const key of PS_FIELDS) this[key][c] *= 1 - f;
    return out;
  }

  totalMass() { let s = 0; for (let c = 0; c < this.g.N; c++) if (this.g.act[c]) s += this.m[c]; return s; }

  valves(dt, pos) {
    const p = this.p, g = this.g, fuel = this.fuel;
    const B = p.bore;
    const Lin = this.valveLift(pos, 360 - p.ivo, 540 + p.ivc, p.ivLift);
    const Lex = this.valveLift(pos, 180 - p.evo, 360 + p.evc, p.evLift);
    this.liftIn = Lin; this.liftEx = Lex;
    if (Lin <= 0 && Lex <= 0) return;
    const M = this.totalMass();
    let Tm = 0, Rm = 0, cpm = 0;
    for (let c = 0; c < g.N; c++) if (g.act[c]) { Tm += this.m[c] * this.T[c]; Rm += this.m[c] * this.R[c]; cpm += this.m[c] * this.cp[c]; }
    Tm /= M; Rm /= M; cpm /= M;
    const gm = cpm / (cpm - Rm);
    const p0 = this.p0;
    const limit = (pr) => 0.3 * Math.abs(pr - p0) * this.Vtot / (gm * Rm * Tm);
    const RF = PS_RU / fuel.M;
    if (Lin > 0) {
      const D = p.ivDiam * B, A = Math.min(2 * Math.PI * D * Lin, 2 * Math.PI * D * D / 4 * 0.9);
      const Cd = 0.65;
      if (p.pIntake >= p0) {
        // Fresh charge: air + fuel vapour (port injection); liquid enters as droplets.
        const pfi = p.injection === 'puerto';
        const f = pfi ? 1 / (p.lambda * fuel.AFR) : 0;
        const vf = fuel.liquid ? p.pfiVapor : 1;
        const YF = f * vf / (1 + f * vf);
        const R = (1 - YF) * PS_AIR.R + YF * RF, cp = (1 - YF) * PS_AIR.cp + YF * fuel.cpV;
        let dm = psOrifice(Cd, A, p.pIntake, p.tIntake, R, cp / (cp - R), p0) * dt;
        dm = Math.min(dm, limit(p.pIntake));
        const rhoIn = p.pIntake / (R * p.tIntake);
        const U = dm / dt / (rhoIn * A);
        this.injectThroughValve(g.inCols, dm, U, true, YF, R, cp, f, vf);
      } else {
        // Backflow into the intake port.
        let dm = psOrifice(Cd, A, p0, Tm, Rm, gm, p.pIntake) * dt;
        dm = Math.min(dm, limit(p.pIntake));
        this.extractThroughValve(g.inCols, dm, this.intakeBuf);
      }
    }
    if (Lex > 0) {
      const D = p.evDiam * B, A = Math.min(2 * Math.PI * D * Lex, 2 * Math.PI * D * D / 4 * 0.9);
      const Cd = 0.65;
      if (p0 > p.pExhaust) {
        let dm = psOrifice(Cd, A, p0, Tm, Rm, gm, p.pExhaust) * dt;
        dm = Math.min(dm, limit(p.pExhaust));
        const out = { m: 0, mF: 0, mP: 0, E: 0 };
        this.extractThroughValve(g.exCols, dm, out);
        if (out.m > 0) {
          const er = this.exhaustRes, a = out.m / (out.m + 2e-4);
          const YF = out.mF / out.m, YP = out.mP / out.m;
          const cpo = (1 - YF - YP) * PS_AIR.cp + YP * PS_PROD.cp + YF * fuel.cpV;
          er.YF += a * (YF - er.YF); er.YP += a * (YP - er.YP); er.T += a * (0.9 * out.E / (out.m * cpo) - er.T);
          this.stats.exhaustOut = (this.stats.exhaustOut || 0) + out.m;
          this.stats.hcOut = (this.stats.hcOut || 0) + out.mF;
        }
      } else {
        // Exhaust gas flows back in.
        const er = this.exhaustRes;
        const R = (1 - er.YF - er.YP) * PS_AIR.R + er.YP * PS_PROD.R + er.YF * RF;
        const cp = (1 - er.YF - er.YP) * PS_AIR.cp + er.YP * PS_PROD.cp + er.YF * fuel.cpV;
        let dm = psOrifice(Cd, A, p.pExhaust, er.T, R, cp / (cp - R), p0) * dt;
        dm = Math.min(dm, limit(p.pExhaust));
        const U = dm / dt / (p.pExhaust / (R * er.T) * A);
        this.spreadInflow(g.exCols, dm, er.YF, er.YP, er.T, cp, U, -1);
      }
    }
  }

  // Distributes an inflow over the roof cells of the given columns, as a jet.
  spreadInflow(cols, dm, YF, YP, T, cp, U, dir) {
    const g = this.g, p = this.p;
    let At = 0;
    for (const i of cols) At += g.Az[i];
    const ja = p.jetAngle * PS_DEG;
    const ux = dir * U * Math.sin(ja), uy = U * Math.cos(ja), vt = U * p.swirl;
    const k = 0.012 * U * U;
    for (const i of cols) {
      const c = this.roofCell(i);
      if (c < 0) continue;
      this.addGas(c, dm * g.Az[i] / At, YF, YP, T, cp, ux, uy, vt, k);
    }
  }

  injectThroughValve(cols, dm, U, fresh, YFfresh, R, cp, f, vf) {
    const buf = this.intakeBuf, fuel = this.fuel, p = this.p;
    // Gas pushed back into the port earlier comes in first.
    if (buf.m > 0) {
      const take = Math.min(buf.m, dm), s = take / buf.m;
      const YF = buf.mF / buf.m, YP = buf.mP / buf.m;
      const cpb = (1 - YF - YP) * PS_AIR.cp + YP * PS_PROD.cp + YF * fuel.cpV;
      const T = buf.E / (buf.m * cpb);
      this.spreadInflow(cols, take, YF, YP, T, cpb, U, 1);
      buf.m -= take; buf.mF -= buf.mF * s; buf.mP -= buf.mP * s; buf.E -= buf.E * s;
      if (buf.m < 1e-12) buf.m = buf.mF = buf.mP = buf.E = 0;
      dm -= take;
    }
    if (dm <= 0) return;
    this.spreadInflow(cols, dm, YFfresh, 0, p.tIntake, cp, U, 1);
    const air = dm * (1 - YFfresh);
    this.stats.fuelIn += dm * YFfresh;
    // Port injected liquid rides in with the charge.
    if (p.injection === 'puerto' && fuel.liquid && vf < 1) {
      const liq = air * f * (1 - vf);
      this.stats.fuelIn += liq;
      this.spawnLiquid(liq, 'valve', U);
    }
  }

  extractThroughValve(cols, dm, sink) {
    const g = this.g;
    let At = 0;
    for (const i of cols) At += g.Az[i];
    for (const i of cols) {
      const c = this.roofCell(i);
      if (c < 0) continue;
      const o = this.removeGas(c, dm * g.Az[i] / At);
      sink.m += o.m; sink.mF += o.mF; sink.mP += o.mP; sink.E += o.E;
    }
  }

  // ------------------------------------------------------------ liquid fuel

  injection(dt, posOld, posNew) {
    const p = this.p, fuel = this.fuel;
    if (p.injection !== 'directa') return;
    const soiPos = (720 - p.diSoi) % 720;
    const crossed = (posOld <= soiPos && posNew > soiPos) || (posOld > posNew && (soiPos >= posOld || soiPos < posNew));
    if (crossed) {
      const air = Math.max(this.lastTrappedAir, 1e-6);
      const mInj = air / (p.lambda * fuel.AFR);
      const dp = Math.max(p.diPressure - this.p0, 1e6);
      const rho = fuel.liquid ? fuel.rhoL : 50;
      const rate = 0.7 * 1.46e-7 * Math.sqrt(2 * rho * dp) * (fuel.liquid ? 1 : 0.05);
      this.inj = { left: mInj, rate, U: 0.7 * Math.sqrt(2 * dp / rho), total: mInj };
      this.stats.fuelIn += mInj;
      this.stats.soiCa = this.ca;
    }
    if (this.inj && this.inj.left > 0) {
      const dm = Math.min(this.inj.left, this.inj.rate * dt);
      this.inj.left -= dm;
      if (fuel.liquid) this.spawnLiquid(dm, 'injector', this.inj.U);
      else {
        // Gas injection: fuel vapour jet from the injector.
        const c = this.roofCell(0);
        if (c >= 0) this.addGas(c, dm, 1, 0, 300, fuel.cpV, 0, Math.min(this.inj.U, 300), 0, 10);
      }
      if (this.inj.left <= 0) this.stats.eoiCa = this.ca;
    }
  }

  spawnLiquid(mass, where, U) {
    if (mass <= 0) return;
    const p = this.p, fuel = this.fuel, g = this.g;
    if (!this.parcelMass) {
      const air = Math.max(this.lastTrappedAir, 1e-6);
      this.parcelMass = air / (p.lambda * fuel.AFR) / 700;
    }
    this.injAcc += mass;
    while (this.injAcc >= this.parcelMass && this.parcels.length < 3000) {
      const pm = this.parcelMass;
      this.injAcc -= pm;
      let d, r, y, u, w, vt, T;
      const rr = Math.random();
      if (where === 'injector') {
        d = p.diDrop * Math.pow(-Math.log(1 - 0.98 * rr), 1 / 3.0) * 0.9;
        const cone = 0.5 * p.diCone * PS_DEG * (0.85 + 0.3 * Math.random());
        const Ui = U * (0.8 + 0.2 * Math.random());
        r = 0.0005 + 0.0005 * Math.random();
        const c0 = this.roofCell(0);
        y = this.yTop[c0 - 0] + 0.0004;
        u = Ui * Math.sin(cone); w = Ui * Math.cos(cone); vt = 0; T = 330;
      } else {
        d = p.pfiDrop * Math.pow(-Math.log(1 - 0.98 * rr), 1 / 3.0) * 0.9;
        const i = g.inCols[(Math.random() * g.inCols.length) | 0];
        r = g.rw[i] + g.dr * Math.random();
        const c0 = this.roofCell(i);
        y = this.yTop[c0 - i * g.NZ] + 0.0003;
        const ja = p.jetAngle * PS_DEG;
        u = 0.5 * U * Math.sin(ja); w = 0.5 * U * Math.cos(ja); vt = 0.5 * U * p.swirl; T = p.tIntake;
      }
      d = Math.max(3e-6, d);
      const md = fuel.rhoL * Math.PI * d * d * d / 6;
      this.parcels.push({ r, y, u, w, vt, d, T, n: pm / md, alive: true });
    }
  }

  droplets(dt) {
    const fuel = this.fuel;
    if (!fuel.liquid || !this.parcels.length) return;
    const p0 = this.p0, g = this.g, NZ = g.NZ, RF = PS_RU / fuel.M;
    const keep = [];
    for (const q of this.parcels) {
      let c = this.cellAt(q.r, q.y);
      if (c < 0) {
        // The piston or a wall moved over the drop: it becomes film on the nearest cell.
        c = this.nearestCell(q.r, q.y);
        this.film[c] += q.n * fuel.rhoL * Math.PI * q.d ** 3 / 6;
        this.filmWall[c] = q.y < this.yTop[c % NZ] ? 0 : 2;
        continue;
      }
      const Tg = this.T[c], rhoG = this.rho[c], mu = psViscosity(Tg), kg = psConductivity(0.5 * (Tg + q.T));
      const i = (c / NZ) | 0;
      const ug = this.u[c], wg = this.w[c], vg = this.v[c] * q.r / g.rc[i];
      const rel = Math.hypot(q.u - ug, q.w - wg, q.vt - vg);
      const Re = rhoG * rel * q.d / mu;
      const Cd = Re < 1000 ? 24 / Math.max(Re, 1e-6) * (1 + Math.pow(Re, 2 / 3) / 6) : 0.424;
      const md = fuel.rhoL * Math.PI * q.d ** 3 / 6;
      const tauV = Re > 1e-6 ? (4 * fuel.rhoL * q.d) / (3 * rhoG * Cd * rel) : fuel.rhoL * q.d * q.d / (18 * mu);
      const e = Math.exp(-dt / tauV);
      const nu = ug + (q.u - ug) * e, nw = wg + (q.w - wg) * e, nv = vg + (q.vt - vg) * e;
      const dmu = q.n * md;
      this.mu[c] -= dmu * (nu - q.u); this.mw[c] -= dmu * (nw - q.w);
      this.mL[c] -= dmu * (nv - q.vt) * q.r;
      // Evaporation (Spalding mass number, Ranz-Marshall).
      const pSat = PS_ATM * Math.exp(fuel.Lb * fuel.M / PS_RU * (1 / fuel.Tb - 1 / Math.min(q.T, fuel.Tc - 1)));
      const Xs = Math.min(0.995, pSat / p0);
      const Mg = PS_RU / this.R[c];
      const Ys = Xs * fuel.M / (Xs * fuel.M + (1 - Xs) * Mg);
      const YF = this.mF[c] / this.m[c];
      const BM = Math.max(0, (Ys - YF) / (1 - Ys));
      const Dv = 1.0e-5 * Math.pow(Tg / 300, 1.75) * (1e5 / p0);
      const Sc = mu / (rhoG * Dv), Pr = 0.7;
      const sq = Math.sqrt(Re);
      const Sh = 2 + 0.552 * sq * Math.cbrt(Sc), Nu = 2 + 0.552 * sq * Math.cbrt(Pr);
      let mdot = Math.PI * q.d * rhoG * Dv * Sh * Math.log1p(BM);
      const L = fuel.Lb * Math.pow(Math.max(0.01, (fuel.Tc - q.T) / (fuel.Tc - fuel.Tb)), 0.38);
      const zc = fuel.cpV * mdot / (Math.PI * q.d * kg * Nu);
      const corr = zc > 1e-8 ? zc / Math.expm1(zc) : 1;
      const Q = Math.PI * q.d * kg * Nu * (Tg - q.T) * corr;
      let dmd = Math.min(md, mdot * dt);
      let dT = (Q * dt - dmd * L) / (md * fuel.cL);
      // keep the droplet temperature update stable for small drops
      const tauT = md * fuel.cL / Math.max(Math.PI * q.d * kg * Nu, 1e-30);
      if (dt > 0.5 * tauT) dT *= 0.5 * tauT / dt * 2 * (1 - Math.exp(-dt / tauT));
      const mNew = md - dmd;
      const evap = dmd * q.n;
      if (evap > 0) {
        this.m[c] += evap; this.mF[c] += evap; this.E[c] += evap * fuel.cpV * q.T;
        this.mu[c] += evap * nu; this.mw[c] += evap * nw; this.mL[c] += evap * nv * q.r;
        this.mk[c] += evap * this.k[c];
      }
      this.E[c] -= Q * dt * q.n;
      q.T = Math.min(fuel.Tc - 2, q.T + dT);
      q.u = nu; q.w = nw; q.vt = nv;
      if (mNew <= md * 1e-4 || q.d < 1.5e-6) {
        // Last bit of a tiny drop flashes to vapour.
        const rest = mNew * q.n;
        this.m[c] += rest; this.mF[c] += rest; this.E[c] += rest * fuel.cpV * q.T;
        q.alive = false;
        continue;
      }
      q.d = Math.cbrt(6 * mNew / (Math.PI * fuel.rhoL));
      // Move; swirl motion of a ring is folded into r through conservation of
      // angular momentum (centrifugal drift).
      q.u += dt * q.vt * q.vt / Math.max(q.r, 1e-4);
      const nr = q.r + q.u * dt, ny = q.y + q.w * dt;
      let cr = nr;
      if (cr < 0) { cr = -cr; q.u = -q.u; }
      const cn = this.cellAt(cr, ny);
      if (cn < 0) {
        // Wall impingement: the drop joins the film of the cell it came from.
        this.film[c] += q.n * fuel.rhoL * Math.PI * q.d ** 3 / 6;
        this.filmWall[c] = ny < this.yTop[c % NZ] ? 0 : (cr >= g.R - 1e-9 ? 1 : 2);
        q.alive = false;
        continue;
      }
      q.vt *= q.r / Math.max(cr, 1e-5);
      q.r = cr; q.y = ny;
      keep.push(q);
    }
    this.parcels = keep.filter(q => q.alive);
  }

  filmEvaporation(dt) {
    const fuel = this.fuel, g = this.g, p = this.p;
    if (!fuel.liquid) return;
    const wallT = [p.tHead, p.tLiner, p.tPiston];
    let total = 0;
    for (let c = 0; c < g.N; c++) {
      if (!g.act[c] || this.film[c] <= 0) continue;
      const Tw = wallT[this.filmWall[c]];
      const pSat = PS_ATM * Math.exp(fuel.Lb * fuel.M / PS_RU * (1 / fuel.Tb - 1 / Math.min(Tw, fuel.Tc - 1)));
      const Xs = Math.min(0.995, pSat / this.p0);
      const Mg = PS_RU / this.R[c];
      const Ys = Xs * fuel.M / (Xs * fuel.M + (1 - Xs) * Mg);
      const YF = this.mF[c] / this.m[c];
      const BM = Math.max(0, (Ys - YF) / (1 - Ys));
      const speed = Math.hypot(this.u[c], this.w[c]) + 1.5 * Math.sqrt(this.k[c]);
      const hm = 0.006 * this.rho[c] * speed;           // kg/m2s
      const Afilm = Math.min(this.g.Az[(c / g.NZ) | 0], this.film[c] / (fuel.rhoL * 2e-5));
      // Diffusion into the gas, plus boiling when the wall is above the saturation temperature.
      const Tsat = 1 / (1 / fuel.Tb - Math.log(this.p0 / PS_ATM) * PS_RU / (fuel.Lb * fuel.M));
      const boil = Tw > Tsat ? 4000 * (Tw - Tsat) / fuel.Lb : 0;
      const dm = Math.min(this.film[c], (hm * Math.log1p(BM) + boil) * Afilm * dt);
      if (dm > 0) {
        this.film[c] -= dm;
        this.m[c] += dm; this.mF[c] += dm; this.E[c] += dm * fuel.cpV * Tw;
      }
      total += this.film[c];
    }
    this.filmTotal = total;
  }

  // ------------------------------------------------------------ combustion

  combustion(dt, posOld, posNew) {
    const p = this.p, g = this.g, N = g.N, NZ = g.NZ, nr = g.nr, act = g.act, fuel = this.fuel, st = this.stats;
    const p0 = this.p0, dz = this.dz;
    const q = this.q;
    for (let c = 0; c < N; c++) q[c] = 0;
    const burnCell = (c, f) => {
      // Burn the fraction f of the reactants left in cell c.
      const m = this.m[c], mF = this.mF[c], mA = m - mF - this.mP[c];
      const dF = f * Math.min(mF, mA / fuel.AFR);
      if (dF <= 0) return 0;
      this.mF[c] -= dF; this.mP[c] += dF * (1 + fuel.AFR);
      this.E[c] += dF * fuel.LHV;
      q[c] += dF * fuel.LHV;
      return dF;
    };
    // Spark: window of 1.5 ms after the spark angle.
    const sparkWin = 1.5e-3 * p.rpm * 6;
    const cyc = posNew >= 360 ? posNew - 720 : posNew, cycOld = posOld >= 360 ? posOld - 720 : posOld;
    const inWin = cyc >= -p.spark && cyc <= -p.spark + sparkWin;
    this.sparking = inWin && !this.ignited;
    let kernel = -1;
    if (inWin && !this.ignited) {
      const c = this.sparkCell();
      if (c >= 0) {
        const YF = this.mF[c] / this.m[c], YA = 1 - YF - this.mP[c] / this.m[c];
        const phi = YA > 0 ? YF / YA * fuel.AFR : 0;
        const Sl = psLaminarSpeed(fuel, phi, this.T[c], p0, this.mP[c] / this.m[c]);
        if (Sl > 0.05) {
          this.ignited = true; st.ignitionCa = this.ca;
          this.mc[c] = this.m[c]; this.c[c] = 1;
          kernel = c;
          st.uSpark = Math.sqrt(2 / 3 * this.k[c]);
        }
      }
    }
    if (kernel >= 0) burnCell(kernel, 1);
    if (!this.ignited && cyc > -p.spark + sparkWin && cycOld <= -p.spark + sparkWin) { st.misfire = true; }

    // Unburned gas temperature (mass weighted over cells the flame has not reached).
    let mu = 0, Tu = 0, Vb = 0;
    for (let c = 0; c < N; c++) if (act[c]) {
      if (this.c[c] < 0.05) { mu += this.m[c]; Tu += this.m[c] * this.T[c]; }
      Vb += this.c[c] * this.Vt[c] || this.c[c] * this.V[c];
    }
    Tu = mu > 0 ? Tu / mu : 800;
    this.Tu = Tu;
    const rf = Math.cbrt(3 * Vb / (2 * Math.PI));
    const fk = 1 - Math.exp(-rf / 0.004);
    let SmaxDl = 0;
    const cOld = this.tmpC || (this.tmpC = new Float64Array(N));
    cOld.set(this.c);
    const RF = PS_RU / fuel.M;
    let qStep = 0, burnedStep = 0;
    if (this.ignited || this.knockActive) {
      for (let i = 0; i < nr; i++) for (let j = 0; j < NZ; j++) {
        const c = i * NZ + j;
        if (!act[c]) continue;
        const ci = cOld[c];
        if (ci >= 0.999) continue;
        const m = this.m[c];
        // Residual fraction and equivalence ratio of the fresh gas, frozen when the
        // front arrives (afterwards the burned part's leftovers would skew them).
        if (ci < 0.02) {
          const mA0 = m - this.mF[c] - this.mP[c];
          this.mxr[c] = this.mP[c];
          this.mphi[c] = m * (mA0 > 0 ? this.mF[c] / mA0 * fuel.AFR : 9);
        }
        const xr = this.mxr[c] / m;
        // |grad c| for a front moving toward lower c.
        let gx = 0, gy = 0;
        if (i > 0) { const d = (ci - cOld[c - NZ]) / g.dr; if (d < 0) gx = d * d; }
        if (i + 1 < nr && act[c + NZ]) { const d = (cOld[c + NZ] - ci) / g.dr; if (d > 0) gx = Math.max(gx, d * d); }
        if (j > 0 && act[c - 1]) { const d = (ci - cOld[c - 1]) / (0.5 * (dz[j] + dz[j - 1])); if (d < 0) gy = d * d; }
        if (j + 1 < NZ && act[c + 1]) { const d = (cOld[c + 1] - ci) / (0.5 * (dz[j] + dz[j + 1])); if (d > 0) gy = Math.max(gy, d * d); }
        let grad = Math.sqrt(gx + gy);
        // Once the front has crossed the cell centre the rest of the cell burns out.
        if (ci > 0.5) grad = Math.max(grad, (1 - ci) / Math.min(g.dr, dz[j]));
        if (grad <= 0) continue;
        const phi = this.mphi[c] / m;
        const T = ci < 0.3 ? this.T[c] : Tu;
        const Sl = psLaminarSpeed(fuel, phi, T, p0, xr);
        if (Sl <= 0) continue;
        // Quenching in thin gaps: Peclet number of about 3 per wall.
        const alpha = psConductivity(T) / (this.rho[c] * this.cp[c]);
        const dq = 3 * alpha / Sl;
        const H = this.Hcol[i];
        const fq = Math.min(1, Math.max(0, (H - 2 * dq) / (2 * dq)));
        if (fq <= 0) continue;
        const up = Math.sqrt(2 / 3 * this.k[c]);
        const St = p.flameFactor * fq * Sl * (1 + 2.0 * fk * Math.pow(up / Sl, 0.7));
        const dl = Math.min(g.dr, 0.5 * (dz[j] + (j + 1 < NZ ? dz[j + 1] : dz[j])));
        SmaxDl = Math.max(SmaxDl, St / Math.min(dl, dz[j]));
        const dc = Math.min(1 - ci, St * grad * dt);
        if (dc <= 0) continue;
        const f = dc / (1 - ci);
        burnedStep += burnCell(c, f);
        this.mc[c] += m * dc;
      }
    }
    // Spark kernel cell finishes quickly.
    // Post flame oxidation of leftovers that meet air in hot burned gas.
    for (let c = 0; c < N; c++) {
      if (!act[c] || this.c[c] < 0.5 || this.T[c] < 1300) continue;
      const k = this.k[c], L = this.tmpL ? this.tmpL[c] : 1e-3;
      const eps = 0.1643 * Math.pow(k, 1.5) / Math.max(L, 5e-5);
      const rate = 4 * eps / k + 300;
      burnedStep += burnCell(c, Math.min(1, rate * dt));
    }
    // End gas autoignition (knock).
    const pos = posNew;
    for (let c = 0; c < N; c++) {
      if (!act[c]) continue;
      const mA = this.m[c] - this.mF[c] - this.mP[c];
      const phiK = mA > 0 ? this.mF[c] / mA * fuel.AFR : 9;
      if (this.c[c] < 0.05 && phiK > fuel.phiMin && phiK < fuel.phiMax && this.T[c] > 600) {
        const tau = psIgnitionDelay(p.octane, p0, this.T[c]);
        this.mik[c] += this.m[c] * dt / tau;
        if (this.mik[c] >= this.m[c]) {
          if (!st.knockCa) {
            st.knockCa = this.ca;
            let Tb = 0, mb = 0;
            for (let d = 0; d < N; d++) if (act[d] && this.c[d] > 0.5) { Tb += this.m[d] * this.T[d]; mb += this.m[d]; }
            const Tk = mb > 0 ? Tb / mb : this.T[c];
            st.knockFreq = psKnockFrequency(p.bore, Tk, 1.27, PS_PROD.R);
          }
          // Volumetric explosion: the whole pocket burns within ~20 us.
          const ci = this.c[c];
          const dc = (1 - ci) * Math.min(1, dt / 2e-5);
          const b = burnCell(c, dc / (1 - ci));
          st.knockFrac += b;
          burnedStep += b;
          this.mc[c] += this.m[c] * dc;
        }
      }
    }
    for (let c = 0; c < N; c++) if (act[c]) {
      qStep += q[c];
      const qv = q[c] / (dt * Math.max(this.V[c], 1e-15));
      this.flame[c] = Math.max(qv, this.flame[c] * Math.exp(-dt / 3e-4));
    }
    st.qComb += qStep; st.fuelBurned += burnedStep;
    this.qStep = qStep;
    this.frontLimit = SmaxDl;
  }

  sparkCell() {
    const g = this.g, c0 = this.roofCell(0);
    if (c0 < 0) return -1;
    const y = this.yTop[c0] + Math.min(0.003, 0.35 * this.Hcol[0]);
    return this.cellAt(0.3 * g.dr, y);
  }

  // ------------------------------------------------------------ diagnostics

  bookkeeping(posOld, posNew, dt) {
    const p = this.p, g = this.g, st = this.stats, fuel = this.fuel;
    // time step limits for the next step
    const lims = [1e-4];
    if (this.maxFluxRate > 0) lims.push(0.35 / this.maxFluxRate);
    if (this.frontLimit > 0) lims.push(0.3 / this.frontLimit);
    let umax = 0;
    for (const q of this.parcels) umax = Math.max(umax, Math.abs(q.u) / g.dr, Math.abs(q.w) / Math.max(this.dz[g.nd] || 1e-3, 2e-4));
    if (umax > 0) lims.push(0.8 / umax);
    this.dtLimit = Math.min(...lims);

    const degPerS = p.rpm * 6;
    if (this.p0 > st.pMax) { st.pMax = this.p0; st.caPMax = this.ca; }
    // squish velocity: mean radial velocity through the face at the recess edge
    let usq = 0;
    if (g.type !== 'plana') {
      const iEdge = Math.max(0, Math.min(g.nr - 2, Math.round(g.rRec / g.dr) - 1));
      let F = 0, A = 0;
      for (let j = g.nd; j < g.nd + g.ng; j++) { const c = iEdge * g.NZ + j; F += this.Fe[c]; A += 2 * Math.PI * g.rw[iEdge + 1] * this.dz[j]; }
      usq = A > 0 ? -F / A : 0;
    }
    this.usq = usq;
    if (usq > st.sqMax) st.sqMax = usq;
    const cyc = posNew >= 360 ? posNew - 720 : posNew;   // -360..360
    // MFB from burned fuel against fuel available at ignition
    const total = this.fuelAtSpark || 0;
    const mfb = total > 0 ? st.fuelBurned / total : 0;
    this.mfb = mfb;
    for (const [k, lv] of [['mfb10', 0.1], ['mfb50', 0.5], ['mfb90', 0.9]]) if (st[k] === null && mfb >= lv) st[k] = cyc;
    let film = this.filmTotal || 0, liq = 0;
    for (const q of this.parcels) liq += q.n * fuel.rhoL * Math.PI * q.d ** 3 / 6;
    liq += this.injAcc;               // injected but not yet released as a parcel
    this.liquidTotal = liq;
    st.filmMax = Math.max(st.filmMax, film); st.liquidMax = Math.max(st.liquidMax, liq);
    // Events
    const cross = (a) => (posOld < a && posNew >= a) || (posOld > posNew && (a >= posOld || a < posNew));
    if (cross(540 + p.ivc)) {
      st.massTrapped = this.totalMass();
      let mP = 0, mF = 0;
      for (let c = 0; c < g.N; c++) if (g.act[c]) { mP += this.mP[c]; mF += this.mF[c]; }
      st.airTrapped = st.massTrapped - mP - mF;
      st.residual = mP / st.massTrapped;
      this.lastTrappedAir = st.airTrapped;
      for (let c = 0; c < g.N; c++) this.mik[c] = 0;
    }
    if (cross((720 - p.spark) % 720)) {
      let mF = 0;
      for (let c = 0; c < g.N; c++) if (g.act[c]) mF += this.mF[c];
      this.fuelAtSpark = mF + liq + film;
      const fuelTot = mF + liq + film;
      st.lambdaTrapped = fuelTot > 0 ? st.airTrapped / (fuelTot * fuel.AFR) : Infinity;
      this.fuelAtSpark = fuelTot;
    }
    if (cross(180 - p.evo)) {
      let mF = 0;
      for (let c = 0; c < g.N; c++) if (g.act[c]) mF += this.mF[c];
      st.hcEvo = this.fuelAtSpark > 0 ? (mF + liq + film) / this.fuelAtSpark : 0;
    }
    // Trace (every ~0.5 deg)
    const tr = this.trace;
    const last = tr.ca.length ? tr.ca[tr.ca.length - 1] : -1e9;
    if (cyc - last >= 0.5 || cyc < last) {
      let Tm = 0, M = 0, U = 0;
      for (let c = 0; c < g.N; c++) if (g.act[c]) { Tm += this.m[c] * this.T[c]; M += this.m[c]; U += this.m[c] * this.k[c]; }
      tr.ca.push(cyc); tr.p.push(this.p0); tr.V.push(this.Vtot);
      tr.hrr.push(this.qStep / (dt * degPerS)); tr.mfb.push(mfb); tr.T.push(Tm / M);
      tr.u.push(Math.sqrt(2 / 3 * U / M)); tr.usq.push(usq);
      let vap = 0;
      for (let c = 0; c < g.N; c++) if (g.act[c]) vap += this.mF[c];
      tr.vap.push(vap); tr.liq.push(liq); tr.film.push(film); tr.burned.push(st.fuelBurned);
    }
    // Cycle boundary at gas exchange TDC
    if (cross(360)) this.endCycle();
  }

  endCycle() {
    const st = this.stats, p = this.p, fuel = this.fuel, g = this.g;
    st.imep = st.work / g.Vd;
    st.imepGross = st.workGross / g.Vd;
    st.efficiency = st.fuelIn > 0 ? st.workGross / (st.fuelIn * fuel.LHV) : 0;
    st.combustionEff = st.fuelIn > 0 ? st.qComb / (st.fuelIn * fuel.LHV) : 0;
    const rhoRef = p.pIntake / (PS_AIR.R * p.tIntake);
    st.volEff = st.airTrapped / (rhoRef * g.Vd);
    st.cycle = this.cycleCount;
    st.complete = this.trace.ca.length > 1000 && this.cycleCount > 0;
    this.prevCycle = { stats: st, trace: this.trace };
    if (st.complete || this.cycleCount > 0) this.history.push(st);
    if (this.history.length > 50) this.history.shift();
    this.cycleCount++;
    this.stats = this.newCycleStats();
    this.trace = this.newTrace();
    // New cycle: clear the flame state.
    for (let c = 0; c < g.N; c++) { this.mc[c] = 0; this.flame[c] = 0; this.mik[c] = 0; }
    this.ignited = false; this.knockActive = false; this.inj = null;
    this.fuelAtSpark = 0;
  }

  // Advance by deg crank degrees or until the time budget (ms) runs out.
  advance(deg, budgetMs) {
    const t0 = (typeof performance !== 'undefined' ? performance : Date).now();
    const target = this.ca + deg;
    let steps = 0;
    while (this.ca < target) {
      this.step(); steps++;
      if (budgetMs && (steps & 3) === 0 && (typeof performance !== 'undefined' ? performance : Date).now() - t0 > budgetMs) break;
    }
    return steps;
  }

  cycleAngle() { const pos = this.cyclePos(this.ca); return pos >= 360 ? pos - 720 : pos; }

  stroke() {
    const pos = this.cyclePos(this.ca);
    return pos < 180 ? 'Expansión' : pos < 360 ? 'Escape' : pos < 540 ? 'Admisión' : 'Compresión';
  }
}

if (typeof module !== 'undefined') {
  module.exports = {
    PS_FUELS, PS_DEFAULTS, PS_RESOLUTIONS, PS_AIR, PS_PROD, PS_ATM, PS_RU, PSCylinder, psBuildGeometry, psPiston,
    psOrifice, psLaminarSpeed, psIgnitionDelay, psKnockFrequency, psSquishEstimate, psBandCholesky, psBandSolve, psTridiag,
  };
}
