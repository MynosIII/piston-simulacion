// Behavioural checks for the single cylinder model. Run: node tests/physics.test.js
const assert = require('assert');
const P = require('../Piston-Sim-Physics.js');

const check = (name, fn) => {
  const t0 = Date.now();
  fn();
  console.log(`ok  ${name} (${((Date.now() - t0) / 1000).toFixed(1)} s)`);
};

const runTo = (sim, cyc) => {   // run until the cycle angle cyc (-360..360) is next reached
  const pos = ((cyc % 720) + 720) % 720;
  let target = sim.ca - sim.cyclePos(sim.ca) + pos;
  if (target <= sim.ca) target += 720;
  while (sim.ca < target) sim.step();
};

const internalEnergy = sim => {
  let U = 0;
  for (let c = 0; c < sim.g.N; c++) if (sim.g.act[c]) U += sim.m[c] * (sim.cp[c] - sim.R[c]) * sim.T[c];
  return U;
};

check('banded Cholesky solves a banded SPD system', () => {
  const n = 60, bw = 7, W = bw + 1, a = new Float64Array(n * W), full = [];
  for (let i = 0; i < n; i++) full.push(new Float64Array(n));
  for (let i = 0; i < n; i++) for (let d = 1; d <= bw && d <= i; d++) {
    const v = Math.random() - 0.5;
    full[i][i - d] = full[i - d][i] = v;
  }
  for (let i = 0; i < n; i++) { let s = 1; for (let j = 0; j < n; j++) if (j !== i) s += Math.abs(full[i][j]); full[i][i] = s; }
  for (let i = 0; i < n; i++) for (let d = 0; d <= bw && d <= i; d++) a[i * W + d] = full[i][i - d];
  const x = Float64Array.from({ length: n }, () => Math.random()), b = new Float64Array(n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) b[i] += full[i][j] * x[j];
  P.psBandCholesky(a, n, bw);
  P.psBandSolve(a, n, bw, b);
  for (let i = 0; i < n; i++) assert(Math.abs(b[i] - x[i]) < 1e-10, `x[${i}]`);
});

check('chamber volume matches the compression ratio', () => {
  for (const chamber of ['cupula', 'bowl', 'plana']) {
    const sim = new P.PSCylinder({ chamber });
    const g = sim.g, Vtdc = sim.geometryAt(0), Vbdc = sim.geometryAt(180);
    assert(Math.abs((Vbdc - Vtdc) / g.Vd - 1) < 1e-9, `${chamber}: swept volume`);
    assert(Math.abs(Vbdc / Vtdc / sim.p.cr - 1) < 1e-9, `${chamber}: ratio ${Vbdc / Vtdc}`);
  }
});

check('motored compression is isentropic and conserves mass', () => {
  const sim = new P.PSCylinder({ lambda: 1e9, heatTransfer: 0, spark: -100 });
  const M0 = sim.totalMass(), p0 = sim.p0, V0 = sim.Vtot;
  let cp = 0, R = 0;
  for (let c = 0; c < sim.g.N; c++) if (sim.g.act[c]) { cp += sim.m[c] * sim.cp[c]; R += sim.m[c] * sim.R[c]; }
  const gamma = cp / (cp - R);
  runTo(sim, 0);
  const expect = p0 * Math.pow(V0 / sim.Vtot, gamma);
  assert(Math.abs(sim.p0 / expect - 1) < 0.01, `p at TDC ${sim.p0} vs ${expect}`);
  assert(Math.abs(sim.totalMass() / M0 - 1) < 1e-9, 'mass');
});

check('first law over the closed part of the cycle', () => {
  const sim = new P.PSCylinder({ heatTransfer: 0, injection: 'puerto' });
  sim.derive();
  const U0 = internalEnergy(sim);
  let W = 0;
  const evo = 180 - sim.p.evo - 2;
  let pPrev = sim.p0, Vprev = sim.Vtot;
  while (sim.cycleAngle() < evo || sim.cycleAngle() > 300) {
    sim.step();
    W += 0.5 * (pPrev + sim.p0) * (sim.Vtot - Vprev);
    pPrev = sim.p0; Vprev = sim.Vtot;
  }
  const U1 = internalEnergy(sim), Q = sim.stats.qComb;
  assert(Q > 0, 'combustion released heat');
  const err = (U1 - U0 + W - Q) / Q;
  assert(Math.abs(err) < 0.03, `energy balance error ${(err * 100).toFixed(2)} %`);
});

check('laminar flame speed and ignition delay correlations', () => {
  const gas = P.PS_FUELS.gasolina;
  const s0 = P.psLaminarSpeed(gas, 1, 298, P.PS_ATM, 0);
  assert(s0 > 0.3 && s0 < 0.45, `S_L at 1 atm ${s0}`);
  assert(P.psLaminarSpeed(gas, 1, 700, 20 * P.PS_ATM, 0) > s0, 'faster when hot');
  assert(P.psLaminarSpeed(gas, 1, 700, 20 * P.PS_ATM, 0.2) < P.psLaminarSpeed(gas, 1, 700, 20 * P.PS_ATM, 0), 'slower with residuals');
  assert(P.psLaminarSpeed(gas, 0.3, 700, 20 * P.PS_ATM, 0) === 0, 'outside flammability');
  assert(P.psIgnitionDelay(95, 40 * P.PS_ATM, 900) < P.psIgnitionDelay(95, 40 * P.PS_ATM, 800), 'delay shrinks with T');
  assert(P.psIgnitionDelay(100, 40 * P.PS_ATM, 900) > P.psIgnitionDelay(85, 40 * P.PS_ATM, 900), 'delay grows with octane');
  const f = P.psKnockFrequency(0.086, 2500, 1.27, 290.5);
  assert(f > 5500 && f < 7500, `knock frequency ${f}`);
});

check('squish velocity is close to the Heywood estimate', () => {
  const sim = new P.PSCylinder({ lambda: 1e9, spark: -100 });
  runTo(sim, -10);
  const g = sim.g;
  const Vrec = g.Vc - g.Ac * g.hTdc;
  const pist = P.psPiston(sim.p, sim.ca, sim.omega);
  const est = P.psSquishEstimate(g.B, 2 * g.rRec, Vrec, g.hTdc + pist.x, Math.abs(pist.v));
  assert(sim.usq > 0.5 * est && sim.usq < 1.5 * est, `squish ${sim.usq.toFixed(2)} m/s vs estimate ${est.toFixed(2)} m/s`);
});

check('fuel is conserved through injection, evaporation and wall film', () => {
  const sim = new P.PSCylinder({ injection: 'directa', diSoi: 120, spark: 0 });
  let vap0 = 0;
  for (let c = 0; c < sim.g.N; c++) if (sim.g.act[c]) vap0 += sim.mF[c];
  runTo(sim, -60);
  assert(sim.parcels.length > 0 || sim.filmTotal > 0, 'liquid present after injection');
  runTo(sim, -5);
  let vap = 0;
  for (let c = 0; c < sim.g.N; c++) if (sim.g.act[c]) vap += sim.mF[c];
  const liq = sim.parcels.reduce((s, q) => s + q.n * sim.fuel.rhoL * Math.PI * q.d ** 3 / 6, 0);
  const total = vap + liq + sim.injAcc + sim.filmTotal + sim.stats.fuelBurned;
  const expect = vap0 + sim.stats.fuelIn;
  assert(Math.abs(total / expect - 1) < 1e-6, `fuel ${total} vs ${expect}`);
});

check('a full firing cycle gives sensible engine numbers', () => {
  const sim = new P.PSCylinder({});
  while (sim.history.length < 2) { sim.step(); assert(isFinite(sim.p0), 'pressure finite'); }
  const st = sim.history[1];
  assert(st.imep > 8e5 && st.imep < 13e5, `IMEP ${st.imep}`);
  assert(st.mfb50 > -5 && st.mfb50 < 15, `MFB50 ${st.mfb50}`);
  assert(st.combustionEff > 0.9, `combustion efficiency ${st.combustionEff}`);
  assert(st.efficiency > 0.3 && st.efficiency < 0.45, `indicated efficiency ${st.efficiency}`);
  assert(st.residual > 0.01 && st.residual < 0.12, `residual ${st.residual}`);
  assert(st.volEff > 0.75 && st.volEff < 1.05, `volumetric efficiency ${st.volEff}`);
  assert(!st.knockCa, 'no knock at default settings');
});

check('advanced spark with low octane knocks', () => {
  const sim = new P.PSCylinder({ cr: 12.5, spark: 38, octane: 87 });
  while (sim.history.length < 1) sim.step();
  assert(sim.history[0].knockCa, 'knock detected');
});
