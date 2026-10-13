/* Piston-Sim-UI.js
 * Page for Piston-Sim-Physics.js: the cylinder section, the plots, the
 * readouts and the controls. Plain DOM and canvas, no libraries.
 */
'use strict';

const PSUI = (() => {
  const $ = (s, el = document) => el.querySelector(s);
  const h = (tag, attrs, ...kids) => {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === 'class') el.className = v;
      else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
      else if (v !== false && v != null) el.setAttribute(k, v === true ? '' : v);
    }
    for (const kid of kids.flat()) if (kid != null) el.append(kid.nodeType ? kid : document.createTextNode(kid));
    return el;
  };
  const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
  const lerp = (a, b, t) => a + (b - a) * t;
  const mix = (c1, c2, t) => [lerp(c1[0], c2[0], t), lerp(c1[1], c2[1], t), lerp(c1[2], c2[2], t)];

  const PRESETS = [
    { id: 'puerto', name: 'Inyección en el puerto', meta: 'λ 1 · mezcla homogénea', params: {} },
    { id: 'directa', name: 'Inyección directa', meta: 'inyecta en la admisión', params: { injection: 'directa', diSoi: 300 } },
    { id: 'estratificada', name: 'Carga estratificada', meta: 'λ 1.4 · inyección tardía', params: { injection: 'directa', diSoi: 125, lambda: 1.4, pIntake: 1.0e5, spark: 16 } },
    { id: 'knock', name: 'Golpe (knock)', meta: 'RC 12.5 · avance 38°', params: { cr: 12.5, spark: 38, octane: 87 } },
    { id: 'pobre', name: 'Mezcla pobre', meta: 'λ 1.5 · avance 30°', params: { lambda: 1.5, spark: 30 } },
    { id: 'parcial', name: 'Carga parcial', meta: '0.45 bar · 1500 rpm', params: { pIntake: 0.45e5, rpm: 1500, spark: 32 } },
    { id: 'squish', name: 'Sin squish', meta: 'cámara plana', params: { chamber: 'plana', spark: 24 } },
    { id: 'etanol', name: 'Etanol', meta: 'RC 12', params: { fuel: 'etanol', cr: 12 } },
    { id: 'gnc', name: 'GNC (metano)', meta: 'gas, sin gotas', params: { fuel: 'metano', spark: 20 } },
  ];

  const FIELDS = {
    temperatura: { label: 'Temperatura', unit: 'K', min: 300, max: 2800 },
    velocidad: { label: 'Velocidad', unit: 'm/s', min: 0, max: 30 },
    turbulencia: { label: "Turbulencia u'", unit: 'm/s', min: 0, max: 12 },
    dosado: { label: 'Dosado φ', unit: '', min: 0, max: 2 },
    swirl: { label: 'Swirl', unit: 'm/s', min: -15, max: 15 },
  };

  const pct = v => `${Math.round(v * 100)} %`;
  const GROUPS = [
    {
      title: 'Operación', hint: 'Cambian en vivo.', items: [
        { key: 'rpm', label: 'Régimen', min: 800, max: 7000, step: 50, fmt: v => `${v} rpm` },
        { key: 'pIntake', label: 'Presión de admisión', min: 0.3e5, max: 1.8e5, step: 0.01e5, fmt: v => `${(v / 1e5).toFixed(2)} bar` },
        { key: 'lambda', label: 'Lambda (λ)', min: 0.6, max: 2.2, step: 0.01, fmt: v => v.toFixed(2) },
        { key: 'spark', label: 'Avance de encendido', min: 0, max: 55, step: 1, fmt: v => `${v}° APMS` },
        { key: 'octane', label: 'Octanaje', min: 70, max: 125, step: 1, fmt: v => `${v}` },
        { key: 'tIntake', label: 'Temp. de admisión', min: 260, max: 380, step: 1, fmt: v => `${Math.round(v - 273.15)} °C` },
        { key: 'pExhaust', label: 'Presión de escape', min: 0.9e5, max: 2.0e5, step: 0.01e5, fmt: v => `${(v / 1e5).toFixed(2)} bar` },
      ],
    },
    {
      title: 'Inyección', hint: 'El combustible líquido se ve como gotas blancas y película celeste.', items: [
        { key: 'injection', label: 'Tipo', options: [['puerto', 'En el puerto'], ['directa', 'Directa']] },
        { key: 'diSoi', label: 'Inicio (directa)', min: 20, max: 360, step: 1, fmt: v => `${v}° APMS` },
        { key: 'diPressure', label: 'Presión (directa)', min: 30e5, max: 350e5, step: 5e5, fmt: v => `${Math.round(v / 1e5)} bar` },
        { key: 'diDrop', label: 'Gota media (directa)', min: 5e-6, max: 60e-6, step: 1e-6, fmt: v => `${Math.round(v * 1e6)} µm` },
        { key: 'diCone', label: 'Ángulo del cono', min: 20, max: 140, step: 1, fmt: v => `${v}°` },
        { key: 'pfiVapor', label: 'Vaporizado en el puerto', min: 0, max: 1, step: 0.01, fmt: pct },
        { key: 'pfiDrop', label: 'Gota media (puerto)', min: 20e-6, max: 250e-6, step: 5e-6, fmt: v => `${Math.round(v * 1e6)} µm` },
      ],
    },
    {
      title: 'Válvulas y flujo', hint: 'Ángulos respecto del PMS de cruce (admisión) y del PMI (escape).', items: [
        { key: 'ivo', label: 'Abre admisión', min: -30, max: 60, step: 1, fmt: v => `${v}° APMS` },
        { key: 'ivc', label: 'Cierra admisión', min: 0, max: 100, step: 1, fmt: v => `${v}° DPMI` },
        { key: 'evo', label: 'Abre escape', min: 0, max: 100, step: 1, fmt: v => `${v}° APMI` },
        { key: 'evc', label: 'Cierra escape', min: -30, max: 60, step: 1, fmt: v => `${v}° DPMS` },
        { key: 'ivLift', label: 'Alzada admisión', min: 0.003, max: 0.014, step: 0.0001, fmt: v => `${(v * 1000).toFixed(1)} mm` },
        { key: 'swirl', label: 'Swirl del puerto', min: 0, max: 1.5, step: 0.01, fmt: v => v.toFixed(2) },
        { key: 'jetAngle', label: 'Ángulo del chorro', min: 0, max: 75, step: 1, fmt: v => `${v}°` },
      ],
    },
    {
      title: 'Combustión y paredes', items: [
        { key: 'flameFactor', label: 'Factor de llama', min: 0.3, max: 2.5, step: 0.01, fmt: v => `× ${v.toFixed(2)}` },
        { key: 'heatTransfer', label: 'Transferencia de calor', min: 0, max: 2.5, step: 0.01, fmt: v => `× ${v.toFixed(2)}` },
        { key: 'tHead', label: 'Temp. culata', min: 320, max: 650, step: 5, fmt: v => `${Math.round(v - 273.15)} °C` },
        { key: 'tLiner', label: 'Temp. camisa', min: 320, max: 600, step: 5, fmt: v => `${Math.round(v - 273.15)} °C` },
        { key: 'tPiston', label: 'Temp. pistón', min: 320, max: 700, step: 5, fmt: v => `${Math.round(v - 273.15)} °C` },
      ],
    },
    {
      title: 'Geometría y combustible', hint: 'Estos cambios reinician la simulación.', reset: true, items: [
        { key: 'fuel', label: 'Combustible', options: Object.entries(PS_FUELS).map(([k, f]) => [k, f.label]) },
        { key: 'bore', label: 'Diámetro', min: 0.06, max: 0.11, step: 0.0005, fmt: v => `${(v * 1000).toFixed(1)} mm` },
        { key: 'stroke', label: 'Carrera', min: 0.05, max: 0.11, step: 0.0005, fmt: v => `${(v * 1000).toFixed(1)} mm` },
        { key: 'rod', label: 'Biela', min: 0.1, max: 0.2, step: 0.001, fmt: v => `${(v * 1000).toFixed(0)} mm` },
        { key: 'cr', label: 'Relación de compresión', min: 7, max: 16, step: 0.1, fmt: v => `${v.toFixed(1)} : 1` },
        { key: 'chamber', label: 'Cámara', options: [['cupula', 'Cúpula en la culata'], ['bowl', 'Bowl en el pistón'], ['plana', 'Plana (sin squish)']] },
        { key: 'squishClear', label: 'Huelgo de squish', min: 0.0005, max: 0.004, step: 0.0001, fmt: v => `${(v * 1000).toFixed(1)} mm` },
        { key: 'squishFrac', label: 'Área de squish', min: 0.05, max: 0.75, step: 0.01, fmt: pct },
        { key: 'resolution', label: 'Malla', options: [['rapida', 'Rápida'], ['normal', 'Normal'], ['fina', 'Fina']] },
      ],
    },
  ];

  // ------------------------------------------------------------ colours

  const C_AIR = [38, 62, 96], C_MIX = [46, 128, 255], C_RICH = [96, 60, 240], C_GAS = [128, 132, 138];
  const C_FLAME = [255, 138, 31], C_HOT = [255, 214, 110];
  function compositionColor(sim, c, out, o, qRef) {
    const m = sim.m[c], YF = sim.mF[c] / m, YP = sim.mP[c] / m, YA = Math.max(0, 1 - YF - YP);
    const fresh = YA + YF;
    const phi = YA > 1e-6 ? YF / YA * sim.fuel.AFR : 3;
    let col = phi < 1 ? mix(C_AIR, C_MIX, Math.sqrt(clamp(phi, 0, 1))) : mix(C_MIX, C_RICH, clamp((phi - 1) / 1.2, 0, 1));
    col = mix(col, C_GAS, clamp(1 - fresh, 0, 1));
    // Reaction zone: heat release rate against this frame's strongest cells.
    const q = sim.flame[c] / qRef;
    if (q > 0.05) {
      const a = clamp(Math.log10(q / 0.05) / 1.3, 0, 1);
      col = mix(col, mix(C_FLAME, C_HOT, clamp(a * 1.6 - 0.8, 0, 1)), Math.min(1, a * 1.2));
    }
    out[o] = col[0]; out[o + 1] = col[1]; out[o + 2] = col[2];
  }

  // Perceptual dark-to-bright ramp (approximates "inferno").
  const RAMP = [[8, 6, 24], [62, 12, 102], [140, 34, 104], [212, 72, 66], [248, 146, 22], [252, 232, 140]];
  const DIVERGE = [[40, 90, 220], [70, 74, 82], [230, 120, 40]];
  function ramp(t, stops) {
    t = clamp(t, 0, 1) * (stops.length - 1);
    const k = Math.min(stops.length - 2, Math.floor(t));
    return mix(stops[k], stops[k + 1], t - k);
  }
  function fieldValue(sim, field, c) {
    switch (field) {
      case 'temperatura': return sim.T[c];
      case 'velocidad': return Math.hypot(sim.u[c], sim.w[c]);
      case 'turbulencia': return Math.sqrt(2 / 3 * sim.k[c]);
      case 'dosado': { const m = sim.m[c], YF = sim.mF[c] / m, YA = 1 - YF - sim.mP[c] / m; return YA > 1e-6 ? YF / YA * sim.fuel.AFR : 2; }
      case 'swirl': return sim.v[c];
      default: return 0;
    }
  }
  function fieldColor(sim, field, c, out, o) {
    const f = FIELDS[field], v = fieldValue(sim, field, c);
    const t = (v - f.min) / (f.max - f.min);
    const col = field === 'swirl' || field === 'dosado' ? ramp(t, DIVERGE) : ramp(t, RAMP);
    out[o] = col[0]; out[o + 1] = col[1]; out[o + 2] = col[2];
  }

  // ------------------------------------------------------------ app

  class App {
    constructor() {
      this.preset = 'puerto';
      this.params = Object.assign({}, PS_DEFAULTS);
      this.sim = new PSCylinder(this.params);
      this.running = true;
      this.speed = 40;           // crank degrees per second of wall time
      this.view = { field: 'temperatura', arrows: true, drops: true, front: true, zoom: true };
      this.colL = null; this.colR = null;
      this.last = performance.now();
      this.simRate = 0;
      this.build();
      this.resize();
      addEventListener('resize', () => this.resize());
      requestAnimationFrame(t => this.frame(t));
    }

    build() {
      const bar = $('#bar');
      this.playBtn = h('button', { class: 'btn primary', onclick: () => this.toggle(), 'aria-pressed': 'true' }, 'Pausa');
      const step = h('button', { class: 'btn', onclick: () => { this.running = false; this.updatePlay(); this.sim.advance(1, 200); this.draw(); }, title: 'Avanzar 1° (→)' }, '+1°');
      const reset = h('button', { class: 'btn', onclick: () => this.reset() }, 'Reiniciar');
      const speed = h('select', { class: 'sel', 'aria-label': 'Velocidad', onchange: e => { this.speed = +e.target.value; } },
        [[10, '10 °/s'], [40, '40 °/s'], [120, '120 °/s'], [360, '360 °/s'], [1e6, 'Máxima']].map(([v, t]) => h('option', { value: v, selected: v === this.speed }, t)));
      this.angleOut = h('span', { class: 'angle' });
      bar.append(h('div', { class: 'brand' }, h('strong', {}, 'Pistón Lab'), h('span', {}, 'un cilindro, en detalle')),
        h('div', { class: 'transport' }, this.playBtn, step, reset, speed), this.angleOut);
      const list = $('#presets');
      this.presetBtns = PRESETS.map(pr => h('button', {
        class: 'preset', role: 'tab', 'aria-selected': pr.id === this.preset ? 'true' : 'false',
        onclick: () => this.loadPreset(pr.id),
      }, h('span', { class: 'preset-name' }, pr.name), h('span', { class: 'preset-meta' }, pr.meta)));
      list.append(...this.presetBtns);

      // view toolbar
      const tools = $('#view-tools');
      const fieldSel = h('select', { class: 'sel', 'aria-label': 'Campo de la derecha', onchange: e => { this.view.field = e.target.value; this.drawLegend(); } },
        Object.entries(FIELDS).map(([k, f]) => h('option', { value: k, selected: k === this.view.field }, f.label)));
      const chk = (key, label) => h('label', { class: 'chk' }, h('input', { type: 'checkbox', checked: this.view[key], onchange: e => { this.view[key] = e.target.checked; } }), label);
      const zoom = h('label', { class: 'chk' }, h('input', { type: 'checkbox', checked: this.view.zoom, onchange: e => { this.view.zoom = e.target.checked; this.resize(); this.draw(); } }), 'Zoom en la cámara');
      tools.append(h('span', { class: 'muted' }, 'Izquierda: composición · Derecha:'), fieldSel, chk('arrows', 'Flechas'), chk('drops', 'Gotas'), chk('front', 'Frente'), zoom);

      this.cv = $('#cyl'); this.ctx = this.cv.getContext('2d');
      this.off = document.createElement('canvas'); this.offCtx = this.off.getContext('2d');
      this.plots = { p: $('#plot-p'), pv: $('#plot-pv'), hrr: $('#plot-hrr'), fuel: $('#plot-fuel') };
      this.legend = $('#legend');
      this.live = $('#live'); this.cycle = $('#cycle');
      this.buildControls();
      this.drawLegend();
      addEventListener('keydown', e => {
        if (e.target.closest && e.target.closest('input,select,textarea')) return;
        if (e.code === 'Space') { e.preventDefault(); this.toggle(); }
        if (e.code === 'ArrowRight') { e.preventDefault(); this.running = false; this.updatePlay(); this.sim.advance(1, 200); this.draw(); }
      });
    }

    buildControls() {
      const root = $('#controls-groups');
      root.innerHTML = '';
      this.inputs = {};
      for (const grp of GROUPS) {
        const box = h('section', { class: 'group' }, h('h2', {}, grp.title), grp.hint ? h('p', { class: 'hint' }, grp.hint) : null);
        for (const it of grp.items) {
          const id = `in-${it.key}`;
          let input, out = null;
          const val = this.params[it.key] == null ? this.sim.p[it.key] : this.params[it.key];
          if (it.options) {
            input = h('select', { id }, it.options.map(([v, t]) => h('option', { value: v, selected: String(v) === String(val) }, t)));
            input.addEventListener('change', () => this.change(it, input.value, grp.reset));
          } else {
            out = h('output', { for: id }, it.fmt(val));
            input = h('input', { id, type: 'range', min: it.min, max: it.max, step: it.step, value: val });
            input.addEventListener('input', () => { out.textContent = it.fmt(+input.value); if (!grp.reset) this.change(it, +input.value, false); });
            input.addEventListener('change', () => { if (grp.reset) this.change(it, +input.value, true); });
          }
          this.inputs[it.key] = { input, out, it };
          box.append(h('div', { class: 'row' }, h('label', { for: id }, it.label), input, out || h('span')));
        }
        root.append(box);
      }
    }

    syncControls() {
      for (const [key, { input, out, it }] of Object.entries(this.inputs)) {
        const v = this.sim.p[key];
        input.value = v;
        if (out) out.textContent = it.fmt(+v);
      }
    }

    change(it, value, needsReset) {
      this.params[it.key] = value;
      if (it.key === 'fuel') this.params.octane = null;
      if (needsReset) this.reset();
      else this.sim.setParam(it.key, value);
    }

    loadPreset(id) {
      const pr = PRESETS.find(x => x.id === id);
      this.preset = id;
      this.params = Object.assign({}, PS_DEFAULTS, pr.params);
      this.presetBtns.forEach((b, k) => b.setAttribute('aria-selected', PRESETS[k].id === id ? 'true' : 'false'));
      this.reset();
    }

    reset() {
      this.sim = new PSCylinder(Object.assign({}, this.params));
      this.colL = this.colR = null;
      this.syncControls();
      this.resize();
      this.draw();
    }

    toggle() { this.running = !this.running; this.updatePlay(); }
    updatePlay() { this.playBtn.textContent = this.running ? 'Pausa' : 'Seguir'; this.playBtn.setAttribute('aria-pressed', String(this.running)); }

    frame(t) {
      const dtw = Math.min(0.1, (t - this.last) / 1000);
      this.last = t;
      if (this.running) {
        const ca0 = this.sim.ca;
        this.sim.advance(Math.max(0.05, this.speed * dtw), 14);
        const rate = (this.sim.ca - ca0) / Math.max(dtw, 1e-3);
        this.simRate = lerp(this.simRate, rate, 0.05);
      }
      this.draw();
      requestAnimationFrame(tt => this.frame(tt));
    }

    // ---------------------------------------------------------- layout

    resize() {
      const cv = this.cv, sim = this.sim, g = sim.g, p = sim.p;
      const wrap = cv.parentElement;
      const cssW = Math.max(280, Math.min(wrap.clientWidth, 760));
      // Physical window: head to a little below the piston at BDC, or just the chamber.
      const head = this.view.zoom ? 0.008 : 0.02, wall = 0.009;
      const yMin = -(g.nd * g.dzRec) - head;
      const yMax = this.view.zoom ? g.hTdc + g.nb * g.dzRec + Math.max(0.3 * p.stroke, 0.012) : g.hTdc + p.stroke + 0.62 * p.bore;
      const wPhys = p.bore + 2 * wall + 0.012, hPhys = yMax - yMin;
      const s = Math.min(cssW / wPhys, (this.view.zoom ? 520 : 680) / hPhys);
      const cssH = Math.round(hPhys * s);
      const dpr = Math.min(2, window.devicePixelRatio || 1);
      cv.style.width = `${Math.round(wPhys * s)}px`; cv.style.height = `${cssH}px`;
      cv.width = Math.round(wPhys * s * dpr); cv.height = Math.round(cssH * dpr);
      this.L = { s, dpr, yMin, yMax, wall, head, cx: wPhys * s / 2, W: wPhys * s, H: cssH };
      for (const pc of Object.values(this.plots)) {
        const w = pc.parentElement.clientWidth || 300, hh = pc.dataset.h ? +pc.dataset.h : 150;
        pc.style.width = `${w}px`; pc.style.height = `${hh}px`;
        pc.width = Math.round(w * dpr); pc.height = Math.round(hh * dpr);
      }
    }

    X(r) { return this.L.cx + r * this.L.s; }
    Y(y) { return (y - this.L.yMin) * this.L.s; }

    // ---------------------------------------------------------- drawing

    draw() {
      const sim = this.sim;
      this.drawCylinder();
      this.drawPlots();
      this.drawReadouts();
      const pos = sim.cycleAngle();
      const fmtA = pos >= 0 ? `${pos.toFixed(0)}° DPMS` : `${(-pos).toFixed(0)}° APMS`;
      const rate = this.running ? ` · ${this.simRate < 1 ? this.simRate.toFixed(1) : Math.round(this.simRate)} °/s` : ' · en pausa';
      this.angleOut.textContent = `${sim.stroke()} · ${fmtA}${rate}`;
    }

    drawCylinder() {
      const { ctx, sim, L } = this, g = sim.g, p = sim.p, s = L.s, dpr = L.dpr;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, L.W, L.H);
      const R = g.R, X = r => this.X(r), Y = y => this.Y(y);
      const metal = '#3a4047', metalDark = '#2a2f35', metalEdge = '#59616a';
      // Liner
      const yLinerEnd = g.hTdc + p.stroke + 0.62 * p.bore;
      ctx.fillStyle = metal;
      ctx.fillRect(X(-R - L.wall), Y(0), L.wall * s, (yLinerEnd) * s);
      ctx.fillRect(X(R), Y(0), L.wall * s, (yLinerEnd) * s);
      // Coolant passages
      ctx.fillStyle = '#16303f';
      ctx.fillRect(X(-R - L.wall * 0.7), Y(0.004), L.wall * 0.35 * s, (p.stroke * 0.7) * s);
      ctx.fillRect(X(R + L.wall * 0.35), Y(0.004), L.wall * 0.35 * s, (p.stroke * 0.7) * s);
      // Head block
      const yHeadTop = L.yMin;
      ctx.fillStyle = metal;
      ctx.fillRect(X(-R - L.wall), Y(yHeadTop), (2 * R + 2 * L.wall) * s, (-yHeadTop) * s);
      // Piston
      const yc = sim.h, ph = 0.6 * p.bore;
      const grad = ctx.createLinearGradient(X(-R), 0, X(R), 0);
      grad.addColorStop(0, '#4d555d'); grad.addColorStop(0.5, '#7b848d'); grad.addColorStop(1, '#4d555d');
      ctx.fillStyle = grad;
      ctx.fillRect(X(-R + 0.0002), Y(yc), (2 * R - 0.0004) * s, ph * s);
      ctx.fillStyle = metalDark;
      for (let k = 0; k < 3; k++) ctx.fillRect(X(-R + 0.0002), Y(yc + 0.003 + k * 0.0035), (2 * R - 0.0004) * s, Math.max(1, 0.0012 * s));
      // Pin and rod
      const pinY = yc + 0.36 * p.bore, pinR = 0.11 * p.bore;
      const a = p.stroke / 2, th = sim.ca * Math.PI / 180;
      const rodAng = Math.asin(a * Math.sin(th) / p.rod);
      ctx.save();
      ctx.beginPath(); ctx.rect(0, Y(yc + 0.12 * p.bore), L.W, L.H); ctx.clip();
      ctx.strokeStyle = '#8a939c'; ctx.lineWidth = 0.22 * p.bore * s; ctx.lineCap = 'round';
      ctx.beginPath(); ctx.moveTo(X(0), Y(pinY)); ctx.lineTo(X(0) + Math.sin(rodAng) * p.rod * s, Y(pinY) + Math.cos(rodAng) * p.rod * s); ctx.stroke();
      ctx.restore();
      ctx.fillStyle = '#9aa3ac'; ctx.beginPath(); ctx.arc(X(0), Y(pinY), pinR * s, 0, 2 * Math.PI); ctx.fill();
      ctx.fillStyle = metalDark; ctx.beginPath(); ctx.arc(X(0), Y(pinY), pinR * 0.55 * s, 0, 2 * Math.PI); ctx.fill();

      // Gas
      this.drawGas();

      // Wall film
      if (sim.fuel.liquid) {
        ctx.fillStyle = 'rgba(150, 225, 255, 0.95)';
        for (let c = 0; c < g.N; c++) {
          const fm = sim.film[c];
          if (!g.act[c] || fm < 2e-9) continue;
          const i = (c / g.NZ) | 0, j = c % g.NZ;
          const th2 = clamp(1 + Math.log10(fm / 2e-9) * 1.2, 1, 4);
          const y0 = sim.yTop[j], y1 = y0 + sim.dz[j];
          const wtype = sim.filmWall[c];
          for (const sg of [-1, 1]) {
            if (wtype === 1 || (wtype === 2 && i + 1 < g.nr && !g.act[c + g.NZ] && sim.g.rowType[j] !== 1)) {
              const xr = sg > 0 ? X(g.rw[i + 1]) - th2 : X(-g.rw[i + 1]);
              ctx.fillRect(xr, Y(y0), th2, sim.dz[j] * s);
            } else {
              const yy = wtype === 0 ? Y(y0) : Y(y1) - th2;
              const xa = sg > 0 ? X(g.rw[i]) : X(-g.rw[i + 1]);
              ctx.fillRect(xa, yy, g.dr * s, th2);
            }
          }
        }
      }
      // Droplets
      if (this.view.drops && sim.parcels.length) {
        ctx.fillStyle = 'rgba(240, 248, 255, 0.9)';
        for (const q of sim.parcels) {
          const rad = clamp(q.d * 1e6 * 0.045, 0.7, 4);
          for (const sg of [-1, 1]) {
            ctx.beginPath(); ctx.arc(X(sg * q.r), Y(q.y), rad, 0, 2 * Math.PI); ctx.fill();
          }
        }
      }
      // Flame front: edges between burned (c >= 0.5) and fresh cells
      if (this.view.front && (sim.ignited || sim.stats.knockCa)) {
        ctx.strokeStyle = 'rgba(255, 236, 160, 0.9)'; ctx.lineWidth = 1.5;
        ctx.beginPath();
        for (let i = 0; i < g.nr; i++) for (let j = 0; j < g.NZ; j++) {
          const c = i * g.NZ + j;
          if (!g.act[c]) continue;
          const b = sim.c[c] >= 0.5;
          if (i + 1 < g.nr && g.act[c + g.NZ] && (sim.c[c + g.NZ] >= 0.5) !== b && sim.mF[b ? c + g.NZ : c] > 1e-4 * sim.m[c]) {
            for (const sg of [-1, 1]) { ctx.moveTo(X(sg * g.rw[i + 1]), Y(sim.yTop[j])); ctx.lineTo(X(sg * g.rw[i + 1]), Y(sim.yTop[j] + sim.dz[j])); }
          }
          if (j + 1 < g.NZ && g.act[c + 1] && (sim.c[c + 1] >= 0.5) !== b && sim.mF[b ? c + 1 : c] > 1e-4 * sim.m[c]) {
            const y = Y(sim.yTop[j + 1]);
            for (const sg of [-1, 1]) { ctx.moveTo(X(sg * g.rw[i]), y); ctx.lineTo(X(sg * g.rw[i + 1]), y); }
          }
        }
        ctx.stroke();
      }
      // Velocity arrows
      if (this.view.arrows) this.drawArrows();

      // Valves
      this.drawValve(g.inR, sim.liftIn || 0, '#6c8fb3', p.ivDiam);
      this.drawValve(g.exR, sim.liftEx || 0, '#b37a6c', p.evDiam);
      // Spark plug
      const c0 = sim.roofCell(0);
      const roofY = sim.yTop[c0];
      const sparkY = roofY + Math.min(0.003, 0.35 * sim.Hcol[0]);
      ctx.fillStyle = '#c9ced3';
      ctx.fillRect(X(-0.0009), Y(roofY) - 2, 0.0018 * s, (sparkY - roofY) * s + 2);
      ctx.fillRect(X(-0.0026), Y(sparkY) - 1, 0.0052 * s * 0.5, Math.max(2, 0.0006 * s));
      ctx.fillStyle = '#e8e2d0'; ctx.fillRect(X(-0.002), Y(yHeadTop), 0.004 * s, (roofY - yHeadTop) * s);
      if (sim.sparking || (sim.ignited && sim.stats.ignitionCa != null && sim.ca - sim.stats.ignitionCa < 3)) {
        const gx = X(0), gy = Y(sparkY + 0.0006);
        const rg = ctx.createRadialGradient(gx, gy, 0, gx, gy, 0.006 * s);
        rg.addColorStop(0, 'rgba(255,255,255,1)'); rg.addColorStop(0.3, 'rgba(160,190,255,0.7)'); rg.addColorStop(1, 'rgba(120,150,255,0)');
        ctx.fillStyle = rg; ctx.beginPath(); ctx.arc(gx, gy, 0.006 * s, 0, 2 * Math.PI); ctx.fill();
      }
      // Injector
      if (p.injection === 'directa') {
        ctx.fillStyle = '#9fb0bf';
        ctx.fillRect(X(0.0028), Y(yHeadTop), 0.0022 * s, (roofY - yHeadTop) * s + 1);
        ctx.fillStyle = sim.inj && sim.inj.left > 0 ? '#d6f1ff' : '#5f6f7d';
        ctx.fillRect(X(0.0030), Y(roofY) - 1, 0.0018 * s, 3);
      }
      // Labels
      ctx.font = '600 11px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.fillStyle = '#c7ced4';
      if (!this.view.zoom || true) ctx.textBaseline = 'alphabetic';
      ctx.fillText('ADM', X(g.inR), Y(yHeadTop) + 13); ctx.fillText('ADM', X(-g.inR), Y(yHeadTop) + 13);
      ctx.fillText('ESC', X(g.exR), Y(yHeadTop) + 13); ctx.fillText('ESC', X(-g.exR), Y(yHeadTop) + 13);
      // Crank inset
      this.drawCrankInset();
    }

    drawValve(r, lift, color, diamFrac) {
      const { ctx, sim, L } = this, g = sim.g, s = L.s;
      const X = rr => this.X(rr), Y = y => this.Y(y);
      const i = Math.min(g.nr - 1, Math.floor(r / g.dr));
      const roof = sim.yTop[sim.roofCell(i) % g.NZ];
      const halfW = Math.min(0.12 * g.R, 0.25 * diamFrac * g.B);
      for (const sg of [-1, 1]) {
        const x0 = X(sg * r);
        ctx.fillStyle = color;
        // stem
        ctx.fillRect(x0 - 0.0011 * s, Y(L.yMin), 0.0022 * s, (roof + lift - L.yMin) * s);
        // head
        ctx.beginPath();
        ctx.moveTo(x0 - halfW * s, Y(roof + lift));
        ctx.lineTo(x0 + halfW * s, Y(roof + lift));
        ctx.lineTo(x0 + 0.0014 * s, Y(roof + lift - 0.0025));
        ctx.lineTo(x0 - 0.0014 * s, Y(roof + lift - 0.0025));
        ctx.closePath(); ctx.fill();
        if (lift > 1e-4) {
          ctx.fillStyle = 'rgba(20, 24, 28, 0.9)';
          ctx.fillRect(x0 - halfW * s, Y(roof) - 1, 2 * halfW * s, Math.max(1, lift * s * 0.15));
        }
      }
    }

    drawCrankInset() {
      const { ctx, sim, L } = this;
      const r = 26, cx = 12 + r, cy = L.H - 12 - r;
      ctx.fillStyle = 'rgba(14,16,18,0.85)'; ctx.beginPath(); ctx.arc(cx, cy, r + 6, 0, 2 * Math.PI); ctx.fill();
      ctx.strokeStyle = '#59616a'; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(cx, cy, r, 0, 2 * Math.PI); ctx.stroke();
      const th = sim.ca * Math.PI / 180;
      ctx.strokeStyle = '#ff7a1a'; ctx.lineWidth = 3; ctx.lineCap = 'round';
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx + Math.sin(th) * r * 0.8, cy - Math.cos(th) * r * 0.8); ctx.stroke();
      ctx.fillStyle = '#c7ced4'; ctx.font = '10px system-ui, sans-serif'; ctx.textAlign = 'center';
      ctx.fillText('PMS', cx, cy - r - 8 < 8 ? 8 : cy - r - 9);
    }

    drawGas() {
      const { sim, L } = this, g = sim.g, N = g.N, NZ = g.NZ;
      if (!this.colL || this.colL.length !== N * 3) { this.colL = new Float32Array(N * 3); this.colR = new Float32Array(N * 3); }
      let qMax = 0;
      for (let c = 0; c < N; c++) if (g.act[c] && sim.flame[c] > qMax) qMax = sim.flame[c];
      const qRef = Math.max(5e9, 0.6 * qMax);
      for (let c = 0; c < N; c++) {
        if (!g.act[c]) continue;
        compositionColor(sim, c, this.colL, c * 3, qRef);
        fieldColor(sim, this.view.field, c, this.colR, c * 3);
        if (sim.mik[c] > 0.7 * sim.m[c] && sim.c[c] < 0.5 && sim.mF[c] > 1e-3 * sim.m[c]) {
          // end gas close to autoignition glows red
          const t = clamp((sim.mik[c] / sim.m[c] - 0.7) / 0.3, 0, 1);
          for (const arr of [this.colL]) { arr[c * 3] = lerp(arr[c * 3], 255, t); arr[c * 3 + 1] = lerp(arr[c * 3 + 1], 40, t); arr[c * 3 + 2] = lerp(arr[c * 3 + 2], 40, t); }
        }
      }
      // Render at half resolution then scale with smoothing.
      const s = L.s, k = 0.5;
      const yTop = sim.yTop[0], yBot = sim.yTop[NZ];
      const wpx = Math.max(2, Math.round(2 * g.R * s * k)), hpx = Math.max(2, Math.round((yBot - yTop) * s * k));
      if (this.off.width !== wpx || this.off.height !== hpx) { this.off.width = wpx; this.off.height = hpx; this.img = this.offCtx.createImageData(wpx, hpx); }
      const img = this.img, data = img.data;
      // Row centres
      const ycs = this.ycs && this.ycs.length === NZ ? this.ycs : (this.ycs = new Float64Array(NZ));
      for (let j = 0; j < NZ; j++) ycs[j] = sim.yTop[j] + 0.5 * sim.dz[j];
      const half = wpx / 2;
      const colI = this.colI && this.colI.length === wpx ? this.colI : (this.colI = new Int32Array(wpx));
      const colF = this.colF && this.colF.length === wpx ? this.colF : (this.colF = new Float32Array(wpx));
      const colC = this.colC && this.colC.length === wpx ? this.colC : (this.colC = new Int32Array(wpx));
      for (let x = 0; x < wpx; x++) {
        const r = Math.abs((x + 0.5 - half) / (s * k));
        const u = r / g.dr - 0.5;
        const i0 = Math.floor(u);
        colI[x] = i0; colF[x] = u - i0; colC[x] = Math.min(g.nr - 1, Math.floor(r / g.dr));
      }
      let jRow = 0;
      for (let y = 0; y < hpx; y++) {
        const yy = yTop + (y + 0.5) / (s * k);
        while (jRow < NZ - 1 && sim.yTop[jRow + 1] <= yy) jRow++;
        let j0 = 0;
        while (j0 < NZ - 1 && ycs[j0 + 1] <= yy) j0++;
        let fy = (yy - ycs[j0]) / ((j0 + 1 < NZ ? ycs[j0 + 1] : ycs[j0] + 1) - ycs[j0]);
        if (yy < ycs[0]) fy = 0;
        fy = clamp(fy, 0, 1);
        for (let x = 0; x < wpx; x++) {
          const o = (y * wpx + x) * 4;
          const cc = colC[x] * NZ + jRow;
          if (!g.act[cc]) { data[o + 3] = 0; continue; }
          const arr = x < half ? this.colL : this.colR;
          let i0 = colI[x], fx = colF[x];
          if (i0 < 0) { i0 = 0; fx = 0; }
          let rr = 0, gg = 0, bb = 0, ww = 0;
          for (let di = 0; di < 2; di++) for (let dj = 0; dj < 2; dj++) {
            const ii = i0 + di, jj = j0 + dj;
            if (ii >= g.nr || jj >= NZ) continue;
            const c = ii * NZ + jj;
            if (!g.act[c]) continue;
            const wgt = (di ? fx : 1 - fx) * (dj ? fy : 1 - fy) + 1e-6;
            rr += wgt * arr[c * 3]; gg += wgt * arr[c * 3 + 1]; bb += wgt * arr[c * 3 + 2]; ww += wgt;
          }
          if (ww <= 0) { const c = cc; rr = arr[c * 3]; gg = arr[c * 3 + 1]; bb = arr[c * 3 + 2]; ww = 1; }
          data[o] = rr / ww; data[o + 1] = gg / ww; data[o + 2] = bb / ww; data[o + 3] = 255;
        }
      }
      this.offCtx.putImageData(img, 0, 0);
      const ctx = this.ctx;
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(this.off, this.X(-g.R), this.Y(yTop), 2 * g.R * s, (yBot - yTop) * s);
      // Divider on the axis
      ctx.strokeStyle = 'rgba(255,255,255,0.18)'; ctx.setLineDash([3, 4]); ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(this.X(0), this.Y(yTop)); ctx.lineTo(this.X(0), this.Y(yBot)); ctx.stroke(); ctx.setLineDash([]);
    }

    drawArrows() {
      const { ctx, sim, L } = this, g = sim.g, s = L.s;
      const step = 17;
      const yTop = sim.yTop[0], yBot = sim.yTop[g.NZ];
      ctx.strokeStyle = 'rgba(255,255,255,0.75)'; ctx.fillStyle = 'rgba(255,255,255,0.75)'; ctx.lineWidth = 1;
      ctx.beginPath();
      for (let py = this.Y(yTop) + step / 2; py < this.Y(yBot); py += step) {
        const y = (py / s) + L.yMin;
        for (let px = L.cx + step / 2; px < this.X(g.R); px += step) {
          const r = (px - L.cx) / s;
          const c = sim.cellAt(r, y);
          if (c < 0) continue;
          const u = sim.u[c], w = sim.w[c], sp = Math.hypot(u, w);
          if (sp < 0.3) continue;
          const len = Math.min(step * 0.9, 2.5 + sp * 0.8), ux = u / sp, uy = w / sp;
          for (const sg of [-1, 1]) {
            const x0 = L.cx + sg * (px - L.cx), dx = sg * ux;
            const x1 = x0 + dx * len, y1 = py + uy * len;
            ctx.moveTo(x0, py); ctx.lineTo(x1, y1);
            ctx.moveTo(x1, y1); ctx.lineTo(x1 - dx * 3 - uy * 2, y1 - uy * 3 + dx * 2);
            ctx.moveTo(x1, y1); ctx.lineTo(x1 - dx * 3 + uy * 2, y1 - uy * 3 - dx * 2);
          }
        }
      }
      ctx.stroke();
    }

    drawLegend() {
      const f = FIELDS[this.view.field];
      const stops = this.view.field === 'swirl' || this.view.field === 'dosado' ? DIVERGE : RAMP;
      const css = stops.map((c, k) => `rgb(${c.map(Math.round).join(',')}) ${Math.round(k / (stops.length - 1) * 100)}%`).join(',');
      this.legend.innerHTML = '';
      const sw = (col, label) => h('span', { class: 'sw' }, h('i', { style: `background:${col}` }), label);
      this.legend.append(
        h('div', { class: 'legend-row' }, h('strong', {}, 'Composición'),
          sw('rgb(46,128,255)', 'mezcla'), sw('rgb(96,60,240)', 'rica'), sw('rgb(38,62,96)', 'aire'),
          sw('rgb(128,132,138)', 'gases quemados'), sw('rgb(255,138,31)', 'llama'), sw('rgb(240,248,255)', 'gotas'), sw('rgb(150,225,255)', 'película')),
        h('div', { class: 'legend-row' }, h('strong', {}, f.label),
          h('span', { class: 'bar-scale' }, h('i', { style: `background:linear-gradient(90deg,${css})` }),
            h('span', {}, `${f.min}`), h('span', {}, `${f.max} ${f.unit}`))));
    }

    // ---------------------------------------------------------- plots

    plotFrame(cv, xr, yr, opts) {
      const ctx = cv.getContext('2d'), dpr = this.L.dpr;
      const W = cv.width / dpr, H = cv.height / dpr;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, W, H);
      const m = { l: 40, r: 8, t: 18, b: 22 };
      const pw = W - m.l - m.r, ph = H - m.t - m.b;
      const tx = opts.logX ? v => m.l + (Math.log(v) - Math.log(xr[0])) / (Math.log(xr[1]) - Math.log(xr[0])) * pw : v => m.l + (v - xr[0]) / (xr[1] - xr[0]) * pw;
      const ty = opts.logY ? v => m.t + ph - (Math.log(Math.max(v, yr[0])) - Math.log(yr[0])) / (Math.log(yr[1]) - Math.log(yr[0])) * ph : v => m.t + ph - (v - yr[0]) / (yr[1] - yr[0]) * ph;
      ctx.strokeStyle = '#2c3237'; ctx.lineWidth = 1;
      ctx.font = '10px system-ui, sans-serif'; ctx.fillStyle = '#8d969c';
      ctx.textAlign = 'center';
      for (const xv of opts.xTicks || []) { const x = tx(xv); ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke(); ctx.fillText(opts.xFmt ? opts.xFmt(xv) : xv, x, H - 7); }
      ctx.textAlign = 'right';
      for (const yv of opts.yTicks || []) { const y = ty(yv); ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(m.l + pw, y); ctx.stroke(); ctx.fillText(opts.yFmt ? opts.yFmt(yv) : yv, m.l - 4, y + 3); }
      ctx.textAlign = 'left'; ctx.fillStyle = '#e6e8e9'; ctx.font = '600 11px system-ui, sans-serif';
      ctx.fillText(opts.title, m.l, 12);
      return { ctx, tx, ty, m, pw, ph, W, H };
    }

    line(fr, xs, ys, color, width, from, dash) {
      const { ctx, tx, ty } = fr;
      ctx.strokeStyle = color; ctx.lineWidth = width || 1.5; ctx.setLineDash(dash || []);
      ctx.beginPath();
      let pen = false, prevX = null;
      for (let k = from || 0; k < xs.length; k++) {
        const x = tx(xs[k]), y = ty(ys[k]);
        if (!isFinite(x) || !isFinite(y)) { pen = false; continue; }
        if (prevX !== null && xs[k] < prevX && !fr.allowBack) pen = false;
        if (pen) ctx.lineTo(x, y); else ctx.moveTo(x, y);
        pen = true; prevX = xs[k];
      }
      ctx.stroke(); ctx.setLineDash([]);
    }

    drawPlots() {
      const sim = this.sim, tr = sim.trace, prev = sim.prevCycle ? sim.prevCycle.trace : null, p = sim.p;
      const bar = v => v / 1e5;
      // Pressure against crank angle
      let pmax = 10;
      for (const t of [tr, prev]) if (t) for (const v of t.p) pmax = Math.max(pmax, bar(v));
      pmax = Math.ceil(pmax / 10) * 10;
      {
        const fr = this.plotFrame(this.plots.p, [-360, 360], [0, pmax], {
          title: 'Presión en el cilindro (bar)', xTicks: [-360, -180, 0, 180, 360], yTicks: [0, pmax / 2, pmax],
          xFmt: v => `${v}°`, yFmt: v => Math.round(v),
        });
        // valve events and spark
        const { ctx, tx, m, ph } = fr;
        const band = (a0, a1, col) => { ctx.fillStyle = col; ctx.fillRect(tx(a0), m.t, tx(a1) - tx(a0), ph); };
        band(-360, -180 + p.ivc, 'rgba(108,143,179,0.10)'); band(360 - p.ivo, 360, 'rgba(108,143,179,0.10)');
        band(180 - p.evo, 360, 'rgba(179,122,108,0.10)'); band(-360, -360 + p.evc, 'rgba(179,122,108,0.10)');
        ctx.strokeStyle = 'rgba(160,190,255,0.8)'; ctx.setLineDash([2, 3]);
        ctx.beginPath(); ctx.moveTo(tx(-p.spark), m.t); ctx.lineTo(tx(-p.spark), m.t + ph); ctx.stroke(); ctx.setLineDash([]);
        if (prev) this.line(fr, prev.ca, prev.p.map(bar), 'rgba(255,122,26,0.3)', 1.2);
        this.line(fr, tr.ca, tr.p.map(bar), '#ff7a1a', 1.8);
        const x = tx(sim.cycleAngle());
        ctx.strokeStyle = 'rgba(255,255,255,0.4)'; ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke();
      }
      // p-V, log-log
      {
        const g = sim.g, vmin = g.Vc * 0.9, vmax = (g.Vc + g.Vd) * 1.1;
        const fr = this.plotFrame(this.plots.pv, [vmin, vmax], [0.2, Math.max(pmax, 20)], {
          title: 'Diagrama p–V (log)', logX: true, logY: true,
          xTicks: [g.Vc, g.Vc + g.Vd], xFmt: v => `${(v * 1e6).toFixed(0)} cm³`, yTicks: [0.5, 1, 5, 10, 50].filter(v => v < Math.max(pmax, 20)), yFmt: v => v,
        });
        fr.allowBack = true;
        if (prev) this.line(fr, prev.V, prev.p.map(bar), 'rgba(63,180,255,0.3)', 1.2);
        this.line(fr, tr.V, tr.p.map(bar), '#3fb4ff', 1.8);
        const { ctx, tx, ty } = fr;
        ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(tx(sim.Vtot), ty(bar(sim.p0)), 3, 0, 2 * Math.PI); ctx.fill();
      }
      // Heat release and burned fraction around TDC
      {
        let hmax = 10;
        for (const t of [tr, prev]) if (t) for (let k = 0; k < t.ca.length; k++) if (t.ca[k] > -60 && t.ca[k] < 120) hmax = Math.max(hmax, t.hrr[k]);
        hmax = Math.ceil(hmax / 10) * 10;
        const fr = this.plotFrame(this.plots.hrr, [-60, 120], [0, hmax], {
          title: 'Calor liberado (J/°) y fracción quemada', xTicks: [-60, -30, 0, 30, 60, 90, 120], yTicks: [0, hmax / 2, hmax],
          xFmt: v => `${v}°`, yFmt: v => Math.round(v),
        });
        if (prev) this.line(fr, prev.ca, prev.hrr, 'rgba(255,138,31,0.3)', 1.2);
        this.line(fr, tr.ca, tr.hrr, '#ff8a1f', 1.8);
        const ty2 = v => fr.m.t + fr.ph - v * fr.ph;
        const fr2 = Object.assign({}, fr, { ty: ty2 });
        if (prev) this.line(fr2, prev.ca, prev.mfb, 'rgba(230,232,233,0.3)', 1.2, 0, [3, 3]);
        this.line(fr2, tr.ca, tr.mfb, '#e6e8e9', 1.4, 0, [3, 3]);
      }
      // Fuel state
      {
        const t = tr;
        let fmax = 1e-9;
        for (const T of [tr, prev]) if (T) for (let k = 0; k < T.ca.length; k++) fmax = Math.max(fmax, T.vap[k] + T.liq[k] + T.film[k] + T.burned[k]);
        const mg = v => v * 1e6;
        const top = Math.max(1, Math.ceil(mg(fmax) / 5) * 5);
        const fr = this.plotFrame(this.plots.fuel, [-360, 360], [0, top], {
          title: 'Combustible (mg): vapor, líquido, película, quemado', xTicks: [-360, -180, 0, 180, 360], yTicks: [0, top / 2, top],
          xFmt: v => `${v}°`, yFmt: v => +v.toFixed(1),
        });
        this.line(fr, t.ca, t.vap.map(mg), '#2e80ff', 1.6);
        this.line(fr, t.ca, t.liq.map(mg), '#f0f8ff', 1.4);
        this.line(fr, t.ca, t.film.map(mg), '#96e1ff', 1.4, 0, [4, 2]);
        this.line(fr, t.ca, t.burned.map(mg), '#ff8a1f', 1.4);
      }
    }

    // ---------------------------------------------------------- readouts

    drawReadouts() {
      const sim = this.sim, d = sim.diagnostics();
      this.frameCount = (this.frameCount || 0) + 1;
      if (this.frameCount % 6) return;
      const f1 = v => v.toFixed(1), mg = v => `${(v * 1e6).toFixed(2)} mg`;
      const lam = isFinite(d.lambda) ? d.lambda.toFixed(2) : '—';
      const rows = [
        ['Presión', `${(d.p / 1e5).toFixed(2)} bar`],
        ['Temperatura media', `${Math.round(d.T)} K`],
        ['Gas fresco / quemado', `${Math.round(d.Tu)} / ${d.Tb ? Math.round(d.Tb) + ' K' : '—'}`],
        ['Masa en el cilindro', mg(d.M)],
        ['λ en el cilindro', lam],
        ['Gases quemados', `${(d.residual * 100).toFixed(1)} % de la masa`],
        ['Fracción quemada', `${(d.mfb * 100).toFixed(0)} %`],
        ["Turbulencia u'", `${f1(d.u)} m/s`],
        ['Velocidad de squish', `${f1(d.usq)} m/s`],
        ['Relación de swirl', d.swirl.toFixed(2)],
        ['Vapor / líquido / película', `${(d.vapor * 1e6).toFixed(1)} / ${(d.liquid * 1e6).toFixed(1)} / ${(d.film * 1e6).toFixed(1)} mg`],
        ['Alzada adm. / esc.', `${(d.liftIn * 1000).toFixed(1)} / ${(d.liftEx * 1000).toFixed(1)} mm`],
      ];
      this.live.innerHTML = '';
      for (const [k, v] of rows) this.live.append(h('dt', {}, k), h('dd', {}, v));
      const pc = sim.prevCycle;
      this.cycle.innerHTML = '';
      if (!pc) { this.cycle.append(h('p', { class: 'hint' }, 'Los resultados aparecen al terminar el primer ciclo.')); return; }
      const st = pc.stats;
      const ang = v => (v == null ? '—' : `${v >= 0 ? v.toFixed(1) + '° DPMS' : (-v).toFixed(1) + '° APMS'}`);
      const caPos = v => { const x = ((v % 720) + 720) % 720; return x >= 360 ? x - 720 : x; };
      const knock = st.knockCa ? `Sí, a ${ang(caPos(st.knockCa))} · ${(st.knockFreq / 1000).toFixed(1)} kHz · ${(st.knockFrac * 1e6).toFixed(2)} mg` : 'No';
      const res = [
        ['PMI neta / bruta', `${(st.imep / 1e5).toFixed(2)} / ${(st.imepGross / 1e5).toFixed(2)} bar`],
        ['Presión máxima', `${(st.pMax / 1e5).toFixed(1)} bar a ${ang(caPos(st.caPMax))}`],
        ['Quemado 10 / 50 / 90 %', [st.mfb10, st.mfb50, st.mfb90].map(v => (v == null ? '—' : `${v.toFixed(1)}°`)).join(' / ')],
        ['Duración 10–90 %', st.mfb10 != null && st.mfb90 != null ? `${(st.mfb90 - st.mfb10).toFixed(1)}°` : '—'],
        ['Rendimiento indicado', `${(st.efficiency * 100).toFixed(1)} %`],
        ['Eficiencia de combustión', `${(st.combustionEff * 100).toFixed(1)} %`],
        ['Calor a las paredes', st.qComb > 0 ? `${(st.qWall / st.qComb * 100).toFixed(1)} % del liberado` : '—'],
        ['Llenado (η volumétrico)', `${(st.volEff * 100).toFixed(1)} %`],
        ['λ atrapado', isFinite(st.lambdaTrapped) ? st.lambdaTrapped.toFixed(2) : '—'],
        ['Combustible sin quemar al abrir escape', `${(st.hcEvo * 100).toFixed(1)} %`],
        ["u' en la chispa", `${st.uSpark.toFixed(1)} m/s`],
        ['Squish máximo', `${st.sqMax.toFixed(1)} m/s`],
        ['Película máxima', `${(st.filmMax * 1e6).toFixed(2)} mg`],
        ['Golpe (knock)', knock],
        ['Encendido', st.misfire ? 'Falló: mezcla fuera de rango en la bujía' : 'Normal'],
      ];
      for (const [k, v] of res) this.cycle.append(h('dt', {}, k), h('dd', { class: k === 'Golpe (knock)' && st.knockCa ? 'warn' : k === 'Encendido' && st.misfire ? 'warn' : '' }, v));
    }
  }

  return { App, PRESETS, GROUPS, FIELDS };
})();

if (typeof window !== 'undefined') window.addEventListener('DOMContentLoaded', () => { window.pistonLab = new PSUI.App(); });
