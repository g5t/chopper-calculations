/* The BIFROST source-intensity figure, drawn from pre-computed frames.
 *
 * scripts/bifrost_intensity_figure.py writes the source image and one frame per setting
 * into docs/assets/bifrost-intensity/; this draws whichever frame the controls select, as
 * SVG, so the page needs nothing but the browser. The static picture the page carries is
 * replaced once this has drawn, and stays for a reader without JavaScript.
 */
(function () {
  'use strict';
  const NS = 'http://www.w3.org/2000/svg';

  function el(tag, attrs, parent, text) {
    const node = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs || {})) node.setAttribute(k, v);
    if (text !== undefined) node.textContent = text;
    if (parent) parent.appendChild(node);
    return node;
  }

  function niceTicks(lo, hi, count) {
    const raw = (hi - lo) / Math.max(count, 1);
    const power = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * power).find(s => s >= raw) || 10 * power;
    const ticks = [];
    for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9 * step; t += step) {
      ticks.push(Math.abs(t) < 1e-12 * step ? 0 : t);
    }
    return ticks;
  }

  function label(value) {
    const a = Math.abs(value);
    if (a !== 0 && (a < 1e-3 || a >= 1e5)) return value.toExponential(1);
    return String(+value.toPrecision(6));
  }

  /* One set of axes in `box`, with ticks along the bottom and the left. */
  function axes(svg, box, xlim, ylim, xlabel, ylabel) {
    const g = el('g', {}, svg);
    const x = v => box.x + (v - xlim[0]) / (xlim[1] - xlim[0]) * box.w;
    const y = v => box.y + box.h - (v - ylim[0]) / (ylim[1] - ylim[0]) * box.h;
    const clip = 'clip' + Math.random().toString(36).slice(2);
    el('rect', {x: box.x, y: box.y, width: box.w, height: box.h},
       el('clipPath', {id: clip}, el('defs', {}, svg)));
    const plot = el('g', {'clip-path': `url(#${clip})`}, svg);
    for (const t of niceTicks(xlim[0], xlim[1], Math.max(3, Math.round(box.w / 80)))) {
      el('line', {x1: x(t), x2: x(t), y1: box.y + box.h, y2: box.y + box.h + 5,
                  stroke: 'currentColor'}, g);
      el('text', {x: x(t), y: box.y + box.h + 22, 'text-anchor': 'middle',
                  'font-size': 16, fill: 'currentColor'}, g, label(t));
    }
    for (const t of niceTicks(ylim[0], ylim[1], Math.max(3, Math.round(box.h / 45)))) {
      el('line', {x1: box.x - 5, x2: box.x, y1: y(t), y2: y(t), stroke: 'currentColor'}, g);
      el('text', {x: box.x - 8, y: y(t) + 4, 'text-anchor': 'end', 'font-size': 16,
                  fill: 'currentColor'}, g, label(t));
    }
    el('text', {x: box.x + box.w / 2, y: box.y + box.h + 46, 'text-anchor': 'middle',
                'font-size': 18, fill: 'currentColor'}, g, xlabel);
    el('text', {x: box.x - 64, y: box.y + box.h / 2, 'text-anchor': 'middle', 'font-size': 18,
                fill: 'currentColor',
                transform: `rotate(-90 ${box.x - 64} ${box.y + box.h / 2})`}, g, ylabel);
    const frame = () => el('rect', {x: box.x, y: box.y, width: box.w, height: box.h,
                                    fill: 'none', stroke: 'currentColor'}, svg);
    return {x, y, plot, g, frame};
  }

  /* Energy ticks at their wavelengths, along one side of a wavelength axis. */
  function energyAxis(ax, box, ticks, lim, side) {
    let last = null;   /* skip a tick that would crowd the one before it */
    for (const [wavelength, text] of ticks) {
      if (wavelength < lim[0] || wavelength > lim[1]) continue;
      const at = side === 'right' ? ax.y(wavelength) : ax.x(wavelength);
      if (last !== null && Math.abs(at - last) < (side === 'right' ? 20 : 34)) continue;
      last = at;
      if (side === 'right') {
        const yy = ax.y(wavelength);
        el('line', {x1: box.x + box.w, x2: box.x + box.w + 5, y1: yy, y2: yy,
                    stroke: 'currentColor'}, ax.g);
        el('text', {x: box.x + box.w + 8, y: yy + 4, 'font-size': 16, fill: 'currentColor'},
           ax.g, text);
      } else {
        const xx = ax.x(wavelength);
        el('line', {x1: xx, x2: xx, y1: box.y - 5, y2: box.y, stroke: 'currentColor'}, ax.g);
        el('text', {x: xx, y: box.y - 9, 'text-anchor': 'middle', 'font-size': 16,
                    fill: 'currentColor'}, ax.g, text);
      }
    }
    if (side === 'right') {
      const cx = box.x + box.w + 50, cy = box.y + box.h / 2;
      el('text', {x: cx, y: cy, 'text-anchor': 'middle', 'font-size': 18, fill: 'currentColor',
                  transform: `rotate(90 ${cx} ${cy})`}, ax.g, 'energy [meV]');
    } else {
      el('text', {x: box.x + box.w / 2, y: box.y - 32, 'text-anchor': 'middle',
                  'font-size': 18, fill: 'currentColor'}, ax.g, 'energy [meV]');
    }
  }

  function stairs(ax, start, step, values, style) {
    let d = `M${ax.x(start)},${ax.y(0)}`;
    values.forEach((v, i) => {
      d += `L${ax.x(start + i * step)},${ax.y(v)}L${ax.x(start + (i + 1) * step)},${ax.y(v)}`;
    });
    d += `L${ax.x(start + values.length * step)},${ax.y(0)}`;
    el('path', Object.assign({d, fill: 'none'}, style), ax.plot);
  }

  function legend(ax, box, entries) {
    entries.forEach(([text, style], i) => {
      const yy = box.y + 18 + 20 * i, x0 = box.x + box.w - 225;
      el('line', Object.assign({x1: x0, x2: x0 + 22, y1: yy, y2: yy}, style), ax.g);
      el('text', {x: x0 + 28, y: yy + 4, 'font-size': 15, fill: 'currentColor'}, ax.g, text);
    });
  }

  const SUPPORT = {stroke: '#999', 'stroke-width': 1.3, 'stroke-dasharray': '5 3'};

  function draw(svg, data, base, key) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    const f = data.frames[key];
    const left = {x: 90, y: 40, w: 420, h: 440};
    const ax = axes(svg, left, data.time, data.wavelength,
                    'emission time at the source [ms]', 'wavelength [Å]');
    el('image', {href: base + 'source.png', x: left.x, y: left.y, width: left.w,
                 height: left.h, preserveAspectRatio: 'none'}, ax.plot);
    for (const polygon of f.outline) {
      const points = polygon.map(([t, w]) => `${ax.x(t)},${ax.y(w)}`).join(' ');
      el('polygon', {points, fill: 'none', stroke: 'white', 'stroke-width': 3}, ax.plot);
      el('polygon', {points, fill: 'none', stroke: 'crimson', 'stroke-width': 1.5}, ax.plot);
    }
    ax.frame();
    energyAxis(ax, left, data.energy_ticks, data.wavelength, 'right');
    el('text', {x: left.x + left.w / 2, y: left.y - 12, 'text-anchor': 'middle',
                'font-size': 19, fill: 'currentColor'}, svg,
       'source, and the region the train transmits');

    if (!f.arrival) return;
    const panels = [
      [{x: 690, y: 70, w: 390, h: 160}, f.arrival, 'arrival time at the sample [ms]',
       'intensity [per ms]', '#1f77b4'],
      [{x: 690, y: 345, w: 390, h: 135}, f.wavelength, 'wavelength [Å]',
       'intensity [per Å]', '#ff7f0e'],
    ];
    panels.forEach(([box, d, xlabel, ylabel, colour], i) => {
      const xlim = [d.start, d.start + d.step * d.support.length];
      const top = Math.max(...d.support) * 1.08 || 1;
      const a = axes(svg, box, xlim, [0, top], xlabel, ylabel);
      stairs(a, d.start, d.step, d.support, SUPPORT);
      stairs(a, d.start, d.step, d.weighted, {stroke: colour, 'stroke-width': 2});
      a.frame();
      legend(a, box, [['whole region open', SUPPORT],
                      ['weighted by transmission', {stroke: colour, 'stroke-width': 2}]]);
      if (i === 1) energyAxis(a, box, data.energy_ticks, xlim, 'top');
    });
    el('text', {x: 885, y: 40, 'text-anchor': 'middle', 'font-size': 19, fill: 'currentColor'},
       svg, `${(100 * f.fraction).toFixed(2)}% of the source rectangle reaches the sample`);
  }

  async function init() {
    const root = document.getElementById('bifrost-intensity');
    if (!root || root.dataset.ready) return;
    root.dataset.ready = '1';
    const fallback = document.querySelector('img[src$="bifrost-intensity/static.png"]');
    if (!fallback) return;
    const base = fallback.src.replace(/static\.png$/, '');
    let data;
    try {
      data = await (await fetch(base + 'frames.json')).json();
    } catch (error) {
      return;   /* leave the static picture */
    }

    const controls = document.createElement('div');
    controls.style.cssText = 'display:flex;flex-wrap:wrap;gap:0.4em 1.6em;margin:0.5em 0;';
    const state = {};
    const fmt = {energies: v => `${v} meV`, shaping_times: v => `${v} ms`,
                 pulse_delays: v => `${v} µs`};
    const names = {energies: 'lowest energy', shaping_times: 'pulse-shaping time',
                   pulse_delays: 'pulse delay'};
    const index = (list, value) => Math.max(0, list.indexOf(value));
    const svg = el('svg', {viewBox: '0 0 1110 570', width: '100%', role: 'img',
                           'aria-label': 'BIFROST source intensity through its chopper train'});
    const render = () => {
      const key = [data.energies[state.energies], data.shaping_times[state.shaping_times],
                   data.pulse_delays[state.pulse_delays], data.psc_speeds[state.psc_speeds]]
        .map(v => String(+v)).join('|');
      draw(svg, data, base, key);
    };
    ['energies', 'shaping_times', 'pulse_delays'].forEach((k, i) => {
      const list = data[k];
      state[k] = index(list, data.default[i]);
      const wrap = document.createElement('label');
      wrap.style.cssText = 'display:flex;align-items:center;gap:0.5em;';
      const value = document.createElement('span');
      value.style.minWidth = '5em';
      const input = document.createElement('input');
      Object.assign(input, {type: 'range', min: 0, max: list.length - 1, step: 1,
                            value: state[k]});
      const show = () => { value.textContent = fmt[k](list[state[k]]); };
      input.addEventListener('input', () => { state[k] = +input.value; show(); render(); });
      show();
      wrap.append(names[k], input, value);
      controls.appendChild(wrap);
    });
    state.psc_speeds = index(data.psc_speeds, data.default[3]);
    const psc = document.createElement('fieldset');
    psc.style.cssText = 'border:none;padding:0;margin:0;display:flex;align-items:center;gap:0.6em;';
    psc.append('pulse-shaping disks');
    /* Buttons rather than radio inputs: the site's theme hides native radios. */
    const buttons = data.psc_speeds.map((speed, i) => {
      const button = document.createElement('button');
      button.type = 'button';
      button.textContent = `${speed} Hz ${i === 0 ? '(designed)' : '(substitute)'}`;
      button.style.cssText = 'padding:0.15em 0.7em;border:1px solid currentColor;'
        + 'border-radius:0.3em;cursor:pointer;background:transparent;color:inherit;';
      button.addEventListener('click', () => { state.psc_speeds = i; mark(); render(); });
      psc.appendChild(button);
      return button;
    });
    const mark = () => buttons.forEach((button, i) => {
      const on = i === state.psc_speeds;
      button.setAttribute('aria-pressed', String(on));
      button.style.fontWeight = on ? 'bold' : 'normal';
      button.style.opacity = on ? '1' : '0.6';
    });
    mark();
    controls.appendChild(psc);
    root.append(controls, svg);
    render();
    fallback.style.display = 'none';
  }

  if (typeof document$ !== 'undefined') {
    document$.subscribe(init);
  } else if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
