// Small perceptual colour maps as [r, g, b] in 0..1. Inferno for physical
// scalars (heat rate, dynamic pressure, density), viridis for ensembles.
const INFERNO = [
  [0.001, 0.000, 0.014], [0.088, 0.036, 0.220], [0.258, 0.039, 0.406], [0.416, 0.090, 0.433],
  [0.578, 0.148, 0.404], [0.735, 0.215, 0.330], [0.866, 0.316, 0.226], [0.954, 0.462, 0.100],
  [0.988, 0.645, 0.040], [0.974, 0.835, 0.223], [0.988, 0.998, 0.645],
];
const VIRIDIS = [
  [0.267, 0.005, 0.329], [0.283, 0.141, 0.458], [0.254, 0.265, 0.530], [0.207, 0.372, 0.553],
  [0.164, 0.471, 0.558], [0.128, 0.567, 0.551], [0.135, 0.659, 0.518], [0.267, 0.749, 0.441],
  [0.478, 0.821, 0.318], [0.741, 0.873, 0.150], [0.993, 0.906, 0.144],
];

function sample(lut, t) {
  const x = Math.min(1, Math.max(0, Number.isFinite(t) ? t : 0)) * (lut.length - 1);
  const i = Math.min(lut.length - 2, Math.floor(x)), f = x - i;
  const a = lut[i], b = lut[i + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}

export function inferno(t) { return sample(INFERNO, t); }
export function viridis(t) { return sample(VIRIDIS, t); }
export const INFERNO_CSS = 'linear-gradient(to right, rgb(0,0,4), rgb(66,10,104), rgb(147,38,103), rgb(221,81,58), rgb(252,165,10), rgb(252,255,164))';
export const VIRIDIS_CSS = 'linear-gradient(to right, rgb(68,1,84), rgb(59,82,139), rgb(33,145,140), rgb(94,201,98), rgb(253,231,37))';
