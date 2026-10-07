// Plain-language descriptions for the explorer. Every number comes from the
// calendar state or from the solver's Earth control; the words only read them.

export const REGIMES = [
  {
    key: "day",
    name: "Daytime",
    colour: "#ffc53d",
    short: "Day",
    blurb: "The Sun is above the horizon.",
  },
  {
    key: "twilight",
    name: "Golden twilight",
    colour: "#ff8a5c",
    short: "Twilight",
    blurb:
      "The Sun has set, and sunlight scattered high in the deep air still lights the ground well enough to work outside without lamps.",
  },
  {
    key: "earthlit",
    name: "Earthlit night",
    colour: "#38c6f4",
    short: "Earthlit night",
    blurb: "Earth outshines what remains of the twilight.",
  },
  {
    key: "night",
    name: "Night",
    colour: "#7b6cf6",
    short: "Night",
    blurb:
      "A dim glow of twilight remains, below the level at which Earth's civil twilight ends.",
  },
  {
    key: "dark",
    name: "Deep night",
    colour: "#3a3f7a",
    short: "Deep night",
    blurb: "Neither the Sun nor Earth lights the sky. The stars are out.",
  },
];

export const EVENT_INFO = {
  sunrise: { name: "Sunrise", icon: "sunrise", colour: "#ffc53d" },
  sunset: { name: "Sunset", icon: "sunset", colour: "#ff9f43" },
  dusk: { name: "Twilight ends", icon: "dusk", colour: "#b07cff" },
  dawn: { name: "Twilight begins", icon: "dawn", colour: "#ff8a5c" },
  earthrise: { name: "Earthrise", icon: "earthrise", colour: "#38c6f4" },
  earthset: { name: "Earthset", icon: "earthset", colour: "#2a8fc4" },
  "full-earth": { name: "Full Earth", icon: "full", colour: "#5ad8ff" },
  "new-earth": { name: "New Earth", icon: "new", colour: "#6b84a8" },
  "eclipse-begins": { name: "Earth eclipses the Sun", icon: "eclipse", colour: "#ff5d8f" },
  "eclipse-ends": { name: "Eclipse ends", icon: "eclipse", colour: "#ff5d8f" },
};

export const PRACTICAL_DUSK_LUX = 2.98;
// Maximum full-Moon illuminance on Earth, about 0.3 lux (Kyba, Mohar & Posch 2017).
export const FULL_MOON_LUX = 0.3;

const DIRS = [
  "north",
  "north-northeast",
  "northeast",
  "east-northeast",
  "east",
  "east-southeast",
  "southeast",
  "south-southeast",
  "south",
  "south-southwest",
  "southwest",
  "west-southwest",
  "west",
  "west-northwest",
  "northwest",
  "north-northwest",
];
export const direction = (az) => DIRS[Math.round(az / 22.5) % 16];
export const compass = (az) =>
  ["N", "NE", "E", "SE", "S", "SW", "W", "NW"][Math.round(az / 45) % 8];

const fmt = (v, d = 0) =>
  v.toLocaleString("en-US", { maximumFractionDigits: d, minimumFractionDigits: 0 });

export function lux(v) {
  if (!(v > 0)) return "0";
  if (v < 0.001) return v.toExponential(1).replace("e", "×10^");
  if (v < 0.1) return fmt(v, 3);
  if (v < 10) return fmt(v, 2);
  if (v < 100) return fmt(v, 1);
  return fmt(Math.round(v / (v < 10000 ? 1 : 100)) * (v < 10000 ? 1 : 100));
}

export function duration(ms, parts = 2) {
  const a = Math.abs(ms);
  const units = [
    ["day", 86400000],
    ["hour", 3600000],
    ["minute", 60000],
  ];
  const out = [];
  let rest = a;
  for (const [name, size] of units) {
    const n = Math.floor(rest / size);
    if (n > 0 || (out.length && out.length < parts)) {
      if (n > 0) out.push(`${n} ${name}${n === 1 ? "" : "s"}`);
      rest -= n * size;
    }
    if (out.length >= parts) break;
  }
  return out.length ? out.join(" ") : "less than a minute";
}
export function shortDuration(ms) {
  const a = Math.abs(ms),
    d = Math.floor(a / 86400000),
    h = Math.floor((a % 86400000) / 3600000),
    m = Math.round((a % 3600000) / 60000);
  if (d > 0) return `${d}d ${h}h`;
  if (h > 0) return `${h}h ${m}m`;
  return `${m}m`;
}
export const relative = (ms) =>
  ms >= 0 ? `in ${duration(ms)}` : `${duration(ms)} ago`;

// Earth's phase as seen from the Moon.
export function earthPhase(lit, waxing) {
  if (lit >= 0.97) return "Full Earth";
  if (lit <= 0.03) return "New Earth";
  if (Math.abs(lit - 0.5) < 0.06) return waxing ? "First-quarter Earth" : "Last-quarter Earth";
  if (lit > 0.5) return waxing ? "Waxing gibbous Earth" : "Waning gibbous Earth";
  return waxing ? "Waxing crescent Earth" : "Waning crescent Earth";
}

// Interpolate the Sun elevation on Earth that gives the same horizontal light,
// from the solver's Earth control (same atmosphere code, Earth's own column).
export function earthSunFor(luxValue, earthHeader) {
  const s = earthHeader.suns,
    L = earthHeader.horizontal_lux;
  if (!(luxValue > 0)) return null;
  if (luxValue > L[L.length - 1]) return { above: true, elevation: 90 };
  for (let i = s.length - 2; i >= 0; i--) {
    if (L[i] <= luxValue && L[i] > 0) {
      const a = Math.log(Math.max(L[i], 1e-30)),
        b = Math.log(L[i + 1]);
      const t = (Math.log(luxValue) - a) / (b - a);
      return { elevation: s[i] + t * (s[i + 1] - s[i]) };
    }
  }
  return null;
}

export function earthComparison(luxValue, earthHeader) {
  if (!(luxValue > 0)) return "No measurable light";
  if (luxValue < 0.4) {
    const ratio = luxValue / FULL_MOON_LUX;
    if (ratio >= 0.6) return "About as bright as a full-Moon night on Earth";
    if (ratio >= 0.02)
      return `About ${fmt(ratio * 100, 0)}% of the light of a full Moon on Earth`;
    return "Darker than any moonlit night on Earth";
  }
  const match = earthSunFor(luxValue, earthHeader);
  if (!match) return "";
  const e = match.elevation;
  if (e >= 60) return "As bright as a sunny midday on Earth";
  if (e >= 25) return `Like a sunny day on Earth with the Sun ${fmt(e, 0)}° high`;
  if (e >= 5) return `Like Earth with the Sun ${fmt(e, 0)}° high, morning or late afternoon`;
  if (e >= 0) return `Like Earth with the Sun just ${fmt(e, 1)}° above the horizon`;
  if (e >= -6) return `Like Earth's civil twilight, with the Sun ${fmt(-e, 1)}° below the horizon`;
  return `Like late twilight on Earth, with the Sun ${fmt(-e, 1)}° down`;
}

export function moonlightRatio(luxValue) {
  const r = luxValue / FULL_MOON_LUX;
  if (r >= 1.5) return `${fmt(r, r < 10 ? 1 : 0)}× the light of a full Moon on Earth`;
  return null;
}

export function altitudeWords(el, az) {
  const a = Math.abs(el);
  if (el >= 85) return "almost straight overhead";
  if (el >= 0) return `${fmt(a, 0)}° up in the ${direction(az)}`;
  return `${fmt(a, 0)}° below the ${direction(az)} horizon`;
}

export const fmtNumber = fmt;
