// Everyday words for the public sky page. Every number comes from the calendar
// state, its events or the solver's Earth control sky; this module only reads them.

export const PRACTICAL_DUSK_LUX = 2.98;
// The brightest full Moon on Earth gives about 0.3 lux (Kyba, Mohar & Posch 2017).
export const FULL_MOON_LUX = 0.3;

export const CONDITIONS = [
  { key: "day", name: "Daytime", colour: "#ffd166" },
  { key: "twilight", name: "Golden twilight", colour: "#ff9f5a" },
  { key: "earthlit", name: "Earthlit night", colour: "#6fd3ff" },
  { key: "night", name: "Night", colour: "#a99bff" },
  { key: "dark", name: "Dark night", colour: "#8790c8" },
];

export function conditionOf(s) {
  if (s.sun_elevation_deg + s.sun_radius_deg > 0) return 0;
  if (s.earth > s.solar && s.total >= 0.05) return 2;
  if (s.total >= PRACTICAL_DUSK_LUX) return 1;
  if (s.total >= 0.05) return 3;
  return 4;
}

export const MOMENTS = {
  sunrise: "Sunrise",
  sunset: "Sunset",
  dusk: "Twilight fades",
  dawn: "First light",
  earthrise: "Earthrise",
  earthset: "Earthset",
  "full-earth": "Full Earth",
  "new-earth": "New Earth",
  "eclipse-begins": "Earth eclipses the Sun",
};

const DIRS = ["north", "northeast", "east", "southeast", "south", "southwest", "west", "northwest"];
export const direction = (az) => DIRS[Math.round(az / 45) % 8];

export function placeInSky(el, az) {
  if (el >= 70) return "high overhead";
  if (el >= 35) return `high in the ${direction(az)}`;
  if (el >= 10) return `in the ${direction(az)}`;
  return `low in the ${direction(az)}`;
}

export function earthPhaseWords(lit) {
  if (lit >= 0.97) return "full";
  if (lit >= 0.8) return "nearly full";
  if (lit >= 0.6) return "more than half lit";
  if (lit >= 0.4) return "half lit";
  if (lit >= 0.15) return "a crescent";
  if (lit >= 0.03) return "a thin crescent";
  return "new and dark";
}

export function fromNow(ms) {
  const a = Math.abs(ms),
    hours = a / 3600000,
    days = hours / 24;
  let text;
  if (days >= 1.5) text = `${Math.round(days)} days`;
  else if (hours >= 20) text = "a day";
  else if (hours >= 1.5) text = `${Math.round(hours)} hours`;
  else if (hours >= 0.75) text = "an hour";
  else {
    const minutes = Math.max(1, Math.round(hours * 60));
    text = `${minutes} minute${minutes === 1 ? "" : "s"}`;
  }
  return ms >= 0 ? `in ${text}` : `${text} ago`;
}
const away = (ms) => fromNow(ms).replace(/^in /, "");

// The Sun height on Earth that gives the same light on the ground, from the
// solver's Earth control sky.
function earthSunFor(lux, control) {
  const s = control.suns,
    L = control.horizontal_lux;
  if (lux >= L[L.length - 1]) return 90;
  for (let i = s.length - 2; i >= 0; i--)
    if (L[i] > 0 && L[i] <= lux) {
      const t = (Math.log(lux) - Math.log(L[i])) / (Math.log(L[i + 1]) - Math.log(L[i]));
      return s[i] + t * (s[i + 1] - s[i]);
    }
  return null;
}

export function brightness(lux, condition, control) {
  if (!(lux > 0)) return "As dark as a starlit night";
  if (condition === 2 || lux < 0.4) {
    const r = lux / FULL_MOON_LUX;
    if (r >= 1.5) return `About ${Math.round(r)} times brighter than a full Moon on Earth`;
    if (r >= 0.6) return "As bright as a full-Moon night on Earth";
    return "Darker than a moonlit night on Earth";
  }
  const e = earthSunFor(lux, control);
  if (e === null) return "";
  if (e >= 30) return "As bright as a sunny day on Earth";
  if (e >= 10) return "As bright as a sunny morning on Earth";
  if (e >= 0) return "As bright as Earth just after sunrise";
  if (e >= -6) return "As bright as Earth just after sunset";
  return "As bright as late dusk on Earth";
}

// One sentence about the Sun's part of the story, from the state and the
// nearest events; Earth and the brightness have their own lines.
export function sentence({ state: s, condition, rising, lastSunset, lastSunrise, nextSunrise, nextSunset, now, reach }) {
  const rise = nextSunrise ? ` Sunrise is ${away(nextSunrise - now)} away.` : "";
  const set = lastSunset ? `The Sun set ${fromNow(lastSunset - now)}` : "The Sun is below the horizon";
  if (condition === 0) {
    const where = placeInSky(s.sun_elevation_deg, s.sun_azimuth_deg);
    if (lastSunrise && now - lastSunrise < 24 * 3600000)
      return `The Sun rose ${fromNow(lastSunrise - now)} and stays up for about two weeks.`;
    if (nextSunset && nextSunset - now < 24 * 3600000) return `The Sun is ${where} and sets ${fromNow(nextSunset - now)}.`;
    return `The Sun is ${where}${nextSunset ? ` and sets ${fromNow(nextSunset - now)}` : ""}. A day here lasts about two weeks.`;
  }
  if (condition === 1) {
    if (rising && nextSunrise) return `Morning twilight. The sky glows brighter by the hour, and the Sun rises ${fromNow(nextSunrise - now)}.`;
    return `${set}, but the deep air still glows with its light.${rise}`;
  }
  if (condition === 2) return `${set}, and Earth now lights the night.${rise}`;
  if (condition === 3) return `${set}, and only a faint glow of twilight remains.${rise}`;
  return `${reach === "never" ? "Earth never rises here, and the" : "The"} sky is dark and full of stars.${rise}`;
}

export function earthLine(s, reach, nextEarthrise, now) {
  if (reach === "never") return "Earth never rises here";
  const phase = earthPhaseWords(s.earth_lit_fraction);
  if (s.earth_elevation_deg + s.earth_radius_deg <= 0)
    return `Earth is below the horizon${nextEarthrise ? `, rising ${fromNow(nextEarthrise - now)}` : ""}`;
  return `Earth is ${phase}, ${placeInSky(s.earth_elevation_deg, s.earth_azimuth_deg)}`;
}
