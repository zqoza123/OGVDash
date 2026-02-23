import { Suspense, lazy, useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  Github,
  Globe,
  Linkedin,
  Play,
  Pause,
  RotateCcw,
  Plus,
  Trash2,
  Download,
  ImageDown
} from "lucide-react";

/** @typedef {{ S:number, K:number, T:number, sigma:number, r:number, q:number, optionType:'call'|'put' }} OptionParams */
/** @typedef {'price'|'delta'|'gamma'|'theta'|'vega'|'rho'|'vanna'|'volga'|'charm'|'color'|'speed'|'dualDelta'|'dualGamma'} GreekName */
/** @typedef {{ id:string, instrumentType:'option'|'stock', optionType:'call'|'put', strike:number, qty:number, entryPrice:number, T:number, multiplier:number }} StrategyLeg */
/** @typedef {{ x:number[], y:number[], z:number[][] }} SurfaceData */
/** @typedef {{ id:string, label:string, instrumentType:string, optionType:string, strike:number, T:number, qty:number, delta:number, gamma:number, theta:number, vega:number, vanna:number, rho:number, volga:number }} PortfolioGreekRow */

const YEAR_DAYS = 365;
const T_EPS = 1e-4;
const SIGMA_EPS = 1e-4;
const MIN_POS = 1e-12;
const SQRT_TWO_PI = Math.sqrt(2 * Math.PI);

export const DEFAULT_PARAMS = {
  S: 100,
  K: 100,
  T: 0.5,
  sigma: 0.3,
  r: 0.05,
  q: 0,
  optionType: "call"
};

const TAB_LIST = [
  { id: "surface", label: "Greek Surfaces" },
  { id: "theta", label: "Theta Decay" },
  { id: "strategy", label: "Strategy P&L" },
  { id: "portfolio", label: "Portfolio Greeks" },
  { id: "learn", label: "Learn" }
];

const GREEK_GROUPS = [
  {
    label: "First Order",
    options: ["price", "delta", "vega", "theta", "rho"]
  },
  {
    label: "Second Order",
    options: ["gamma", "vanna", "volga", "charm"]
  },
  {
    label: "Third Order",
    options: ["color", "speed"]
  },
  {
    label: "Dual Greeks",
    options: ["dualDelta", "dualGamma"]
  }
];

const GREEK_LABELS = {
  price: "Price",
  delta: "Delta",
  gamma: "Gamma",
  theta: "Theta",
  vega: "Vega",
  rho: "Rho",
  vanna: "Vanna",
  volga: "Volga",
  charm: "Charm",
  color: "Color",
  speed: "Speed",
  dualDelta: "Dual Delta",
  dualGamma: "Dual Gamma"
};

export const PRESET_STRATEGIES = {
  "Long Call": [{ instrumentType: "option", optionType: "call", strikeOffset: 0, qty: 1 }],
  "Long Put": [{ instrumentType: "option", optionType: "put", strikeOffset: 0, qty: 1 }],
  "Covered Call": [
    { instrumentType: "stock", qty: 1 },
    { instrumentType: "option", optionType: "call", strikeOffset: 5, qty: -1 }
  ],
  "Protective Put": [
    { instrumentType: "stock", qty: 1 },
    { instrumentType: "option", optionType: "put", strikeOffset: -5, qty: 1 }
  ],
  Straddle: [
    { instrumentType: "option", optionType: "call", strikeOffset: 0, qty: 1 },
    { instrumentType: "option", optionType: "put", strikeOffset: 0, qty: 1 }
  ],
  Strangle: [
    { instrumentType: "option", optionType: "call", strikeOffset: 10, qty: 1 },
    { instrumentType: "option", optionType: "put", strikeOffset: -10, qty: 1 }
  ],
  "Bull Call Spread": [
    { instrumentType: "option", optionType: "call", strikeOffset: -5, qty: 1 },
    { instrumentType: "option", optionType: "call", strikeOffset: 5, qty: -1 }
  ],
  "Bear Put Spread": [
    { instrumentType: "option", optionType: "put", strikeOffset: 5, qty: 1 },
    { instrumentType: "option", optionType: "put", strikeOffset: -5, qty: -1 }
  ],
  "Iron Condor": [
    { instrumentType: "option", optionType: "put", strikeOffset: -15, qty: 1 },
    { instrumentType: "option", optionType: "put", strikeOffset: -5, qty: -1 },
    { instrumentType: "option", optionType: "call", strikeOffset: 5, qty: -1 },
    { instrumentType: "option", optionType: "call", strikeOffset: 15, qty: 1 }
  ],
  "Iron Butterfly": [
    { instrumentType: "option", optionType: "put", strikeOffset: -10, qty: 1 },
    { instrumentType: "option", optionType: "put", strikeOffset: 0, qty: -1 },
    { instrumentType: "option", optionType: "call", strikeOffset: 0, qty: -1 },
    { instrumentType: "option", optionType: "call", strikeOffset: 10, qty: 1 }
  ],
  "Calendar Spread": [
    {
      instrumentType: "option",
      optionType: "call",
      strikeOffset: 0,
      qty: -1,
      expiryMult: 0.5
    },
    {
      instrumentType: "option",
      optionType: "call",
      strikeOffset: 0,
      qty: 1,
      expiryMult: 1.0
    }
  ]
};

const LEARN_ROWS = [
  { greek: "Delta", symbol: "Δ", measure: "dV/dS", intuition: "Option value move per $1 spot move" },
  { greek: "Gamma", symbol: "Γ", measure: "d²V/dS²", intuition: "How fast delta changes" },
  { greek: "Theta", symbol: "Θ", measure: "dV/dt", intuition: "Time decay for long options" },
  { greek: "Vega", symbol: "ν", measure: "dV/dσ", intuition: "Sensitivity to implied vol" },
  { greek: "Rho", symbol: "ρ", measure: "dV/dr", intuition: "Sensitivity to rates" },
  { greek: "Vanna", symbol: "", measure: "d²V/dSdσ", intuition: "Delta-Vol interaction" },
  { greek: "Volga", symbol: "", measure: "d²V/dσ²", intuition: "Vega convexity" },
  { greek: "Charm", symbol: "", measure: "d²V/dSdt", intuition: "Delta decay over time" },
  { greek: "Color", symbol: "", measure: "d³V/dS²dt", intuition: "Gamma decay over time" },
  { greek: "Dual Delta", symbol: "", measure: "dV/dK", intuition: "Strike sensitivity" },
  { greek: "Dual Gamma", symbol: "", measure: "d²V/dK²", intuition: "Change of strike sensitivity" }
];

const CHART_COLORS = ["#1d4ed8", "#0f766e", "#b45309", "#7e22ce", "#be123c", "#065f46", "#1d4ed8"];

const Plot = lazy(() => import("react-plotly.js"));
const PortfolioGreekBarChart = lazy(() => import("./components/PortfolioGreekBarChart"));

const clamp = (v, min, max) => Math.min(max, Math.max(min, v));

const isFiniteNumber = (value) => Number.isFinite(value) && !Number.isNaN(value);

const linspace = (start, end, count) => {
  if (count <= 1) {
    return [start];
  }
  const step = (end - start) / (count - 1);
  const values = new Array(count);
  for (let i = 0; i < count; i += 1) {
    values[i] = start + step * i;
  }
  return values;
};

const parseNumeric = (value, fallback) => {
  const parsed = Number(value);
  return isFiniteNumber(parsed) ? parsed : fallback;
};

const formatCurrency = (value, decimals = 2) =>
  `$${Number(value || 0).toLocaleString(undefined, {
    minimumFractionDigits: decimals,
    maximumFractionDigits: decimals
  })}`;

const formatSigned = (value, decimals = 4) => {
  const num = Number(value || 0);
  const prefix = num > 0 ? "+" : "";
  return `${prefix}${num.toFixed(decimals)}`;
};

const formatPercent = (value, decimals = 2) => `${(value * 100).toFixed(decimals)}%`;

const optionSign = (optionType) => (optionType === "call" ? 1 : -1);

const sanitizeParams = (params) => {
  const S = Math.max(parseNumeric(params.S, DEFAULT_PARAMS.S), MIN_POS);
  const K = Math.max(parseNumeric(params.K, DEFAULT_PARAMS.K), MIN_POS);
  const T = Math.max(parseNumeric(params.T, DEFAULT_PARAMS.T), 0);
  const sigma = Math.max(parseNumeric(params.sigma, DEFAULT_PARAMS.sigma), 0);
  const r = parseNumeric(params.r, DEFAULT_PARAMS.r);
  const q = parseNumeric(params.q, DEFAULT_PARAMS.q);
  const optionType = params.optionType === "put" ? "put" : "call";
  return { S, K, T, sigma, r, q, optionType };
};

export function normalPDF(x) {
  return Math.exp(-0.5 * x * x) / SQRT_TWO_PI;
}

export function normalCDF(x) {
  if (x <= -10) return 0;
  if (x >= 10) return 1;

  const a1 = 0.254829592;
  const a2 = -0.284496736;
  const a3 = 1.421413741;
  const a4 = -1.453152027;
  const a5 = 1.061405429;
  const p = 0.3275911;

  const sign = x < 0 ? -1 : 1;
  const absX = Math.abs(x) / Math.sqrt(2);
  const t = 1 / (1 + p * absX);
  const y = 1 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.exp(-absX * absX));
  return 0.5 * (1 + sign * y);
}

const computeD1D2 = (params) => {
  const safe = sanitizeParams(params);
  const TEff = Math.max(safe.T, T_EPS);
  const sigmaEff = Math.max(safe.sigma, SIGMA_EPS);
  const sqrtT = Math.sqrt(TEff);
  const logSK = Math.log(clamp(safe.S / safe.K, 1e-8, 1e8));
  const carry = safe.r - safe.q;
  const d1 = (logSK + (carry + 0.5 * sigmaEff * sigmaEff) * TEff) / (sigmaEff * sqrtT);
  const d2 = d1 - sigmaEff * sqrtT;
  return { ...safe, TEff, sigmaEff, sqrtT, d1, d2, carry };
};

const deterministicPrice = (params) => {
  const safe = sanitizeParams(params);
  const sign = optionSign(safe.optionType);
  if (safe.T <= T_EPS) {
    return Math.max(sign * (safe.S - safe.K), 0);
  }

  const forward = safe.S * Math.exp((safe.r - safe.q) * safe.T);
  const payoff = Math.max(sign * (forward - safe.K), 0);
  return payoff * Math.exp(-safe.r * safe.T);
};

export function blackScholesPrice(params) {
  const safe = sanitizeParams(params);
  if (safe.T <= T_EPS || safe.sigma <= SIGMA_EPS) {
    return Math.max(deterministicPrice(safe), 0);
  }

  const { S, K, T, r, q, optionType } = safe;
  const { d1, d2 } = computeD1D2(safe);
  const discQ = Math.exp(-q * T);
  const discR = Math.exp(-r * T);
  let price;

  if (optionType === "call") {
    price = S * discQ * normalCDF(d1) - K * discR * normalCDF(d2);
  } else {
    price = K * discR * normalCDF(-d2) - S * discQ * normalCDF(-d1);
  }

  if (!isFiniteNumber(price)) {
    return 0;
  }

  return Math.max(price, 0);
}

const finiteDiffGreek = (params, greekName) => {
  const safe = sanitizeParams(params);
  const hS = Math.max(0.01, safe.S * 1e-3);
  const hT = Math.max(1 / 3650, safe.T * 1e-3);
  const hVol = Math.max(1e-4, safe.sigma * 1e-3);
  const hRate = 1e-4;
  const hK = Math.max(0.01, safe.K * 1e-3);

  const withS = (S) => blackScholesPrice({ ...safe, S: Math.max(S, MIN_POS) });
  const withT = (T) => blackScholesPrice({ ...safe, T: Math.max(T, 0) });
  const withVol = (sigma) => blackScholesPrice({ ...safe, sigma: Math.max(sigma, SIGMA_EPS) });
  const withRate = (r) => blackScholesPrice({ ...safe, r });
  const withK = (K) => blackScholesPrice({ ...safe, K: Math.max(K, MIN_POS) });

  const base = blackScholesPrice(safe);

  switch (greekName) {
    case "price":
      return base;
    case "delta":
      return (withS(safe.S + hS) - withS(safe.S - hS)) / (2 * hS);
    case "gamma":
      return (withS(safe.S + hS) - 2 * base + withS(safe.S - hS)) / (hS * hS);
    case "theta":
      return (withT(safe.T - hT) - withT(safe.T + hT)) / (2 * hT);
    case "vega":
      return (withVol(safe.sigma + hVol) - withVol(safe.sigma - hVol)) / (2 * hVol);
    case "rho":
      return (withRate(safe.r + hRate) - withRate(safe.r - hRate)) / (2 * hRate);
    case "dualDelta":
      return (withK(safe.K + hK) - withK(safe.K - hK)) / (2 * hK);
    case "dualGamma":
      return (withK(safe.K + hK) - 2 * base + withK(safe.K - hK)) / (hK * hK);
    case "charm": {
      const deltaTMinus = finiteDiffGreek({ ...safe, T: Math.max(safe.T - hT, 0) }, "delta");
      const deltaTPlus = finiteDiffGreek({ ...safe, T: safe.T + hT }, "delta");
      return (deltaTMinus - deltaTPlus) / (2 * hT);
    }
    case "color": {
      const gammaTMinus = finiteDiffGreek({ ...safe, T: Math.max(safe.T - hT, 0) }, "gamma");
      const gammaTPlus = finiteDiffGreek({ ...safe, T: safe.T + hT }, "gamma");
      return (gammaTMinus - gammaTPlus) / (2 * hT);
    }
    case "speed": {
      const gammaPlus = finiteDiffGreek({ ...safe, S: safe.S + hS }, "gamma");
      const gammaMinus = finiteDiffGreek({ ...safe, S: safe.S - hS }, "gamma");
      return (gammaPlus - gammaMinus) / (2 * hS);
    }
    case "vanna": {
      const deltaVolPlus = finiteDiffGreek({ ...safe, sigma: safe.sigma + hVol }, "delta");
      const deltaVolMinus = finiteDiffGreek({ ...safe, sigma: Math.max(safe.sigma - hVol, SIGMA_EPS) }, "delta");
      return (deltaVolPlus - deltaVolMinus) / (2 * hVol);
    }
    case "volga": {
      const vegaPlus = finiteDiffGreek({ ...safe, sigma: safe.sigma + hVol }, "vega");
      const vegaMinus = finiteDiffGreek({ ...safe, sigma: Math.max(safe.sigma - hVol, SIGMA_EPS) }, "vega");
      return (vegaPlus - vegaMinus) / (2 * hVol);
    }
    default:
      return 0;
  }
};

export function computeGreek(params, greekName) {
  const safe = sanitizeParams(params);
  const unstable = safe.T < 2 * T_EPS || safe.sigma < 2 * SIGMA_EPS;
  const fallbackGreeks = new Set(["charm", "color", "speed"]);

  if (greekName === "price") {
    return blackScholesPrice(safe);
  }

  if (safe.T <= T_EPS || safe.sigma <= SIGMA_EPS) {
    return finiteDiffGreek(safe, greekName);
  }

  const { S, K, T, r, q, optionType } = safe;
  const { d1, d2, sqrtT, sigmaEff, TEff, carry } = computeD1D2(safe);
  const discQ = Math.exp(-q * T);
  const discR = Math.exp(-r * T);
  const n1 = normalPDF(d1);
  const nd1 = normalCDF(d1);
  const nd2 = normalCDF(d2);

  let analytic;

  switch (greekName) {
    case "delta":
      analytic = optionType === "call" ? discQ * nd1 : discQ * (nd1 - 1);
      break;
    case "gamma":
      analytic = (discQ * n1) / (S * sigmaEff * sqrtT);
      break;
    case "theta": {
      const firstTerm = -(S * discQ * n1 * sigmaEff) / (2 * sqrtT);
      if (optionType === "call") {
        analytic = firstTerm - r * K * discR * nd2 + q * S * discQ * nd1;
      } else {
        analytic = firstTerm + r * K * discR * normalCDF(-d2) - q * S * discQ * normalCDF(-d1);
      }
      break;
    }
    case "vega":
      analytic = S * discQ * n1 * sqrtT;
      break;
    case "rho":
      analytic = optionType === "call" ? K * T * discR * nd2 : -K * T * discR * normalCDF(-d2);
      break;
    case "vanna":
      analytic = -(discQ * n1 * d2) / sigmaEff;
      break;
    case "volga":
      analytic = (S * discQ * n1 * sqrtT * d1 * d2) / sigmaEff;
      break;
    case "charm": {
      const common = (discQ * n1 * (2 * carry * TEff - d2 * sigmaEff * sqrtT)) / (2 * TEff * sigmaEff * sqrtT);
      if (optionType === "call") {
        analytic = q * discQ * nd1 - common;
      } else {
        analytic = -q * discQ * normalCDF(-d1) - common;
      }
      break;
    }
    case "color": {
      const bracket =
        2 * q * TEff +
        1 +
        (d1 * (2 * carry * TEff - d2 * sigmaEff * sqrtT)) / (sigmaEff * sqrtT);
      analytic = (-(discQ * n1) / (2 * S * TEff * sigmaEff * sqrtT)) * bracket;
      break;
    }
    case "speed": {
      const gamma = (discQ * n1) / (S * sigmaEff * sqrtT);
      analytic = -(gamma / S) * (d1 / (sigmaEff * sqrtT) + 1);
      break;
    }
    case "dualDelta":
      analytic = optionType === "call" ? -discR * nd2 : discR * normalCDF(-d2);
      break;
    case "dualGamma":
      analytic = (discR * normalPDF(d2)) / (K * sigmaEff * sqrtT);
      break;
    default:
      analytic = 0;
  }

  if (!isFiniteNumber(analytic) || (unstable && fallbackGreeks.has(greekName))) {
    return finiteDiffGreek(safe, greekName);
  }

  return analytic;
}

export function computeSurfaceData(params, greekName, resolution = 50) {
  const safe = sanitizeParams(params);
  const strikes = linspace(Math.max(1, safe.S * 0.5), safe.S * 1.5, resolution);
  const expiries = linspace(1 / YEAR_DAYS, 2, resolution);
  const z = new Array(expiries.length);

  for (let yi = 0; yi < expiries.length; yi += 1) {
    const T = expiries[yi];
    const row = new Array(strikes.length);
    for (let xi = 0; xi < strikes.length; xi += 1) {
      const K = strikes[xi];
      row[xi] = computeGreek({ ...safe, K, T }, greekName);
    }
    z[yi] = row;
  }

  return { x: strikes, y: expiries, z };
}

const createLegId = () => Math.random().toString(36).slice(2, 10);

const computeLegEntryPrice = (leg, params) => {
  if (leg.instrumentType === "stock") {
    return params.S;
  }
  return blackScholesPrice({
    ...params,
    optionType: leg.optionType,
    K: leg.strike,
    T: leg.T
  });
};

export function instantiatePreset(name, params) {
  const templates = PRESET_STRATEGIES[name] || PRESET_STRATEGIES["Long Call"];
  return templates.map((template) => {
    const leg = {
      id: createLegId(),
      instrumentType: template.instrumentType,
      optionType: template.optionType || "call",
      strike: template.instrumentType === "option" ? Math.max(1, params.S + (template.strikeOffset || 0)) : params.S,
      qty: template.qty,
      T: Math.max(1 / YEAR_DAYS, params.T * (template.expiryMult || 1)),
      multiplier: 100,
      entryPrice: 0
    };
    return {
      ...leg,
      entryPrice: computeLegEntryPrice(leg, params)
    };
  });
}

const getLegRemainingT = (leg, evalDte, maxDte) => {
  const elapsedYears = Math.max(0, (maxDte - evalDte) / YEAR_DAYS);
  return Math.max(leg.T - elapsedYears, 0);
};

const computeLegValue = (leg, stateParams, remainingT) => {
  if (leg.instrumentType === "stock") {
    return stateParams.S;
  }
  return blackScholesPrice({
    ...stateParams,
    optionType: leg.optionType,
    K: leg.strike,
    T: remainingT
  });
};

export function computeStrategyPnL(legs, params, evalDte) {
  if (!Array.isArray(legs) || legs.length === 0) {
    return 0;
  }

  const safe = sanitizeParams(params);
  const maxDte = Math.max(...legs.map((leg) => Math.round(leg.T * YEAR_DAYS)), 1);
  const targetDte = clamp(Math.round(evalDte ?? maxDte), 0, maxDte);

  let total = 0;
  for (const leg of legs) {
    const remainingT = getLegRemainingT(leg, targetDte, maxDte);
    const value = computeLegValue(leg, safe, remainingT);
    total += leg.qty * (leg.multiplier || 100) * (value - leg.entryPrice);
  }

  return total;
}

export function computePnLHeatmap(legs, params, resolution = 80, evalDte) {
  const safe = sanitizeParams(params);
  const spotRange = linspace(Math.max(1, safe.S * 0.7), safe.S * 1.3, resolution);
  const volRange = linspace(0.05, Math.max(0.06, Math.min(1.5, safe.sigma * 2.5)), resolution);
  const z = new Array(volRange.length);

  for (let yi = 0; yi < volRange.length; yi += 1) {
    const sigma = volRange[yi];
    const row = new Array(spotRange.length);
    for (let xi = 0; xi < spotRange.length; xi += 1) {
      const S = spotRange[xi];
      row[xi] = computeStrategyPnL(legs, { ...safe, S, sigma }, evalDte);
    }
    z[yi] = row;
  }

  return { x: spotRange, y: volRange, z };
}

export function findBreakevens(legs, params, sampleCount = 400) {
  const safe = sanitizeParams(params);
  const x = linspace(Math.max(1, safe.S * 0.5), safe.S * 1.5, sampleCount);
  const y = x.map((spot) => computeStrategyPnL(legs, { ...safe, S: spot }, 0));

  const points = [];
  for (let i = 1; i < x.length; i += 1) {
    const y1 = y[i - 1];
    const y2 = y[i];
    if (y1 === 0) points.push(x[i - 1]);
    if (y1 * y2 < 0) {
      const ratio = Math.abs(y1) / (Math.abs(y1) + Math.abs(y2));
      points.push(x[i - 1] + (x[i] - x[i - 1]) * ratio);
    }
  }

  return Array.from(new Set(points.map((point) => Number(point.toFixed(4)))));
}

export function computePortfolioRows(legs, params) {
  const safe = sanitizeParams(params);
  /** @type {PortfolioGreekRow[]} */
  const rows = [];
  const totals = {
    delta: 0,
    gamma: 0,
    theta: 0,
    vega: 0,
    vanna: 0,
    rho: 0,
    volga: 0
  };

  legs.forEach((leg, index) => {
    const multiplier = leg.multiplier || 100;
    let greekValues;

    if (leg.instrumentType === "stock") {
      greekValues = {
        delta: 1,
        gamma: 0,
        theta: 0,
        vega: 0,
        vanna: 0,
        rho: 0,
        volga: 0
      };
    } else {
      const legParams = {
        ...safe,
        optionType: leg.optionType,
        K: leg.strike,
        T: leg.T
      };
      greekValues = {
        delta: computeGreek(legParams, "delta"),
        gamma: computeGreek(legParams, "gamma"),
        theta: computeGreek(legParams, "theta"),
        vega: computeGreek(legParams, "vega"),
        vanna: computeGreek(legParams, "vanna"),
        rho: computeGreek(legParams, "rho"),
        volga: computeGreek(legParams, "volga")
      };
    }

    const row = {
      id: leg.id,
      label: `Leg ${index + 1}`,
      instrumentType: leg.instrumentType,
      optionType: leg.optionType,
      strike: leg.strike,
      T: leg.T,
      qty: leg.qty,
      delta: leg.qty * multiplier * greekValues.delta,
      gamma: leg.qty * multiplier * greekValues.gamma,
      theta: leg.qty * multiplier * greekValues.theta,
      vega: leg.qty * multiplier * greekValues.vega,
      vanna: leg.qty * multiplier * greekValues.vanna,
      rho: leg.qty * multiplier * greekValues.rho,
      volga: leg.qty * multiplier * greekValues.volga
    };

    totals.delta += row.delta;
    totals.gamma += row.gamma;
    totals.theta += row.theta;
    totals.vega += row.vega;
    totals.vanna += row.vanna;
    totals.rho += row.rho;
    totals.volga += row.volga;
    rows.push(row);
  });

  return { rows, totals };
}

const useRafSetter = (setter) => {
  const frameRef = useRef(0);
  const queuedRef = useRef(null);

  useEffect(() => {
    return () => {
      if (frameRef.current) cancelAnimationFrame(frameRef.current);
    };
  }, []);

  return useCallback(
    (nextState) => {
      queuedRef.current = nextState;
      if (frameRef.current) return;

      frameRef.current = requestAnimationFrame(() => {
        frameRef.current = 0;
        setter((prev) => {
          if (typeof queuedRef.current === "function") {
            return queuedRef.current(prev);
          }
          return queuedRef.current;
        });
      });
    },
    [setter]
  );
};

const downloadFile = (filename, text, contentType) => {
  const blob = new Blob([text], { type: contentType });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
};

const exportSurfaceCsv = (surfaceData, greekName) => {
  const headers = ["expiry\\strike", ...surfaceData.x.map((v) => v.toFixed(4))];
  const rows = surfaceData.y.map((expiry, rowIndex) => [
    expiry.toFixed(6),
    ...surfaceData.z[rowIndex].map((v) => (isFiniteNumber(v) ? String(v) : "0"))
  ]);
  const csv = [headers.join(","), ...rows.map((row) => row.join(","))].join("\n");
  downloadFile(`surface-${greekName}.csv`, csv, "text/csv;charset=utf-8;");
};

const exportHeatmapCsv = (heatmapData) => {
  const headers = ["vol\\spot", ...heatmapData.x.map((v) => v.toFixed(4))];
  const rows = heatmapData.y.map((vol, rowIndex) => [
    vol.toFixed(6),
    ...heatmapData.z[rowIndex].map((v) => (isFiniteNumber(v) ? String(v) : "0"))
  ]);
  const csv = [headers.join(","), ...rows.map((row) => row.join(","))].join("\n");
  downloadFile("strategy-heatmap.csv", csv, "text/csv;charset=utf-8;");
};

const exportPlotImage = async (graphDiv, filename) => {
  if (!graphDiv) return;
  try {
    const plotlyModule = await import("plotly.js-dist-min");
    const plotly = plotlyModule?.default || plotlyModule;
    if (!plotly?.toImage) return;
    const dataUrl = await plotly.toImage(graphDiv, {
      format: "png",
      height: 720,
      width: 1200
    });
    const link = document.createElement("a");
    link.href = dataUrl;
    link.download = filename;
    link.click();
  } catch {
    // Ignore export failure when the plotting bundle is unavailable.
  }
};

const numericInputBaseClass =
  "w-full rounded border border-slate-500 bg-white px-2 py-1 text-sm text-slate-900 dark:bg-slate-700 dark:text-slate-100";

function SliderField({ label, value, min, max, step, onChange, numberStep, suffix, scale = 1 }) {
  const displayValue = value * scale;
  const [inputText, setInputText] = useState(displayValue.toString());

  useEffect(() => {
    setInputText(displayValue.toString());
  }, [displayValue]);

  const handleSlider = (evt) => {
    onChange(parseNumeric(evt.target.value, displayValue) / scale);
  };

  const commitNumberInput = () => {
    const parsed = Number(inputText);
    if (!Number.isFinite(parsed)) {
      setInputText(displayValue.toString());
      return;
    }
    const clamped = clamp(parsed, min * scale, max * scale);
    onChange(clamped / scale);
    setInputText(clamped.toString());
  };

  const handleNumberChange = (evt) => {
    setInputText(evt.target.value);
  };

  const handleNumberKeyDown = (evt) => {
    if (evt.key === "Enter") {
      evt.currentTarget.blur();
    }
  };

  return (
    <div className="space-y-1">
      <div className="flex items-center justify-between text-sm">
        <label>{label}</label>
        <span className="font-mono text-xs">
          {displayValue.toFixed(2)}
          {suffix || ""}
        </span>
      </div>
      <div className="grid grid-cols-[1fr_110px] gap-2">
        <input
          type="range"
          min={min * scale}
          max={max * scale}
          step={(step || 0.01) * scale}
          value={displayValue}
          onChange={handleSlider}
        />
        <input
          type="number"
          className={numericInputBaseClass}
          value={inputText}
          min={min * scale}
          max={max * scale}
          step={(numberStep || step || 0.01) * scale}
          onChange={handleNumberChange}
          onBlur={commitNumberInput}
          onKeyDown={handleNumberKeyDown}
        />
      </div>
    </div>
  );
}

function ParameterPanel({ params, onChangeParam, onReset, quickGreeks }) {
  return (
    <div className="space-y-4 rounded border border-slate-500 bg-slate-100 p-3 dark:bg-slate-800">
      <h2 className="text-sm font-semibold uppercase tracking-wide">Parameter Panel</h2>
      <SliderField label="Spot Price ($)" value={params.S} min={10} max={500} step={0.5} onChange={(v) => onChangeParam("S", clamp(v, 10, 500))} />
      <SliderField
        label="Implied Volatility"
        value={params.sigma}
        min={0.05}
        max={1.5}
        step={0.005}
        numberStep={0.001}
        scale={100}
        suffix="%"
        onChange={(v) => onChangeParam("sigma", clamp(v, 0.05, 1.5))}
      />
      <SliderField
        label="Risk-Free Rate"
        value={params.r}
        min={0}
        max={0.15}
        step={0.001}
        scale={100}
        suffix="%"
        onChange={(v) => onChangeParam("r", clamp(v, 0, 0.15))}
      />
      <SliderField
        label="Dividend Yield"
        value={params.q}
        min={0}
        max={0.1}
        step={0.001}
        scale={100}
        suffix="%"
        onChange={(v) => onChangeParam("q", clamp(v, 0, 0.1))}
      />
      <SliderField
        label="Strike ($)"
        value={params.K}
        min={10}
        max={500}
        step={0.5}
        onChange={(v) => onChangeParam("K", clamp(v, 10, 500))}
      />
      <SliderField
        label="Expiry (years)"
        value={params.T}
        min={1 / YEAR_DAYS}
        max={2}
        step={0.0025}
        numberStep={0.01}
        onChange={(v) => onChangeParam("T", clamp(v, 1 / YEAR_DAYS, 2))}
      />

      <div className="space-y-2 text-sm">
        <div>Option Type</div>
        <div className="flex gap-3">
          <label className="inline-flex items-center gap-2">
            <input
              type="radio"
              checked={params.optionType === "call"}
              onChange={() => onChangeParam("optionType", "call")}
            />
            Call
          </label>
          <label className="inline-flex items-center gap-2">
            <input
              type="radio"
              checked={params.optionType === "put"}
              onChange={() => onChangeParam("optionType", "put")}
            />
            Put
          </label>
        </div>
      </div>

      <button
        type="button"
        className="w-full rounded border border-slate-600 bg-slate-200 px-3 py-2 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
        onClick={onReset}
      >
        Reset Defaults
      </button>

      <div className="rounded border border-slate-500 bg-white p-2 text-xs dark:bg-slate-700">
        <div className="mb-2 font-semibold">Quick Greeks</div>
        <div className="grid grid-cols-2 gap-x-3 gap-y-1 font-mono">
          <div>Price</div>
          <div className="text-right">{formatCurrency(quickGreeks.price)}</div>
          <div>Delta</div>
          <div className="text-right">{quickGreeks.delta.toFixed(4)}</div>
          <div>Gamma</div>
          <div className="text-right">{quickGreeks.gamma.toFixed(6)}</div>
          <div>Theta</div>
          <div className="text-right">{quickGreeks.theta.toFixed(4)}</div>
          <div>Vega</div>
          <div className="text-right">{quickGreeks.vega.toFixed(4)}</div>
          <div>Rho</div>
          <div className="text-right">{quickGreeks.rho.toFixed(4)}</div>
          <div>Vanna</div>
          <div className="text-right">{quickGreeks.vanna.toFixed(4)}</div>
          <div>Volga</div>
          <div className="text-right">{quickGreeks.volga.toFixed(4)}</div>
        </div>
      </div>
    </div>
  );
}

function GreekSurfaceView({
  params,
  selectedGreek,
  setSelectedGreek,
  resolution,
  setResolution,
  surfaceData,
  onExportCsv,
  onExportImage,
  onGraphReady
}) {
  const plotData = useMemo(
    () => [
      {
        type: "surface",
        x: surfaceData.x,
        y: surfaceData.y,
        z: surfaceData.z,
        colorscale: "RdBu",
        colorbar: { title: GREEK_LABELS[selectedGreek], thickness: 14 },
        contours: {
          z: {
            show: true,
            usecolormap: true,
            highlightcolor: "#111827",
            project: { z: true }
          }
        },
        hovertemplate: `Strike: $%{x:.2f}<br>Expiry: %{y:.2f}y<br>${GREEK_LABELS[selectedGreek]}: %{z:.6f}<extra></extra>`
      }
    ],
    [surfaceData, selectedGreek]
  );

  return (
    <div className="space-y-3">
      <div className="flex flex-wrap items-center gap-2">
        <label className="text-sm">Greek</label>
        <select
          className={`${numericInputBaseClass} max-w-56`}
          value={selectedGreek}
          onChange={(evt) => setSelectedGreek(evt.target.value)}
        >
          {GREEK_GROUPS.map((group) => (
            <optgroup key={group.label} label={group.label}>
              {group.options.map((name) => (
                <option key={name} value={name}>
                  {GREEK_LABELS[name]}
                </option>
              ))}
            </optgroup>
          ))}
        </select>
        <label className="text-sm">Resolution</label>
        <select
          className={`${numericInputBaseClass} max-w-24`}
          value={resolution}
          onChange={(evt) => setResolution(Number(evt.target.value))}
        >
          {[30, 40, 50].map((res) => (
            <option key={res} value={res}>
              {res}
            </option>
          ))}
        </select>
        <button
          type="button"
          className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
          onClick={onExportCsv}
        >
          <Download className="h-4 w-4" /> CSV
        </button>
        <button
          type="button"
          className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
          onClick={onExportImage}
        >
          <ImageDown className="h-4 w-4" /> PNG
        </button>
      </div>
      <div className="h-[560px] rounded border border-slate-500 bg-white p-1 dark:bg-slate-900">
        <Suspense fallback={<div className="grid h-full place-items-center text-sm text-slate-500">Loading surface chart...</div>}>
          <Plot
            data={plotData}
            onInitialized={(figure, graphDiv) => onGraphReady(graphDiv, figure)}
            onUpdate={(figure, graphDiv) => onGraphReady(graphDiv, figure)}
            layout={{
              margin: { l: 0, r: 0, t: 30, b: 0 },
              paper_bgcolor: "rgba(0,0,0,0)",
              plot_bgcolor: "rgba(0,0,0,0)",
              scene: {
                xaxis: { title: "Strike Price ($)" },
                yaxis: { title: "Time to Expiry (years)" },
                zaxis: { title: GREEK_LABELS[selectedGreek] },
                camera: { eye: { x: 1.5, y: 1.4, z: 1.2 } }
              },
              title: `${GREEK_LABELS[selectedGreek]} Surface (${params.optionType.toUpperCase()})`
            }}
            config={{ responsive: true, displaylogo: false }}
            style={{ width: "100%", height: "100%" }}
          />
        </Suspense>
      </div>
    </div>
  );
}

function ThetaDecayAnimator({ params }) {
  const [animState, setAnimState] = useState("IDLE");
  const [speed, setSpeed] = useState(1);
  const initialDte = Math.max(1, Math.round(params.T * YEAR_DAYS));
  const [currentDte, setCurrentDte] = useState(initialDte);

  useEffect(() => {
    setAnimState("IDLE");
    setCurrentDte(initialDte);
  }, [initialDte, params.S, params.K, params.sigma, params.r, params.q, params.optionType]);

  useEffect(() => {
    if (animState !== "RUNNING") return;

    const interval = setInterval(() => {
      setCurrentDte((prev) => {
        const next = prev - speed;
        if (next <= 0) {
          setAnimState("FINISHED");
          return 0;
        }
        return next;
      });
    }, 100);

    return () => clearInterval(interval);
  }, [animState, speed]);

  const curve = useMemo(() => {
    const rows = [];
    const startValue = blackScholesPrice(params);
    for (let dte = initialDte; dte >= 0; dte -= 1) {
      const value = blackScholesPrice({ ...params, T: dte / YEAR_DAYS });
      rows.push({ dte, value, decay: startValue - value });
    }
    return rows;
  }, [initialDte, params]);

  const currentValue = useMemo(
    () => blackScholesPrice({ ...params, T: Math.max(currentDte, 0) / YEAR_DAYS }),
    [currentDte, params]
  );

  const thetaChart = useMemo(() => {
    const width = 980;
    const height = 360;
    const padding = { top: 14, right: 16, bottom: 30, left: 24 };

    if (curve.length === 0) {
      return {
        width,
        height,
        linePath: "",
        areaPath: "",
        markerX: padding.left,
        axisY: height - padding.bottom
      };
    }

    const usableWidth = width - padding.left - padding.right;
    const usableHeight = height - padding.top - padding.bottom;
    const minValue = Math.min(...curve.map((point) => point.value));
    const maxValue = Math.max(...curve.map((point) => point.value));
    const valueRange = Math.max(maxValue - minValue, 1e-9);

    const pointToSvg = (point) => {
      const x =
        padding.left +
        ((initialDte - point.dte) / Math.max(initialDte, 1)) * usableWidth;
      const y =
        padding.top +
        (1 - (point.value - minValue) / valueRange) * usableHeight;
      return { x, y };
    };

    const points = curve.map(pointToSvg);
    const linePath = points
      .map((point, index) => `${index === 0 ? "M" : "L"} ${point.x.toFixed(2)} ${point.y.toFixed(2)}`)
      .join(" ");

    const axisY = height - padding.bottom;
    const areaPath = `${linePath} L ${points[points.length - 1].x.toFixed(2)} ${axisY.toFixed(2)} L ${points[0].x.toFixed(2)} ${axisY.toFixed(2)} Z`;

    const markerX =
      padding.left +
      ((initialDte - Math.max(0, Math.min(currentDte, initialDte))) / Math.max(initialDte, 1)) *
        usableWidth;

    return { width, height, linePath, areaPath, markerX, axisY };
  }, [curve, currentDte, initialDte]);

  return (
    <div className="space-y-3">
      <div className="grid gap-2 rounded border border-slate-500 bg-slate-100 p-3 dark:bg-slate-800 md:grid-cols-[1fr_auto]">
        <div>
          <div className="text-sm font-semibold">Theta Decay Animator</div>
          <div className="text-2xl font-bold">{currentDte} days remaining</div>
          <div className="text-xs text-slate-600 dark:text-slate-300">Current value: {formatCurrency(currentValue)}</div>
        </div>
        <div className="flex flex-wrap items-center gap-2">
          <button
            type="button"
            className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
            onClick={() => setAnimState((prev) => (prev === "RUNNING" ? "PAUSED" : "RUNNING"))}
            disabled={animState === "FINISHED"}
          >
            {animState === "RUNNING" ? <Pause className="h-4 w-4" /> : <Play className="h-4 w-4" />}
            {animState === "RUNNING" ? "Pause" : "Play"}
          </button>
          <button
            type="button"
            className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
            onClick={() => {
              setAnimState("IDLE");
              setCurrentDte(initialDte);
            }}
          >
            <RotateCcw className="h-4 w-4" /> Reset
          </button>
          <select
            value={speed}
            onChange={(evt) => setSpeed(Number(evt.target.value))}
            className={`${numericInputBaseClass} w-24`}
          >
            {[1, 2, 5, 10].map((s) => (
              <option key={s} value={s}>
                {s}x
              </option>
            ))}
          </select>
        </div>
      </div>
      <div className="h-[500px] rounded border border-slate-500 bg-white p-2 dark:bg-slate-900">
        <svg viewBox={`0 0 ${thetaChart.width} ${thetaChart.height}`} className="h-full w-full">
          <rect x="0" y="0" width={thetaChart.width} height={thetaChart.height} fill="transparent" />
          <line
            x1="24"
            y1={thetaChart.axisY}
            x2={thetaChart.width - 16}
            y2={thetaChart.axisY}
            stroke="#64748b"
            strokeWidth="1"
          />
          <path d={thetaChart.areaPath} fill="rgba(148, 163, 184, 0.24)" />
          <path d={thetaChart.linePath} fill="none" stroke="#475569" strokeWidth="2.2" />
          <line x1={thetaChart.markerX} y1="14" x2={thetaChart.markerX} y2={thetaChart.axisY} stroke="#0f172a" strokeWidth="1.5" />
          <text x={thetaChart.markerX + 4} y="24" fill="#334155" fontSize="12">
            Current
          </text>
        </svg>
      </div>
      <div className="grid gap-2 rounded border border-slate-500 bg-slate-100 p-3 text-xs dark:bg-slate-800 md:grid-cols-3">
        <div>Start value: {formatCurrency(curve[0]?.value ?? 0)}</div>
        <div>End value: {formatCurrency(curve[curve.length - 1]?.value ?? 0)}</div>
        <div>Total decay: {formatCurrency((curve[0]?.value ?? 0) - (curve[curve.length - 1]?.value ?? 0))}</div>
      </div>
    </div>
  );
}

function StrategyBuilder({
  params,
  legs,
  setLegs,
  preset,
  setPreset,
  evalDte,
  setEvalDte,
  heatmap,
  breakevens,
  onLoadPreset,
  onMarkEntries,
  onExportCsv,
  onExportImage,
  onHeatmapGraphReady
}) {
  const maxDte = Math.max(1, ...legs.map((leg) => Math.round(leg.T * YEAR_DAYS)));
  const flat = heatmap.z.flat().filter((value) => isFiniteNumber(value));
  const maxProfit = flat.length ? Math.max(...flat) : 0;
  const maxLoss = flat.length ? Math.min(...flat) : 0;

  const updateLeg = (id, field, value) => {
    setLegs((prev) =>
      prev.map((leg) => {
        if (leg.id !== id) return leg;
        if (field === "instrumentType") {
          const nextType = value;
          return {
            ...leg,
            instrumentType: nextType,
            optionType: nextType === "stock" ? "call" : leg.optionType,
            strike: nextType === "stock" ? params.S : leg.strike,
            T: nextType === "stock" ? params.T : leg.T
          };
        }
        return { ...leg, [field]: value };
      })
    );
  };

  const addLeg = () => {
    const leg = {
      id: createLegId(),
      instrumentType: "option",
      optionType: "call",
      strike: params.K,
      qty: 1,
      T: params.T,
      multiplier: 100,
      entryPrice: blackScholesPrice({ ...params, K: params.K, optionType: "call", T: params.T })
    };
    setLegs((prev) => [...prev, leg]);
  };

  const removeLeg = (id) => {
    setLegs((prev) => prev.filter((leg) => leg.id !== id));
  };

  const heatmapData = [
    {
      type: "heatmap",
      x: heatmap.x,
      y: heatmap.y,
      z: heatmap.z,
      zmid: 0,
      colorscale: [
        [0, "#b91c1c"],
        [0.45, "#fca5a5"],
        [0.5, "#ffffff"],
        [0.55, "#86efac"],
        [1, "#166534"]
      ],
      hovertemplate: "Spot: $%{x:.2f}<br>Vol: %{y:.1%}<br>P&L: $%{z:,.0f}<extra></extra>"
    }
  ];

  return (
    <div className="space-y-3">
      <div className="grid gap-3 rounded border border-slate-500 bg-slate-100 p-3 dark:bg-slate-800 xl:grid-cols-[1fr_auto]">
        <div className="flex flex-wrap items-center gap-2">
          <label className="text-sm">Preset</label>
          <select className={`${numericInputBaseClass} max-w-64`} value={preset} onChange={(evt) => setPreset(evt.target.value)}>
            {Object.keys(PRESET_STRATEGIES).map((name) => (
              <option key={name} value={name}>
                {name}
              </option>
            ))}
          </select>
          <button
            type="button"
            className="rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
            onClick={onLoadPreset}
          >
            Load Preset
          </button>
          <button
            type="button"
            className="rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
            onClick={onMarkEntries}
          >
            Mark Entries @ Current
          </button>
          <button
            type="button"
            className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
            onClick={addLeg}
          >
            <Plus className="h-4 w-4" /> Add Leg
          </button>
        </div>
        <div className="flex flex-wrap items-center gap-2">
          <button
            type="button"
            className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
            onClick={onExportCsv}
          >
            <Download className="h-4 w-4" /> CSV
          </button>
          <button
            type="button"
            className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
            onClick={onExportImage}
          >
            <ImageDown className="h-4 w-4" /> PNG
          </button>
        </div>
      </div>

      <div className="overflow-x-auto rounded border border-slate-500 bg-white dark:bg-slate-900">
        <table className="min-w-full border-collapse text-xs">
          <thead className="bg-slate-200 dark:bg-slate-700">
            <tr>
              <th className="border border-slate-400 p-1 text-left">Instrument</th>
              <th className="border border-slate-400 p-1 text-left">Type</th>
              <th className="border border-slate-400 p-1 text-left">Strike</th>
              <th className="border border-slate-400 p-1 text-left">T (y)</th>
              <th className="border border-slate-400 p-1 text-left">Qty</th>
              <th className="border border-slate-400 p-1 text-left">Entry</th>
              <th className="border border-slate-400 p-1 text-left">Action</th>
            </tr>
          </thead>
          <tbody>
            {legs.map((leg) => (
              <tr key={leg.id}>
                <td className="border border-slate-300 p-1">
                  <select
                    className={numericInputBaseClass}
                    value={leg.instrumentType}
                    onChange={(evt) => updateLeg(leg.id, "instrumentType", evt.target.value)}
                  >
                    <option value="option">Option</option>
                    <option value="stock">Stock</option>
                  </select>
                </td>
                <td className="border border-slate-300 p-1">
                  <select
                    className={numericInputBaseClass}
                    value={leg.optionType}
                    disabled={leg.instrumentType === "stock"}
                    onChange={(evt) => updateLeg(leg.id, "optionType", evt.target.value)}
                  >
                    <option value="call">Call</option>
                    <option value="put">Put</option>
                  </select>
                </td>
                <td className="border border-slate-300 p-1">
                  <input
                    type="number"
                    className={numericInputBaseClass}
                    value={leg.strike}
                    disabled={leg.instrumentType === "stock"}
                    onChange={(evt) => updateLeg(leg.id, "strike", Math.max(1, parseNumeric(evt.target.value, leg.strike)))}
                  />
                </td>
                <td className="border border-slate-300 p-1">
                  <input
                    type="number"
                    className={numericInputBaseClass}
                    value={leg.T}
                    disabled={leg.instrumentType === "stock"}
                    step="0.01"
                    onChange={(evt) =>
                      updateLeg(leg.id, "T", Math.max(1 / YEAR_DAYS, parseNumeric(evt.target.value, leg.T)))
                    }
                  />
                </td>
                <td className="border border-slate-300 p-1">
                  <input
                    type="number"
                    className={numericInputBaseClass}
                    value={leg.qty}
                    step="1"
                    onChange={(evt) => updateLeg(leg.id, "qty", parseNumeric(evt.target.value, leg.qty))}
                  />
                </td>
                <td className="border border-slate-300 p-1">
                  <input
                    type="number"
                    className={numericInputBaseClass}
                    value={leg.entryPrice}
                    step="0.01"
                    onChange={(evt) => updateLeg(leg.id, "entryPrice", parseNumeric(evt.target.value, leg.entryPrice))}
                  />
                </td>
                <td className="border border-slate-300 p-1">
                  <button
                    type="button"
                    className="inline-flex items-center rounded border border-slate-500 p-1 hover:bg-slate-200 dark:hover:bg-slate-700"
                    onClick={() => removeLeg(leg.id)}
                    disabled={legs.length <= 1}
                  >
                    <Trash2 className="h-4 w-4" />
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className="rounded border border-slate-500 bg-slate-100 p-3 dark:bg-slate-800">
        <div className="mb-2 flex items-center justify-between text-sm">
          <label>Evaluation DTE (days): {evalDte}</label>
          <span>Max DTE: {maxDte}</span>
        </div>
        <input type="range" min="0" max={maxDte} step="1" value={evalDte} onChange={(evt) => setEvalDte(Number(evt.target.value))} className="w-full" />
      </div>

      <div className="grid gap-2 rounded border border-slate-500 bg-slate-100 p-3 text-sm dark:bg-slate-800 md:grid-cols-3">
        <div>Max Profit (grid): {formatCurrency(maxProfit, 0)}</div>
        <div>Max Loss (grid): {formatCurrency(maxLoss, 0)}</div>
        <div>Breakeven(s): {breakevens.length ? breakevens.map((b) => formatCurrency(b)).join(", ") : "None"}</div>
      </div>

      <div className="h-[540px] rounded border border-slate-500 bg-white p-1 dark:bg-slate-900">
        <Suspense fallback={<div className="grid h-full place-items-center text-sm text-slate-500">Loading heatmap...</div>}>
          <Plot
            data={heatmapData}
            onInitialized={(figure, graphDiv) => onHeatmapGraphReady(graphDiv, figure)}
            onUpdate={(figure, graphDiv) => onHeatmapGraphReady(graphDiv, figure)}
            layout={{
              title: "Strategy P&L Heatmap",
              margin: { l: 50, r: 20, t: 40, b: 50 },
              xaxis: { title: "Spot Price ($)" },
              yaxis: { title: "Implied Volatility", tickformat: ".0%" },
              shapes: breakevens.map((x) => ({
                type: "line",
                x0: x,
                x1: x,
                y0: heatmap.y[0],
                y1: heatmap.y[heatmap.y.length - 1],
                line: { color: "#111827", width: 1, dash: "dot" }
              })),
              paper_bgcolor: "rgba(0,0,0,0)",
              plot_bgcolor: "rgba(0,0,0,0)"
            }}
            config={{ responsive: true, displaylogo: false }}
            style={{ width: "100%", height: "100%" }}
          />
        </Suspense>
      </div>
    </div>
  );
}

function PortfolioGreeksTable({ params, legs, setLegs, rows, totals, onLoadStrategy }) {
  const updateLeg = (id, field, value) => {
    setLegs((prev) =>
      prev.map((leg) => {
        if (leg.id !== id) return leg;
        if (field === "instrumentType") {
          const nextType = value;
          return {
            ...leg,
            instrumentType: nextType,
            optionType: nextType === "stock" ? "call" : leg.optionType,
            strike: nextType === "stock" ? params.S : leg.strike,
            T: nextType === "stock" ? params.T : leg.T
          };
        }
        return { ...leg, [field]: value };
      })
    );
  };

  const addLeg = () => {
    setLegs((prev) => [
      ...prev,
      {
        id: createLegId(),
        instrumentType: "option",
        optionType: "call",
        strike: params.K,
        qty: 1,
        T: params.T,
        multiplier: 100,
        entryPrice: 0
      }
    ]);
  };

  const removeLeg = (id) => {
    setLegs((prev) => prev.filter((leg) => leg.id !== id));
  };

  const chartData = useMemo(() => {
    const greekKeys = ["delta", "gamma", "theta", "vega", "vanna", "rho", "volga"];
    return greekKeys.map((key) => {
      const row = { greek: GREEK_LABELS[key] };
      rows.forEach((item, index) => {
        row[`leg${index}`] = item[key];
      });
      return row;
    });
  }, [rows]);

  return (
    <div className="space-y-3">
      <div className="flex flex-wrap items-center gap-2 rounded border border-slate-500 bg-slate-100 p-3 dark:bg-slate-800">
        <button
          type="button"
          className="inline-flex items-center gap-1 rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
          onClick={addLeg}
        >
          <Plus className="h-4 w-4" /> Add Leg
        </button>
        <button
          type="button"
          className="rounded border border-slate-600 bg-slate-200 px-2 py-1 text-sm hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600"
          onClick={onLoadStrategy}
        >
          Load Strategy Legs
        </button>
      </div>

      <div className="overflow-x-auto rounded border border-slate-500 bg-white dark:bg-slate-900">
        <table className="min-w-full border-collapse text-xs">
          <thead className="bg-slate-200 dark:bg-slate-700">
            <tr>
              <th className="border border-slate-400 p-1">Leg</th>
              <th className="border border-slate-400 p-1">Instrument</th>
              <th className="border border-slate-400 p-1">Type</th>
              <th className="border border-slate-400 p-1">Strike</th>
              <th className="border border-slate-400 p-1">T</th>
              <th className="border border-slate-400 p-1">Qty</th>
              <th className="border border-slate-400 p-1">Delta</th>
              <th className="border border-slate-400 p-1">Gamma</th>
              <th className="border border-slate-400 p-1">Theta</th>
              <th className="border border-slate-400 p-1">Vega</th>
              <th className="border border-slate-400 p-1">Vanna</th>
              <th className="border border-slate-400 p-1">Action</th>
            </tr>
          </thead>
          <tbody>
            {legs.map((leg, index) => {
              const row = rows[index];
              return (
                <tr key={leg.id}>
                  <td className="border border-slate-300 p-1">Leg {index + 1}</td>
                  <td className="border border-slate-300 p-1">
                    <select
                      className={numericInputBaseClass}
                      value={leg.instrumentType}
                      onChange={(evt) => updateLeg(leg.id, "instrumentType", evt.target.value)}
                    >
                      <option value="option">Option</option>
                      <option value="stock">Stock</option>
                    </select>
                  </td>
                  <td className="border border-slate-300 p-1">
                    <select
                      className={numericInputBaseClass}
                      value={leg.optionType}
                      disabled={leg.instrumentType === "stock"}
                      onChange={(evt) => updateLeg(leg.id, "optionType", evt.target.value)}
                    >
                      <option value="call">Call</option>
                      <option value="put">Put</option>
                    </select>
                  </td>
                  <td className="border border-slate-300 p-1">
                    <input
                      type="number"
                      className={numericInputBaseClass}
                      value={leg.strike}
                      disabled={leg.instrumentType === "stock"}
                      onChange={(evt) => updateLeg(leg.id, "strike", Math.max(1, parseNumeric(evt.target.value, leg.strike)))}
                    />
                  </td>
                  <td className="border border-slate-300 p-1">
                    <input
                      type="number"
                      className={numericInputBaseClass}
                      value={leg.T}
                      step="0.01"
                      disabled={leg.instrumentType === "stock"}
                      onChange={(evt) =>
                        updateLeg(leg.id, "T", Math.max(1 / YEAR_DAYS, parseNumeric(evt.target.value, leg.T)))
                      }
                    />
                  </td>
                  <td className="border border-slate-300 p-1">
                    <input
                      type="number"
                      className={numericInputBaseClass}
                      value={leg.qty}
                      onChange={(evt) => updateLeg(leg.id, "qty", parseNumeric(evt.target.value, leg.qty))}
                    />
                  </td>
                  <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(row?.delta || 0, 2)}</td>
                  <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(row?.gamma || 0, 4)}</td>
                  <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(row?.theta || 0, 2)}</td>
                  <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(row?.vega || 0, 2)}</td>
                  <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(row?.vanna || 0, 2)}</td>
                  <td className="border border-slate-300 p-1 text-center">
                    <button
                      type="button"
                      className="inline-flex items-center rounded border border-slate-500 p-1 hover:bg-slate-200 dark:hover:bg-slate-700"
                      onClick={() => removeLeg(leg.id)}
                      disabled={legs.length <= 1}
                    >
                      <Trash2 className="h-4 w-4" />
                    </button>
                  </td>
                </tr>
              );
            })}
            <tr className="bg-slate-100 font-semibold dark:bg-slate-700">
              <td className="border border-slate-300 p-1" colSpan={6}>
                TOTAL
              </td>
              <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(totals.delta, 2)}</td>
              <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(totals.gamma, 4)}</td>
              <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(totals.theta, 2)}</td>
              <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(totals.vega, 2)}</td>
              <td className="border border-slate-300 p-1 text-right font-mono">{formatSigned(totals.vanna, 2)}</td>
              <td className="border border-slate-300 p-1" />
            </tr>
          </tbody>
        </table>
      </div>

      <div className="rounded border border-slate-500 bg-slate-100 p-3 text-sm dark:bg-slate-800">
        Net exposure: delta {formatSigned(totals.delta, 2)}, gamma {formatSigned(totals.gamma, 4)}, theta {formatSigned(totals.theta, 2)}.
      </div>

      <div className="h-[420px] rounded border border-slate-500 bg-white p-2 dark:bg-slate-900">
        <Suspense fallback={<div className="grid h-full place-items-center text-sm text-slate-500">Loading decomposition chart...</div>}>
          <PortfolioGreekBarChart chartData={chartData} rows={rows} colors={CHART_COLORS} />
        </Suspense>
      </div>
    </div>
  );
}

function LearnTab() {
  return (
    <div className="space-y-3">
      <div className="rounded border border-slate-500 bg-slate-100 p-3 text-sm dark:bg-slate-800">
        This tab summarizes Greeks for study and trading interpretation. The dashboard uses European Black-Scholes-Merton with continuous dividend yield.
      </div>
      <div className="overflow-x-auto rounded border border-slate-500 bg-white dark:bg-slate-900">
        <table className="min-w-full border-collapse text-sm">
          <thead className="bg-slate-200 dark:bg-slate-700">
            <tr>
              <th className="border border-slate-400 p-2 text-left">Greek</th>
              <th className="border border-slate-400 p-2 text-left">Symbol</th>
              <th className="border border-slate-400 p-2 text-left">Measures</th>
              <th className="border border-slate-400 p-2 text-left">Intuition</th>
            </tr>
          </thead>
          <tbody>
            {LEARN_ROWS.map((row) => (
              <tr key={row.greek}>
                <td className="border border-slate-300 p-2">{row.greek}</td>
                <td className="border border-slate-300 p-2">{row.symbol || "-"}</td>
                <td className="border border-slate-300 p-2">{row.measure}</td>
                <td className="border border-slate-300 p-2">{row.intuition}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export default function OptionsGreeksDashboard() {
  const [activeTab, setActiveTab] = useState("surface");
  const [params, setParams] = useState(DEFAULT_PARAMS);
  const setParamsRaf = useRafSetter(setParams);

  const [selectedGreek, setSelectedGreek] = useState("delta");
  const [surfaceResolution, setSurfaceResolution] = useState(40);

  const [strategyPreset, setStrategyPreset] = useState("Long Call");
  const [strategyLegs, setStrategyLegs] = useState(() => instantiatePreset("Long Call", DEFAULT_PARAMS));
  const [strategyEvalDte, setStrategyEvalDte] = useState(Math.round(DEFAULT_PARAMS.T * YEAR_DAYS));

  const [portfolioLegs, setPortfolioLegs] = useState(() => instantiatePreset("Long Call", DEFAULT_PARAMS));

  const [surfaceGraphDiv, setSurfaceGraphDiv] = useState(null);
  const [heatmapGraphDiv, setHeatmapGraphDiv] = useState(null);

  const warnings = useMemo(() => {
    const notes = [];
    if (params.sigma > 1.2) notes.push("Very high volatility can amplify numerical noise.");
    if (params.T < 1 / YEAR_DAYS) notes.push("Very short expiry can cause steep Greek behavior.");
    if (params.S <= 0 || params.K <= 0) notes.push("Spot and strike must be positive.");
    return notes;
  }, [params]);

  const updateParam = useCallback(
    (field, value) => {
      setParamsRaf((prev) => {
        const next = { ...prev, [field]: value };
        return sanitizeParams(next);
      });
    },
    [setParamsRaf]
  );

  const resetParams = useCallback(() => {
    setParams(DEFAULT_PARAMS);
  }, []);

  const quickGreeks = useMemo(() => {
    const single = sanitizeParams(params);
    return {
      price: blackScholesPrice(single),
      delta: computeGreek(single, "delta"),
      gamma: computeGreek(single, "gamma"),
      theta: computeGreek(single, "theta"),
      vega: computeGreek(single, "vega"),
      rho: computeGreek(single, "rho"),
      vanna: computeGreek(single, "vanna"),
      volga: computeGreek(single, "volga")
    };
  }, [params]);

  const surfaceData = useMemo(() => {
    return computeSurfaceData(params, selectedGreek, surfaceResolution);
  }, [params, selectedGreek, surfaceResolution]);

  const loadPreset = useCallback(() => {
    const legs = instantiatePreset(strategyPreset, params);
    setStrategyLegs(legs);
    const dte = Math.max(1, ...legs.map((leg) => Math.round(leg.T * YEAR_DAYS)));
    setStrategyEvalDte(dte);
  }, [params, strategyPreset]);

  const markStrategyEntries = useCallback(() => {
    setStrategyLegs((prev) =>
      prev.map((leg) => ({
        ...leg,
        entryPrice: computeLegEntryPrice(leg, params)
      }))
    );
  }, [params]);

  const strategyMaxDte = useMemo(
    () => Math.max(1, ...strategyLegs.map((leg) => Math.round(leg.T * YEAR_DAYS))),
    [strategyLegs]
  );

  useEffect(() => {
    setStrategyEvalDte((prev) => clamp(prev, 0, strategyMaxDte));
  }, [strategyMaxDte]);

  const strategyHeatmap = useMemo(
    () => computePnLHeatmap(strategyLegs, params, 80, strategyEvalDte),
    [strategyLegs, params, strategyEvalDte]
  );

  const breakevens = useMemo(() => findBreakevens(strategyLegs, params), [strategyLegs, params]);

  const portfolioData = useMemo(() => computePortfolioRows(portfolioLegs, params), [portfolioLegs, params]);

  const themeClass = "dark min-h-screen bg-zinc-900 text-zinc-100";

  return (
    <div className={themeClass}>
      <div className="mx-auto max-w-[1700px] p-3">
        <header className="mb-3 flex flex-wrap items-center justify-between rounded border border-slate-500 bg-slate-100 p-3 dark:bg-slate-800">
          <div className="space-y-1">
            <h1 className="text-lg font-bold uppercase tracking-wide">OGVDash</h1>
            <p className="text-xs text-slate-600 dark:text-slate-300">Real-time interactive options Greeks surface and heatmap tooling software.</p>
            <div className="flex flex-wrap items-center gap-3 text-xs text-slate-700 dark:text-slate-300">
              <span>Made by David Park</span>
              <a
                href="https://parkdavid.com"
                target="_blank"
                rel="noreferrer"
                className="inline-flex items-center gap-1 hover:underline"
              >
                <Globe className="h-3.5 w-3.5" />
                parkdavid.com
              </a>
              <a
                href="https://github.com/zqoza123"
                target="_blank"
                rel="noreferrer"
                className="inline-flex items-center gap-1 hover:underline"
              >
                <Github className="h-3.5 w-3.5" />
                GitHub
              </a>
              <a
                href="https://www.linkedin.com/in/jpddd/"
                target="_blank"
                rel="noreferrer"
                className="inline-flex items-center gap-1 hover:underline"
              >
                <Linkedin className="h-3.5 w-3.5" />
                LinkedIn
              </a>
            </div>
          </div>
        </header>

        {warnings.length > 0 && (
          <div className="mb-3 rounded border border-amber-600 bg-amber-100 p-2 text-sm text-amber-900">
            {warnings.map((note) => (
              <div key={note}>{note}</div>
            ))}
          </div>
        )}

        <div className="mb-3 flex flex-wrap gap-1 rounded border border-slate-500 bg-slate-100 p-2 dark:bg-slate-800">
          {TAB_LIST.map((tab) => (
            <button
              key={tab.id}
              type="button"
              className={`rounded border px-3 py-1 text-sm ${
                activeTab === tab.id
                  ? "border-slate-700 bg-slate-700 text-white"
                  : "border-slate-500 bg-slate-200 text-slate-900 hover:bg-slate-300 dark:bg-slate-700 dark:text-slate-100 dark:hover:bg-slate-600"
              }`}
              onClick={() => setActiveTab(tab.id)}
            >
              {tab.label}
            </button>
          ))}
        </div>

        <div className="grid gap-3 xl:grid-cols-[340px_1fr]">
          <ParameterPanel params={params} onChangeParam={updateParam} onReset={resetParams} quickGreeks={quickGreeks} />

          <main className="rounded border border-slate-500 bg-slate-100 p-3 dark:bg-slate-800">
            {activeTab === "surface" && (
              <GreekSurfaceView
                params={params}
                selectedGreek={selectedGreek}
                setSelectedGreek={setSelectedGreek}
                resolution={surfaceResolution}
                setResolution={setSurfaceResolution}
                surfaceData={surfaceData}
                onExportCsv={() => exportSurfaceCsv(surfaceData, selectedGreek)}
                onExportImage={() => exportPlotImage(surfaceGraphDiv, `surface-${selectedGreek}.png`)}
                onGraphReady={(graphDiv) => setSurfaceGraphDiv(graphDiv)}
              />
            )}

            {activeTab === "theta" && <ThetaDecayAnimator params={params} />}

            {activeTab === "strategy" && (
              <StrategyBuilder
                params={params}
                legs={strategyLegs}
                setLegs={setStrategyLegs}
                preset={strategyPreset}
                setPreset={setStrategyPreset}
                evalDte={strategyEvalDte}
                setEvalDte={setStrategyEvalDte}
                heatmap={strategyHeatmap}
                breakevens={breakevens}
                onLoadPreset={loadPreset}
                onMarkEntries={markStrategyEntries}
                onExportCsv={() => exportHeatmapCsv(strategyHeatmap)}
                onExportImage={() => exportPlotImage(heatmapGraphDiv, "strategy-heatmap.png")}
                onHeatmapGraphReady={(graphDiv) => setHeatmapGraphDiv(graphDiv)}
              />
            )}

            {activeTab === "portfolio" && (
              <PortfolioGreeksTable
                params={params}
                legs={portfolioLegs}
                setLegs={setPortfolioLegs}
                rows={portfolioData.rows}
                totals={portfolioData.totals}
                onLoadStrategy={() =>
                  setPortfolioLegs(strategyLegs.map((leg) => ({ ...leg, id: createLegId() })))
                }
              />
            )}

            {activeTab === "learn" && <LearnTab />}
          </main>
        </div>
      </div>
    </div>
  );
}
