import { describe, expect, it, vi } from "vitest";

vi.mock("react-plotly.js", () => ({
  default: () => null
}));

vi.mock("plotly.js-dist-min", () => ({
  default: {
    toImage: vi.fn()
  }
}));

import {
  blackScholesPrice,
  computeGreek,
  computePortfolioRows,
  computeStrategyPnL,
  instantiatePreset,
  DEFAULT_PARAMS
} from "./OptionsGreeksDashboard";

describe("Black-Scholes math", () => {
  it("matches baseline call price", () => {
    const call = blackScholesPrice({
      S: 100,
      K: 100,
      T: 1,
      sigma: 0.3,
      r: 0.05,
      q: 0,
      optionType: "call"
    });

    expect(call).toBeCloseTo(14.23, 2);
  });

  it("satisfies put-call parity", () => {
    const params = {
      S: 100,
      K: 100,
      T: 1,
      sigma: 0.3,
      r: 0.05,
      q: 0
    };

    const call = blackScholesPrice({ ...params, optionType: "call" });
    const put = blackScholesPrice({ ...params, optionType: "put" });
    const lhs = put;
    const rhs = call - params.S * Math.exp(-params.q * params.T) + params.K * Math.exp(-params.r * params.T);

    expect(lhs).toBeCloseTo(rhs, 4);
  });

  it("has sane first/second order greek signs around ATM", () => {
    const params = {
      S: 100,
      K: 100,
      T: 1,
      sigma: 0.3,
      r: 0.05,
      q: 0,
      optionType: "call"
    };

    expect(computeGreek(params, "delta")).toBeGreaterThan(0.5);
    expect(computeGreek(params, "gamma")).toBeGreaterThan(0);
    expect(computeGreek(params, "vega")).toBeGreaterThan(0);
    expect(computeGreek(params, "theta")).toBeLessThan(0);
  });
});

describe("Strategy P&L", () => {
  it("matches long call expiry payoff", () => {
    const leg = {
      id: "a",
      instrumentType: "option",
      optionType: "call",
      strike: 100,
      qty: 1,
      entryPrice: 10,
      T: 0.5,
      multiplier: 100
    };

    const pnl = computeStrategyPnL(
      [leg],
      {
        ...DEFAULT_PARAMS,
        S: 120,
        sigma: 0.3,
        r: 0.05,
        q: 0
      },
      0
    );

    expect(pnl).toBeCloseTo(1000, 5);
  });

  it("shows iron condor center better than far wings at expiry", () => {
    const base = {
      ...DEFAULT_PARAMS,
      S: 100,
      K: 100,
      T: 0.5,
      sigma: 0.25,
      optionType: "call"
    };

    const legs = instantiatePreset("Iron Condor", base);
    const center = computeStrategyPnL(legs, { ...base, S: 100 }, 0);
    const leftWing = computeStrategyPnL(legs, { ...base, S: 60 }, 0);
    const rightWing = computeStrategyPnL(legs, { ...base, S: 140 }, 0);

    expect(center).toBeGreaterThan(leftWing);
    expect(center).toBeGreaterThan(rightWing);
  });
});

describe("Portfolio aggregation", () => {
  it("is linear across legs", () => {
    const leg = {
      id: "l1",
      instrumentType: "option",
      optionType: "call",
      strike: 100,
      qty: 1,
      entryPrice: 0,
      T: 0.5,
      multiplier: 100
    };

    const { totals: one } = computePortfolioRows([leg], DEFAULT_PARAMS);
    const { totals: two } = computePortfolioRows(
      [leg, { ...leg, id: "l2" }],
      DEFAULT_PARAMS
    );

    expect(two.delta).toBeCloseTo(one.delta * 2, 8);
    expect(two.gamma).toBeCloseTo(one.gamma * 2, 8);
    expect(two.theta).toBeCloseTo(one.theta * 2, 8);
    expect(two.vega).toBeCloseTo(one.vega * 2, 8);
  });
});
