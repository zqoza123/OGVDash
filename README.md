# OGVDash

A simple, browser-based options Greeks dashboard.

I built this to make options behavior easier to *see* instead of just reading numbers in a table. You can move inputs like spot, vol, rates, and time, then instantly see how Greeks and strategy P&L shift. 
This project was originally inspired by Andy Nguyen's 2008 QuantNet post challenging users to replicate a euro-options calculator. But I've decided to modernize while adding a bit of my own twists based on my own experience with options and greeks. : )

## What it does
- Shows Greek surfaces across strike and expiry
- Runs a theta decay view so you can watch time value fade
- Builds common options strategies and shows P&L heatmaps
- Breaks portfolio Greeks down by leg

## Math (quick version)
Under the hood it uses Black-Scholes-Merton for European options.

Core pieces are:
- `d1 = [ln(S/K) + (r - q + 0.5*sigma^2)*T] / (sigma*sqrt(T))`
- `d2 = d1 - sigma*sqrt(T)`

From those, option price + Greeks (delta, gamma, theta, vega, rho, plus higher-order Greeks) are calculated directly. For edge cases like tiny `T` or tiny `sigma`, it falls back to safer finite-difference style calculations so the app stays stable.

## Run it
```bash
npm install
npm run dev
```

Then open the local URL Vite prints (usually `http://localhost:5173`).
