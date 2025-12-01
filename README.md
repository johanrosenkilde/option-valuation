# Employee Stock Option Valuation: Black-Scholes and Binomial Models

This script computes the value of employee stock options (ESOs), which are a
form of compensation granting employees the right to purchase company stock at
a fixed price (the strike price) after a vesting period. ESOs differ from
standard exchange-traded options in several ways: they are typically not
transferable, may be subject to vesting and forfeiture, and are often
exercised early due to job changes or liquidity needs.

Two main approaches are implemented:

1. Black-Scholes Model:
   - A closed-form solution for pricing European-style options (exercisable
     only at expiry).
   - Assumes constant volatility, no early exercise, and no vesting/forfeiture
     features.
   - Useful as a baseline, but does not capture the complexities of ESOs.

2. Binomial Model (with Hull-White adjustments):
   - A flexible, tree-based approach that can model American-style options
     (exercisable at any time), vesting schedules, early exercise, and
     employee exit/forfeiture.
   - The model simulates the evolution of the stock price over discrete time
     steps, allowing for the possibility of early exercise and the impact of
     employee turnover (exit_rate) and vesting (Vs).
   - The Hull-White extension further accommodates the probability of employee
     exit and the forced exercise or forfeiture of unvested options.

Key Features for Employee Stock Options:
   - Vesting periods: Options become exercisable only after a specified period
     (Vs).
   - Forfeiture: If the employee leaves before vesting, unvested options are
     lost.
   - Early exercise: Employees may exercise options before expiry, assumed to
     occur when the stock price reaches a certain multiple of the strike price
     (exercise_factor).
   - Multiple scenarios: The model supports weighted summing over any number
     of scenarios (e.g., early/late exit), each with its own probability. This
     reflects common practice when valuating the stocks in 409A.

This script expects input parameters to be provided in a YAML file, specified
as a command-line argument.

Example YAML format:

```yaml
s0: 1.95                # Current stock price (float)
K: 1.10                 # Strike price (float)
exit_rate: 0.15         # Annual probability the employee exits the company,
                          forfeiting further vests (float, e.g. 0.15 for 15%)
exercise_factor: 3      # Factor over strike price at which the employee exercises the option. (float)
scenarios:
  late_exit:
    T: 10                   # Time to expiration in years (float)
    Vs: [1,2,3,4]           # List of vesting periods in years (list of floats)
    r: 0.039                # Risk-free interest rate (float)
    sigma: 0.72             # Volatility (float)
    prob: 0.25              # Probability of scenario (float)
  early_exit:
    T: 1.55                 # This scenario the company exits early, so remaining 
                            # options vest immediately.
    Vs: [1, 1.5, 1.5, 1.5]
    r: 0.045
    sigma: 0.16
    prob: 0.75
```

All fields are required. Any number of scenarios can be provided under
`scenarios`, with arbitrary names. The sum of all scenario 'prob' fields must
be 1.0.
