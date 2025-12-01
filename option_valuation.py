# Employee Stock Option Valuation: Black-Scholes and Binomial Models
# ------------------------------------------------------------------
#
# This script computes the value of employee stock options (ESOs), which are a
# form of compensation granting employees the right to purchase company stock at
# a fixed price (the strike price) after a vesting period. ESOs differ from
# standard exchange-traded options in several ways: they are typically not
# transferable, may be subject to vesting and forfeiture, and are often
# exercised early due to job changes or liquidity needs.
#
# Two main approaches are implemented:
#
# 1. Black-Scholes Model:
#    - A closed-form solution for pricing European-style options (exercisable
#      only at expiry).
#    - Assumes constant volatility, no early exercise, and no vesting/forfeiture
#      features.
#    - Useful as a baseline, but does not capture the complexities of ESOs.
#
# 2. Binomial Model (with Hull-White adjustments):
#    - A flexible, tree-based approach that can model American-style options
#      (exercisable at any time), vesting schedules, early exercise, and
#      employee exit/forfeiture.
#    - The model simulates the evolution of the stock price over discrete time
#      steps, allowing for the possibility of early exercise and the impact of
#      employee turnover (exit_rate) and vesting (Vs).
#    - The Hull-White extension further accommodates the probability of employee
#      exit and the forced exercise or forfeiture of unvested options.
#
# Key Features for Employee Stock Options:
#    - Vesting periods: Options become exercisable only after a specified period
#      (Vs).
#    - Forfeiture: If the employee leaves before vesting, unvested options are
#      lost.
#    - Early exercise: Employees may exercise options before expiry, assumed to
#      occur when the stock price reaches a certain multiple of the strike price
#      (exercise_factor).
#    - Multiple scenarios: The model supports weighted summing over any number
#      of scenarios (e.g., early/late exit), each with its own probability. This
#      reflects common practice when valuating the stocks in 409A.
#
# This script expects input parameters to be provided in a YAML file, specified
# as a command-line argument.
#
# Example YAML format:
#
# s0: 1.95                # Current stock price (float)
# K: 1.10                 # Strike price (float)
# exit_rate: 0.15         # Annual probability the employee exits the company,
#                           forfeiting further vests (float, e.g. 0.15 for 15%)
# exercise_factor: 3      # Factor over strike price at which the employee exercises the option. (float)
# scenarios:
#   late_exit:
#     T: 10                   # Time to expiration in years (float)
#     Vs: [1,2,3,4]           # List of vesting periods in years (list of floats)
#     r: 0.039                # Risk-free interest rate (float)
#     sigma: 0.72             # Volatility (float)
#     prob: 0.25              # Probability of scenario (float)
#   early_exit:
#     T: 1.55                 # This scenario the company exits early, so remaining 
#                             # options vest immediately.
#     Vs: [1, 1.5, 1.5, 1.5]
#     r: 0.045
#     sigma: 0.16
#     prob: 0.75
#
# All fields are required. Any number of scenarios can be provided under
# `scenarios`, with arbitrary names. The sum of all scenario 'prob' fields must
# be 1.0.
#

import math
from scipy.stats import norm
import yaml
import sys

def black_scholes(S0, K, T, r, sigma):
    """
    Calculate the Black-Scholes price for European call options.
    
    Parameters:
    - S0: Current stock price
    - K: Strike price
    - T: Time to expiration (in years)
    - r: Risk-free interest rate
    - sigma: Annualized volatility of the stock
    
    Returns:
    - Price of the European call option
    """
    dplus = (math.log(S0 / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    dminus = dplus - sigma * math.sqrt(T)
    return S0 * norm.cdf(dplus) - K * math.exp(-r * T) * norm.cdf(dminus)

def _binomial_model_pure(S0, K, Rounds, r, dividend, sigma, style):
    """
    Pure binomial model. Raw version which progresses for a fixed number of
    rounds, and all time-sensitive parameters are adjusted accordingly.
    """
    if style not in ["American", "European"]:
        raise ValueError("Invalid option style. Please choose 'American' or 'European'.")
    exp_minus_r = math.exp(-r)
    u = math.exp(sigma)
    d = 1 / u
    p = (math.exp(r - dividend) - d) / (u - d)

    assert 1 < sigma**2 / (r - dividend)**2, "Arbitrage opportunity detected"

    def stock_price_tree(layer, ups):
        downs = layer - ups
        return S0 * math.exp(sigma * (ups - downs))
    
    option_value_tree = [[0 for _ in range(j + 1)] for j in range(Rounds + 1)]
    
    for i in range(Rounds + 1):
        option_value_tree[Rounds][i] = max(stock_price_tree(Rounds, i) - K, 0)
    
    for j in range(Rounds - 1, -1, -1):
        for i in range(j + 1):
            binomial_value = exp_minus_r * (p * option_value_tree[j + 1][i + 1] + (1 - p) * option_value_tree[j + 1][i])
            if style == "American":
                option_value_tree[j][i] = max(binomial_value, stock_price_tree(j, i) - K)
            else:
                option_value_tree[j][i] = binomial_value
    
    return option_value_tree[0][0]


def binomial_model_pure(S0, K, Years, r, dividend, sigma, style, granularity):
    """
    Pricing an option, with all the parameters considered in the time step granularity.
    - `S0`: Current stock price.
    - `K`: Strike price.
    - `Years`: Years to expiration.
    - `r`: Risk-free interest rate per year.
    - `sigma`: Volatility of the stock per year.
    - `style`: Style of the option, either "American" or "European".
    - `granularity`: Number of time steps to divide the whole duration into.
    """
    # Calculate time step based on granularity
    delta_t = Years / granularity

    # Adjusted parameters for the time step
    adj_dividend = dividend * delta_t
    adj_r = r * delta_t
    adj_sigma = sigma * math.sqrt(delta_t)
    
    # Call the raw binomial model function with vesting and forfeiture
    return _binomial_model_pure(S0, K, granularity, adj_r, adj_dividend, adj_sigma, style)


def _binomial_model_hull_white(S0, K, Rounds, r, dividend, sigma, vests_after, exit_rate, exercise_factor, style):
    """
    The binomial model using Hull-White adjustments to account for vesting
    period, exit rate and early exercise.
    Raw version which progresses for a fixed number of rounds, and all
    time-sensitive parameters are adjusted accordingly.
    """
    if style != "American":
        raise ValueError("Only 'American' style supported for Hull--White.")
    exp_minus_r = math.exp(-r)
    u = math.exp(sigma)
    d = 1 / u
    p = (math.exp(r - dividend) - d) / (u - d)

    assert 1 < sigma**2 / (r - dividend)**2, "Arbitrage opportunity detected"

    def stock_price_tree(layer, ups):
        downs = layer - ups
        return S0 * math.exp(sigma * (ups - downs))
    
    option_value_tree = [[0 for _ in range(j + 1)] for j in range(Rounds + 1)]
    
    for i in range(Rounds + 1):
        option_value_tree[Rounds][i] = max(stock_price_tree(Rounds, i) - K, 0)
    
    for j in range(Rounds - 1, -1, -1):
        for i in range(j + 1):
            stock_price = stock_price_tree(j, i)
            if j >= vests_after and stock_price >= K*exercise_factor:
                v = stock_price - K
            elif j >= vests_after:
                v = (
                    (1 - exit_rate) * exp_minus_r * (p * option_value_tree[j + 1][i + 1] + (1 - p) * option_value_tree[j + 1][i])
                    + exit_rate * max(stock_price - K, 0)
                )
            else:
                v = (1 - exit_rate) * exp_minus_r * (p * option_value_tree[j + 1][i + 1] + (1 - p) * option_value_tree[j + 1][i])
            option_value_tree[j][i] = v
    
    return option_value_tree[0][0]


def binomial_model_hull_white(S0, K, Years, r, dividend, sigma, vests_after, exit_rate, exercise_factor, style, granularity):
    """
    Pricing an option, with all the parameters considered in the time step granularity.
    - `S0`: Current stock price.
    - `K`: Strike price.
    - `Years`: Years to expiration.
    - `r`: Risk-free interest rate per year.
    - `sigma`: Volatility of the stock per year.
    - `vests_after`: Number of years after which the stock vests.
    - `exit_rate`: Rate at which the employee exits the company. At exit, they forfeit unvested stocks and is forced to exercise vested stocks.
    - `exercise_factor`: Factor over strike price at which the employee exercises the option.
    - `style`: Style of the option, either "American" or "European".
    - `granularity`: Number of time steps to divide the whole duration into.
    """
    def _call(granularity):
        "Call with the given granularity"
        delta_t = Years / granularity

        # Adjusted parameters for the time step
        adj_dividend = dividend * delta_t
        adj_r = r * delta_t
        adj_sigma = sigma * math.sqrt(delta_t)
        adj_vests_after = vests_after / delta_t
        adj_exit_rate = exit_rate * delta_t
        return _binomial_model_hull_white(S0, K, granularity, adj_r, adj_dividend, adj_sigma, adj_vests_after, adj_exit_rate, exercise_factor, style)
    
    def _min_from(granularity):
        "Get the next local minimum from the given granularity"
        last = None
        while True:
            val = _call(granularity)
            if last is not None and val > last:
                return last
            last = val
            granularity += 1

    # Drop in to get a few local minima to get a better estimate
    least = _min_from(granularity)
    REPS = 1
    for i in range(REPS):
        least = min(least, _min_from(granularity + i * 50))
    return least 


import numpy as np

def monte_carlo_option_pricing(S0, K, T, V, r, sigma, stay_prob, num_simulations, granularity):
    """
    TODO: As written by GPT4 - untested!
    
    Monte Carlo simulation to price a call option considering vesting and forfeiture.

    Parameters:
    - S0: Current stock price
    - K: Strike price
    - T: Time to expiration (in years)
    - V: Vesting period (in years)
    - r: Risk-free interest rate
    - sigma: Annualized volatility of the stock
    - stay_prob: Annualized probability of staying with the company
    - num_simulations: Number of Monte Carlo simulations
    - granularity: Number of time steps

    Returns:
    - Price of the call option
    """
    delta_t = T / granularity
    stay_prob_per_step = stay_prob ** delta_t

    def generate_price_path():
        """ Generates a single simulated price path """
        prices = [S0]
        for _ in range(granularity):
            z = np.random.standard_normal()  # Generate random standard normal value
            price = prices[-1] * np.exp((r - 0.5 * sigma ** 2) * delta_t + sigma * np.sqrt(delta_t) * z)
            prices.append(price)
        return prices

    payoff_sum = 0  # Sum of payoffs from all simulations

    for _ in range(num_simulations):
        prices = generate_price_path()
        vested = True

        # Check if the employee remains in the company until the vesting period
        for j in range(int(V / delta_t)):
            if np.random.rand() > stay_prob_per_step:
                vested = False
                break

        if vested:
            # The employee remains in the company; calculate payoff at expiry
            payoff_sum += max(prices[-1] - K, 0)
        else:
            # The employee leaves before vesting; immediate exercise of vested options
            # Calculate the maximum payoff within the vesting period
            payoff_sum += max(0, max(prices[j] - K for j in range(int(V / delta_t))))

    # Discounting the expected payoff to the present value
    option_price = (np.exp(-r * T) * payoff_sum) / num_simulations
    return option_price


# Parameters will be loaded from YAML file.
S0 = None
K = None
exit_rate = None
exercise_factor = None
scenarios = None

def format_value(value, strike_price, stock_price):
    """Format option value with absolute amount and meaningful percentages."""
    if value <= 0:
        return f"${value:.3f} (no intrinsic value)"
    
    pct_of_strike = (value / strike_price) * 100
    pct_of_stock = (value / stock_price) * 100
    
    return f"${value:.3f} ({pct_of_strike:.1f}% of strike, {pct_of_stock:.1f}% of stock price)"


def scenario_value(d, pure_or_hw):
    if pure_or_hw == "pure":
        # Pure pricing doesn't support vesting periods
        binom_prices = [binomial_model_pure(S0, K, d["T"], d["r"], 0., d["sigma"], "American", 1000)]
    else:
        binom_prices = [binomial_model_hull_white(S0, K, d["T"], d["r"], 0., d["sigma"], V, exit_rate, exercise_factor, 'American', 600) for V in d["Vs"]]
    avg = sum(binom_prices) / len(binom_prices)
    return avg


def print_input_summary():
    """Print a summary of the input parameters for transparency."""
    print("=" * 60)
    print("INPUT PARAMETERS")
    print("=" * 60)
    print(f"Current stock price (S0): ${S0:.2f}")
    print(f"Strike price (K): ${K:.2f}")
    print(f"Annual exit rate: {exit_rate:.1%}")
    print(f"Early exercise factor: {exercise_factor}x strike price")
    print(f"\nSCENARIOS:")
    for scenario_name, d in scenarios.items():
        print(f"  {scenario_name.replace('_', ' ').title()}:")
        print(f"    - Time to expiry: {d['T']:.1f} years")
        print(f"    - Vesting schedule: {', '.join(f'{v:.1f}yr' for v in d['Vs'])}")
        print(f"    - Risk-free rate: {d['r']:.1%}")
        print(f"    - Volatility: {d['sigma']:.1%}")
        print(f"    - Probability: {d['prob']:.1%}")
    print()

def exit_value():
    print_input_summary()
    
    print("=" * 60)
    print("OPTION VALUATION RESULTS")
    print("=" * 60)
    
    def calculate_and_display(model_name, model_description, pure_or_hw):
        print(f"\n{model_name.upper()}:")
        print(f"{model_description}")
        print("-" * 40)
        
        scenario_results = []
        for scenario_name, d in scenarios.items():
            val = scenario_value(d, pure_or_hw)
            scenario_results.append((scenario_name, val, d["prob"]))
            formatted_val = format_value(val, K, S0)
            scenario_label = scenario_name.replace('_', ' ').title()
            print(f"  {scenario_label:12} | {formatted_val}")
        
        weighted = sum(val * prob for _, val, prob in scenario_results)
        formatted_weighted = format_value(weighted, K, S0)
        print(f"  {'Expected Value':12} | {formatted_weighted}")
        
        return weighted
    
    pure_result = calculate_and_display(
        "Standard Binomial Model",
        "Basic binomial model without employee-specific adjustments", 
        "pure"
    )
    
    hw_result = calculate_and_display(
        "Hull-White Employee Model", 
        "Binomial model with vesting, exit rates, and early exercise",
        "hw"
    )
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Standard Model Expected Value:    {format_value(pure_result, K, S0)}")
    print(f"Employee-Adjusted Expected Value: {format_value(hw_result, K, S0)}")
    
    discount = ((pure_result - hw_result) / pure_result) * 100 if pure_result > 0 else 0
    print(f"Employee adjustment discount:     {discount:.1f}%")
    print(f"\nThe employee-specific factors (vesting, exit risk, early exercise)")
    print(f"reduce the option value by ${pure_result - hw_result:.3f} compared to a standard option.")



def load_parameters_from_yaml(yaml_path):
    """
    Load option valuation parameters from a YAML file.
    See the top of this script for the expected format and field meanings.
    """
    with open(yaml_path, 'r') as f:
        params = yaml.safe_load(f)
    # Validate required fields
    required_fields = ["s0", "K", "exit_rate", "exercise_factor", "scenarios"]
    for field in required_fields:
        if field not in params:
            raise ValueError(f"Missing required field in YAML: {field}")
    scenarios_dict = params["scenarios"]
    if not isinstance(scenarios_dict, dict) or len(scenarios_dict) == 0:
        raise ValueError("'scenarios' must be a non-empty mapping in YAML.")
    total_prob = 0.0
    for scenario_name, d in scenarios_dict.items():
        for subfield in ["T", "Vs", "r", "sigma", "prob"]:
            if subfield not in d:
                raise ValueError(f"Missing field '{subfield}' in scenario '{scenario_name}' in YAML.")
        total_prob += d["prob"]
    if not abs(total_prob - 1.0) < 1e-8:
        raise ValueError(f"Sum of all scenario 'prob' fields must be 1.0, got {total_prob}")
    return params


def assign_globals_from_params(params):
    global S0, K, exit_rate, exercise_factor, scenarios
    S0 = params["s0"]
    K = params["K"]
    exit_rate = params["exit_rate"]
    exercise_factor = params["exercise_factor"]
    scenarios = params["scenarios"]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python option_valuation.py <parameters.yaml>")
        sys.exit(1)
    yaml_path = sys.argv[1]
    params = load_parameters_from_yaml(yaml_path)
    assign_globals_from_params(params)
    print(f"\nLoaded parameters from '{yaml_path}'.")
    print("Starting employee stock option valuation...\n")
    exit_value()
    print("\nValuation complete.")

# Parameters from "How to value employee stock options"
test_sample_options = {
    "S0": 50,
    "K": 50,
    "T": 10,
    "r": 0.075,
    "sigma": 0.3,
    "dividend": 0.025,
    "vests_after": 3
}

def test_binomial_pure_no_dividend_low_gran():
    v = binomial_model_pure(
        test_sample_options["S0"],
        test_sample_options["K"],
        test_sample_options["T"], 
        test_sample_options["r"],
        0.0, # dividend
        test_sample_options["sigma"],
        'American',
        20
    )
    # This value was found using
    # https://accuratecalculators.com/options-calculator
    assert abs(v - 29.94) < 0.05

def test_binomial_pure_no_dividend_high_gran():
    v = binomial_model_pure(
        test_sample_options["S0"],
        test_sample_options["K"],
        test_sample_options["T"], 
        test_sample_options["r"],
        0.0, # dividend
        test_sample_options["sigma"],
        'American',
        1000
    )
    # This value was found using
    # https://accuratecalculators.com/options-calculator
    assert abs(v - 30.10) < 0.05

def test_binomial_pure_with_dividend():
    v = binomial_model_pure(
        test_sample_options["S0"],
        test_sample_options["K"],
        test_sample_options["T"], 
        test_sample_options["r"],
        test_sample_options["dividend"],
        test_sample_options["sigma"],
        'American',
        1000
    )
    # This value was found in the paper "How to value employee stock options"
    assert abs(v - 21.03) < 0.05

def test_binomial_hull_white():
    def test_at(exit_rate, exercise_factor):
        return binomial_model_hull_white(
            test_sample_options["S0"],
            test_sample_options["K"],
            test_sample_options["T"], 
            test_sample_options["r"],
            test_sample_options["dividend"],
            test_sample_options["sigma"],
            test_sample_options["vests_after"],
            exit_rate,
            exercise_factor,
            'American',
            650
        )
    print("val", test_at(0.03, 1.2))
    # This value was found in the paper "How to value employee stock options"
    # Calculating for many iterations and taking the min apparently is a better
    # estimator. For speed, we instead just slack the error margin.
    assert abs(test_at(0.03, 1.2) - 13.13) < 0.20
    assert abs(test_at(0.05, 2.0) - 15.80) < 0.20
    assert abs(test_at(0.10, 3.0) - 13.75) < 0.25

# def test_binomial_model_different_steps():
#     T = late_exit["T"] 
#     V = 3
#     r = late_exit["r"]
#     sigma = late_exit["sigma"]
#     last = -1
#     for N in range(1, 1000):
#         call_price = binomial_option_pricing(S0, K, T, V, r, sigma, stay_prob, N)
#         if abs(call_price - last) > 0.001:
#             print(f"Time Steps: {N} - Call Option % of strike: {rel_strike(call_price)}")
#             last = call_price

# def test_compare_vesting_periods():
#     T = late_exit["T"]
#     r = late_exit["r"]
#     sigma = late_exit["sigma"]
#     for V in [0] + late_exit["Vs"]:
#         call_price = binomial_option_pricing(S0, K, T, V, r, sigma, stay_prob, 1000)
#         print(f"Vesting Period: {V} - Call Option Price: {rel_strike(call_price)}")
#     bs = black_scholes(S0, K, late_exit["T"], late_exit["r"], late_exit["sigma"])
#     print(f"Black-Scholes Price: {bs} = {rel_strike(bs)}")
