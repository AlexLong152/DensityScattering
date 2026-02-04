import numpy as np

# 1N results
results_1N = {
    "3H": {
        "FT_SpV": [1.481, 0.046],
        "FT_SmV": [-0.047, 0.005],
        "E": [-0.995, 0.030],
        "FL_SpV": [1.481, 0.046],
        "FL_SmV": [-0.035, 0.005],
        "L": [-1.047, 0.034],
    },
    "3He": {
        "FT_SpV": [-0.044, 0.005],
        "FT_SmV": [1.475, 0.045],
        "E": [1.746, 0.053],
        "FL_SpV": [-0.032, 0.005],
        "FL_SmV": [1.475, 0.045],
        "L": [-1.921, 0.060],
    },
    "6Li": {
        "FT_SpV": [0.468, 0.014],
        "FT_SmV": [0.468, 0.014],
        "E": [0.360, 0.011],
        "FL_SpV": [0.470, 0.014],
        "FL_SmV": [0.470, 0.014],
        "L": [-1.979, 0.060],
    },
}

# 2N results
results_2N = {
    "3He": {
        "FT": {"O2": [-30.65, 0.67], "O4": [-30.04, 0.70]},
        "FL": {"O2": [-24.19, 0.52], "O4": [-24.30, 0.53]},
        "E": {"O2": [-4.086, 0.089], "O4": [-4.005, 0.093]},
        "L": {"O2": [-3.225, 0.069], "O4": [-3.240, 0.071]},
    },
    "3H": {
        "FT": {"O2": [-30.89, 0.68], "O4": [-30.25, 0.70]},
        "FL": {"O2": [-24.40, 0.53], "O4": [-24.50, 0.54]},
        "E": {"O2": [-4.119, 0.090], "O4": [-4.034, 0.094]},
        "L": {"O2": [-3.253, 0.070], "O4": [-3.267, 0.072]},
    },
    "6Li": {
        "FT": {"O2": [-16.84, 0.35], "O4": [-16.54, 0.37]},
        "FL": {"O2": [-13.11, 0.27], "O4": [-13.10, 0.28]},
        "E": {"O2": [-2.298, 0.048], "O4": [-2.256, 0.050]},
        "L": {"O2": [-1.789, 0.037], "O4": [-1.787, 0.039]},
    },
}


def add_with_error(val1, val2):
    """Add two [value, error] arrays with error propagation."""
    central = val1[0] + val2[0]
    error = np.sqrt(val1[1] ** 2 + val2[1] ** 2)
    return [central, error]


# Compute total E0+ and L0+ (1N + 2N)
results_total = {}
for nucleus in ["3H", "3He", "6Li"]:
    results_total[nucleus] = {}
    for order in ["O2", "O4"]:
        results_total[nucleus][order] = {
            "E": add_with_error(
                results_1N[nucleus]["E"], results_2N[nucleus]["E"][order]
            ),
            "L": add_with_error(
                results_1N[nucleus]["L"], results_2N[nucleus]["L"][order]
            ),
        }

# Print results
print("Total E0+ and L0+ (1N + 2N):")
print("-" * 60)
for nucleus in ["3H", "3He", "6Li"]:
    print(f"\n{nucleus}:")
    for order in ["O4"]:
        E = results_total[nucleus][order]["E"]
        L = results_total[nucleus][order]["L"]
        print(
            f"  {order}: E0+ = {E[0]:.3f}({E[1] * 1000:.0f})  L0+ = {L[0]:.3f}({L[1] * 1000:.0f})"
        )
