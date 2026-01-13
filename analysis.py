import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    import impedance_rs  # This will be the compiled Rust module
except ImportError:
    print("Error: impedance_rs module not found. Please build and install it using 'maturin develop --release'.")
    sys.exit(1)

# ==========================================
# 1. DATA LOADING (Kept in Python)
# ==========================================
def load_chi_data(filepath):
    """Reads CHI760e CSV format, skipping headers robustly."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        start_idx = 0
        for i, line in enumerate(lines):
            if "Freq/Hz" in line:
                start_idx = i + 2 
                break
                
        df = pd.read_csv(filepath, skiprows=start_idx, header=None)
        df = df.dropna(how='all')
        
        # Ensure numpy arrays are float type
        freq = df[0].values.astype(float)
        z_real = df[1].values.astype(float)
        z_imag = df[2].values.astype(float)
        phase = df[4].values.astype(float)
        
        # Construct complex impedance array
        # Explicitly cast to complex128 to avoid type ambiguity
        # Use numpy arrays directly to avoid pandas ExtensionArray issues
        z_real_np = np.array(z_real, dtype=np.float64)
        z_imag_np = np.array(z_imag, dtype=np.float64)
        # z_meas = z_real_np + 1j * z_imag_np
        return freq, z_real_np, z_imag_np, phase
        
    except Exception as e:
        print(f"Error parsing CSV: {e}")
        return None, None, None, None
# ==========================================
# 2. MAIN EXECUTION
# ==========================================
def main():
    filepath = "./data/example.csv"
    print(f"Loading {filepath}...")
    freq, z_real_np, z_imag_np, phase = load_chi_data(filepath)
    
    if freq is None or z_real_np is None or z_imag_np is None or phase is None:
        print("Failed to load data.")
        return

    # circuit_str = "R0-p(CPE1,R1-W1)"
    circuit_str = "R0-p(CPE1,R1-Gw1)"
    print(f"Fitting circuit: {circuit_str}")

    try:
        # Call the Rust fitting function
        # fit_circuit(circuit_str, frequencies, z_real, z_imag, phase_deg)
        # Note: Rust expects Vec<f64>, numpy arrays convert automatically or via .tolist()
        # PyO3 usually handles numpy array to Vec<f64> conversion if typed correctly, 
        # but explicit tolist() is safer if not using numpy-specific PyO3 features.
        
        fitted_params, param_names, param_units, fit_real, fit_imag, fit_mag, fit_phase = impedance_rs.fit_circuit( # type: ignore
            circuit_str, 
            freq.tolist(), 
            z_real_np.tolist(), 
            z_imag_np.tolist(), 
            phase.tolist()
        )

        print("\nFitting Successful!")
        print("Fitted Parameters:")
        for name, value, unit in zip(param_names, fitted_params, param_units):
            print(f"  {name}: {value:.4e} [{unit}]")

        # ==========================================
        # 3. PLOTTING
        # ==========================================
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Nyquist Plot
        ax1.plot(z_real_np, -z_imag_np, 'o', label='Data')
        ax1.plot(fit_real, [-y for y in fit_imag], '-', label='Fit', linewidth=2)
        ax1.set_xlabel("Z' (Ohm)")
        ax1.set_ylabel("-Z'' (Ohm)")
        ax1.set_title("Nyquist Plot")
        ax1.legend()
        ax1.grid(True)
        ax1.set_aspect('equal', 'box')

        # Bode Plot (Magnitude)
        ax2_mag = ax2
        ax2_phase = ax2.twinx()

        ax2_mag.loglog(freq, np.sqrt(z_real_np**2 + z_imag_np**2), 'o', label='|Z| Data')
        ax2_mag.loglog(freq, fit_mag, '-', label='|Z| Fit')
        
        ax2_phase.semilogx(freq, phase, 's', color='orange', label='Phase Data', alpha=0.5)
        ax2_phase.semilogx(freq, fit_phase, '--', color='red', label='Phase Fit')

        ax2_mag.set_xlabel("Frequency (Hz)")
        ax2_mag.set_ylabel("|Z| (Ohm)")
        ax2_phase.set_ylabel("Phase (deg)")
        ax2.set_title("Bode Plot")
        
        # Combine legends
        lines1, labels1 = ax2_mag.get_legend_handles_labels()
        lines2, labels2 = ax2_phase.get_legend_handles_labels()
        ax2_mag.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        ax2_mag.grid(True)

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"An error occurred during fitting: {e}")


if __name__ == "__main__":
    main()
