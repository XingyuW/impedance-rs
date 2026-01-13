use pyo3::prelude::*;
use nalgebra::DVector;
use levenberg_marquardt::LevenbergMarquardt;
use std::f64::consts::PI;

mod elements;
mod circuits;
mod fitting;

use circuits::{parse_circuit_string, Impedance};
use fitting::{guess_parameters, ImpedanceFitter, lin_kk_solver, transform_forward, transform_backward};

#[pyfunction]
fn fit_circuit(
    circuit_str: String,
    frequencies: Vec<f64>,
    z_real: Vec<f64>,
    z_imag: Vec<f64>,
    phase_deg: Vec<f64>,
) -> PyResult<(Vec<f64>, Vec<String>, Vec<String>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    // 1. Parse
    let circuit = parse_circuit_string(&circuit_str)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;

    // 2. Guess
    let initial_params_physical = guess_parameters(&circuit, &frequencies, &z_real, &z_imag, &phase_deg);
    let constraints = circuit.get_constraints();

    // Transform to internal parameters
    let initial_params_internal: Vec<f64> = initial_params_physical.iter().zip(constraints.iter())
        .map(|(&p, &c)| transform_forward(p, c))
        .collect();

    // 3. Optimize
    let fitter = ImpedanceFitter {
        circuit: circuit.clone(),
        frequencies: frequencies.clone(),
        z_real_data: z_real,
        z_imag_data: z_imag,
        params: DVector::from_vec(initial_params_internal),
        constraints: constraints.clone(),
    };

    let (result, _report) = LevenbergMarquardt::new().minimize(fitter);
    
    let fitted_params_internal = result.params.as_slice();
    
    // Transform back to physical parameters
    let fitted_params_physical: Vec<f64> = fitted_params_internal.iter().zip(constraints.iter())
        .map(|(&p, &c)| transform_backward(p, c))
        .collect();

    // Get names and units
    let param_names = circuit.get_param_names();
    let param_units = circuit.get_param_units();

    // 4. Simulate & 5. Derived Metrics
    let mut fitted_real = Vec::with_capacity(frequencies.len());
    let mut fitted_imag = Vec::with_capacity(frequencies.len());
    let mut fitted_mag = Vec::with_capacity(frequencies.len());
    let mut fitted_phase = Vec::with_capacity(frequencies.len());

    for &f in &frequencies {
        let omega = 2.0 * PI * f;
        let z = circuit.calculate(omega, &fitted_params_physical);
        
        fitted_real.push(z.re);
        fitted_imag.push(z.im);
        fitted_mag.push(z.norm());
        fitted_phase.push(z.im.atan2(z.re).to_degrees());
    }

    Ok((
        fitted_params_physical,
        param_names,
        param_units,
        fitted_real,
        fitted_imag,
        fitted_mag,
        fitted_phase,
    ))
}

#[pyfunction]
fn ignore_below_x(frequencies: Vec<f64>, z_real: Vec<f64>, z_imag: Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut f_out = Vec::new();
    let mut z_real_out = Vec::new();
    let mut z_imag_out = Vec::new();

    for i in 0..frequencies.len() {
        if z_imag[i] < 0.0 {
            f_out.push(frequencies[i]);
            z_real_out.push(z_real[i]);
            z_imag_out.push(z_imag[i]);
        }
    }
    (f_out, z_real_out, z_imag_out)
}

#[pyfunction]
#[pyo3(signature = (frequencies, z_real, z_imag, freqmin=None, freqmax=None))]
fn crop_frequencies(
    frequencies: Vec<f64>,
    z_real: Vec<f64>,
    z_imag: Vec<f64>,
    freqmin: Option<f64>,
    freqmax: Option<f64>,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut f_out = Vec::new();
    let mut z_real_out = Vec::new();
    let mut z_imag_out = Vec::new();

    let min_f = freqmin.unwrap_or(0.0);
    let max_f = freqmax.unwrap_or(f64::INFINITY);

    for i in 0..frequencies.len() {
        if frequencies[i] >= min_f && frequencies[i] <= max_f {
            f_out.push(frequencies[i]);
            z_real_out.push(z_real[i]);
            z_imag_out.push(z_imag[i]);
        }
    }
    (f_out, z_real_out, z_imag_out)
}

#[pyfunction]
fn lin_kk(
    frequencies: Vec<f64>,
    z_real: Vec<f64>,
    z_imag: Vec<f64>,
    c: f64,
    max_m: usize,
) -> PyResult<(usize, f64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    let result = lin_kk_solver(&frequencies, &z_real, &z_imag, c, max_m);
    Ok(result)
}

#[pymodule]
fn impedance_rs(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fit_circuit, m)?)?;
    m.add_function(wrap_pyfunction!(ignore_below_x, m)?)?;
    m.add_function(wrap_pyfunction!(crop_frequencies, m)?)?;
    m.add_function(wrap_pyfunction!(lin_kk, m)?)?;
    Ok(())
}
