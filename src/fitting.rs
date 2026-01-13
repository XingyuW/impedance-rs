use nalgebra::{DMatrix, DVector, Dyn, Owned};
use levenberg_marquardt::LeastSquaresProblem;
use std::f64::consts::PI;
use num_complex::Complex64;
use crate::circuits::{CircuitNode, Impedance};
use crate::elements::{ElementType, Constraint};

struct GuessState {
    rs: f64,
    r_ct: f64,
    q_cpe: f64,
    sigma_gw: f64,
    alpha_gw: f64,
    r_count: usize,
}

/// Guesses initial parameters for the circuit based on impedance data.
pub fn guess_parameters(
    node: &CircuitNode,
    frequencies: &[f64],
    z_real: &[f64],
    z_imag: &[f64],
    _phase_deg: &[f64],
) -> Vec<f64> {
    let total_params = node.count_total_params();
    let mut params = vec![0.0; total_params];
    
    // 1. Identify HF and LF indices
    let (idx_hf, _) = frequencies.iter().enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();
    let (idx_lf, _) = frequencies.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();

    // 2. Rs Guess: Real part at High Freq
    let rs_guess = z_real[idx_hf].abs();

    // 3. Smooth Imaginary part to find dips
    let n = z_imag.len();
    let mut im_smooth = z_imag.to_vec();
    if n >= 3 {
        for i in 1..n-1 {
            im_smooth[i] = (z_imag[i-1] + z_imag[i] + z_imag[i+1]) / 3.0;
        }
    }

    // 4. Find Dip (Local Maxima in Im(Z) - closest to zero)
    // We look for peaks in the smoothed imaginary data
    let mut peaks = Vec::new();
    for i in 1..n-1 {
        if im_smooth[i] > im_smooth[i-1] && im_smooth[i] > im_smooth[i+1] {
            peaks.push(i);
        }
    }

    let rct_guess;
    let idx_dip;

    if !peaks.is_empty() {
        // Find the highest peak value (closest to 0)
        let idx_in_peaks = peaks.iter()
            .max_by(|&&a, &&b| im_smooth[a].partial_cmp(&im_smooth[b]).unwrap())
            .unwrap();
        idx_dip = *idx_in_peaks;
        rct_guess = (z_real[idx_dip] - rs_guess).abs();
    } else {
        // No dip, use arc top (minimum of imaginary part, i.e., most negative)
        let (idx_peak, _) = z_imag.iter().enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        let z_peak_real = z_real[idx_peak];
        rct_guess = (z_peak_real - rs_guess).abs() * 2.0;
        idx_dip = if idx_hf < idx_lf { n - 1 } else { 0 }; // End of arc search range
    }

    let final_rct_guess = if rct_guess < 1e-9 { 1000.0 } else { rct_guess };

    // 5. Q (CPE) Guess
    // Search for peak in the arc region (between HF and Dip)
    let start = idx_hf.min(idx_dip);
    let end = idx_hf.max(idx_dip);
    let range_start = if start == end { 0 } else { start };
    let range_end = if start == end { n } else { end + 1 }; // Exclusive

    let slice_im = &z_imag[range_start..range_end];
    let (idx_local, _) = slice_im.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap_or((0, &0.0));
    
    let idx_top = range_start + idx_local;
    let w_peak = 2.0 * PI * frequencies[idx_top];
    let q_guess = 1.0 / (w_peak * final_rct_guess);

    // 6. Sigma & Alpha (Warburg) Guess
    // Use LF point
    let z_lf = Complex64::new(z_real[idx_lf], z_imag[idx_lf]);
    let phase_lf_rad = z_lf.im.atan2(z_lf.re);
    let phase_lf_deg = phase_lf_rad.to_degrees();
    
    let mut alpha_guess = phase_lf_deg.abs() / 90.0;
    if alpha_guess > 0.6 { alpha_guess = 0.5; }
    if alpha_guess < 0.2 { alpha_guess = 0.25; }

    let w_low = 2.0 * PI * frequencies[idx_lf];
    let z_mag_low = z_lf.norm();
    // |Z_Gw| = sigma * w^-alpha => sigma = |Z| * w^alpha
    let sigma_guess = z_mag_low * w_low.powf(alpha_guess);

    let mut state = GuessState {
        rs: rs_guess,
        r_ct: final_rct_guess,
        q_cpe: q_guess,
        sigma_gw: sigma_guess,
        alpha_gw: alpha_guess,
        r_count: 0,
    };

    fill_guesses(node, &mut params, &mut state);
    
    params
}

fn fill_guesses(
    node: &CircuitNode,
    params: &mut [f64],
    state: &mut GuessState,
) {
    match node {
        CircuitNode::Series(nodes) | CircuitNode::Parallel(nodes) => {
            for n in nodes {
                fill_guesses(n, params, state);
            }
        }
        CircuitNode::Element(etype, idx, _) => {
            match etype {
                ElementType::R => {
                    if state.r_count == 0 {
                        params[*idx] = state.rs;
                    } else {
                        params[*idx] = state.r_ct;
                    }
                    state.r_count += 1;
                }
                ElementType::C => {
                    params[*idx] = 1e-6; 
                }
                ElementType::L => {
                    params[*idx] = 1e-6;
                }
                ElementType::W => {
                    params[*idx] = state.sigma_gw;
                }
                ElementType::CPE => {
                    params[*idx] = state.q_cpe; 
                    params[*idx + 1] = 0.9; // Alpha usually close to 1 for CPE
                }
                ElementType::Wo | ElementType::Ws => {
                    params[*idx] = state.r_ct; 
                    params[*idx + 1] = 1.0; 
                }
                ElementType::La => {
                    params[*idx] = 1e-6; 
                    params[*idx + 1] = 1.0; 
                }
                ElementType::Gw => {
                    params[*idx] = state.sigma_gw;
                    params[*idx + 1] = state.alpha_gw;
                }
                ElementType::G | ElementType::Gs => {
                    params[*idx] = state.r_ct; 
                    params[*idx + 1] = 1.0; 
                    if let ElementType::Gs = etype {
                        params[*idx + 2] = 0.5;
                    }
                }
                ElementType::K | ElementType::Zarc | ElementType::TLMQ => {
                    params[*idx] = state.r_ct; 
                    params[*idx + 1] = 1.0; 
                    if let ElementType::TLMQ = etype {
                        params[*idx + 1] = state.q_cpe;
                    }
                    if etype.param_count() > 2 {
                        params[*idx + 2] = 0.8;
                    }
                }
                ElementType::T => {
                    params[*idx] = state.r_ct; 
                    params[*idx + 1] = state.r_ct; 
                    params[*idx + 2] = 1.0; 
                    params[*idx + 3] = 1.0; 
                }
            }
        }
    }
}

/// Struct for Levenberg-Marquardt optimization.
pub struct ImpedanceFitter {
    pub circuit: CircuitNode,
    pub frequencies: Vec<f64>,
    pub z_real_data: Vec<f64>,
    pub z_imag_data: Vec<f64>,
    pub params: DVector<f64>,
    pub constraints: Vec<Constraint>,
}

// Helper functions for parameter transformation
pub fn transform_forward(physical: f64, constraint: Constraint) -> f64 {
    match constraint {
        Constraint::Positive => {
            if physical <= 0.0 { -23.0 } else { physical.ln() } // Handle non-positive input gracefully
        },
        Constraint::ZeroOne => {
            if physical <= 0.0 { -23.0 } else if physical >= 1.0 { 23.0 } else { (physical / (1.0 - physical)).ln() }
        },
        Constraint::None => physical,
    }
}

pub fn transform_backward(internal: f64, constraint: Constraint) -> f64 {
    match constraint {
        Constraint::Positive => internal.exp(),
        Constraint::ZeroOne => 1.0 / (1.0 + (-internal).exp()),
        Constraint::None => internal,
    }
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for ImpedanceFitter {
    type ParameterStorage = Owned<f64, Dyn>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;

    fn set_params(&mut self, params: &DVector<f64>) {
        self.params = params.clone();
    }

    fn params(&self) -> DVector<f64> {
        self.params.clone()
    }

    fn residuals(&self) -> Option<DVector<f64>> {
        let n_points = self.frequencies.len();
        let mut residuals = DVector::zeros(2 * n_points);
        
        // Transform internal params to physical params
        let physical_params: Vec<f64> = self.params.iter().zip(self.constraints.iter())
            .map(|(&p, &c)| transform_backward(p, c))
            .collect();

        for i in 0..n_points {
            let omega = 2.0 * PI * self.frequencies[i];
            let z_model = self.circuit.calculate(omega, &physical_params);
            
            residuals[i] = z_model.re - self.z_real_data[i];
            residuals[i + n_points] = z_model.im - self.z_imag_data[i];
        }
        Some(residuals)
    }

    fn jacobian(&self) -> Option<DMatrix<f64>> {
        let n_points = self.frequencies.len();
        let n_params = self.params.len();
        let mut jacobian = DMatrix::zeros(2 * n_points, n_params);
        let epsilon = 1e-8;

        // Current physical parameters
        let physical_params: Vec<f64> = self.params.iter().zip(self.constraints.iter())
            .map(|(&p, &c)| transform_backward(p, c))
            .collect();
        
        let mut base_residuals = DVector::zeros(2 * n_points);
        for i in 0..n_points {
            let omega = 2.0 * PI * self.frequencies[i];
            let z = self.circuit.calculate(omega, &physical_params);
            base_residuals[i] = z.re - self.z_real_data[i];
            base_residuals[i + n_points] = z.im - self.z_imag_data[i];
        }

        let mut p_curr_internal = self.params.clone();

        for j in 0..n_params {
            let original_val = p_curr_internal[j];
            p_curr_internal[j] += epsilon;
            
            // Calculate perturbed physical params
            let perturbed_physical_val = transform_backward(p_curr_internal[j], self.constraints[j]);
            
            let mut p_perturbed_physical = physical_params.clone();
            p_perturbed_physical[j] = perturbed_physical_val;

            for i in 0..n_points {
                let omega = 2.0 * PI * self.frequencies[i];
                let z = self.circuit.calculate(omega, &p_perturbed_physical);
                
                let res_re = z.re - self.z_real_data[i];
                let res_im = z.im - self.z_imag_data[i];

                jacobian[(i, j)] = (res_re - base_residuals[i]) / epsilon;
                jacobian[(i + n_points, j)] = (res_im - base_residuals[i + n_points]) / epsilon;
            }

            p_curr_internal[j] = original_val; 
        }

        Some(jacobian)
    }
}

/// Solves the Linear Kramers-Kronig validation.
pub fn lin_kk_solver(
    frequencies: &[f64],
    z_real: &[f64],
    z_imag: &[f64],
    c: f64,
    max_m: usize,
) -> (usize, f64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = frequencies.len();
    
    let mut best_m = 0;
    let mut best_mu = 1.0;
    let mut best_z_fit_re = vec![0.0; n];
    let mut best_z_fit_im = vec![0.0; n];
    let mut best_res_real = vec![0.0; n];
    let mut best_res_imag = vec![0.0; n];

    for m in 3..=max_m {
        let min_f = frequencies.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_f = frequencies.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let tau_min = 1.0 / (2.0 * PI * max_f);
        let tau_max = 1.0 / (2.0 * PI * min_f);
        
        let mut taus = Vec::with_capacity(m);
        if m > 1 {
            for k in 0..m {
                let log_tau = tau_min.ln() + (k as f64 / (m as f64 - 1.0)) * (tau_max / tau_min).ln();
                taus.push(log_tau.exp());
            }
        } else {
            taus.push(tau_min);
        }

        let n_params = m + 1;
        let mut a_mat = DMatrix::<f64>::zeros(2 * n, n_params);
        let mut b_vec = DVector::<f64>::zeros(2 * n);

        for i in 0..n {
            let w = 2.0 * PI * frequencies[i];
            
            a_mat[(i, 0)] = 1.0; 
            a_mat[(i + n, 0)] = 0.0; 
            
            b_vec[i] = z_real[i];
            b_vec[i + n] = z_imag[i];

            for k in 0..m {
                let tau = taus[k];
                let denom = 1.0 + (w * tau).powi(2);
                
                a_mat[(i, k + 1)] = 1.0 / denom;
                a_mat[(i + n, k + 1)] = -(w * tau) / denom;
            }
        }

        let svd = a_mat.clone().svd(true, true);
        let x = svd.solve(&b_vec, 1e-9).unwrap_or_else(|_| DVector::zeros(n_params));
        
        let r_k = x.rows(1, m);
        let mut sum_pos = 0.0;
        let mut sum_neg = 0.0;
        for &val in r_k.iter() {
            if val >= 0.0 {
                sum_pos += val.abs();
            } else {
                sum_neg += val.abs();
            }
        }
        
        let mu = if sum_neg > 1e-12 {
            1.0 - sum_pos / sum_neg
        } else {
            -1.0 
        };
        
        let ax = &a_mat * &x;
        let mut z_fit_re = Vec::with_capacity(n);
        let mut z_fit_im = Vec::with_capacity(n);
        for i in 0..n {
            z_fit_re.push(ax[i]);
            z_fit_im.push(ax[i + n]);
        }
        
        best_m = m;
        best_mu = mu;
        best_res_real = z_real.iter().zip(z_fit_re.iter()).map(|(a, b)| a - b).collect();
        best_res_imag = z_imag.iter().zip(z_fit_im.iter()).map(|(a, b)| a - b).collect();
        best_z_fit_re = z_fit_re;
        best_z_fit_im = z_fit_im;

        if mu < c {
            break;
        }
    }

    (best_m, best_mu, best_z_fit_re, best_z_fit_im, best_res_real, best_res_imag)
}
