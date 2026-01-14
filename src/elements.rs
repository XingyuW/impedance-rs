use num_complex::Complex64;
use std::f64::consts::PI;

/// Represents the constraint on a parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Constraint {
    /// Parameter must be strictly positive (p > 0).
    /// Transformation: p = exp(x)
    Positive,
    /// Parameter must be between 0 and 1 (0 < p < 1).
    /// Transformation: p = 1 / (1 + exp(-x))
    ZeroOne,
    /// No constraint.
    #[allow(dead_code)]
    None,
}

/// Represents the type of an electrical element in the equivalent circuit.
///
/// Each variant corresponds to a specific electrochemical component with a defined impedance behavior.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementType {
    /// Resistor (R)
    ///
    /// Represents an ideal resistor.
    /// Impedance: Z = R
    /// Parameters: [R] (Ohms)
    R,

    /// Capacitor (C)
    ///
    /// Represents an ideal capacitor.
    /// Impedance: Z = 1 / (j * omega * C)
    /// Parameters: [C] (Farads)
    C,

    /// Inductor (L)
    ///
    /// Represents an ideal inductor.
    /// Impedance: Z = j * omega * L
    /// Parameters: [L] (Henries)
    L,

    /// Warburg Element (W) - Semi-infinite
    ///
    /// Represents semi-infinite linear diffusion.
    /// Impedance: Z = sigma * (1 - j) / sqrt(omega)
    /// Parameters: [sigma] (Ohm * sec^-1/2)
    W,

    /// Constant Phase Element (CPE)
    ///
    /// Represents a non-ideal capacitor or distributed time constant.
    /// Impedance: Z = 1 / (Q * (j * omega)^alpha)
    /// Parameters: [Q, alpha] (Ohm^-1 * sec^alpha, dimensionless 0 < alpha <= 1)
    Cpe,

    /// Warburg Open (Wo) - Finite-space
    ///
    /// Represents finite-length diffusion with a reflective boundary (open circuit).
    /// Impedance: Z = Z0 * coth(sqrt(j * omega * tau)) / sqrt(j * omega * tau)
    /// Parameters: [Z0, tau] (Ohms, seconds)
    Wo,

    /// Warburg Short (Ws) - Finite-length
    ///
    /// Represents finite-length diffusion with a transmissive boundary (short circuit).
    /// Impedance: Z = Z0 * tanh(sqrt(j * omega * tau)) / sqrt(j * omega * tau)
    /// Parameters: [Z0, tau] (Ohms, seconds)
    Ws,

    /// Modified Inductance (La)
    ///
    /// Represents an inductance with a non-integer power dependence.
    /// Impedance: Z = L * (j * omega)^alpha
    /// Parameters: [L, alpha] (Henries * sec^(alpha-1), dimensionless)
    La,

    /// Generalized Warburg Element (Gw)
    ///
    /// Represents a generalized Warburg element.
    /// Impedance: Z = sigma * (j * omega)^(-alpha)
    /// Parameters: [sigma, alpha] (Ohm * sec^alpha, dimensionless)
    Gw,

    /// Gerischer Element (G)
    ///
    /// Represents a chemical reaction coupled with diffusion (semi-infinite).
    /// Impedance: Z = R_G / sqrt(1 + j * omega * t_G)
    /// Parameters: [R_G, t_G] (Ohms, seconds)
    G,

    /// Finite-length Gerischer Element (Gs)
    ///
    /// Represents a chemical reaction coupled with diffusion in a finite length.
    /// Impedance: Z = R_G / (sqrt(1 + j * omega * t_G) * tanh(phi * sqrt(1 + j * omega * t_G)))
    /// Parameters: [R_G, t_G, phi] (Ohms, seconds, dimensionless)
    Gs,

    /// K Element (K)
    ///
    /// An RC element used in Lin-KK validation models.
    /// Impedance: Z = R / (1 + j * omega * tau_k)
    /// Parameters: [R, tau_k] (Ohms, seconds)
    K,

    /// Zarc Element (Zarc)
    ///
    /// Represents a Cole-Cole relaxation (RQ element rewritten).
    /// Impedance: Z = R / (1 + (j * omega * tau_k)^gamma)
    /// Parameters: [R, tau_k, gamma] (Ohms, seconds, dimensionless)
    Zarc,

    /// Transmission Line Model (TLMQ)
    ///
    /// Simplified transmission line model for porous electrodes.
    /// Impedance: Z = sqrt(Rion * Zs) * coth(sqrt(Rion / Zs)), where Zs = 1 / (Qs * (j * omega)^gamma)
    /// Parameters: [Rion, Qs, gamma] (Ohms, Ohm^-1 * sec^gamma, dimensionless)
    Tlmq,

    /// Porous Electrode Model (T)
    ///
    /// Macrohomogeneous porous electrode model (Paasch et al.).
    /// Impedance: Z = A * coth(beta)/beta + B / (beta * sinh(beta)), where beta = sqrt(a + j * omega * b)
    /// Parameters: [A, B, a, b] (Ohm, Ohm, dimensionless, seconds)
    T,
}

impl ElementType {
    /// Returns the number of parameters required for this element type.
    pub fn param_count(&self) -> usize {
        match self {
            ElementType::R => 1,
            ElementType::C => 1,
            ElementType::L => 1,
            ElementType::W => 1,
            ElementType::Cpe => 2,
            ElementType::Wo => 2,
            ElementType::Ws => 2,
            ElementType::La => 2,
            ElementType::Gw => 2,
            ElementType::G => 2,
            ElementType::Gs => 3,
            ElementType::K => 2,
            ElementType::Zarc => 3,
            ElementType::Tlmq => 3,
            ElementType::T => 4,
        }
    }

    /// Returns the constraints for the parameters of this element.
    pub fn constraints(&self) -> Vec<Constraint> {
        match self {
            ElementType::R => vec![Constraint::Positive],
            ElementType::C => vec![Constraint::Positive],
            ElementType::L => vec![Constraint::Positive],
            ElementType::W => vec![Constraint::Positive],
            ElementType::Cpe => vec![Constraint::Positive, Constraint::ZeroOne],
            ElementType::Wo => vec![Constraint::Positive, Constraint::Positive],
            ElementType::Ws => vec![Constraint::Positive, Constraint::Positive],
            ElementType::La => vec![Constraint::Positive, Constraint::Positive], // alpha usually > 0
            ElementType::Gw => vec![Constraint::Positive, Constraint::ZeroOne],
            ElementType::G => vec![Constraint::Positive, Constraint::Positive],
            ElementType::Gs => vec![
                Constraint::Positive,
                Constraint::Positive,
                Constraint::Positive,
            ],
            ElementType::K => vec![Constraint::Positive, Constraint::Positive],
            ElementType::Zarc => vec![
                Constraint::Positive,
                Constraint::Positive,
                Constraint::ZeroOne,
            ],
            ElementType::Tlmq => vec![
                Constraint::Positive,
                Constraint::Positive,
                Constraint::ZeroOne,
            ],
            ElementType::T => vec![
                Constraint::Positive,
                Constraint::Positive,
                Constraint::Positive,
                Constraint::Positive,
            ],
        }
    }

    /// Returns the names of the parameters for this element.
    pub fn param_names(&self) -> Vec<&'static str> {
        match self {
            ElementType::R => vec!["R"],
            ElementType::C => vec!["C"],
            ElementType::L => vec!["L"],
            ElementType::W => vec!["sigma"],
            ElementType::Cpe => vec!["Q", "alpha"],
            ElementType::Wo => vec!["Z0", "tau"],
            ElementType::Ws => vec!["Z0", "tau"],
            ElementType::La => vec!["L", "alpha"],
            ElementType::Gw => vec!["sigma", "alpha"],
            ElementType::G => vec!["R_G", "t_G"],
            ElementType::Gs => vec!["R_G", "t_G", "phi"],
            ElementType::K => vec!["R", "tau_k"],
            ElementType::Zarc => vec!["R", "tau_k", "gamma"],
            ElementType::Tlmq => vec!["Rion", "Qs", "gamma"],
            ElementType::T => vec!["A", "B", "a", "b"],
        }
    }

    /// Returns the units of the parameters for this element.
    pub fn param_units(&self) -> Vec<&'static str> {
        match self {
            ElementType::R => vec!["Ohm"],
            ElementType::C => vec!["F"],
            ElementType::L => vec!["H"],
            ElementType::W => vec!["Ohm s^-1/2"],
            ElementType::Cpe => vec!["Ohm^-1 s^alpha", ""],
            ElementType::Wo => vec!["Ohm", "s"],
            ElementType::Ws => vec!["Ohm", "s"],
            ElementType::La => vec!["H s^(alpha-1)", ""],
            ElementType::Gw => vec!["Ohm s^alpha", ""],
            ElementType::G => vec!["Ohm", "s"],
            ElementType::Gs => vec!["Ohm", "s", ""],
            ElementType::K => vec!["Ohm", "s"],
            ElementType::Zarc => vec!["Ohm", "s", ""],
            ElementType::Tlmq => vec!["Ohm", "Ohm^-1 s^gamma", ""],
            ElementType::T => vec!["Ohm", "Ohm", "", "s"],
        }
    }

    /// Calculates the complex impedance of the element at a given angular frequency.
    ///
    /// # Arguments
    /// * `omega` - Angular frequency (2 * PI * f) in rad/s.
    /// * `params` - Slice of parameters for this element.
    ///
    /// # Returns
    /// * `Complex64` - The complex impedance Z.
    pub fn calculate(&self, omega: f64, params: &[f64]) -> Complex64 {
        let p = params;
        match self {
            ElementType::R => Complex64::new(p[0], 0.0),
            ElementType::C => {
                if omega > 1e-9 {
                    Complex64::new(0.0, -1.0 / (omega * p[0]))
                } else {
                    Complex64::new(1e12, 0.0) // DC block
                }
            }
            ElementType::L => Complex64::new(0.0, omega * p[0]),
            ElementType::W => {
                // Warburg: Z = sigma * (1-j) / sqrt(omega)
                if omega > 1e-9 {
                    let s = p[0] / omega.sqrt();
                    Complex64::new(s, -s)
                } else {
                    Complex64::new(1e6, -1e6) // Limit at DC
                }
            }
            ElementType::Cpe => {
                // Z = 1 / (Q * (j*omega)^alpha)
                let q = p[0];
                let alpha = p[1];
                if omega > 1e-9 {
                    let term = omega.powf(alpha);
                    let phase = PI / 2.0 * alpha;
                    let denom_real = q * term * phase.cos();
                    let denom_imag = q * term * phase.sin();
                    let denom = Complex64::new(denom_real, denom_imag);
                    1.0 / denom
                } else {
                    Complex64::new(1e12, 0.0)
                }
            }
            ElementType::Wo => {
                // Z = Z0 * coth(sqrt(j*w*tau)) / sqrt(j*w*tau)
                let z0 = p[0];
                let tau = p[1];
                if omega > 1e-9 {
                    let sqrt_jwt = (Complex64::i() * omega * tau).sqrt();
                    let tanh = sqrt_jwt.tanh();
                    if tanh.norm_sqr() > 1e-16 {
                        z0 / (sqrt_jwt * tanh)
                    } else {
                        Complex64::new(1e12, 0.0)
                    }
                } else {
                    Complex64::new(1e12, 0.0)
                }
            }
            ElementType::Ws => {
                // Z = Z0 * tanh(sqrt(j*w*tau)) / sqrt(j*w*tau)
                let z0 = p[0];
                let tau = p[1];
                if omega > 1e-9 {
                    let sqrt_jwt = (Complex64::i() * omega * tau).sqrt();
                    if sqrt_jwt.norm_sqr() > 1e-16 {
                        z0 * sqrt_jwt.tanh() / sqrt_jwt
                    } else {
                        Complex64::new(z0, 0.0)
                    }
                } else {
                    Complex64::new(z0, 0.0)
                }
            }
            ElementType::La => {
                // Z = L * (j*w)^alpha
                let l = p[0];
                let alpha = p[1];
                let jw_alpha = (Complex64::i() * omega).powf(alpha);
                Complex64::new(l, 0.0) * jw_alpha
            }
            ElementType::Gw => {
                // Z = sigma * (j * omega)^(-alpha)
                let sigma = p[0];
                let alpha = p[1];
                if omega > 1e-9 {
                    let jw = Complex64::i() * omega;
                    let term = jw.powf(-alpha);
                    Complex64::new(sigma, 0.0) * term
                } else {
                    Complex64::new(1e12, 0.0)
                }
            }
            ElementType::G => {
                // Z = R_G / sqrt(1 + j*w*t_G)
                let r_g = p[0];
                let t_g = p[1];
                let denom = (Complex64::new(1.0, 0.0) + Complex64::i() * omega * t_g).sqrt();
                Complex64::new(r_g, 0.0) / denom
            }
            ElementType::Gs => {
                // Z = R_G / (sqrt(1 + j*w*t_G) * tanh(phi * sqrt(1 + j*w*t_G)))
                let r_g = p[0];
                let t_g = p[1];
                let phi = p[2];
                let sqrt_term = (Complex64::new(1.0, 0.0) + Complex64::i() * omega * t_g).sqrt();
                let tanh_term = (Complex64::new(phi, 0.0) * sqrt_term).tanh();
                if tanh_term.norm_sqr() > 1e-16 {
                    Complex64::new(r_g, 0.0) / (sqrt_term * tanh_term)
                } else {
                    Complex64::new(1e12, 0.0)
                }
            }
            ElementType::K => {
                // Z = R / (1 + j*w*tau_k)
                let r = p[0];
                let tau_k = p[1];
                let denom = Complex64::new(1.0, 0.0) + Complex64::i() * omega * tau_k;
                Complex64::new(r, 0.0) / denom
            }
            ElementType::Zarc => {
                // Z = R / (1 + (j*w*tau_k)^gamma)
                let r = p[0];
                let tau_k = p[1];
                let gamma = p[2];
                let term = (Complex64::i() * omega * tau_k).powf(gamma);
                let denom = Complex64::new(1.0, 0.0) + term;
                Complex64::new(r, 0.0) / denom
            }
            ElementType::Tlmq => {
                // Z = sqrt(Rion * Zs) * coth(sqrt(Rion / Zs))
                // Zs = 1 / (Qs * (j*w)^gamma)
                let r_ion = p[0];
                let qs = p[1];
                let gamma = p[2];

                if omega > 1e-9 {
                    let jw_gamma = (Complex64::i() * omega).powf(gamma);
                    let z_s = 1.0 / (Complex64::new(qs, 0.0) * jw_gamma);

                    let r_ion_c = Complex64::new(r_ion, 0.0);
                    let sqrt_prod = (r_ion_c * z_s).sqrt();
                    let sqrt_quot = (r_ion_c / z_s).sqrt();
                    let tanh_quot = sqrt_quot.tanh();

                    if tanh_quot.norm_sqr() > 1e-16 {
                        sqrt_prod / tanh_quot
                    } else {
                        Complex64::new(1e12, 0.0)
                    }
                } else {
                    Complex64::new(1e12, 0.0)
                }
            }
            ElementType::T => {
                // Z = A * coth(beta)/beta + B / (beta * sinh(beta))
                // beta = sqrt(a + j*w*b)
                let a_big = p[0];
                let b_big = p[1];
                let a_small = p[2];
                let b_small = p[3];

                let beta = (Complex64::new(a_small, 0.0) + Complex64::i() * omega * b_small).sqrt();

                if beta.norm_sqr() > 1e-16 {
                    let tanh_beta = beta.tanh();
                    let sinh_beta = beta.sinh();

                    let term1 = if tanh_beta.norm_sqr() > 1e-16 {
                        Complex64::new(a_big, 0.0) / (beta * tanh_beta)
                    } else {
                        Complex64::new(1e12, 0.0)
                    };

                    let term2 = if sinh_beta.norm_sqr() > 1e-16 {
                        Complex64::new(b_big, 0.0) / (beta * sinh_beta)
                    } else {
                        Complex64::new(0.0, 0.0) // Or huge?
                    };

                    term1 + term2
                } else {
                    Complex64::new(1e12, 0.0)
                }
            }
        }
    }
}
