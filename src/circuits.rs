use num_complex::Complex64;
use nom::{
    branch::alt,
    bytes::complete::{tag},
    character::complete::{digit1},
    combinator::{map},
    multi::separated_list1,
    sequence::{delimited, pair},
    IResult,
    Parser,
};
use crate::elements::ElementType;
use crate::elements::Constraint;

/// Trait for calculating impedance.
pub trait Impedance {
    /// Calculates the impedance at a given frequency.
    fn calculate(&self, omega: f64, params: &[f64]) -> Complex64;
    /// Returns the number of parameters required.
    #[allow(dead_code)]
    fn param_count(&self) -> usize;
}

/// Represents a node in the circuit AST (Abstract Syntax Tree).
///
/// A circuit is composed of elements connected in series or parallel.
#[derive(Debug, Clone, PartialEq)]
pub enum CircuitNode {
    /// A single electrochemical element (e.g., R, C, CPE).
    /// (Type, Global Param Index, Label Suffix)
    Element(ElementType, usize, String), 
    /// Elements connected in series (Z_total = Z1 + Z2 + ...).
    Series(Vec<CircuitNode>),
    /// Elements connected in parallel (1/Z_total = 1/Z1 + 1/Z2 + ...).
    Parallel(Vec<CircuitNode>),
}

impl CircuitNode {
    /// Assigns parameter indices to elements in the circuit.
    /// This is useful if we want to map flat parameter arrays to specific elements.
    pub fn assign_indices(&mut self, current_index: &mut usize) {
        match self {
            CircuitNode::Series(nodes) | CircuitNode::Parallel(nodes) => {
                for node in nodes {
                    node.assign_indices(current_index);
                }
            }
            CircuitNode::Element(etype, idx, _) => {
                *idx = *current_index;
                *current_index += etype.param_count();
            }
        }
    }

    /// Returns the total number of parameters in the circuit.
    pub fn count_total_params(&self) -> usize {
        match self {
            CircuitNode::Series(nodes) | CircuitNode::Parallel(nodes) => {
                nodes.iter().map(|n| n.count_total_params()).sum()
            }
            CircuitNode::Element(etype, _, _) => etype.param_count(),
        }
    }

    /// Returns the constraints for all parameters in the circuit.
    pub fn get_constraints(&self) -> Vec<Constraint> {
        match self {
            CircuitNode::Series(nodes) | CircuitNode::Parallel(nodes) => {
                nodes.iter().flat_map(|n| n.get_constraints()).collect()
            }
            CircuitNode::Element(etype, _, _) => etype.constraints(),
        }
    }

    /// Returns the names of all parameters in the circuit.
    pub fn get_param_names(&self) -> Vec<String> {
        match self {
            CircuitNode::Series(nodes) | CircuitNode::Parallel(nodes) => {
                nodes.iter().flat_map(|n| n.get_param_names()).collect()
            }
            CircuitNode::Element(etype, _idx, label) => {
                etype.param_names().iter().map(|&s| format!("{}_{}", s, label)).collect()
            }
        }
    }

    /// Returns the units of all parameters in the circuit.
    pub fn get_param_units(&self) -> Vec<String> {
        match self {
            CircuitNode::Series(nodes) | CircuitNode::Parallel(nodes) => {
                nodes.iter().flat_map(|n| n.get_param_units()).collect()
            }
            CircuitNode::Element(etype, _, _) => {
                etype.param_units().iter().map(|&s| s.to_string()).collect()
            }
        }
    }

    /// Internal evaluation using a parameter iterator.
    ///
    /// This ensures parameters are assigned to elements in the correct order (depth-first).
    fn eval(&self, omega: f64, params: &mut impl Iterator<Item = f64>) -> Complex64 {
        match self {
            CircuitNode::Element(el, _, _) => {
                let count = el.param_count();
                let element_params: Vec<f64> = params.take(count).collect();
                // If we run out of params, use 0.0 (should not happen if validated)
                if element_params.len() < count {
                    return Complex64::new(0.0, 0.0);
                }
                el.calculate(omega, &element_params)
            }
            CircuitNode::Series(nodes) => {
                nodes.iter()
                    .map(|n| n.eval(omega, params))
                    .fold(Complex64::new(0.0, 0.0), |acc, z| acc + z)
            }
            CircuitNode::Parallel(nodes) => {
                let admittance_sum = nodes.iter()
                    .map(|n| {
                        let z = n.eval(omega, params);
                        if z.norm_sqr() > 1e-18 { 1.0 / z } else { Complex64::new(1e12, 0.0) } // Avoid div by zero
                    })
                    .fold(Complex64::new(0.0, 0.0), |acc, y| acc + y);
                
                if admittance_sum.norm_sqr() > 1e-18 {
                    1.0 / admittance_sum
                } else {
                    Complex64::new(1e12, 0.0)
                }
            }
        }
    }
}

impl Impedance for CircuitNode {
    fn calculate(&self, omega: f64, params: &[f64]) -> Complex64 {
        let mut iter = params.iter().cloned();
        self.eval(omega, &mut iter)
    }

    fn param_count(&self) -> usize {
        self.count_total_params()
    }
}

// -----------------------------------------------------------------------------
// Parser (nom)
// -----------------------------------------------------------------------------

fn parse_element_type(input: &str) -> IResult<&str, ElementType> {
    alt((
        map(tag("CPE"), |_| ElementType::CPE),
        map(tag("Wo"), |_| ElementType::Wo),
        map(tag("Ws"), |_| ElementType::Ws),
        map(tag("La"), |_| ElementType::La),
        map(tag("Gw"), |_| ElementType::Gw),
        map(tag("Gs"), |_| ElementType::Gs),
        map(tag("G"), |_| ElementType::G),
        map(tag("K"), |_| ElementType::K),
        map(tag("Zarc"), |_| ElementType::Zarc),
        map(tag("TLMQ"), |_| ElementType::TLMQ),
        map(tag("T"), |_| ElementType::T),
        map(tag("R"), |_| ElementType::R),
        map(tag("C"), |_| ElementType::C),
        map(tag("L"), |_| ElementType::L),
        map(tag("W"), |_| ElementType::W),
    )).parse(input)
}

fn parse_element(input: &str) -> IResult<&str, CircuitNode> {
    map(
        pair(parse_element_type, digit1),
        |(etype, idx_str)| CircuitNode::Element(etype, 0, idx_str.to_string()), // Index assigned later, label stored
    ).parse(input)
}

fn parse_parallel(input: &str) -> IResult<&str, CircuitNode> {
    map(
        delimited(
            tag("p("),
            separated_list1(tag(","), parse_node),
            tag(")"),
        ),
        CircuitNode::Parallel,
    ).parse(input)
}

fn parse_atom(input: &str) -> IResult<&str, CircuitNode> {
    alt((parse_parallel, parse_element)).parse(input)
}

fn parse_node(input: &str) -> IResult<&str, CircuitNode> {
    map(
        separated_list1(tag("-"), parse_atom),
        |nodes| {
            if nodes.len() == 1 {
                nodes[0].clone()
            } else {
                CircuitNode::Series(nodes)
            }
        },
    ).parse(input)
}

/// Parses a circuit string (e.g., "R0-p(C1,R2)") into a CircuitNode AST.
pub fn parse_circuit_string(input: &str) -> Result<CircuitNode, String> {
    let (remaining, mut node) = parse_node(input).map_err(|e| e.to_string())?;
    if !remaining.is_empty() {
        return Err(format!("Unparsed input: {}", remaining));
    }
    let mut idx = 0;
    node.assign_indices(&mut idx);
    Ok(node)
}
